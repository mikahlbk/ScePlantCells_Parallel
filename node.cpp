//node.cpp
//=========================
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
//=========================
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//=========================

//========================================
/** class Node Functions **/
Node::Node(Coord loc) {
	my_loc = loc;
	new_force = Coord();
	int vtk_id;
}
void Node::set_Damping(double new_damping){
	this-> damping = new_damping;
	return;
}
void Node::update_Location() {
	//cout << "New Force after" << new_force << endl;
	//cout << "Location" << my_loc << endl;
	//cout << "damping" << damping << endl;
	//cout << "dt" << dt << endl;
	my_loc += new_force*dt*damping;
	//cout << "updated" << my_loc << endl;
	return;
}
void Node::update_VTK_Id(int id) {
	vtk_id = id;
	return;
}
Node::~Node() {}
//========================================
/**class Cyt Node Functions**/
//constructor
Cyt_Node::Cyt_Node(Coord loc,shared_ptr<Cell> my_cell) : Node(loc) {
	this->my_cell = my_cell;
	return;
}

void Cyt_Node::update_Cell(shared_ptr<Cell> cell){
	this->my_cell = cell;
	return;
}
void Cyt_Node::new_location(Coord location) {
	this->my_loc = location;
	return;
}
void Cyt_Node::calc_Forces(int Ti) {
	//for cytoplasm, just need morse potential for int-int and int-membr
	Coord Fii = calc_Morse_II(Ti);
	//cout << "II" << Fii << endl;
	Coord Fmi = calc_Morse_MI(my_cell->get_Left_Corner(),Ti);
	//cout << "MI" << Fmi << endl;
   	new_force = Fmi + Fii;
	//cout << "New Force before" << new_force << endl;
	return;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell

Coord Cyt_Node::calc_Morse_II(int Ti) {
	//calc force for II
	Coord Fii; //initialized to zero
	if (my_cell==NULL) {
		cout << "Error: Trying to access NULL Pointer. Aborting!" << endl;
		exit(1);
	}
	
	vector<shared_ptr<Cyt_Node>>cyts;
	my_cell->get_Cyt_Nodes_Vec(cyts);
	shared_ptr<Cyt_Node> me= shared_from_this();
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fii) schedule(static,1)
		for (unsigned int j = 0; j < cyts.size(); j++) {
			//don't calculate yourself
			if (cyts.at(j) != me) {
				//calc morse between this node and node j
				Fii +=  me->morse_Equation(cyts.at(j), Ti);
			}
		}
	}

	return Fii;
}

Coord Cyt_Node::calc_Morse_MI(shared_ptr<Wall_Node> orig, int Ti) {
	//calc force for MI
	Coord Fmi;
	if (my_cell==NULL) {
		cout << "Error: Trying to access NULL Pointer. Aborting!" << endl;
		exit(1);
	}
	
	vector<shared_ptr<Wall_Node>> walls;
	this->get_My_Cell()->get_Wall_Nodes_Vec(walls);
	shared_ptr<Cyt_Node> me = shared_from_this();
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fmi) schedule(static,1)
		for(unsigned int i = 0; i < walls.size(); i++) {
			Fmi +=me-> morse_Equation(walls.at(i),Ti);
			//update curr_wall
		}
	}

	return Fmi;
}

Coord Cyt_Node::morse_Equation(shared_ptr<Cyt_Node> cyt, int Ti) {
	if (cyt == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	}
	
	//use Int-Int variables
    	Coord Fii;
    	Coord diff_vect = cyt->get_Location() - my_loc;
    	double diff_len = diff_vect.length();
	double attract = (U_II/xsi_II)*exp(diff_len*(-1)/xsi_II);
    	double repel = (W_II/gamma_II)*exp(diff_len*(-1)/gamma_II);
	Fii = diff_vect*((-attract + repel)/diff_len);
	return Fii;
}

Coord Cyt_Node::morse_Equation(shared_ptr<Wall_Node> wall, int Ti) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	} 

	//use Mem-Int variables
	Coord Fmi;
	Coord diff_vect = wall->get_Location() - my_loc; 
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
   	double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);
     	Fmi = diff_vect*((-attract + repel)/diff_len);
	return Fmi;
}

Cyt_Node::~Cyt_Node() {
	//dont need to do anything
	//because using smart pointer
}


//===============================================================================
/** class Wall Node Functions **/

// Constructors-----------------
Wall_Node::Wall_Node(Coord loc,shared_ptr<Cell> my_cell) : Node(loc) {
	//functions that use this must set
	//left
	//right
	this->my_cell = my_cell;
	//equilibrium length
	//angle
	//klinear
	//kbend
	//equi angle
	//cross prod
	added = 0;
	this->cyt_force = Coord(0,0);
	this->is_connected = 0;
	//this->closest = NULL;
	//this->closest_len = 100;
	//adhesion pairs vec
}

Wall_Node::Wall_Node(Coord loc,shared_ptr<Cell> my_cell, shared_ptr<Wall_Node> left, shared_ptr<Wall_Node> right) : Node(loc)   {
	//functions that use this must set
	this->left = left;
    	this->right = right;
	this-> my_cell = my_cell;
	//equilibrium length
	update_Angle();
	//klinear
	//kbend
	//equi angle
	//cross prod
	added = 1;
	this->cyt_force = Coord(0,0);
	this->is_connected = 0;
	//this->closest = NULL;
	//this->closest_len = 100;
	//adhesion pairs vec
}

Wall_Node::~Wall_Node() {
	//unnecessary since using
	//smartpointers
}

//  Getters and Setters-----------------------------------------------------------
void Wall_Node::set_Left_Neighbor(shared_ptr<Wall_Node> new_Left) {
	this->left = new_Left;
	return;
}
void Wall_Node::set_Right_Neighbor(shared_ptr<Wall_Node> new_Right) {
	this->right = new_Right;
	return;
}
void Wall_Node::update_Cell(shared_ptr<Cell> new_cell) {
	this->my_cell = new_cell;
	return;
}
void Wall_Node::set_membr_len(double length){
	this->membr_equ_len = length;
	return;
}
void Wall_Node::update_Angle() {
	Coord left_vect = get_Left_Neighbor()->get_Location() - get_Location();
	Coord right_vect = get_Right_Neighbor()->get_Location() - get_Location();

	
	double left_len = left_vect.length();
	double right_len = right_vect.length();

	double costheta = left_vect.dot(right_vect) / (left_len * right_len);
	double theta = acos( min( max(costheta,-1.0), 1.0) );

	double crossProd = left_vect.cross(right_vect);

	if (crossProd < 0.0) {
		theta = 2 * pi - theta;
	}
	
	//update protected member variables
	my_angle = theta;
	//cout << "Angle: " <<  theta << endl;
	cross_Prod = crossProd;
	
	return;
}
void Wall_Node::set_K_LINEAR(double k_lin){
	this->K_LINEAR = k_lin;
	return;
}
void Wall_Node::set_K_BEND(double k_bend) {
	this->K_BEND = k_bend;
	return;
}
void Wall_Node::update_Equi_Angle(double new_theta) {
	this->equi_angle = new_theta;
	return;
}
void Wall_Node::set_added(int update){
	this->added = update;
	return;
}
//==========================================================
//Adhesion functions

shared_ptr<Wall_Node> Wall_Node::find_Closest_Node(Coord this_location,shared_ptr<Cell> closest_neighbor, vector<shared_ptr<Wall_Node>> unavailable_nodes){
	//for closest neighbor cell get wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	closest_neighbor->get_Wall_Nodes_Vec(walls);
	//find the closest node
	shared_ptr<Wall_Node> closest_node;
	double curr_dist;
	double smallest_distance = 100;
	bool unavailable = false;
	for(unsigned int i= 0; i<walls.size(); i++) {
		//cout << "computing distance between node of interest and curr wall" << endl;
		curr_dist = (walls.at(i)->get_Location() - this_location).length();
		//cout << "check within adhesion range" << endl;
		if(curr_dist < ADHThresh){
			//cout << "if within adhesion range is it smaller than current possible" << endl;
			if(curr_dist < smallest_distance) {
				//cout << "it is the smallest but will we keep it" << endl;
				for(unsigned int j=0; j < unavailable_nodes.size();j++) {
					if(walls.at(i) == unavailable_nodes.at(j)){
						unavailable = true;
					}
				}
				if(unavailable){
					//dont connect
				}
				else{
					//cout << "closest found in find closest" << endl;
					smallest_distance = curr_dist;
					closest_node = walls.at(i);
				}
			}
		}
	}
	return closest_node;
}
void Wall_Node::make_connection(vector<shared_ptr<Cell>> neighbors) {
	shared_ptr<Wall_Node> this_ptr=shared_from_this();
	Coord this_location = this_ptr->get_Location();
	//neighbor cells passed in
	//find closest neighbor
	shared_ptr<Cell> closest_neighbor;
	Coord curr_center;
	double curr_distance_to_center;
	double smallest_distance_to_center = 100;
	vector<shared_ptr<Wall_Node>> unavailable_nodes;
	//cout << "Finding closest neighbor" << endl;
	for(unsigned int i= 0; i < neighbors.size() ; i++) {
		//cout << "getting center of current neighbor" << endl;
		curr_center = neighbors.at(i)->get_Cell_Center();
		//cout << "location from current node to center" << endl;
		curr_distance_to_center = (curr_center - this_ptr->get_Location()).length();
		//cout << "if statement" << endl;
		if(curr_distance_to_center < smallest_distance_to_center) {
			smallest_distance_to_center = curr_distance_to_center;
			closest_neighbor = neighbors.at(i);
		}
	}
	//cout << "Closest neighbor " << closest_neighbor << endl;
	//cout<< "out of if" << endl;
	//find the closest node on that neighbor cell within adhesion range
	shared_ptr<Wall_Node> possible_connection;
	//cout << "find closest node" << endl;
	possible_connection = find_Closest_Node(this_location,closest_neighbor,unavailable_nodes);
	//cout << "found closest node first " << possible_connection <<  endl;
	int counter = 0;
	//Determine connection
	//is this node already connected
	bool not_connected = true;
	//cout << "determining connection" << endl;
	do{
		if(possible_connection != NULL) {
		//cout << "back at top" << endl;
		vector<pair<double, shared_ptr<Wall_Node>>> adh_pairs = possible_connection->get_adhesion_vec();
		//cout << "adhesion vector size prior to other cell: " << adh_pairs.size() << endl;
		double curr_dist_adhesion = (possible_connection->get_Location() -this_ptr->get_Location()).length();
	
		if(adh_pairs.size() == 0){
			//no ----> connect
			//cout << "Connected, curr dist adhesion" << curr_dist_adhesion <<endl;
			possible_connection->adh_push_back(curr_dist_adhesion,this_ptr);
			possible_connection->set_is_connected(1);
			this_ptr->adh_push_back(curr_dist_adhesion,possible_connection);
			this_ptr->set_is_connected(1);
			not_connected = false;
		}
		else{
			//yes ---> am i closer?
			//cout << "curr dist adhesion: " << curr_dist_adhesion << endl;
			//cout << "adh curr dist: " << adh_pairs.at(0).first << endl;
			if(curr_dist_adhesion < adh_pairs.at(0).first){
				//yes --->  make connection
				//cout << "Connected" << endl;
				//cout << "current connect" << adh_pairs.at(0).second << endl;
				adh_pairs.at(0).second->clear_adh_vec();
				adh_pairs.at(0).second->set_is_connected(0);
				possible_connection->clear_adh_vec();
				possible_connection->adh_push_back(curr_dist_adhesion,this_ptr);
				possible_connection->set_is_connected(1);
				this_ptr->adh_push_back(curr_dist_adhesion,possible_connection);
				this_ptr->set_is_connected(1);
				not_connected = false;
			}
			else {
				//go through this up to fifth choice then quit
				counter++;
				//cout << "unavailale: " << counter << endl;
				if(counter < 3) {
					//no ---> update unavailable vector
					//cout << "pushed back unavailable" << endl;
					unavailable_nodes.push_back(possible_connection);
					//find next closest
					//cout << "find next closest" << endl;
					possible_connection = find_Closest_Node(this_location,closest_neighbor,unavailable_nodes);
					//cout << "found next closest" << possible_connection << endl;
				}
				else{
					not_connected = false;
				}
			}
		}
		}
		else {
			not_connected = false;
		}
	}while(not_connected);
	//cout << "out of while" << endl;
	return;
}
/*void Wall_Node::set_Closest(shared_ptr<Wall_Node>  closest, double closest_len) {
	this->closest = closest;
	this->closest_len = closest_len;
	return;
}*/
void Wall_Node::set_is_connected(int is_connected_truth){
	this->is_connected = is_connected_truth;
	return;
}
void Wall_Node::clear_adh_vec(){
	this->adhesion_pairs.clear();	
	return;
}
void Wall_Node::adh_push_back(double distance, shared_ptr<Wall_Node> closest){
	this->adhesion_pairs.push_back(make_pair(distance,closest));
	return;
}
void Wall_Node::remove_from_adh_vecs(){
	shared_ptr<Wall_Node> curr_closest;
	shared_ptr<Wall_Node> me = shared_from_this();
	vector<pair<double, shared_ptr<Wall_Node>>> adh_pairs;

	for(unsigned int i=0; i< adhesion_pairs.size();i++) {
		curr_closest = adhesion_pairs[i].second;
		adh_pairs = curr_closest->get_adhesion_vec();
		curr_closest->clear_adh_vec();
		for(unsigned int j = 0; j < adh_pairs.size();j++){
			if(adh_pairs[j].second != me){
				curr_closest->adh_push_back(adh_pairs[j].first, adh_pairs[j].second);
			}
		}
	}
	return;
}
/*void Wall_Node::clear_closest_in_adh_vec(){
	cout << "undo closest" << endl;
	for(unsigned int i = 0; i<adhesion_pairs.size(); i++){
		adhesion_pairs.at(i)->set_Closest(NULL, 100);
	}
	return;
}*/
//===========================================================
// Calc Force Functions -----------------------
//calculates total force on current wall node
void Wall_Node::calc_Forces(int Ti) {
	//holds the current force from all
	//interacting nodes
	Coord sum;
	//start with forces from nodes
	//in the same cell
	//calculates force from internal nodes
	//on current wall node
	sum += calc_Morse_SC(Ti);
	//cout << "SC success" << endl;
	//cout << sum << endl;
	//each wall node wiil keep 
	//track of force from cytoplasm
	//for pressure measurements	
	this->cyt_force = sum;
	//cout << "cyt force success" << endl;
	//linear force from left/right neighbor
	sum += calc_Linear();
	//cout << "linear success" << endl;
	//bending force from left/right neighbor
	sum += calc_Bending();
	//cout << "bending sucess"<< endl;
	//different cell forces on wall nodes 
	//from a different cell
	//morse between neighboring cell
	//adhesion between neighboring cell
	sum += calc_Morse_DC(Ti);
	//cout << "DC Success" << calc_Morse_DC(Ti) << endl;
	new_force = sum;
	return;
}
//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC(int Ti) {
	vector<shared_ptr<Cyt_Node>> cyt_nodes;
	my_cell->get_Cyt_Nodes_Vec(cyt_nodes);
	shared_ptr<Wall_Node> curr_wall= shared_from_this();
	Coord Fmi;
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fmi) schedule(static,1)
		for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
			Fmi += curr_wall->morse_Equation(cyt_nodes.at(i), Ti);
		}
	}

	//cout << "morse_sc: " << Fmi << endl;	
	return Fmi;
}
//linear spring force of neighboring springs
Coord Wall_Node::calc_Linear() {
	Coord F_lin;

//	cout << "calc left" << endl;
	F_lin += linear_Equation(left);

//	cout << "calc right" << endl;
	F_lin += linear_Equation(right);
	return F_lin;
}
//bending force of node
Coord Wall_Node::calc_Bending() {
	Coord F_bend;
	F_bend += bending_Equation_Center();
	//cout << "left" << endl;
	F_bend += bending_Equation_Left();
	//cout << "right" << endl;
	F_bend += bending_Equation_Right();
	//cout << "done" << endl;
	if (cross_Prod < 0.0) {
		F_bend = F_bend*(-1);
	}	
	return F_bend;
}
//morse potential between wall node i and every wall node in neighboring cell
Coord Wall_Node::calc_Morse_DC(int Ti) {
	Coord Fdc;
	vector<shared_ptr<Cell>> cells;
	my_cell->get_Neighbor_Cells(cells);	
	//cout << "Neighbor cells: " << cells.size()<<endl;
	//cout << "getting neighbors" << endl;
	#pragma omp parallel 
	{
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fdc) schedule(static,1) 
		for (unsigned int i = 0; i < cells.size(); i++) {
			Fdc += neighbor_nodes(cells.at(i), Ti);
		}
	}
	//cout << "Fdc" << Fdc << endl;
	//cout << "adhesion" << endl;
	//if(this->closest != NULL){
	//	cout << "closest not null" << endl;
	//	closest->get_Location();
	//	cout << "got location" << endl;
	//	Fdc += this->linear_Equation_ADH(this->closest);
	//	cout << "computed adh successfully" << endl;
	//}
	for(unsigned int i = 0; i < adhesion_pairs.size(); i++){
		Fdc += this->linear_Equation_ADH(adhesion_pairs[i].second);
	}
	return Fdc;
}
//function to get all wall nodes of neighbor cell i
Coord Wall_Node::neighbor_nodes(shared_ptr<Cell> neighbor, int Ti) {
	Coord sum;
	vector<shared_ptr<Wall_Node>> walls;
	neighbor->get_Wall_Nodes_Vec(walls);
	//cout << "Number walls" << walls.size() << endl;
	shared_ptr<Wall_Node> me = shared_from_this();
	#pragma omp parallel
	{
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:sum) schedule(static,1) 
		for(unsigned int j =0; j< walls.size(); j++) {
			sum += me->morse_Equation(walls.at(j), Ti);
		}
	}
	//cout<< "Sum: " << sum << endl;
	return sum;
}			
//===========================================================
// Mathematical force calculations
Coord Wall_Node::morse_Equation(shared_ptr<Cyt_Node> cyt, int Ti) {
	if (cyt == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	}
	
	//use Membr - int variables
	Coord Fmi;
	Coord diff_vect = cyt->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
	double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);

	Fmi = diff_vect*((-attract + repel)/diff_len);
	
	//cout << Fmi << endl;
	return Fmi;
}

Coord Wall_Node::morse_Equation(shared_ptr<Wall_Node> wall, int Ti) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	}

	//use Mem-Mem variables
   	Coord Fmmd;
   	Coord diff_vect = wall->get_Location() - my_loc;
    	double diff_len = diff_vect.length();
      	double attract = (U_MM/xsi_MM)*exp(diff_len*(-1)/xsi_MM);
    	double repel = (W_MM/gamma_MM)*exp(diff_len*(-1)/gamma_MM);
   	
	Fmmd = diff_vect*((-attract + repel)/diff_len);
	
	//cout << Fmmd << endl;
	return Fmmd;
}
Coord Wall_Node::bending_Equation_Center() {
	Coord F_center;
	double self_Constant; 
	
	double eps = 0.0001;

	if (abs(my_angle - pi) < eps) {
		return F_center;
	}
	else {
		self_Constant = K_BEND*(my_angle - equi_angle)/(sqrt(1-pow(cos(my_angle),2)));
	}
	Coord left_vect = left->get_Location() - my_loc;
	Coord right_vect = right->get_Location() - my_loc;
	double left_len = left_vect.length();
	double right_len = right_vect.length();
	Coord term_l1 = (left_vect*(-1))/(left_len*right_len);
	Coord term_l2 = left_vect*cos(my_angle)/pow(left_len,2);
	Coord term_r1 = (right_vect*(-1))/(left_len*right_len);
	Coord term_r2 = right_vect*cos(my_angle)/pow(right_len,2);

	F_center = (term_l1 + term_l2 + term_r1 + term_r2) * self_Constant;
	
	//cout << "Bending center: " << F_center << endl;	
	return F_center;
}
Coord Wall_Node::bending_Equation_Left() {
	Coord F_left;
	double K_BEND = left->get_K_BEND();
	double left_equi_angle = left->get_Equi_Angle();
	double left_angle = left->get_Angle();
	double left_Constant;
	
	double eps = 0.0001;

	if (abs(left_angle - pi) < eps) {
		return F_left;
	}
	else {
		left_Constant = K_BEND*(left_angle - left_equi_angle)/(sqrt(1-pow(cos(left_angle),2)));
	}
		

	Coord left_vect = left->get_Location() - my_loc;
	double left_len = left_vect.length();
	Coord left_left_vect = left->get_Left_Neighbor()->get_Location()-left->get_Location();
	double left_left_len = left_left_vect.length();
	Coord left_left_term1 = left_left_vect/(left_left_len*left_len);
	Coord left_term2 = left_vect*cos(left_angle)/pow(left_len,2);

	F_left = (left_left_term1 + left_term2) * left_Constant;
	
	//cout << "Bending left: " << F_left << endl;
	return F_left;
}

Coord Wall_Node::bending_Equation_Right() {
	Coord F_right;
	double K_BEND = right->get_K_BEND();
	double right_equ_angle = right->get_Equi_Angle();
	double right_angle = right->get_Angle();
	double right_Constant;
	
	double eps = 0.0001;

	if (abs(right_angle - pi) < eps) {
		return F_right;
	}
	else{
		right_Constant = K_BEND*(right_angle-right_equ_angle)/(sqrt(1-pow(cos(right_angle),2)));
		//}
	}

	Coord right_vect = right->get_Location() - my_loc;
	double right_len = right_vect.length();
	Coord right_right_vect = right->get_Right_Neighbor()->get_Location()-right->get_Location();
	double right_right_len = right_right_vect.length();
	Coord right_right_term1 = right_right_vect/(right_right_len*right_len);
	Coord right_term2 = right_vect*cos(right_angle)/pow(right_len,2);

	F_right = (right_right_term1 + right_term2)*right_Constant;
	
	//cout << "Bending right: " << F_right << endl;
	return F_right;
}
Coord Wall_Node::linear_Equation(shared_ptr<Wall_Node> wall) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	}
	
	//use spring constant variables
	Coord F_lin;
	Coord diff_vect = wall->get_Location() - my_loc;
	double diff_len = diff_vect.length();

	F_lin = (diff_vect/diff_len)*K_LINEAR*(diff_len - this->membr_equ_len);
	//cout << "linear" << F_lin << endl;
	return F_lin;	
}

Coord Wall_Node::linear_Equation_ADH(shared_ptr<Wall_Node>& wall) {
	if (wall == NULL) {
		cout << "Error: Trying to access NULL pointer Aborting!" << endl;
		exit(1);
	};
	Coord F_lin;
//	cout << "compute diff vec" << endl;
	Coord wall_loc = wall->get_Location();
//	cout << "wall loc"  << endl;
	Coord loc = my_loc;
//	cout << "my loc " << endl;
	Coord diff_vect = wall_loc - loc;
//	cout << "coord diff is : " << diff_vect << endl;
	double diff_len = diff_vect.length();
	F_lin = (diff_vect/diff_len)*(K_ADH*(diff_len - MembrEquLen_ADH));
	return F_lin;
}
//==========================================================
// End of node.cpp
