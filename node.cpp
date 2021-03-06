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
//determines which nodes on neighbor cells are within adhesion
//threshold and pushes them onto connection vector
//number of connections specified in phys.h
void Wall_Node::make_connection(vector<shared_ptr<Wall_Node>> neighbor_walls) {
	shared_ptr<Wall_Node> this_ptr = shared_from_this();
	Coord this_ptr_loc = this_ptr->get_Location();
	vector<shared_ptr<Wall_Node>> this_ptr_adh_vec;
	this_ptr_adh_vec = this_ptr->get_adh_vec();
	Coord neighbor_node_loc;
	double biggest_dist;
	double curr_distance;
	int counter = 0;
	//cout << "neighbor walls size : " << neighbor_walls.size() << endl;
	for(unsigned int i= 0; i < neighbor_walls.size() ; i++) {
		counter++;
		//cout << "neighbor wall: " << counter << endl;
		neighbor_node_loc = neighbor_walls.at(i)->get_Location();
		curr_distance = (this_ptr_loc - neighbor_node_loc).length();
			if(curr_distance < ADHThresh){
			//	cout << "within adh thresh" << endl;
				if(this_ptr_adh_vec.size() < NUMBER_ADH_CONNECTIONS){
					this_ptr->adh_push_back(neighbor_walls.at(i));
					this_ptr_adh_vec = this_ptr->get_adh_vec();
				}
				else{
					//sort in descending order
					reverse(this_ptr_adh_vec.begin(),this_ptr_adh_vec.end());
					biggest_dist = (this_ptr_loc - this_ptr_adh_vec.at(0)->get_Location()).length();
					
					if(curr_distance < biggest_dist){
						this_ptr->update_adh_vec(neighbor_walls.at(i));
						this_ptr_adh_vec = this_ptr->get_adh_vec();
					}
				}
			}
		//cout << "adhesion vec size: " << this_ptr_adh_vec.size()<< endl;
	}
	return;
}
//ensures that if a node has made a connection
//both nodes apply force to each other
void Wall_Node::one_to_one_check(){
	shared_ptr<Wall_Node> this_ptr = shared_from_this();
	vector<shared_ptr<Wall_Node>> this_adh_vec;
	this_adh_vec = this_ptr->adhesion_vector;
	vector<shared_ptr<Wall_Node>> connection_adh_vec;
	for(unsigned int i = 0; i < this_adh_vec.size(); i++){
		connection_adh_vec = this_adh_vec.at(i)->get_adh_vec();
		if(find(connection_adh_vec.begin(), connection_adh_vec.end(), this_ptr) != connection_adh_vec.end()){
			//do nothing its already in there
		}
		else {
			this_adh_vec.at(i)->adh_push_back(this_ptr);
		}
	}
	return;
}
//clears adhesion vector of current node
void Wall_Node::clear_adhesion_vec(){
	this->adhesion_vector.clear();	
	return;
}
//push back a cell wall node onto adhesion vector 
//of current node
void Wall_Node::adh_push_back(shared_ptr<Wall_Node> neighbor_node){
	this->adhesion_vector.push_back(neighbor_node);
	return;
}
void Wall_Node::update_adh_vec(shared_ptr<Wall_Node> node) {
	reverse(this->adhesion_vector.begin(),this->adhesion_vector.end());
	this->adhesion_vector.at(0) = node;
	return;
}
//removes the current node from adhesion vector of
//cell wall nodes from neighboring cells 
void Wall_Node::remove_from_adh_vecs(){
	shared_ptr<Wall_Node> me = shared_from_this();
	vector<shared_ptr<Wall_Node>> neighbor_connections;
	for(unsigned int i = 0; i < adhesion_vector.size(); i++){
		neighbor_connections = adhesion_vector.at(i)->get_adh_vec();
		adhesion_vector.at(i)->clear_adhesion_vec();
		for(unsigned int j = 0; j<neighbor_connections.size();j++){
			if(neighbor_connections.at(j) != me){
				adhesion_vector.at(i)->adh_push_back(neighbor_connections.at(j));
			}
		}
	}
	return;
}
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
	vector<shared_ptr<Wall_Node>> adhesion_pairs;
	adhesion_pairs = this->get_adh_vec();
	for(unsigned int i = 0; i < adhesion_pairs.size(); i++){
		Fdc += this->linear_Equation_ADH(adhesion_pairs.at(i));	
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
	if(this->get_My_Cell()->get_Layer() == 1){
		F_lin = (diff_vect/diff_len)*(K_ADH_L1*(diff_len - MembrEquLen_ADH));
	}
	else if(this->get_My_Cell()->get_Layer() == 2){
		F_lin = (diff_vect/diff_len)*(K_ADH_L2*(diff_len - MembrEquLen_ADH));
	}
	else{
		F_lin = (diff_vect/diff_len)*(K_ADH*(diff_len - MembrEquLen_ADH));
	}
	return F_lin;
}
//==========================================================
// End of node.cpp
