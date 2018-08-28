//node.cpp
//=========================
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
//=========================
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
//=========================

//========================================
/** class Node Functions **/
Node::Node(Coord loc) {
	my_loc = loc;
	new_force = Coord();
	int vtk_id;
}

Coord Node::get_Location() {
	return my_loc;
}

Coord Node::get_Force() {
	return new_force;
}
void Node::update_VTK_Id(int id) {
	vtk_id = id;
	return;
}
void Node::update_Location(double& new_damping) {
	my_loc += new_force*dt*new_damping;
	return;
}

Node::~Node() {}
//========================================
/**class Cyt Node Functions**/
//constructor
Cyt_Node::Cyt_Node(Coord loc, Cell*  my_cell) : Node(loc) {
	this->my_cell = my_cell;
	return;
}

void Cyt_Node::update_Cell(Cell* cell){
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
	Coord Fmi = calc_Morse_MI(my_cell->get_Wall_Nodes(),Ti);
   	new_force = Fmi + Fii;
	return;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell

Coord Cyt_Node::calc_Morse_II(int Ti) {
	//calc force for II
	Coord Fii; //initialized to zero

	vector<shared_ptr<Cyt_Node>>cyts;
	if (my_cell==NULL) {
		cout << "Error: Trying to access NULL Pointer. Aborting!" << endl;
		exit(1);
	}
	//cout << "gonna make the pointer" << endl;
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
	//cout << "oaralelel" << endl;

	return Fii;
}

Coord Cyt_Node::calc_Morse_MI(shared_ptr<Wall_Node> orig, int Ti) {
	//calc force for MI
	Coord Fmi;
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
	//int div_time = this->get_My_Cell()->get_Cell_Progress_div();
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
//	int div_time = this->get_My_Cell()->get_Cell_Progress_div();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
   	double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);
     	Fmi = diff_vect*((-attract + repel)/diff_len);
	return Fmi;
}

Cyt_Node::~Cyt_Node() {
//	my_cell = NULL;
}


//===============================================================================
/** class Wall Node Functions **/

// Constructors-----------------
Wall_Node::Wall_Node(Coord loc,Cell* my_cell) : Node(loc) {
	//functions that use this must set
	//left
	//right
	this->my_cell = my_cell;
	//set equilibrium length
	//set angle
	//set klinear
	//set equi angle
	//cross prod
	this->cyt_force = Coord(0,0);
	this->closest = NULL;
	this->curr_closest = NULL;
	//this->microfibril_pair = NULL;
	//this->pressure = 0;
	this->closest_len = 100;
	this->curr_slope = 100;
	this->is_added = false;
	//this->is_delete = false;
}

Wall_Node::Wall_Node(Coord loc,Cell* my_cell, shared_ptr<Wall_Node> left, shared_ptr<Wall_Node> right) : Node(loc)   {
	this->left = left;
    	this->right = right;
	this-> my_cell = my_cell;
	//equilibrium length
	update_Angle();
	//set klinear
	//set equi anglei
	//cross prod
	this->cyt_force = Coord(0,0);

	this->closest = NULL;
	this->curr_closest = NULL;
	//this->microfibril_pair= NULL;
	this->curr_slope = 100;
	this->closest_len = 100;
	//this->is_added;
	//this->pressure = 0;
	//this->is_new = true;
	//this->is_delete = false;
}

Wall_Node::~Wall_Node() {
//	my_cell = NULL;	
//	left = NULL;
//	right = NULL;
//	//microfibril_pair = NULL;
//	closest = NULL;
}

//  Getters and Setters-----------------------------------------------------------
void Wall_Node::set_is_new(){
	this->is_added = true;
	return;
}
void Wall_Node::set_Left_Neighbor(shared_ptr<Wall_Node> new_Left) {
	this->left = new_Left;
	return;
}
//void Wall_Node::set_pressure(double& new_press){
//	this->pressure = new_press;
//	return;
//}
void Wall_Node::set_Right_Neighbor(shared_ptr<Wall_Node> new_Right) {
	this->right = new_Right;
	return;
}
//void Wall_Node::set_Delete(int y) {
//	if(y == 0){
//		this->is_delete = false;
//	}
//	if(y == 1){
//		this->is_delete = true;
//	}
//	return;
//}
//void Wall_Node::set_is_new(bool yes) {
//	this->is_new = yes;
//	return;
//}
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
void Wall_Node::set_membr_len(double length){
	this->membr_equ_len = length;
	return;
}
void Wall_Node::update_Cell(Cell* new_cell) {
	this->my_cell = new_cell;
	return;
}
void Wall_Node::update_Equi_Angle(double new_theta) {
	this->equi_angle = new_theta;
	return;
}
void Wall_Node::set_K_LINEAR(double& k_lin){
	this->K_LINEAR = k_lin;
	return;
}
void Wall_Node::clear_adh_vec(){
	this->adhesion_pairs.clear();	
	return;
}
void Wall_Node::add_adh_pair(shared_ptr<Wall_Node> pair) {
	this->adhesion_pairs.push_back(pair);
	return;
}
//void Wall_Node::clear_Closest_Vec() {
//	while ( !closest_vec.empty()) {
//		closest_vec.pop_back();
//	}
//}

//void Wall_Node::set_Closest_Vec(Wall_Node* closest) {
//	this->closest_vec.push_back(closest);
//	return;
//}
void Wall_Node::set_curr_Closest(shared_ptr<Wall_Node> curr_closest){
	this->curr_closest = curr_closest;
	return;
}
void Wall_Node::set_Closest(shared_ptr<Wall_Node>  closest, double closest_len) {
	this->closest = closest;
	this->closest_len = closest_len;
	return;
}
//void Wall_Node::set_microfibril_pair(Wall_Node* pair,double curr_slope){
//	this->microfibril_pair = pair;
//	this->curr_slope = curr_slope;
//	return;
//}
void Wall_Node::print_info(){
		cout << "left" <<this->left<<endl;
		cout << "right" << this->right<<endl;
        	cout << "cell" << this->my_cell <<endl;
		cout << "membrane length" << this-> membr_equ_len << endl;
		cout << "angle" << this->my_angle<<endl;
		cout << "linear constatn" <<this->K_LINEAR<<endl;
		cout << "equilibrium" <<this->equi_angle<<endl;
		cout << "cross" <<this->cross_Prod << endl;
		cout << "cyt" << this->cyt_force<<endl;
		cout << "closest" <<this->closest<<endl;
		cout << "slope" <<this->curr_slope<<endl;
		cout << "closest" <<this->closest_len<<endl;
//		cout << "micro" <<this->microfibril_pair<<endl;
		cout << "added" << this->is_added <<endl;
	return;
}
	

// Calc Force Functions -----------------------
//calculates total force on current wall node
void Wall_Node::calc_Forces(int Ti) {
	Coord sum;
	//same cell forces on wall nodes
	//from cytoplasm nodes
	/*if(Ti>15000){
	Coord sum;
	//same cell forces on wall nodes
	//from cytoplasm nodes
	sum += calc_Morse_SC(Ti);
	cout << "SC success" << endl;
	//linear force from left/right neighbor
	sum += calc_Linear();
	cout << "linear" << endl;
	//bending force from left/right neighbor
	sum += calc_Bending();
//	cout << "bending"<< endl;
	//each wall node wiil keep 
	//track of force from cytoplasm
	//for pressure measurements	
	this->cyt_force = sum;
//	cout << "cyt force" << endl;
	//different cell forces on wall nodes
	//morse between neighboring cell
	//adhesion between neighboring cell
	sum += calc_Morse_DC(Ti);
//	cout << "DC" << calc_Morse_DC(Ti) << endl;
	}
	else{*/
	sum += calc_Morse_SC(Ti);
	cout << "SC success" << endl;
	//linear force from left/right neighbor
	sum += calc_Linear();
	cout << "linear" << endl;
	//bending force from left/right neighbor
	sum += calc_Bending();
	cout << "bending"<< endl;
	//each wall node wiil keep 
	//track of force from cytoplasm
	//for pressure measurements	
	this->cyt_force = sum;
	cout << "cyt force" << endl;
	//different cell forces on wall nodes
	//morse between neighboring cell
	//adhesion between neighboring cell
	sum += calc_Morse_DC(Ti);
	cout << "DC" << calc_Morse_DC(Ti) << endl;
	//}
	new_force = sum;

	return;
}

//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC(int Ti) {
	vector<shared_ptr<Cyt_Node>> cyt_nodes;
	my_cell->get_Cyt_Nodes_Vec(cyt_nodes);
	shared_ptr<Wall_Node> curr_wall= shared_from_this();
	Coord Fmi;
//	cout << "Cyt nodes" << endl;
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fmi) schedule(static,1)
		for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
			Fmi += curr_wall->morse_Equation(cyt_nodes.at(i), Ti);
		}
	}
	//if(this->microfibril_pair != NULL){
		//cout << "closest not null" << endl;
	//	this->microfibril_pair->get_Location();
		//cout << "got location" << endl;
	//	Fmi += linear_Equation_microfibril(this->microfibril_pair);
		//cout << "computed adh successfully" << endl;
	//}
	
//	cout << "morse_sc: " << Fmi << endl;	
	return Fmi;
}


//morse potential between wall node i and every wall node in neighboring cell
Coord Wall_Node::calc_Morse_DC(int Ti) {
	Coord Fdc;
	vector<shared_ptr<Cell>> cells;
	my_cell->get_Neighbor_Cells(cells);	
	//cout << "getting neighbors" << endl;
	shared_ptr<Wall_Node> curr = NULL;
	shared_ptr<Wall_Node> orig = NULL;
	cout << "meighbor nodes" << endl;
	#pragma omp parallel 
	{
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fdc) schedule(static,1) 
		for (unsigned int i = 0; i < cells.size(); i++) {
			Fdc += neighbor_nodes(cells.at(i), Ti);
		}
	}
	cout << "adhesion" << endl;
	if(this->closest != NULL){
	//	cout << "closest not null" << endl;
		closest->get_Location();
	//	cout << "got location" << endl;
		Fdc += linear_Equation_ADH(this->closest);
	//	cout << "computed adh successfully" << endl;
	}
	return Fdc;
}
//function to get all wall nodes of neighbor cell i
Coord Wall_Node::neighbor_nodes(shared_ptr<Cell> neighbor, int Ti) {
	Coord sum;
	vector<shared_ptr<Wall_Node>> walls;
	neighbor->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:sum) schedule(static,1) 
		for(unsigned int j =0; j< walls.size(); j++) {
			//cout << "getting wall nodes" << endl;
			sum += morse_Equation(walls.at(j), Ti);
			//cout << "morse" << endl;
		}
	}
	return sum;
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

//linear spring force of neighboring springs
Coord Wall_Node::calc_Linear() {
	Coord F_lin;

//	cout << "calc left" << endl;
	F_lin += linear_Equation(left);

//	cout << "calc right" << endl;
	F_lin += linear_Equation(right);
	return F_lin;
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
		//if((this->get_My_Cell()->get_Layer() == 1) || (this->get_My_Cell()->get_Layer() ==2)){
		//	self_Constant = K_BEND_L1*(my_angle - equi_angle)/(sqrt(1-pow(cos(my_angle),2)));
		//}
		//else {
			self_Constant = K_BEND*(my_angle - equi_angle)/(sqrt(1-pow(cos(my_angle),2)));
		//}
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
	//double left_k_bend = left->get_Bending_Spring();
	double left_equi_angle = left->get_Equi_Angle();
	double left_angle = left->get_Angle();
	double left_Constant;
	
	double eps = 0.0001;

	if (abs(left_angle - pi) < eps) {
		return F_left;
	}
	else {
		//if((this->get_My_Cell()->get_Layer() == 1) || (this->get_My_Cell()->get_Layer() ==2)){
		//	left_Constant = K_BEND_L1*(left_angle - left_equi_angle)/(sqrt(1-pow(cos(left_angle),2)));
		//}
		//else {
			left_Constant = K_BEND*(left_angle - left_equi_angle)/(sqrt(1-pow(cos(left_angle),2)));
		//}
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
	//double right_k_bend = right->get_Bending_Spring();
	double right_equ_angle = right->get_Equi_Angle();
	double right_angle = right->get_Angle();
	double right_Constant;
	
	double eps = 0.0001;

	if (abs(right_angle - pi) < eps) {
		return F_right;
	}
	else {
		//if((this->get_My_Cell()->get_Layer() == 1) || (this->get_My_Cell()->get_Layer() ==2)) {
		///	right_Constant = K_BEND_L1*(right_angle-right_equ_angle)/(sqrt(1-pow(cos(right_angle),2)));
		//}
		//else{
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
	Coord diff_vect = wall->get_Location() - my_loc;
//	cout << "coord diff is : " << diff_vect << endl;
	double diff_len = diff_vect.length();
	if((wall->get_My_Cell()->get_Layer() == 1) && (this->get_My_Cell()->get_Layer() == 1)) {
		F_lin = (diff_vect/diff_len)*(K_ADH_L1*(diff_len - MembrEquLen_ADH));
	}
	else {
		F_lin = (diff_vect/diff_len)*(K_ADH*(diff_len - MembrEquLen_ADH));
	}

	return F_lin;
}
Coord Wall_Node::linear_Equation_microfibril(shared_ptr<Wall_Node>& wall) {
	//cout << "wall node is: " << wall << endl;
//	Wall_Node* wall = NULL;
//	for(unsigned int i = 0;i < closest_nodes.size();i++){
//		wall = closest_nodes.at(i);
//	if (wall == NULL) {
//		cout << "Problems for days" << endl;
//	};
	Coord F_lin;
//	cout << "compute diff vec" << endl;
	Coord wall_loc = wall->get_Location();
//	cout << "wall loc"  << endl;
	Coord loc = my_loc;
////	cout << "my loc " << endl;
	Coord diff_vect = wall->get_Location() - my_loc;
//	cout << "coord diff is : " << diff_vect << endl;
	double diff_len = diff_vect.length();
	F_lin = (diff_vect/diff_len)*(K_microfibril*(diff_len - MembrEquLen_microfibril));
//	}

	return F_lin;
}
//==========================================================
//Adhesion functions

shared_ptr<Wall_Node> Wall_Node::find_Closest_Node(vector<shared_ptr<Cell>>& neighbors) {
	shared_ptr<Wall_Node> curr = NULL;
	shared_ptr<Wall_Node> orig = NULL;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Cell> curr_cell = NULL;
	shared_ptr<Wall_Node> closest = NULL;
;	double curr_dist = 0;
	double smallest = 100;
	for(int i = 0; i < neighbors.size(); i++) {
		curr_cell = neighbors.at(i);
		//find the closest node on curr_Side
		curr = curr_cell->get_Left_Corner();
		orig = curr;
		do{
			next = curr->get_Left_Neighbor();
			curr_dist = (this->my_loc - curr->get_Location()).length();
			if(curr_dist < ADHThresh) {
				if(curr_dist < smallest) {
					closest = curr;
					smallest = curr_dist;
				}
			}
			curr = next;
		} while (next != orig);
	}
	return closest;
}


void Wall_Node::make_Connection(shared_ptr<Wall_Node> curr_Closest) {
	double curr_dist = 0;
	shared_ptr<Wall_Node> this_ptr=shared_from_this();
	if(curr_Closest != NULL) {
	//	if(curr_Closest->get_curr_Closest() == this){
			this_ptr->closest = curr_Closest;
			curr_Closest->add_adh_pair(this_ptr);
			curr_dist = (this_ptr->get_Location() - curr_Closest->get_Location()).length();
			this_ptr->closest_len = curr_dist;
	//	}	
	}
	return;
}
/*void Wall_Node::find_microfibril_pair_horiz(vector<Wall_Node*> side2){
	double x_1;
	double x_2;
	double y_1;
	double y_2;
	double slope;
	Wall_Node* curr = NULL;
	double smallest;
//
	for(unsigned int i=0;i<side2.size();i++){
		x_1 = this->get_Location().get_X();
		y_1 = this->get_Location().get_Y();
		x_2 = side2.at(i)->get_Location().get_X();
		y_2 = side2.at(i)->get_Location().get_Y();
		slope = (x_2-x_1)/(y_2-y_1);
		//cout << "slope" << slope << endl;
		if(side2.at(i)->get_micro_pair() == NULL) {
			if(slope < this->curr_slope){
				smallest = slope;
				curr = side2.at(i);
			}
		}
	}
	if(curr != NULL) {
		this->curr_slope = smallest;
		this->microfibril_pair = curr;
		curr->set_microfibril_pair(this,smallest);
	}

	return;
}
void Wall_Node::find_microfibril_pair_vert(vector<Wall_Node*> side2){
	double x_1;
	double x_2;
	double y_1;
	double y_2;
	double slope;
	Wall_Node* curr = NULL;
	double smallest;
//
	for(unsigned int i=0;i<side2.size();i++){
		x_1 = this->get_Location().get_X();
		y_1 = this->get_Location().get_Y();
		x_2 = side2.at(i)->get_Location().get_X();
		y_2 = side2.at(i)->get_Location().get_Y();
		slope = (y_2-y_1)/(x_2-x_1);
//	//	cout << "slope" << slope << endl;
		if(side2.at(i)->get_micro_pair() == NULL) {
			if(slope < this->curr_slope){
				smallest = slope;
				curr = side2.at(i);
			}
		}
	}
	if(curr != NULL) {
		this->curr_slope = smallest;
		this->microfibril_pair = curr;
		curr->set_microfibril_pair(this,smallest);
	}
//
	return;
}*/
//==========================================================
// End of node.cpp

