//node.cpp
//=========================
#include <iostream>
#include <vector>
#include <cmath>
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
Cyt_Node::Cyt_Node(Coord loc, Cell* my_cell) : Node(loc) {
	this->my_cell = my_cell;
//	isStationary = false;

}

void Cyt_Node::calc_Forces() {
	//for cytoplasm, just need morse potential for int-int and int-membr
	Coord Fii = calc_Morse_II();
	Coord Fmi = calc_Morse_MI(my_cell->get_Wall_Nodes());
   	new_force = Fmi + Fii;
	
	return;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell

Coord Cyt_Node::calc_Morse_II() {
	//calc force for II
	Coord Fii; //initialized to zero

	vector<Cyt_Node*>cyts;
	if (my_cell==NULL) {
		cout << "Error: Trying to access NULL Pointer. Aborting!" << endl;
		exit(1);
	}
	my_cell->get_Cyt_Nodes_Vec(cyts);
	Cyt_Node* me = this;
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fii) schedule(static,1)
		for (unsigned int j = 0; j < cyts.size(); j++) {
			//don't calculate yourself
			if (cyts.at(j) != me) {
				//calc morse between this node and node j
				Fii +=  me->morse_Equation(cyts.at(j));
			}
		}
	}

	return Fii;
}

Coord Cyt_Node::calc_Morse_MI(Wall_Node* orig) {
	//calc force for IM
	Coord Fmi;
	vector<Wall_Node*> walls;
	this->get_My_Cell()->get_Wall_Nodes_Vec(walls);
	Cyt_Node* me = this;
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fmi) schedule(static,1)
		for(unsigned int i = 0; i < walls.size(); i++) {
			Fmi +=me-> morse_Equation(walls.at(i));
			//update curr_wall
		}
	}

	return Fmi;
}

Coord Cyt_Node::morse_Equation(Cyt_Node* cyt) {
	
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
    
//	Fii = diff_vect * (1*ALPHA*DELTA*(1-exp(-ALPHA*(diff_len -MORSE_EQ)))*(exp(-ALPHA*(diff_len-MORSE_EQ)))*(1.0/diff_len));
	Fii = diff_vect*((-attract + repel)/diff_len);
	return Fii;
}

Coord Cyt_Node::morse_Equation(Wall_Node* wall) {

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
    	//Fmi = diff_vect * (1*ALPHA_MI*DELTA_MI*(1-exp(-ALPHA_MI*(diff_len -MORSE_EQ_MI)))*(exp(-ALPHA_MI*(diff_len-MORSE_EQ_MI)))*(1.0/diff_len));
  	 Fmi = diff_vect*((-attract + repel)/diff_len);
	return Fmi;
}

Cyt_Node::~Cyt_Node() {
	my_cell = NULL;
}


//======================================================
/** class Wall Node Functions **/

// Constructors-----------------
Wall_Node::Wall_Node(Coord loc, Cell* my_cell) : Node(loc) {
	this->my_cell = my_cell;
	this->closest = NULL;
	this->closest_len = 100;
}

Wall_Node::Wall_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right) : Node(loc)   {
	this->left = left;
    	this->right = right;
	this-> my_cell = my_cell;
	this->closest = NULL;
	this->closest_len = 100;
	update_Angle();
}

Wall_Node::~Wall_Node() {
	my_cell = NULL;	
	left = NULL;
	right = NULL;
	closest = NULL;
	//Wall_Node* wall = NULL;
	/*while ( !closest_vec.empty()) {
		wall = closest_vec.at(closest_vec.size() - 1);
		delete wall;
		closest_vec.pop_back();
	}*/
}

//  Getters and Setters--------------------
void Wall_Node::set_Equi_Angle(double angle) {
	equi_angle = angle;
	return;
}

void Wall_Node::set_Left_Neighbor(Wall_Node* new_Left) {
	this->left = new_Left;
	return;
}

void Wall_Node::set_Right_Neighbor(Wall_Node* new_Right) {
	this->right = new_Right;
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
	cross_Prod = crossProd;
	
	return;
}
void Wall_Node::update_Cell(Cell* new_cell) {
	this->my_cell = new_cell;
}
void Wall_Node::update_Equi_Angle(double new_theta) {
	equi_angle = new_theta;

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
void Wall_Node::set_Closest(Wall_Node*  closest, double closest_len) {
	this->closest = closest;
	this->closest_len = closest_len;
	return;
}
// Calc Force Functions -----------------------
void Wall_Node::calc_Forces(int Ti) {
	// Initialize force sum to zero by default constructor
	Coord sum;
	sum += calc_Morse_SC();
//	cout << "SC success" << endl;	
	cyt_force = sum;

	sum += calc_Morse_DC();
//	cout << "DC" << calc_Morse_DC() << endl;
	sum += calc_Linear();
//	cout << "linear" << endl;
	sum += calc_Bending();
//	cout << "bending" << endl;

	// Update new_force variable for location updating
	new_force = sum;

	return;
}

//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC() {
	vector<Cyt_Node*> cyt_nodes;
	my_cell->get_Cyt_Nodes_Vec(cyt_nodes);
	Wall_Node* curr_wall = this;
	Coord Fmi;
	#pragma omp parallel
	{	
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fmi) schedule(static,1)
		for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
			Fmi += curr_wall->morse_Equation(cyt_nodes.at(i));
		}
	}
//	cout << "	morse_sc: " << Fmi << endl;	
	return Fmi;
}


//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_DC() {
	Coord Fdc;
	vector<Cell*> cells;
	my_cell->get_Neighbor_Cells(cells);	
	vector<Wall_Node*> walls;
	//cout << "getting neighbors" << endl;
	Wall_Node* curr = NULL;
	Wall_Node* orig = NULL;
	#pragma omp parallel 
	{
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:Fdc) schedule(static,1) 
		for (unsigned int i = 0; i < cells.size(); i++) {
			//for(unsigned int j = 0; j <300/* cells.at(i)->get_wall_count()*/; j ++) {
				//cells.at(i)->get_Wall_Nodes_Vec(walls);
				//Fdc += morse_Equation(walls.at(j));
			//}
		//}
			
			Fdc += neighbor_nodes(cells.at(i));
	}
	}

	//cout << "made it out of loop" << endl;
	//cout << closest << endl;
	if(this->closest != NULL){
	//	cout << "closest not null" << endl;
		closest->get_Location();
	//	cout << "got location" << endl;
		Fdc += linear_Equation_ADH(this->closest);
	//	cout << "computed adh successfully" << endl;
	}
	//cout << " morse_DC: " << Fdc << endl;
	return Fdc;
}

Coord Wall_Node::neighbor_nodes(Cell* neighbor) {
	Coord sum;
	vector<Wall_Node*> walls;
	neighbor->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:sum) schedule(static,1) 
		for(unsigned int j =0; j< walls.size(); j++) {
			//cout << "getting wall nodes" << endl;
			sum += morse_Equation(walls.at(j));
			//cout << "morse" << endl;
		}
	}

	return sum;
}			

//bending force of node
Coord Wall_Node::calc_Bending() {
	Coord F_bend;

	F_bend += bending_Equation_Center();
	F_bend += bending_Equation_Left();
	F_bend += bending_Equation_Right();
	
	if (cross_Prod < 0.0) {
		F_bend = F_bend*(-1);
	}	
//	cout << "	bending: " << F_bend << endl;
	return F_bend;
}

//spring force of neighboring springs
Coord Wall_Node::calc_Linear() {
	Coord F_lin;

//	cout << "calc left" << endl;
	F_lin += linear_Equation(left);

//	cout << "calc right" << endl;
	F_lin += linear_Equation(right);
	
//	cout << "	linear: " << F_lin << endl;
	return F_lin;
}


//===========================================================
// Mathematical force calculations


Coord Wall_Node::morse_Equation(Cyt_Node* cyt) {
	if (cyt == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
	}
	
	//use Membr - int variables
	Coord Fmi;
	Coord diff_vect = cyt->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
	double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);

	//Fmi = diff_vect * (1*ALPHA_MI*DELTA_MI*(1-exp(-ALPHA_MI*(diff_len -MORSE_EQ_MI)))*(exp(-ALPHA_MI*(diff_len-MORSE_EQ_MI)))*(1.0/diff_len));
  	 Fmi = diff_vect*((-attract + repel)/diff_len);
	
	//cout << Fmi << endl;
	return Fmi;
}

Coord Wall_Node::morse_Equation(Wall_Node* wall) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
	}

	//use Mem-Mem variables
   	Coord Fmmd;
   	Coord diff_vect = wall->get_Location() - my_loc;
    	double diff_len = diff_vect.length();
    	double attract = (U_MM/xsi_MM)*exp(diff_len*(-1)/xsi_MM);
    	double repel = (W_MM/gamma_MM)*exp(diff_len*(-1)/gamma_MM);
    
   	 //Fmmd =  diff_vect*(0);//(-attract + repel)/diff_len);
	 Fmmd = diff_vect*((-attract + repel)/diff_len);
	
	//	cout << Fmmd << endl;
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

Coord Wall_Node::linear_Equation(Wall_Node* wall) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
	}
	
	//use spring constant variables
	Coord F_lin;
	Coord k_linear = this->get_My_Cell()->get_K_LINEAR();
//	Coord K_LINEAR;
//	if(this->get_My_Cell()->get_Layer() == 1) {
//		K_LINEAR = K_LINEAR_WIDE;
//	}
//	else {
//		K_LINEAR = K_LINEAR_LONG;
//	}
	Coord diff_vect = wall->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	Coord scaled_k = k_linear*(diff_len - MembrEquLen);
	F_lin = (diff_vect/diff_len).distribute(scaled_k);

	return F_lin;	
}

Coord Wall_Node::linear_Equation_ADH(Wall_Node*& wall) {
	//cout << "wall node is: " << wall << endl;
	//Wall_Node* wall = NULL;
//	for(unsigned int i = 0;i < closest_nodes.size();i++){
//		wall = closest_nodes.at(i);
	if (wall == NULL) {
		cout << "Problems for days" << endl;
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
//	}

	return F_lin;
}

//==========================================================
//Adhesion functions

Wall_Node* Wall_Node::find_Closest_Node(vector<Cell*>& neighbors) {
	Wall_Node* curr = NULL;
	Wall_Node* orig = NULL;
	Wall_Node* next = NULL;
	Cell* curr_cell = NULL;
	Wall_Node* closest = NULL;
	double curr_dist = 0;
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


void Wall_Node::make_Connection(Wall_Node* curr_Closest) {
	double curr_dist = 0;
	if(curr_Closest != NULL) {
		curr_dist = (this->get_Location() - curr_Closest->get_Location()).length();
		
	//	if(curr_Closest->get_Closest() ==NULL) {
	//		if(this->closest != NULL) {
	//			this->closest->set_Closest(NULL,100);
	//		}
			this->closest = curr_Closest;
			this->closest_len = curr_dist;
	//		curr_Closest->set_Closest(this,curr_dist);
	//	}	
	//	else {
	//		if(curr_Closest->get_Closest()== this) {
	//			this->closest = curr_Closest;
	//			this->closest_len = curr_dist;
	//		}
	//		else if(curr_dist < curr_Closest->get_Closest_Len()) { 
	//			if(this->closest != NULL) {
	//				this->closest->set_Closest(NULL,100);
	//			}
	//			
	//			this->closest = curr_Closest;
	//			this->closest_len = curr_dist;
	//			curr_Closest->get_Closest()->set_Closest(NULL,100);
	//			curr_Closest->set_Closest(this,curr_dist);
	//		
	//		}
	//	}
	
	/*	else if((this->get_Closest() == NULL) && (curr_Closest->get_Closest() != NULL)) {
			if(curr_dist < curr_Closest->get_Closest_Len()) {
				this->closest = curr_Closest;
				this->closest_len = curr_dist;
				curr_Closest->get_Closest()->set_Closest(NULL,100);
				curr_Closest->set_Closest(this,curr_dist);
			}
		}
		else if((this->get_Closest()!= NULL) && (curr_Closest->get_Closest() != NULL)) {
			if(this->closest == curr_Closest->get_Closest()) {
				//do nothin
			}
			else if(curr_dist < this->closest_len) {
				if(curr_dist < curr_Closest->get_Closest_Len()) {
					this->closest->set_Closest(NULL,100);
					this->closest = curr_Closest;
					this->closest_len = curr_dist;
					curr_Closest->get_Closest()->set_Closest(NULL,100);
					curr_Closest->set_Closest(this,curr_dist);
				}
		
			}
			
		}*/
	}
	return;
}

//==========================================================
// End of node.cpp

