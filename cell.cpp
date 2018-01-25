//cell.cpp
//===================
// Forward Declarations

//===================
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "rand.h"
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors
Cell::Cell(Tissue* tissue) {
	my_tissue = tissue;
	//just divided so reset life length
	life_length = 0;
	//cyt nodes and neighbor cells assigned
	//at time of division
	//rank assigned at time of division
	damping = 1;
	//num_cyt_nodes assigned at time of division
	//starts at 0
	//k_linear assigned at time of division
	num_cyt_nodes = 0;
	//layer inherited from parent cell
	//cell center computed at time of division
	//num_wall nodes computed at time of division
	num_wall_nodes = 0;
	//left_Corner assigned at time of division
	//WUS computed based on new cell center
	wuschel = 0;
	//cytokinin computed based on new cell center
//	counter = 0;
//	curr_area = this->calc_Area();
	cytokinin = 0;
	total_signal = 0;
	Cell_Progress = 0;
	Cell_Progress_add_node = 0;
	Cell_Progress_div = 0;
	//	top = closest_node_top();
//	bottom = closest_node_bottom();
	counter_left = 0;
	counter_right = 0;
}


Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer)    {

	this->rank = rank;
//	cout << "rank " << rank << endl; 
	this->my_tissue = tiss;
	num_cyt_nodes = 0;
	this->layer = layer;
//	int init_radius = radius;
	this->cell_center = center;
//	cout << "length: " << cell_center.get_Y() << endl;
	num_wall_nodes = 0;
	life_length = 0;
	Cell_Progress = 0;
	Cell_Progress_add_node = 0;
	this->calc_WUS();
	this->calc_CYT();
	this->calc_Total_Signal();
//  make functions for these*****
	
	double K_LINEAR = -3.3673*(this->cytokinin) + 5.7335*(this->wuschel) + 269.4673;
	double K_LINEAR_X;
	double K_LINEAR_Y;
	
	if(this->layer == 1) {
		K_LINEAR_Y = 600;
		K_LINEAR_X = 150;
	}	
	else {
		K_LINEAR_X = 150;
		K_LINEAR_Y = 600;
	}	

	this->K_LINEAR = Coord(K_LINEAR_X, K_LINEAR_Y);
	
	//assemble the membrane
	int num_Init_Wall_Nodes = Init_Wall_Nodes;
	double angle_increment = (2*pi)/num_Init_Wall_Nodes;
	
	//make all wall nodes
	double curr_X;
	double curr_Y;
	Coord location = this->cell_center;;
	Wall_Node* currW;
	Wall_Node* prevW;
	Wall_Node* orig;
	double curr_theta = 0;
	curr_X = cell_center.get_X() + radius*cos(curr_theta);
	curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
	location = Coord(curr_X,curr_Y);
	//make the first node
	prevW = new Wall_Node(location,this);
	num_wall_nodes++;
	orig = prevW;
	//this will be the "starter" node
	this->left_Corner = orig;

	//make successive nodes
	for(int i = 0; i<num_Init_Wall_Nodes-1; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = cell_center.get_X() + radius*cos(curr_theta);
		curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		currW = new Wall_Node(location, this);
		num_wall_nodes++;
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}
	
	//connect last node to starter node
	currW->set_Left_Neighbor(orig);
	orig->set_Right_Neighbor(currW);

	//insert cytoplasm nodes
	int	num_init_cyt_nodes = Init_Num_Cyt_Nodes;
	double scal_x_offset = 0.8;
	//Coord location;
	Cyt_Node* cyt;
	double x;
	double y;
	for (int i = 0; i < num_init_cyt_nodes; i++) {
		// USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		double rand_radius = (static_cast<double>(rand()) / RAND_MAX)*scal_x_offset*radius;
		double rand_angle = (static_cast<double>(rand()) / RAND_MAX)*2*pi;
		x = cell_center.get_X()+ rand_radius*cos(rand_angle);
		y = cell_center.get_Y()+ rand_radius*sin(rand_angle);
		location = Coord(x,y);

		cyt = new Cyt_Node(location,this);
		cyt_nodes.push_back(cyt);
		num_cyt_nodes++;
	}
	counter = 0;
	//update equilibrium angle
	update_Wall_Equi_Angles();
	//update wall angles
	update_Wall_Angles();
	//set damping for cells that act as anchor points
	if(layer == 6) {
		this->damping = .3;
	}
	else {
		this->damping = 1;
	}
	top = closest_node_top();
	bottom = closest_node_bottom();
	counter_left = 0;
	counter_right = 0;
	Cell_Progress_div = 0;
}


// Destructor
Cell::~Cell() {

	// Delete Cyt Nodes
	Cyt_Node* cyt = NULL;
	while ( !cyt_nodes.empty()) {
		cyt = cyt_nodes.at(cyt_nodes.size() - 1);
		delete cyt;
		cyt_nodes.pop_back();
	}
	// Delete Wall Nodes
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* last = curr->get_Right_Neighbor();
	last->set_Left_Neighbor(NULL);

	while (curr != NULL) {
		next = curr->get_Left_Neighbor();
		delete curr;
		curr = next;
	}

	my_tissue = NULL;
}

//=============================================================
//========================================
// Getters and Setters
//========================================
//=============================================================
int Cell::get_Node_Count() {
	return num_wall_nodes + num_cyt_nodes;
}
void Cell::get_Cyt_Nodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}
void Cell::get_Neighbor_Cells(vector<Cell*>& cells) {
	cells = neigh_cells;
	return;
}
void Cell::set_K_LINEAR(double& x, double& y) {
	K_LINEAR = Coord(x,y);
	return;
}
void Cell::set_Damping(double& new_damping) {
	this->damping = new_damping;
	return;
}
void Cell::set_Rank(const int id) {
	this->rank = id;
	return;
}
void Cell::set_Layer(int layer) {
	this->layer = layer;
	return;
}
void Cell::reset_Cell_Progress(){
	this->Cell_Progress = 0;
	return;
}
void Cell::set_Left_Corner(Wall_Node*& new_left_corner) {
	this->left_Corner = new_left_corner;
	return;
}

void Cell::set_Wall_Count(int& number_nodes) {
	this->num_wall_nodes = number_nodes;
	return;
}

void Cell::update_Cell_Progress_add_node(int& time) {
	this->Cell_Progress_add_node = time;
	return;
}
void Cell::update_Life_Length() {
	life_length++;
//	cout << life_length << endl;
	return;
}
void Cell::set_div_time(int& Ti) {
	this->Cell_Progress_div = Ti;
}
void Cell::calc_WUS() {
	this->wuschel = -0.229814360326137*pow(cell_center.length(),2) + 9.2927997842538854*cell_center.length() + -26.5308519096016;
	if(wuschel < 0) {
		wuschel = 0;
	}
	//cout << "wuschel: " << wuschel << endl;
	return;
}
void Cell::calc_CYT() {
	this->cytokinin = -0.0000008*pow(cell_center.length(),4) + 0.0002445*pow(cell_center.length(),3) + -0.0213*pow(cell_center.length(),2) + .4733*(cell_center.length()) + 17.555;
	if(cytokinin < 0) {
		cytokinin = 0;
	}
	//	cout << "cytokinin: " << cytokinin << endl;
	return;
}
void Cell::calc_Total_Signal() {
	this->total_signal = wuschel + cytokinin;
	if(total_signal < 0) {
		total_signal = 0;
	}
	return;
}

//=============================================================
//=========================================
// Keep Track of neighbor cells and Adhesion springs
//=========================================
//=============================================================
void Cell::update_Neighbor_Cells() {
	//clear prev vector of neigh cells
	neigh_cells.clear();
	//grab all cells from tissue
	vector<Cell*> all_Cells;
	my_tissue->get_Cells(all_Cells);

	
	// Empty variables for holding info about other cells
	Cell* curr = NULL;
	Coord curr_Cent;
	Coord distance;
	
	double prelim_threshold = 18;
	//double sec_threshold = 1;

	// iterate through all cells
	for (unsigned int i = 0; i < all_Cells.size(); i++) {
		curr = all_Cells.at(i);
		if (curr != this) {
			curr_Cent = curr->get_Cell_Center();
			// Check if cell centers are close enough together
			distance = this->cell_center - curr_Cent;
			//cout << "Distance = " << distance << endl;
			if ( distance.length() < prelim_threshold ) {
				neigh_cells.push_back(curr);
				//cout << rank << "has neighbor" << curr->get_Rank() << endl;
			}
			
		}
		//else you're pointing at yourself and shouldnt do anything
	}	
	
//	cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}

void Cell::update_adhesion_springs() {
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* next = NULL;
	do {
		next = curr->get_Left_Neighbor();
		curr->set_Closest(NULL, 100);
		curr = next;
	} while(next != orig);
	
	vector<Cell*>neighbors;
	this->get_Neighbor_Cells(neighbors);
	Wall_Node* curr_Node = NULL;
	Wall_Node* next_Node = NULL;
	Wall_Node* curr_Closest = NULL;
	double curr_len = 0;
//	for(int i = 0; i<num_wall_nodes;i++) {
	curr_Node = left_Corner;
	do {
		next_Node = curr_Node->get_Left_Neighbor();
		curr_Closest = curr_Node->find_Closest_Node(neighbors);
		curr_Node->make_Connection(curr_Closest);
		curr_Node = next_Node;
	} while(next_Node != left_Corner);
	//}
	return;
}

//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
void Cell::calc_New_Forces(int Ti) {
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces();
	}

	//calc forces on wall nodes
	Wall_Node* curr = left_Corner; 
	Wall_Node* orig = curr;
	int counter = 0;
	do {
		counter++;
//		cout << "Wall node number: " << counter << endl;
		curr->calc_Forces(Ti);
	//	cout << "Forces calculated" << endl;
		curr = curr->get_Left_Neighbor();
	
	} while (curr != orig);
	//cout << "out of forces function" << endl;
	return;
}

void Cell::update_Node_Locations() {
	//update cyt nodes
	double new_damping = 0;
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		new_damping = cyt_nodes.at(i)->get_My_Cell()->get_Damping();
		cyt_nodes.at(i)->update_Location(new_damping);

	}

	//update wall nodes
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = left_Corner;
	
//	if(this->top == NULL) {
//		cout << "NULL top" << endl;
//	}
//	double x_coord = this->top->get_Location().get_X();
//	double direction;
//	Wall_Node* neighbor = NULL;
	do {
		new_damping = curr->get_My_Cell()->get_Damping();
	//	cout << "check if needs stationary" << endl;
	//	direction = curr->get_Location().get_X() - x_coord;
	//	if(direction > 0) {
	//		neighbor = curr->get_Left_Neighbor();
	//	}
	//	else if(direction < 0) {
	//		neighbor = curr->get_Right_Neighbor();
	//	}

	//	if((curr->get_Location().get_Y() > neighbor->get_Location().get_Y() - .02) && (curr->get_Location().get_Y() < neighbor->get_Location().get_Y() + .02)){
//		if(curr->get_isStationary()) {
//			//do nothing
//		}
	//	else {
//		if((curr->get_Angle() - pi < .00001) && (curr->get_Angle() -pi > -.00001)) {
//			curr->set_isStationary();
//		}
	//	}
	//	cout << "update locaation" << endl;
		curr->update_Location(new_damping);
		curr = curr->get_Left_Neighbor();
	
	} while(curr != orig);

	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	update_Wall_Angles();
//	cout << "done" << endl;
	return;
}

void Cell::update_Wall_Angles() {
//	cout << "wall angles" << endl;
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = curr;
	do {
		curr->update_Angle();
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);	

	return;
}

void Cell::update_Wall_Equi_Angles() {
//	cout << "equi angles" << endl;
	double new_equi_angle = (num_wall_nodes-2)*pi/num_wall_nodes;
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = left_Corner;
	
	do {
		curr->update_Equi_Angle(new_equi_angle);
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);	
	
	return;
}

void Cell::update_Cell_Center() {
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = curr;
	Coord total_location = Coord();

	do {
		total_location += curr->get_Location();
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);

	cell_center = total_location*(1.0/num_wall_nodes);
	
	return;
}
//=====================================================================
//==========================================
// Growth of Cell
//==========================================
//=====================================================================
void Cell::update_Cell_Progress(int Ti) {
	//update life length of the current cell
	this->update_Life_Length();
	if(Ti%217 == 0) {
		if(Cell_Progress_div > 0) {
			if(Ti-Cell_Progress_div > 1000) {
				this->wall_Node_Check();
			}
		}
		else {
			this->wall_Node_Check();
		}
	}
	//variables needed if division occurs
	Cell* new_Cell= NULL;
	vector<Cell*> cells;
	vector<Cell*> neighbor_cells;
	this->my_tissue->get_Cells(cells);
	int number_cells = cells.size();
	//variables for determining growth rate
	double sigma = (((double) rand()/(double) RAND_MAX))+.004;
//	double rate = (.004 + sigma)*exp(GROWTH_RATE*life_length);
	double rate = 1;	
	this->curr_area = this->calc_Area();
	//update cell progress
	Cell_Progress = Cell_Progress + rate*dt;
//	cout << "Rank: " << this->rank << "and Progress: " << Cell_Progress << " and life length: " << life_length << endl;
	//division check
//	cout << "Sigma"<< sigma << endl;
	if(Ti==5000) { //(this->Cell_Progress >= 1) && (curr_area >= AREA_THRESH)) {
		cout << "Cell Prog" << Cell_Progress << endl;
		new_Cell = this->divide();
		Cell_Progress_div = Ti;
		new_Cell->set_div_time(Ti);
		cout << "division success" << endl;
		if(new_Cell == NULL) {
			cout << "womp womp" << endl;
		}
		new_Cell->set_Rank(number_cells);
		cout << "set rank" << endl;
		cout << "Parent rank: " << this->rank << endl;
		cout << "sister rank: " << new_Cell->get_Rank() << endl;
		my_tissue->update_Num_Cells(new_Cell);
		cout << "added cell" << endl;
		this->update_Neighbor_Cells();
		this->update_adhesion_springs();
		cout << "updated neighbor" << endl;
		this->life_length = 0;
		this->Cell_Progress = 0;
		this->Cell_Progress_add_node = 0;
		this-> counter = 0;
		this->get_Neighbor_Cells(neighbor_cells);
		for(unsigned int i=0; i < neighbor_cells.size(); i++) {
			neighbor_cells.at(i)->update_adhesion_springs();
		}
		cout << "updated adhesion" << endl;
		new_Cell->update_Neighbor_Cells();
		cout << "new cell neighbor update" << endl;
		new_Cell->update_adhesion_springs();
		cout << "new cell adhesion update" << endl;
		new_Cell->get_Neighbor_Cells(neighbor_cells);
		for(unsigned int i = 0; i< neighbor_cells.size(); i++) {
			neighbor_cells.at(i)->update_adhesion_springs();
		}
	}
	//if no division check if internal node should be added
	else {
	
//		if(Cell_Progress - Cell_Progress_add_node > 2.5) { 
    	//	this->add_Cyt_Node();
		
		//	this->counter = counter +1;
			//cout << "added cyt node" << this->counter << endl;
		//	cout << curr_area << " is area" << endl;
		//	cout << "time: " << Ti << endl;
//			Cell_Progress_add_node = Cell_Progress;
		}
	
	//cout << "Cell Prog: " << Cell_Progress_add_node << endl;
	return;
}

double Cell::calc_Area() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	Coord a_i;
	Coord a_j;
	double area = 0;
	double curr_area = 0;
	do {
		next = curr->get_Left_Neighbor();
		a_i = curr->get_Location() - cell_center;
		a_j = next->get_Location() - cell_center;
		curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
		area += curr_area;
		curr = next;
	} while(next != orig);
	//cout << "Area: " << area << endl;
	return area;
}
void Cell::wall_Node_Check() {
	//cout << "adding a wall node" << endl;
		add_Wall_Node();
	return;
}
void Cell::add_Wall_Node() {
	//find node to the right of largest spring
	Wall_Node* right = NULL;
	find_Largest_Length(right);
	Wall_Node* left = NULL;
	Coord location;
	Wall_Node* added_node = NULL;
	if(right != NULL) {
		cout << "wasnt null" << endl;
		left = right->get_Left_Neighbor();
		location  = (right->get_Location() + left->get_Location())*0.5;
		added_node = new Wall_Node(location, this, left, right);
		//cout << "made new node" << endl;
		right->set_Left_Neighbor(added_node);
		left->set_Right_Neighbor(added_node);
		num_wall_nodes++;
		update_Wall_Equi_Angles();
		update_Wall_Angles();
	}
	else {
		//cout << "null" << endl;
	}
	return;
}
void Cell::find_Largest_Length(Wall_Node*& right) {
	Wall_Node* curr = left_Corner;
	Wall_Node* biggest = NULL;
	Wall_Node* orig = curr;
	Coord left_Neighb_loc;
	Coord curr_Loc;
	Coord diff_vect;
	double max_len = 0;
	double len;
	int big_gaps = 0;
	//loop through all possible Cell Wall 'links' to find biggest
	do {
	//	cout << "finding current lengths and comparing" << endl;
		left_Neighb_loc = curr->get_Left_Neighbor()->get_Location();
		curr_Loc = curr->get_Location();
		diff_vect = left_Neighb_loc - curr_Loc;
		len = diff_vect.length();
		if (len > MEMBR_THRESH_LENGTH) {
			big_gaps++;
			if(len > max_len) {
					max_len = len;
					right = curr;
			}
		}
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);

	//cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;
	return;
}
void Cell::add_Cyt_Node() {

	Cyt_Node* cyt = new Cyt_Node(cell_center, this);
	cyt_nodes.push_back(cyt);

	num_cyt_nodes++;
	return;
}

//=====================================================
//==================================
//Functions for Division
//=================================
//=====================================================
double Cell::find_radius() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	double distance = (cell_center-curr->get_Location()).length();
	double min_distance = distance;

	do {
		next = curr->get_Left_Neighbor();
		distance = (cell_center - curr->get_Location()).length();
		if(distance < min_distance) {
			min_distance = distance;
		}
		curr = next;
	} while (next != orig);
	
	return min_distance;
}

void Cell::add_Cyt_Node_Div(double radius) {
	//USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
	// ---distributes more evenly throughout start cell
	double offset = .8;
	Coord location;
	Cyt_Node* cyt;
	double x;
	double y;
	double rand_radius = unifRand()*offset*radius; 
	double rand_angle = unifRand()*2*pi;
	x = cell_center.get_X()+ rand_radius*cos(rand_angle);
	y = cell_center.get_Y()+ rand_radius*sin(rand_angle);
	location = Coord(x,y);
	cout << location << endl;
	cyt = new Cyt_Node(location,this);
	cyt_nodes.push_back(cyt);
	num_cyt_nodes++;
	return;
}

double Cell::compute_Stress_Tensor_X() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	double stress = 0;
	double area = this->calc_Area();
	double force_curr;
	double force_next;
	double force_curr_next;
	double force_next_curr;
	do {
		next = curr->get_Left_Neighbor();
		force_curr = (curr->get_Force().get_X())*(curr->get_Location().get_X());
		force_next = (next->get_Force().get_X())*(next->get_Location().get_X());
		force_curr_next = (curr->get_Force().get_X())*(next->get_Location().get_X());
		force_next_curr = (next->get_Force().get_X())*(curr->get_Location().get_X());
		stress = stress + ((force_curr + force_next)/(3*area)) + ((force_curr_next + force_next_curr)/(6*area));
		curr = next;
	} while (next != orig);
	
	return stress;
}

double Cell::compute_Stress_Tensor_Y() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	double stress = 0;
	double area = this->calc_Area();
	double force_curr;
	double force_next;
	double force_curr_next;
	double force_next_curr;
	do {
		next = curr->get_Left_Neighbor();
		force_curr = (curr->get_Force().get_Y())*(curr->get_Location().get_Y());
		force_next = (next->get_Force().get_Y())*(next->get_Location().get_Y());
		force_curr_next = (curr->get_Force().get_Y())*(next->get_Location().get_Y());
		force_next_curr = (next->get_Force().get_Y())*(curr->get_Location().get_Y());
		stress = stress + ((force_curr + force_next)/(3*area)) + ((force_curr_next + force_next_curr)/(6*area));
		curr = next;
	} while (next != orig);
	
	return stress;
}

double Cell::compute_Stress_Tensor_XY() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	double stress = 0;
	double area = this->calc_Area();
	double term1;
	double term2;
	double term3;
	double term4;
	double term5;
	double term6;
	double term7;
	double term8;
	do {
		next = curr->get_Left_Neighbor();
		term1 = (curr->get_Force().get_X())*(curr->get_Location().get_Y());
		term2 = (curr->get_Force().get_Y())*(curr->get_Location().get_X());
		term3 = (next->get_Force().get_X())*(next->get_Location().get_Y());
		term4 = (next->get_Force().get_Y())*(next->get_Location().get_X());
		term5 = (curr->get_Force().get_X())*(next->get_Location().get_Y());
		term6 = (curr->get_Force().get_Y())*(next->get_Location().get_X());
		term7 = (next->get_Force().get_X())*(curr->get_Location().get_Y());
		term8 = (next->get_Force().get_Y())*(curr->get_Location().get_X());
		stress = stress + ((term1+ term2 + term3 + term4)/(6*area)) + ((term5 + term6 + term7 + term8)/(12*area));
		curr = next;
	} while (next != orig);
	
	return stress;
}
void Cell::stress_Tensor_Eigenvalues(double& a, double& b, double & c, double& d, vector<double>& eigen_Max) {
//	cout << "Made it" << endl;
	double trace = a + d;
//	cout << "Trace" << trace << endl;
	double determinant = a*d - b*c;
//	cout << "Determinant" << determinant <<  endl;
//	cout << sqrt(pow(trace,2)*.25 - determinant) << endl;
	double lambda_1 = trace*.5 + sqrt(pow(trace,2)*.25 - determinant);
//	cout << "Lambda 1" << lambda_1 << endl;
	double lambda_2 = trace*.5 - sqrt(pow(trace,2)*.25 - determinant);
//	cout << "Lambda 2" << lambda_2 << endl;
	if(lambda_1 > lambda_2) {
		if((c!=0) || (b!=0)) {
			if (c !=0) {
				//cout << "first one" << endl;
				eigen_Max.push_back(lambda_1 - d);
				eigen_Max.push_back(c);
				//cout << "pushed" << endl;
			}	
			else if(b !=0) {
				//cout << "first one" << endl;
				eigen_Max.push_back(b);
				eigen_Max.push_back(lambda_1 -a);
				//cout << "pushed" << endl;
			}
		}
		else {
			//cout << "First" << endl;
			eigen_Max.push_back(1);
			eigen_Max.push_back(0);
			//cout << "pushed" << endl;
		}
	}
	else {
		if((c!=0) || (b!=0)) {
			if (c !=0) {
				//cout << "first" << endl;
				eigen_Max.push_back(lambda_2 - d);
				eigen_Max.push_back(c);
				//cout << "pushed" << endl;
			}	
			else if(b !=0) {
				//cout << "First" << endl;
				eigen_Max.push_back(b);
				eigen_Max.push_back(lambda_2 -a);
				//cout << "pushed" << endl;
			}
		}
		else {
			//cout << "First" << endl;
			eigen_Max.push_back(1);
			eigen_Max.push_back(0);
			//cout << "Pushed" << endl;
		}
	}
	//cout << "Eigen values: 1 " << eigen_Max.at(0) << " 2 " << eigen_Max.at(1) << endl; 
	return;
}
//===========================================================
//==================================
// Output Functions
//==================================
//===========================================================

void Cell::print_Data_Output(ofstream& ofs) {
	
	ofs << "This is where data output goes" << endl;

	return;
}

int Cell::update_VTK_Indices(int& id) {
//	cout << "ID before: " << id << endl;
	int rel_cnt = 0;

	Wall_Node* curr_wall = left_Corner;
	do { 
		curr_wall->update_VTK_Id(id);
		id++;
		if(curr_wall->get_Closest()!= NULL) {
			rel_cnt++;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);
	
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			cyt_nodes.at(i)->update_VTK_Id(id);
			id++;
	}
//	cout << "ID after: " << id << endl;
	return rel_cnt;
}
void Cell::print_VTK_Adh(ofstream& ofs) {

	int my_id, nei_id;
	Wall_Node* neighbor = NULL;
	Wall_Node* curr_wall = left_Corner;

	do {
		neighbor = curr_wall->get_Closest();
		if(neighbor != NULL) {
			my_id = curr_wall->get_VTK_Id();
			nei_id = neighbor->get_VTK_Id();
			ofs << 2 << ' ' << my_id << ' ' << nei_id << endl;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while(curr_wall != left_Corner);
	return;
}
void Cell::print_VTK_Points(ofstream& ofs, int& count) {

	Wall_Node* curr_wall = left_Corner;
	Wall_Node* orig = curr_wall;
	//	cout << "knows left corner" << endl;
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();
		count++;
		//cout << "did it  " << count << endl;
	} while (curr_wall != orig);
	
	//cout << "walls worked" << endl;
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		count++;
	}
	//cout << "points worked" << endl;
	return;
}

void Cell::print_VTK_Scalars_Force(ofstream& ofs) {

	Wall_Node* curr_wall = left_Corner;
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.length() << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.length() << endl;
	}

	return;
}

void Cell::print_VTK_Scalars_WUS(ofstream& ofs) {

	//double concentration = 0;
	Wall_Node* curr_wall = left_Corner;
	do {
		double concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		double concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;
	}
	return;
}

void Cell::print_VTK_Scalars_CYT(ofstream& ofs) {

//	double concentration = 0;

	Wall_Node* curr_wall = left_Corner;
	do {
		double concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		double concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;
	}
	return;
}

void Cell::print_VTK_Scalars_Total(ofstream& ofs) {

//	double concentration = 0;

	Wall_Node* curr_wall = left_Corner;
	do {
		double concentration = curr_wall->get_My_Cell()->get_total_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		double concentration = cyt_nodes.at(i)->get_My_Cell()->get_total_concentration();
		ofs << concentration << endl;
	}
	return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

	Wall_Node* curr_wall = left_Corner;
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while(curr_wall != left_Corner);

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;
	}

	return;
}

//===========================================================
//=======================================================
//Functions for finding top/right/left/bottom-most nodes
//=======================================================
//===========================================================
Wall_Node* Cell::closest_node_top() {
	//finds node with highest y value that is
	//in a certain window around x value of center
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double smallest_dist = 500;
//	cout << "smallest dist top" << endl;
	Coord curr_coord;
	Wall_Node* closest = NULL;

	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(x_coord -curr_x,2));
		if(curr_dist < smallest_dist) {
		//	cout << "passed smallest check" << endl;
			if(curr_y > y_coord) {
		//		cout << "passed y coord check" << endl;	
				//cout << "passed smallest check" << endl;
					smallest_dist = curr_dist;
					//cout << "curr " << curr << endl;
					closest = curr;
				}
		}
		curr = next;
	} while (next != orig);
	//cout << "closest" << closest << endl;
	//up = closest;
	//cout << "up in function" << up << endl;
	return closest;
}

Wall_Node* Cell::closest_node_bottom() {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double smallest_dist = 500;
	Coord curr_coord;
	Wall_Node* closest = NULL;

	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(x_coord - curr_x, 2));
		if(curr_dist < smallest_dist) {
			if(curr_y < y_coord) {
					smallest_dist = curr_dist;
					closest = curr;
			}
		}
		curr = next;
	} while (next != orig);
//	cout << "bottom" << endl;
	//down = closest;
	return closest;
}

/*void Cell::closest_node_left(Wall_Node*& left) {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double smallest_dist = 500;
	Coord curr_coord;
	Wall_Node* closest = NULL;

	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(y_coord - curr_y,2));
		if(curr_dist < smallest_dist) {
			//cout << "passed window check" << endl;
			if(curr_x < x_coord){ 
					smallest_dist = curr_dist;
					closest = curr;
			}
		}
		curr = next;
	} while (next != orig);
//	cout << "left" << endl;
	left = closest;
	return;
}

void Cell::closest_node_right(Wall_Node*& right) {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double smallest_dist = 500;
	Coord curr_coord;
	Wall_Node* closest = NULL;
	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(y_coord-curr_y,2));
		if(curr_dist < smallest_dist) {
			if(curr_x > x_coord) {
					smallest_dist = curr_dist;
					closest = curr;
			}
		}
		curr = next;
	} while (next != orig);
//	cout << "right" << endl;
	right = closest;
	return;
}*/






