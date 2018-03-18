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
#include "rand.h"
#include "cell.h"
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors
Cell::Cell(Tissue* tissue) {
	my_tissue = tissue;
	//rank assigned in division function
	//layer inherited from parent	
	//damping calculated in div function
	//just divided so reset life length
	life_length = 0;
	//growth_rate = unifRand(LOW_GROWTH,HIGH_GROWTH);
	//cyt nodes created in division function
	//start at zero
	num_cyt_nodes = 0;
	//wall nodes renumbered in division function
	num_wall_nodes = 0;
	//center calculate in division function
	Cell_Progress = 0;
	Cell_Progress_add_node = 0;
	Cell_Progress_div = 0;
	is_deleted = false;
	//will calculate signals in div function
	wuschel = 0;
	cytokinin = 0;
	total_signal = 0;	
	//spring constants calculated in div function
	//neighbors, adhesion, wall angles set in division function
}


Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer)    {
	this->my_tissue = tiss;
	this->rank = rank;
	this->layer = layer;
	//set damping for cells that act as anchor points
	if(layer == 6) {
		this->damping = .3;
	}
	else {
		this->damping = 1;
	}
	life_length = 0;
	//cyt nodes initialized further down
	num_cyt_nodes = 0;
	//wall nodes initialized further down
	num_wall_nodes = 0;
	this->cell_center = center;
	Cell_Progress = unifRandInt(0,10);
	Cell_Progress_add_node = 0;
	Cell_Progress_div = 0;
	is_deleted = false;
	this->calc_WUS();
	this->calc_CYT();
	this->calc_Total_Signal();
	this->set_growth_rate();
	//double K_LINEAR = -3.3673*(this->cytokinin) + 5.7335*(this->wuschel) + 269.4673;
	double K_LINEAR_X;
	double K_LINEAR_Y;
	
	if((this->layer == 1)||(this->layer == 2)) {
		K_LINEAR_Y = 650;
		K_LINEAR_X = 150;
	}	
	else {
		K_LINEAR_X = 650;
		K_LINEAR_Y = 150;
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
	wall_nodes.push_back(prevW);
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
		wall_nodes.push_back(currW);
		num_wall_nodes++;
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}
	
	//connect last node to starter node
	currW->set_Left_Neighbor(orig);
	orig->set_Right_Neighbor(currW);

	//insert cytoplasm nodes
	int num_init_cyt_nodes = Cell_Progress;
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
	//update equilibrium angle
	update_Wall_Equi_Angles();
	//update wall angles
	update_Wall_Angles();	
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
void Cell::set_Rank(const int id) {
	this->rank = id;
	return;
}
void Cell::set_Layer(int layer) {
	this->layer = layer;
	return;
}
void Cell::set_Damping(double& new_damping) {
	this->damping = new_damping;
	return;
}
void Cell::update_Life_Length() {
	life_length++;
	return;
}
int Cell::get_Node_Count() {
	return num_wall_nodes + num_cyt_nodes;
}
void Cell::get_Wall_Nodes_Vec(vector<Wall_Node*>& walls) {
	walls = wall_nodes;
	return;
}
void Cell::add_wall_node_vec(Wall_Node* curr) {
	wall_nodes.push_back(curr);
	return;
}
void Cell::get_Cyt_Nodes_Vec(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}
void Cell::reset_Cell_Progress(){
	this->Cell_Progress = 0;
	return;
}
void Cell::update_Cell_Progress_add_node(double& add_node_prog) {
	this->Cell_Progress_add_node = add_node_prog;
	return;
}
void Cell::update_Cell_Progress_div(int& div_time) {
	this->Cell_Progress_div = div_time;
	return;
}
void Cell::set_K_LINEAR(double& x, double& y) {
	this->K_LINEAR = Coord(x,y);
	return;
}
void Cell::get_Neighbor_Cells(vector<Cell*>& cells) {
	cells = neigh_cells;
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
void Cell::calc_WUS() {
	this->wuschel = 109.6*exp(-0.03127*cell_center.length());

/*-0.189814360326137*pow(cell_center.length(),2) + 11.2927997842538854*cell_center.length() + -90.5308519096016;
	if(wuschel < 0) {
		wuschel = 0;
	}
*/
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
void Cell::set_growth_rate() {
	//this->growth_rate = 1000;
	if(this->wuschel < 30){
		this->growth_rate = unifRandInt(700, 1000);;
	}
	else if((this->wuschel >= 30) &&(this->wuschel < 50)) {
		this->growth_rate = unifRandInt(1000, 1300);
	}
	else if((this->wuschel >= 50) && (this->wuschel <70)){
		this->growth_rate = unifRandInt(1300, 1600);
	}
	else if ((this->wuschel >= 70) && (this->wuschel < 90)){
		this->growth_rate = unifRandInt(1600,1900);
	}
	else if ((this->wuschel >= 90) && (this->wuschel < 110)){
		this->growth_rate = unifRandInt(1900,2100);
	}	
	else if ((this->wuschel >= 70) && (this->wuschel < 80)){
		this->growth_rate = unifRandInt(2100,2300);
	}
	else if ((this->wuschel >= 80) && (this->wuschel < 90)){
		this->growth_rate = unifRandInt(2300,2400);
	}
	else if ((this->wuschel >= 90) && (this->wuschel < 108)){
		this->growth_rate = unifRandInt(2400,2500);
	}
	else if(this->wuschel >=110) {
		this->growth_rate = unifRandInt(2500,2700);
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
	double prelim_threshold = 18;
	//double sec_threshold = 1;
	Cell* me = this;
	// iterate through all cells
	#pragma omp parallel
	{
		Cell* curr = NULL;
		Coord curr_Cent;
		Coord distance;
		//vector<Cell*>my_neighbs;
		#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < all_Cells.size(); i++) {
			curr = all_Cells.at(i);
			if (curr != me) {
				curr_Cent = curr->get_Cell_Center();
				// Check if cell centers are close enough together
				distance = me->cell_center - curr_Cent;
				//cout << "Distance = " << distance << endl;
				if ( distance.length() < prelim_threshold ) {
					#pragma omp critical
					neigh_cells.push_back(curr);
					//cout << rank << "has neighbor" << curr->get_Rank() << endl;
				}
			
			}
			//else you're pointing at yourself and shouldnt do anything
	
		}
	}	
	
	//cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}

void Cell::update_adhesion_springs() {
	//cout << "Cell level" << endl;
	//Wall_Node* curr = left_Corner;
	//Wall_Node* orig = curr;
	//Wall_Node* next = NULL;
	//do {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{	

		#pragma omp for schedule(static,1)	
		for(unsigned int i=0; i< walls.size();i++) {
			walls.at(i)->set_Closest(NULL, 100);
		//	walls.at(i)->clear_Closest_Vec();	
		}
	}
	//cout << "cleared" << endl;
	vector<Cell*>neighbors;
	this->get_Neighbor_Cells(neighbors);
	Wall_Node* curr_Node = NULL;
        Wall_Node* next_Node = NULL;
        Wall_Node* curr_Closest = NULL;
        double curr_len = 0;
	curr_Node = left_Corner;
	/*do {
		next_Node = curr_Node->get_Left_Neighbor();
                curr_Closest = curr_Node->find_Closest_Node(neighbors);
                curr_Node->make_Connection(curr_Closest);
                curr_Node = next_Node;
       	} while(next_Node != left_Corner);*/
	#pragma omp parallel 
	{
		Wall_Node* curr_Closest = NULL;
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			curr_Closest = walls.at(i)->find_Closest_Node(neighbors);
		//	cout << "found closest" << endl;			
			walls.at(i)->make_Connection(curr_Closest);
		//	cout << "made connection" << endl;
		}	
	}
	return;
}

//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
void Cell::calc_New_Forces(int Ti) {
	#pragma omp parallel for
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces();
	}

	//calc forces on wall nodes
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);

	//Wall_Node* curr = left_Corner; 
	//Wall_Node* orig = curr;
	//int counter = 0;
	//do {
	#pragma omp parallel
	{
		Wall_Node* curr;
		//counter++;
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			curr = walls.at(i);
//			cout << "Wall node number: " << counter << endl;
			curr->calc_Forces(Ti);
			curr->set_Delete(0);
//			cout << "Forces calculated" << endl;
			//curr = curr->get_Left_Neighbor();
	
		} //while (curr != orig);
	//cout << "out of forces function" << endl;
	}
	return;
}

void Cell::update_Node_Locations() {
	//update cyt nodes
	#pragma omp parallel 
	{
		double new_damping;
		#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
			  new_damping = cyt_nodes.at(i)->get_My_Cell()->get_Damping();
			  cyt_nodes.at(i)->update_Location(new_damping);
		}	
	}

	//update wall nodes
	//Wall_Node* curr = left_Corner;
	//Wall_Node* orig = left_Corner;
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	//do {
	#pragma omp parallel 
	{	
		double new_damping;
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i< walls.size();i++) {
			new_damping = walls.at(i)->get_My_Cell()->get_Damping();
	//	cout << "update locaation" << endl;
			walls.at(i)->update_Location(new_damping);
		}
	}
	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	update_Wall_Angles();
//	cout << "done" << endl;
	return;
}

void Cell::update_Wall_Angles() {
//	cout << "wall angles" << endl;
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
//	Wall_Node* curr = left_Corner;
//	Wall_Node* orig = curr;
//	do {
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i< walls.size();i++) {
		//cout<< "updating" <<endl;
		walls.at(i)->update_Angle();
	}
//	cout << "Success" << endl;
	return;
}

void Cell::update_Wall_Equi_Angles() {
//	cout << "equi angles" << endl;
	double new_equi_angle = (num_wall_nodes-2)*pi/num_wall_nodes;
	//Wall_Node* curr = left_Corner;
	//Wall_Node* orig = left_Corner;
	
	//do {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i < walls.size(); i++) {
		walls.at(i)->update_Equi_Angle(new_equi_angle);
	}
	return;
}

void Cell::update_Cell_Center() {
	//Wall_Node* curr = left_Corner;
	//Wall_Node* orig = curr;
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);

	Coord total_location = Coord();
//	#pragma omp parallel
//	{
		Coord curr_loc;
//		#pragma omp for reduction(CoordPlus:total_location)
		for(unsigned int i=0;i<walls.size();i++) {
			//do {
			curr_loc = walls.at(i)->get_Location();
			total_location += curr_loc;
			//curr = curr->get_Left_Neighbor();
		} //while(curr != orig);
//	}
	cell_center = total_location*(1.0/num_wall_nodes);
	
	return;
}
//=====================================================================
//==========================================
// Growth of Cell
//==========================================
//=====================================================================
void Cell::update_Cell_Progress(int& Ti) {
	//update life length of the current cell
	this->update_Life_Length();

	//variables needed if division occurs
	Cell* new_Cell= NULL;
	vector<Cell*> cells;
	vector<Cell*> neighbor_cells;
	this->my_tissue->get_Cells(cells);
	int number_cells = cells.size();
	//variables for determining growth rate
	//if(Ti >= 0) {
		if(Ti%growth_rate == (growth_rate -1)) {
			this->add_Cyt_Node();
			Cell_Progress = Cell_Progress + 1;
		}
	if((this->Cell_Progress > 30)) {//&& (this->calc_Area() > 50)){
	//	if(this->rank == 41) {
		//cout << "Cell Prog" << Cell_Progress << endl;
		new_Cell = this->divide();
		cout << "division success" << endl;
		if(new_Cell == NULL) {
			cout << "womp womp" << endl;
		}
		my_tissue->update_Num_Cells(new_Cell);
		//setting info about new cell
		new_Cell->set_Rank(number_cells);
		//cout << "set rank" << endl;
		//cout << "Parent rank: " << this->rank << endl;
		//cout << "sister rank: " << new_Cell->get_Rank() << endl;
		//layer in division function		
		//damping in division function
		//life length set to 0 in constructor
		//cyt nodes in divison function
		//wall nodes in division function
		//all cell progress set to 0 in constructor
		new_Cell->update_Cell_Progress_div(Ti);
		//cell center in division function
		//cyt and wus in division function
		//k linear in division function
		//left corner in divison function  
		new_Cell->update_Neighbor_Cells();
		//cout << "new cell neighbor update" << endl;
		new_Cell->update_adhesion_springs();
		//cout << "new cell adhesion update" << endl;
		new_Cell->get_Neighbor_Cells(neighbor_cells);
		#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i< neighbor_cells.size(); i++) {
			neighbor_cells.at(i)->update_adhesion_springs();
		}
		this->life_length = 0;	
		this->Cell_Progress =0;		
		this->Cell_Progress_add_node = 0;
		this->Cell_Progress_div = Ti; 
		this->update_Neighbor_Cells();
		this->update_adhesion_springs();
		//cout << "updated neighbor" << endl;
		this->get_Neighbor_Cells(neighbor_cells);
		#pragma omp parallel for schedule(static, 1)
		for(unsigned int i=0; i < neighbor_cells.size(); i++) {
			neighbor_cells.at(i)->update_adhesion_springs();
		}
	//	}
	}
	//if no division check if internal node should be added
	//else {
	
	/*if((Cell_Progress - Cell_Progress_add_node) >= 1) { 
    			this->add_Cyt_Node();
			cout << "Progress" << (Cell_Progress - Cell_Progress_add_node) << endl;  
	
			cout << "added cyt node" << this->rank << endl;
			//cout << curr_area << " is area" << endl;
			cout << "time: " << Ti << endl;
			Cell_Progress_add_node = Cell_Progress;
	}
	}
	*/
		
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
void Cell::add_wall_Node_Check() {
	//cout << "adding a wall node" << endl;
	add_Wall_Node();
	return;
}
void Cell::delete_wall_Node_Check(){
	delete_Wall_Node();
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
	//	cout << "wasnt null" << endl;
		left = right->get_Left_Neighbor();
		location  = (right->get_Location() + left->get_Location())*0.5;
		added_node = new Wall_Node(location, this, left, right);
		wall_nodes.push_back(added_node);
		//cout << "made new node" << endl;
		right->set_Left_Neighbor(added_node);
		left->set_Right_Neighbor(added_node);
		num_wall_nodes++;
		update_Wall_Equi_Angles();
		update_Wall_Angles();
		//vector<Cell*> neighbs;
		//this->get_Neighbor_Cells(neighbs);
		//Wall_Node* curr_Closest = added_node->find_Closest_Node(neighbs);
		//added_node->make_Connection(curr_Closest);
	}
	else {
		//cout << "null" << endl;
	}
	return;
}
void Cell::delete_Wall_Node() {
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	Wall_Node* small = NULL;
	vector<Cell*>neighbors;
	find_Smallest_Length(small);
	if(small !=NULL) {
	//	cout << "deletion" << endl;
		left = small->get_Left_Neighbor();
		right = small->get_Right_Neighbor();
	
	
		if(this->left_Corner == small) {
			this->set_Left_Corner(left);
		}
	
		delete small;
	
		left->set_Right_Neighbor(right);
		right->set_Left_Neighbor(left);
		left->set_Delete(1);
		right->set_Delete(1);
		this->wall_nodes.clear();
		this->num_wall_nodes = 0;
	
		Wall_Node* curr = this->left_Corner;
		Wall_Node* next = NULL;
		Wall_Node* orig = curr;
	
		do {
			this->wall_nodes.push_back(curr);
			next = curr->get_Left_Neighbor();
			this->num_wall_nodes++;
			curr = next;
		}while(next != orig);
	
		update_Wall_Equi_Angles();
		update_Wall_Angles();
		this->is_deleted = true;
		//update_adhesion_springs();
		//this->get_Neighbor_Cells(neighbors);
		//cout<< "neighbors loop" << endl;
		//for(unsigned int i =0; i< neighbors.size();i++){
		//	neighbors.at(i)->update_adhesion_springs();
		//}
		//cout << "updated" << endl;
	}
	return;
}

void Cell::find_Smallest_Length(Wall_Node*& right) {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 100;
	#pragma omp parallel
	{
		Wall_Node* left_neighbor;
		double curr_len = 0;
		#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < walls.size();i++) {
			left_neighbor = walls.at(i)->get_Left_Neighbor();
			curr_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
			if(curr_len < .06){
				if(curr_len < max_len) {
					#pragma omp critical
					max_len = curr_len;
					right = walls.at(i);
				}
			}
		}
	}
	return;
}
void Cell::find_Largest_Length(Wall_Node*& right) {
	vector<Wall_Node*> walls; 
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 0;
//	Wall_Node* right = NULL;
	#pragma omp parallel 
	{
		Wall_Node* left_neighbor;
		double curr_len = 0;
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			left_neighbor = walls.at(i)->get_Left_Neighbor();
			curr_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
			if(curr_len > MEMBR_THRESH_LENGTH) {			
				if(curr_len > max_len) {
					#pragma omp critical
					max_len = curr_len;
					right = walls.at(i);
				}
			}
		}		
	//cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;
	}
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

void Cell::add_Cyt_Node_Div(double radius_x,double radius_y) {
	//USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
	// ---distributes more evenly throughout start cell
	double offset = .8;
	Coord location;
	Cyt_Node* cyt;
	double x;
	double y;
	double rand_radius_x = unifRand(0.0,1.0)*offset*radius_x; 
	double rand_angle = unifRand(0.0,1.0)*2*pi;
	double rand_radius_y = unifRand(0.0,1.0)*offset*radius_y;
	x = cell_center.get_X()+ rand_radius_x*cos(rand_angle);
	y = cell_center.get_Y()+ rand_radius_y*sin(rand_angle);
	location = Coord(x,y);
	//cout << location << endl;
	cyt = new Cyt_Node(location,this);
	cyt_nodes.push_back(cyt);
	num_cyt_nodes++;
	return;
}
void Cell::compute_Main_Strain_Direction(double& x_length, double& y_length) {
	//double x_length = 0;
	//double y_length = 0;
	vector<Wall_Node*>walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{	
		#pragma omp for reduction(+:x_length) reduction(+:y_length) schedule(static,1)
		for(unsigned int i = 0; i< walls.size();i++) {
			x_length += fabs(walls.at(i)->get_Location().get_X() - walls.at(i)->get_Left_Neighbor()->get_Location().get_X());
			y_length += fabs(walls.at(i)->get_Location().get_Y() - walls.at(i)->get_Left_Neighbor()->get_Location().get_Y());
		}
	}
	return;
}
		
void Cell::find_Largest_Length_Div(Wall_Node*& right, Wall_Node*& second_right) {
	vector<Wall_Node*> walls; 
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 0;
//	Wall_Node* right = NULL;
	#pragma omp parallel 
	{
		Wall_Node* left_neighbor;
		double curr_len = 0;
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			left_neighbor = walls.at(i)->get_Left_Neighbor();
			curr_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
			if(curr_len > max_len) {
				#pragma omp critical
				max_len = curr_len;
				right = walls.at(i);
			}
		}
	}
	Wall_Node* start = right;	
	Wall_Node* end = right;
	Wall_Node* left_neighbor = NULL;
	double second_max_len = 0;
	double curr_len = 0;
	for(unsigned int i = 0; i < 150; i++) {
		start = start->get_Left_Neighbor();
		end = end->get_Right_Neighbor();
	}
	do {
		left_neighbor = start->get_Left_Neighbor();
		curr_len = (start->get_Location()-left_neighbor->get_Location()).length();
		if(curr_len > second_max_len) {
			second_max_len = curr_len;
			second_right = start;
		}
		start = start->get_Left_Neighbor();
	} while(left_neighbor != end);
		
		
	//cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;
	return;
}
Wall_Node* Cell::find_closest_node_top() {
	double x_val_center = this->cell_center.get_X();
	double y_val_center = this->cell_center.get_Y();
	Wall_Node* curr = this->left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* closest = NULL;
	double  curr_diff = 0;
	double smallest_diff =100;
	do {
		if(curr->get_Location().get_Y() < y_val_center) {
			curr_diff = fabs(x_val_center - curr->get_Location().get_X()); 
			if(curr_diff < smallest_diff) {
				smallest_diff = curr_diff;
				closest = curr;
			}
		}
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);	
	return closest;
}
Wall_Node* Cell::find_closest_node_bottom() {
	double x_val_center = this->cell_center.get_X();
	double y_val_center = this->cell_center.get_Y();
	Wall_Node* curr = this->left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* closest = NULL;
	double  curr_diff = 0;
	double smallest_diff = 100;
	do {
		if(curr->get_Location().get_Y() > y_val_center) {
			curr_diff = fabs(x_val_center - curr->get_Location().get_X()); 
			if(curr_diff < smallest_diff) {
				smallest_diff = curr_diff;
				closest = curr;
			}
		}
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);	
	return closest;
}
Wall_Node* Cell::find_closest_node_left() {
	double x_val_center = this->cell_center.get_X();
	double y_val_center = this->cell_center.get_Y();
	Wall_Node* curr = this->left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* closest = NULL;
	double  curr_diff = 0;
	double smallest_diff = 100;
	do {
		if(curr->get_Location().get_X() < x_val_center) {
			curr_diff = fabs(y_val_center - curr->get_Location().get_Y()); 
			if(curr_diff < smallest_diff) {
				smallest_diff = curr_diff;
				closest = curr;
			}
		}
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);	
	return closest;
}
Wall_Node* Cell::find_closest_node_right() {
	double x_val_center = this->cell_center.get_X();
	double y_val_center = this->cell_center.get_Y();
	Wall_Node* curr = this->left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* closest = NULL;
	double  curr_diff = 0;
	double smallest_diff = 100;
	do {
		if(curr->get_Location().get_X() > x_val_center) {
			curr_diff = fabs(y_val_center - curr->get_Location().get_Y()); 
			if(curr_diff < smallest_diff) {
				smallest_diff = curr_diff;
				closest = curr;
			}
		}
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);	
	return closest;
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






