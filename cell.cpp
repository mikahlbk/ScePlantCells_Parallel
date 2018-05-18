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
	//growth_rate assigned in div function
	//cyt nodes reassigned in division function
	//start at zero
	num_cyt_nodes = 0;
	//wall nodes renumbered in division function
	num_wall_nodes = 0;
	//center calculate in division function
	Cell_Progress = 0;
	//will calculate signals in div function
	wuschel = 0;
	cytokinin = 0;
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
	this->calc_WUS();
	this->calc_CYT();
	this->set_growth_rate();
	
	if((this->layer == 1)||(this->layer == 2)) {
                 this->growth_direction = Coord(1,0);
         }
         else {
                 this->growth_direction = Coord(0,1);
         }	
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
	
	//reindex wall node vector for vertical growing cells
	if(this->growth_direction == Coord(0,1)){
		for(int i = 0; i<35; i++){
			this->left_Corner = left_Corner->get_Left_Neighbor();
		}
	Wall_Node* curr = left_Corner;
	wall_nodes.clear();
	do{
		wall_nodes.push_back(curr);
		curr = curr->get_Left_Neighbor();
	} while(curr != left_Corner);
	}
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
        double k_lin = 0;
        double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	//give each node its spring constant and membrane length
	for(unsigned int i = 0; i < walls.size();i++) {	
		curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
		curr_len = curr_vec.length();	
		growth_len = this->growth_direction.length();
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		walls.at(i)->set_membr_len(.07);
		}
		//k_lin = 250;
		k_lin = 150 + 500*(1-pow(costheta,2));
		walls.at(i)->set_K_LINEAR(k_lin);
	}
	//insert cytoplasm nodes
	int num_init_cyt_nodes = Init_Num_Cyt_Nodes + Cell_Progress;
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
		wall_nodes.pop_back();
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
	walls = this->wall_nodes;
	return;
}
void Cell::add_wall_node_vec(Wall_Node* curr) {
	this->wall_nodes.push_back(curr);
	this->num_wall_nodes++;
	return;
}
void Cell::get_Cyt_Nodes_Vec(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}
void Cell::reset_Cell_Progress(int cyt_size){
	this->Cell_Progress = cyt_size;
	return;
}
void Cell::update_cyt_node_vec(Cyt_Node* new_node){
	this->cyt_nodes.push_back(new_node);
	this->num_cyt_nodes++;
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
	//wildtype values
	this->wuschel = 149.7*exp(-0.01042*cell_center.length());
	//signaling domain increased by 2
	//this->wuschel = 109.6*exp(-0.012*cell_center.length());
	//this->wuschel = 120.6*exp(-0.01127*cell_center.length());
	
	return;
}

void Cell::calc_WUSwildtype() {
        //wildtype values
        this->wuschel = 88.35*exp(-0.008762*cell_center.length());
     	return;
}
void Cell::calc_WUSBAP12hr() {
	this->wuschel = 113.2*exp(-0.01374*cell_center.length());
	return;
}
void Cell::calc_WUSBAP24hr() {
        this->wuschel = 41.94*exp(-0.003069*cell_center.length());
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
void Cell::set_growth_rate() {
	//this->growth_rate = 2500;
	if(this->wuschel < 30){
		this->growth_rate = unifRandInt(700, 1200);;
	}
	else if((this->wuschel >= 30) &&(this->wuschel < 50)) {
		this->growth_rate = unifRandInt(1200,1700);
	}
	else if((this->wuschel >= 50) && (this->wuschel <70)){
		this->growth_rate = unifRandInt(1700,2200);
	}
	else if ((this->wuschel >= 70) && (this->wuschel < 90)){
		this->growth_rate = unifRandInt(2200,2700);
	}
	else if ((this->wuschel >= 90) && (this->wuschel < 110)){
		this->growth_rate = unifRandInt(2700,3200);
	}	
	else if ((this->wuschel >= 70) && (this->wuschel < 80)){
		this->growth_rate = unifRandInt(3200,3700);
	}
	else if ((this->wuschel >= 80) && (this->wuschel < 90)){
		this->growth_rate = unifRandInt(3700,4200);
	}
	else if ((this->wuschel >= 90) && (this->wuschel < 108)){
		this->growth_rate = unifRandInt(4200,4700);
	}
	else if(this->wuschel >=110) {
		this->growth_rate = unifRandInt(4700,5200);
	}

	return;
}
void Cell::set_growth_direction(Coord gd){
	this->growth_direction = gd;
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
	double prelim_threshold = 25;
	//double sec_threshold = 1;
	Cell* me = this;
	// iterate through all cells
//	#pragma omp parallel
//	{
		Cell* curr = NULL;
		Coord curr_Cent;
		Coord distance;
		//vector<Cell*>my_neighbs;
//		#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < all_Cells.size(); i++) {
			curr = all_Cells.at(i);
			if (curr != me) {
				curr_Cent = curr->get_Cell_Center();
				// Check if cell centers are close enough together
				distance = me->cell_center - curr_Cent;
				//cout << "Distance = " << distance << endl;
				if ( distance.length() < prelim_threshold ) {
//					#pragma omp critical
					neigh_cells.push_back(curr);
					//cout << rank << "has neighbor" << curr->get_Rank() << endl;
				}
			
			}
			//else you're pointing at yourself and shouldnt do anything
	
		}
//	}	
	
	//cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}

/*void Cell::update_adhesion_springs() {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{	

//		#pragma omp for schedule(static,1)	
		for(unsigned int i=0; i< walls.size();i++) {
			walls.at(i)->set_Closest(NULL, 100);
		}
	}
	//cout << "cleared" << endl;
	vector<Cell*>neighbors;
	this->get_Neighbor_Cells(neighbors);
//	#pragma omp parallel 
//	{
		Wall_Node* curr_Closest = NULL;
//		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			curr_Closest = walls.at(i)->find_Closest_Node(neighbors);
		//	cout << "found closest" << endl;			
			walls.at(i)->make_Connection(curr_Closest);
		//	cout << "made connection" << endl;
		}	
//	}
	return;
}
void Cell::update_microfibril_springs() {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{
		for(unsigned int i=0; i < walls.size();i++) {
			walls.at(i)->set_microfibril_pair(NULL,100);
		}
	}
	//divide nodes into left right or top bottom 
	//depending on growth direction
	//if anticlinal split top bottom
	vector<Wall_Node*>top;
	vector<Wall_Node*>bottom;
	vector<Wall_Node*>left;
	vector<Wall_Node*>right;
	vector<Wall_Node*>side1;
	vector<Wall_Node*>side2;
	double len1;
	double len2;
	if(this->get_growth_direction() == Coord(1,0)) {
		//top bottom split
		this->make_top_bottom_vectors(top,bottom);
		len1 = top.size();
		cout << "lengths top" << endl;
		cout << len1 << endl;
		len2 = bottom.size();
		cout << len2 << endl;
		if(len1<len2) {
			side1 = top;
			side2 = bottom;
		}
		else {
			side1 = bottom;
			side2 = top;
		}
	for(unsigned int i = 0; i< side1.size();i++) {
		side1.at(i)->find_microfibril_pair_horiz(side2);
	}
	
	}
	//if periclinal split left right
	else{
		///left right split
		this->make_left_right_vectors(left,right);
		cout << "lengths right" << endl;
		len1 = right.size();
		cout << len1 << endl;
		len2 = left.size();
		cout << len2 << endl;
		if(len1<len2) {
			side1 = right;
			side2 = left;
		}
		else {
			side1 = left;
			side2 = right;
		}
	for(unsigned int i = 0; i< side1.size();i++) {
		side1.at(i)->find_microfibril_pair_vert(side2);
	}
	}
	
	//take side with the least and serach through other side to connect
	//for(unsigned int i = 0; i< side1.size();i++) {
	//	side1.at(i)->find_microfibril_pair(side2);
	//}
	for(unsigned int i = 0; i< walls.size() ; i++){
		cout << walls.at(i)->get_micro_pair() << endl;
	}
	//cout << "cleared" << endl
	return;
}
void Cell::make_top_bottom_vectors(vector<Wall_Node*>&top, vector<Wall_Node*>& bottom) {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	
	for(unsigned int i = 0; i< walls.size();i++){
		curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
		curr_len = curr_vec.length();	
		growth_len = this->growth_direction.length();
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		if((theta <=.785398)||(theta >= 2.35619)){
			if((walls.at(i)->get_Location().get_Y() > this->cell_center.get_Y())) {
				top.push_back(walls.at(i));
			}
			else {
				bottom.push_back(walls.at(i));;
			}
		}
	}
	
	return;
}
void Cell::make_left_right_vectors(vector<Wall_Node*>&left, vector<Wall_Node*>&right) {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	Wall_Node* start = left_Corner;
	//for(int i = 0; i< 20; i++){
	//	start = start->get_Left_Neighbor();
	//}
	Wall_Node* curr = start;
	
	//for(unsigned int i = 0; i< walls.size();i++){
	do{
		curr_vec = curr->get_Left_Neighbor()->get_Location() - curr->get_Location();
		curr_len = curr_vec.length();	
		growth_len = this->growth_direction.length();
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		if((theta <=.785398)||(theta >=2.35619)){
			if((curr->get_Location().get_X() > this->cell_center.get_X())) {
				right.push_back(curr);
			}
			else {
				left.push_back(curr);
			}
		}
	curr = curr->get_Left_Neighbor();
	}while(curr != start);
	
	return;
}*/
//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
void Cell::calc_New_Forces(int Ti) {
	
	#pragma omp parallel for
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces(Ti);
	}

	//calc forces on wall nodes
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);

	#pragma omp parallel
	{
		Wall_Node* curr;
		//counter++;
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			curr = walls.at(i);
		//	cout << "Wall node number: " << counter << endl;
			curr->calc_Forces(Ti);
		}	
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
			//cout << "update locaation" << endl;
			walls.at(i)->update_Location(new_damping);
		}
	}
	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	update_Wall_Angles();
	//cout << "done" << endl;

	//if((this->rank == 37)|| (this->rank == 57)){
	//cout << "update locations" << endl;
	//int counter = 0;
	//for(unsigned int i = 0; i < wall_nodes.size(); i++) {
	//	cout << wall_nodes.at(i)->get_Location() << endl;
	//	counter++;
	//}
	//cout << counter << endl;
	//}

	return;
}

void Cell::update_Wall_Angles() {
//	cout << "wall angles" << endl;
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	
	//#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i< walls.size();i++) {
		//cout<< "updating" <<endl;
		walls.at(i)->update_Angle();
	}
//	cout << "Success" << endl;
	return;
}

void Cell::update_Wall_Equi_Angles() {
//	cout << "equi angles" << endl;
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);
	//#pragma omp parallel for schedule(static,1)
	double k_lin = 0;
        double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	double new_equi_angle = 0; 
	double circle_angle  = (num_wall_nodes-2)*pi/num_wall_nodes;
	/*for(unsigned int i = 0; i < walls.size();i++) {	
		curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
		curr_len = curr_vec.length();	
		growth_len = this->growth_direction.length();
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		new_equi_angle = circle_angle*(1-pow(cos(theta),2)) + pi*pow(cos(theta),2);
		walls.at(i)->update_Equi_Angle(new_equi_angle);
	}*/
	for(unsigned int i=0; i < walls.size(); i++) {
		walls.at(i)->update_Equi_Angle(circle_angle);
	}
	return;
}

void Cell::update_Cell_Center() {
	vector<Wall_Node*> walls;
	this->get_Wall_Nodes_Vec(walls);

	Coord total_location = Coord();
	#pragma omp parallel
	{
		Coord curr_loc;
		#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
		#pragma omp for reduction(+:total_location) schedule(static,1)
		for(unsigned int i=0;i<walls.size();i++) {
			curr_loc = walls.at(i)->get_Location();
			total_location += curr_loc;
		}
	}
	this->cell_center = total_location*((1.0/static_cast<double>(num_wall_nodes)));
	//cout << cell_center<<"center"<< endl;
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
	if(Ti%growth_rate == (growth_rate -1)) {
		this->add_Cyt_Node();
		Cell_Progress++;
	}
	if(Ti==20000){// && (this->calc_Area() > 50)){
		//if((this->rank == 37)){//||(this->rank ==2)||(this->rank ==1)||(this->rank ==0)) {
	//	cout << "Cell Prog" << Cell_Progress << endl;
		//for(unsigned int i =0; i < wall_nodes.size(); i++) {
		//	cout << wall_nodes.at(i)->get_Location() << endl;
		//} 
		new_Cell = this->divide();
		cout << "division success" << endl;
		if(new_Cell == NULL) {
			cout << "womp womp" << endl;
		}
		//cout << "add cell" << endl;
		this->my_tissue->update_Num_Cells(new_Cell);
		//setting info about new cell
	
		new_Cell->set_Rank(number_cells);
		cout << "set rank" << endl;
		cout << "Parent rank: " << this->rank << endl;
		cout << "sister rank: " << new_Cell->get_Rank() << endl;
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
		//new_Cell->update_adhesion_springs();
		//cout << "new cell adhesion update" << endl;
		//new_Cell->get_Neighbor_Cells(neighbor_cells);
		//cout << "pre sister" << endl;
		//for(unsigned int i =0; i < wall_nodes.size(); i++) {
		//	cout << wall_nodes.at(i)->get_Location() << endl;
		//}
		//for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		//	cout << cyt_nodes.at(i)->get_Location() << endl;
		//}
		//#pragma omp parallel for schedule(static,1)
		//for(unsigned int i = 0; i< neighbor_cells.size(); i++) {
		//	neighbor_cells.at(i)->update_adhesion_springs();
		//}
		this->Cell_Progress_div = Ti; 
		
		this->get_Tissue()->update_Neighbor_Cells();
	//	this->get_Tissue()->update_Adhesion(Ti);
	//	cout << "at end of cell" << endl;
	//	cout << "this" << endl;
	//	vector<Wall_Node*>walls;
	//	new_Cell->get_Wall_Nodes_Vec(walls);
	//	int counter = 0;
	//	for(unsigned int i =0; i < wall_nodes.size(); i++) {
	//		cout << wall_nodes.at(i)->get_Location() << endl;
	//		counter++;
	//	}
	//	cout << counter << endl;
	///	counter = 0;
	//	cout << "sister" << endl;
		
	//	for(unsigned int i =0; i < walls.size(); i++) {
	//		cout << walls.at(i)->get_Location() << endl;
	//		counter++;
	//	}
	//	cout << counter << endl;
		//}
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
void Cell::nematic(Coord& avg_vec, double& angle){
	Coord curr_vec;
	int counter = 0;
	for(unsigned int i = 0; i< wall_nodes.size(); i++){
		curr_vec+= wall_nodes.at(i)->get_Location()-cell_center;
		counter++;
	}
	avg_vec = curr_vec/static_cast<double>(counter);
	double curr_len = avg_vec.length();
	Coord horizontal = Coord(1,0);
	double costheta = horizontal.dot(avg_vec)/(curr_len);
	angle = acos( min( max(costheta,-1.0), 1.0) );
	
	return;
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
	double new_length;
	Wall_Node* added_node = NULL;
	if(right != NULL) {
	//	cout << "wasnt null" << endl;
		left = right->get_Left_Neighbor();
		location  = (right->get_Location() + left->get_Location())*0.5;
		added_node = new Wall_Node(location, this, left, right);
		wall_nodes.push_back(added_node);
		cout << "made new node" << endl;
		right->set_Left_Neighbor(added_node);
		left->set_Right_Neighbor(added_node);
		new_length = (right->get_membr_len() + left->get_membr_len())*.5;
		num_wall_nodes++;
		update_Wall_Equi_Angles();
		update_Wall_Angles();
		Coord curr_vec = added_node->get_Left_Neighbor()->get_Location() - added_node->get_Location();
		double curr_len = curr_vec.length();	
		double growth_len = this->growth_direction.length();
		double costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		double theta = acos( min( max(costheta,-1.0), 1.0) );
		double k_lin = 150+ 500*(1-pow(cos(theta),2));
		added_node->set_K_LINEAR(k_lin);
		added_node->set_membr_len(new_length);
		//vector<Cell*> neighbs;
		added_node->set_is_new(true);
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
		//if(rank == 13) {
		//cout << "rank" << rank << endl;
		//cout << "deletion" << endl;
		//}
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
			num_wall_nodes++;
			curr = next;
		}while(next != orig);
		//if(rank ==13){
		//cout << "Wall nodes after delete" << num_wall_nodes << endl;
		//}
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
//finds right neighbor node of smallest length on membrane
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
//finds right neighbor node of largest length on membrane
void Cell::find_Largest_Length(Wall_Node*& right) {
	vector<Wall_Node*> walls; 
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 0;
	Wall_Node* right_side = NULL;
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
					right_side = walls.at(i);
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
	double offset = .5;
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
	this->cyt_nodes.push_back(cyt);
	this->num_cyt_nodes++;
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
	Wall_Node* right_side = NULL;
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
				right_side = walls.at(i);
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
		if(curr->get_Location().get_Y() > y_val_center) {
			curr_diff = fabs(x_val_center - curr->get_Location().get_X()); 
			if(curr_diff < smallest_diff) {
				smallest_diff = curr_diff;
				closest = curr;
			}
		}
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);	
	//cout << "closest" << closest << endl;
	return closest;
}
Wall_Node* Cell::find_closest_node_bottom() {
	double x_val_center = this->cell_center.get_X();
	double y_val_center = this->cell_center.get_Y();
	Wall_Node* curr = this->left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* closest = NULL;
	double  curr_diff = 0;
	double smallest_diff = 100.0;
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
	
	//cout << "closest" << closest << endl;
	return closest;
}
Wall_Node* Cell::find_closest_node_left() {
	double x_val_center = this->cell_center.get_X();
	double y_val_center = this->cell_center.get_Y();
	Wall_Node* curr = this->left_Corner;
	Wall_Node* orig = curr;
	Wall_Node* closest = NULL;
	double  curr_diff = 0;
	double smallest_diff = 100.0;
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
	
	//cout << "closest" << closest << endl;
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
	
	//cout << "closest" << closest << endl;

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
double Cell::average_Pressure(){
	//force 
	double force;
	double average_pressure;
	double area = 0;
	Wall_Node* curr_wall = left_Corner;
	Wall_Node* neighbor = NULL;
	do {
		force += (curr_wall->get_CytForce()).length();
		neighbor = curr_wall->get_Left_Neighbor();
		area += (curr_wall->get_Location()-neighbor->get_Location()).length();
		curr_wall = neighbor;
		
	} while (curr_wall != left_Corner);
	average_pressure = force/area;
	return average_pressure;
}
void Cell::wall_Pressure(){
	int num_nodes = this->get_wall_count();
	int block = num_nodes/20;
	double force = 0;
	double perim = 0;
	double average_pressure;
	for(unsigned int i=0; i< 21; i++){
		force += (wall_nodes.at(i)->get_CytForce()).length();
		perim += (wall_nodes.at(i)->get_Location()-wall_nodes.at(i)->get_Left_Neighbor()->get_Location()).length();
	}
	average_pressure = force/perim;
	for(unsigned int i = 0; i<21;i++){
		 wall_nodes.at(i)->set_pressure(average_pressure);
	}
	for(unsigned int i=21; i< 41; i++){
		force +=  (wall_nodes.at(i)->get_CytForce()).length();
		perim += (wall_nodes.at(i)->get_Location()-wall_nodes.at(i)->get_Left_Neighbor()->get_Location()).length();
	}
	average_pressure = force/perim;
	for(unsigned int i = 21; i<41;i++){
		 wall_nodes.at(i)->set_pressure(average_pressure);
	}
	for(unsigned int i=41; i< 61; i++){
		force += (wall_nodes.at(i)->get_CytForce()).length();
		perim += (wall_nodes.at(i)->get_Location()-wall_nodes.at(i)->get_Left_Neighbor()->get_Location()).length();
	}
	average_pressure = force/perim;
	for(unsigned int i = 41; i<61;i++){
		 wall_nodes.at(i)->set_pressure(average_pressure);
	}
	for(unsigned int i=61; i< 81; i++){
		force += (wall_nodes.at(i)->get_CytForce()).length();
		perim += (wall_nodes.at(i)->get_Location()-wall_nodes.at(i)->get_Left_Neighbor()->get_Location()).length();
	}
	average_pressure = force/perim;
	for(unsigned int i = 61; i<81;i++){
		 wall_nodes.at(i)->set_pressure(average_pressure);
	}
	for(unsigned int i=81; i< wall_nodes.size(); i++){
		force += (wall_nodes.at(i)->get_CytForce()).length();
		perim += (wall_nodes.at(i)->get_Location()-wall_nodes.at(i)->get_Left_Neighbor()->get_Location()).length();
	}
	average_pressure = force/perim;
	for(unsigned int i = 81; i< wall_nodes.size();i++){
		 wall_nodes.at(i)->set_pressure(average_pressure);
	}
	
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
			cout << "updated" << endl;
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
void Cell::print_VTK_Scalars_Wall_Pressure(ofstream& ofs){
	this->wall_Pressure();
	Wall_Node* curr_wall = this->left_Corner;
	

	double pressure;
	do{
		pressure = curr_wall->get_pressure();
		curr_wall = curr_wall->get_Left_Neighbor();
		ofs << pressure << endl;
	}while(curr_wall!= left_Corner);
	for(unsigned int i=0; i<cyt_nodes.size(); i++){
		ofs << 0 << endl;
	}
	return;
}
void Cell::print_VTK_Scalars_Average_Pressure(ofstream& ofs) {
	double pressure = this->average_Pressure();
	for(unsigned int i=0;i<wall_nodes.size();i++){
		ofs << pressure << endl;
	}
	for(unsigned int i=0;i<cyt_nodes.size();i++){
		ofs<< pressure << endl;
	}
	return;
}
void Cell::print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs) {
	double pressure = this->average_Pressure();
	for(unsigned int i=0;i<1;i++){
		ofs << pressure << endl;
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
void Cell::print_VTK_Scalars_WUS_cell(ofstream& ofs) {
	for(unsigned int i = 0; i < 1; i++) {
		
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
void Cell::print_VTK_Scalars_Node(ofstream& ofs) {
	Wall_Node* currW = left_Corner;
	double color;
	do {
		if(currW->get_color()){
			color = 30.0;
		}
		else{
			color = 0.0;
		}
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		color = 0.0;
		ofs << color << endl;
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






