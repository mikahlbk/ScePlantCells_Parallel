//cell.cpp:
//===================
// Forward Declarations
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "rand.h"
#include "cell.h"
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors
// this constructor is used in
// the divsion function
Cell::Cell(Tissue* tissue) {
	my_tissue = tissue;
	//rank assigned in division function
	//layer inherited from parent	
	//boundary cells don't divide
	boundary = 0;
	//stem cells don't divide 
	stem = 0;
	//damping assigned in div function
	//just divided so reset life length
	life_length = 0;
	//cyt nodes reassigned in division function
	//start at zero
	num_cyt_nodes = 0;
	//wall nodes renumbered in division function
	//start at zero
	num_wall_nodes = 0;
	//cell progress decided in div function
	Cell_Progress = 0;
	//center calculate in division function
	//will calculate signals in div function
	wuschel = 0;
	cytokinin = 0;
	//growth_rate assigned in div function
	//growth direction assigned in division
	//neighbors assigned in div function
	//left corner assigned in division
}
//this constructor is used to initialize first set of cells
//calls set_growth_rate which detemrines growth rate based on WUS CONC
Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer, int boundary, int stem)    {
	this->my_tissue = tiss;
	this->rank = rank;
	this->layer = layer;
	//if boundary is equal to one 
	//then the cell will have higher damping
	//which is assigned below
	//boundary conditions are read in from initial text file
	this->boundary = boundary;
	this-> stem = stem;
	//set damping for cells that act as anchor points
	if(this->stem == 1) {
		this->damping = STEM_DAMP;
	}
	else if((this->boundary == 1)){
		this->damping =  BOUNDARY_DAMP;
	}
	else{
		this->damping = REG_DAMP;
	}
	life_length = 0;
	//cyt nodes initialized in tissue constructor which
	//calls the makes nodes function on each new cell
	num_cyt_nodes = 0;
	//wall nodes initialized in tissue constructor which 
	//calls the make nodes function on each new cell
	num_wall_nodes = 0;
	Cell_Progress = unifRandInt(0,10);
	this->cell_center = center;
	//this gets reupdated after singal is assigned
	//in tissue constructor
	if(this->boundary == 1){
		this->growth_direction = Coord(0,0);
	}
	else if(this->stem == 1){
		this->growth_direction = Coord(0,1);
	}
        else if((this->layer == 1)||(this->layer == 2)) {
                 this->growth_direction = Coord(1,0);
        }
        else{
	 	this->growth_direction = Coord(0,1);
	}

	//cout << "layer" << this->layer << endl;
	//cout << "stem" << this->stem << endl;
	//cout << "boundary" << this-> boundary << endl;
	//cout << "gd" << this->growth_direction << endl;
	//cout << "damping" << this->damping << endl;
	//neighbors update function called after initialization
	//left corner set in make nodes function called by tissue constuctor
}
//calls compute membrane equi length for each node
//calls compute linear spring constant for each node
//calls compute bending spring constant for each node
//calls update wall equi angles for each node
//calls update wall angles to get initial angle of each node
void Cell::make_nodes(double radius){
	
	//assemble the membrane
	int num_Init_Wall_Nodes = Init_Wall_Nodes + 2*(Cell_Progress);
	double angle_increment = (2*pi)/num_Init_Wall_Nodes;
	
	//make all wall nodes
	double curr_X;
	double curr_Y;
	Coord location = this->cell_center;;
	double curr_theta = 0;
	curr_X = cell_center.get_X() + radius*cos(curr_theta);
	curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
	location = Coord(curr_X,curr_Y);
	//make the first node
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Wall_Node> prevW = make_shared<Wall_Node>(location,this_cell);
	shared_ptr<Wall_Node> currW = prevW;
	wall_nodes.push_back(prevW);
	num_wall_nodes++;
	shared_ptr<Wall_Node> orig(prevW);
	//this will be the "starter" node
	this->left_Corner = orig;
	//make successive nodes
	for(int i = 0; i<num_Init_Wall_Nodes-1; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = cell_center.get_X() + radius*cos(curr_theta);
		curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		shared_ptr<Wall_Node> new_node =make_shared<Wall_Node>(location,this_cell);
		currW = new_node;
		wall_nodes.push_back(currW);
		num_wall_nodes++;
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}
			
	//connect last node to starter node
	currW->set_Left_Neighbor(orig);
	orig->set_Right_Neighbor(currW);
	this->perimeter = this->get_curr_perimeter();
	//where is where most private member variables are set
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//damping inherited from cell
	double new_damping = this->get_Damping();
        //variable to hold membr_equ_len
	double l_thresh = 0;
	double k_lin = 0;
	double k_bend = 0;
        for(unsigned int i = 0; i < walls.size();i++) {	
		walls.at(i)->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(walls.at(i));
		k_lin = compute_k_lin(walls.at(i));
		k_bend = compute_k_bend(walls.at(i));
		walls.at(i)->set_membr_len(l_thresh);
		walls.at(i)->set_K_LINEAR(k_lin);
		walls.at(i)->set_K_BEND(k_bend);
	}	
	//double new_damping = this->get_Damping();
	//insert cytoplasm nodes
	int num_init_cyt_nodes = Init_Num_Cyt_Nodes + Cell_Progress;
	this->Cell_Progress = num_init_cyt_nodes;
	double scal_x_offset = 0.8;
	//Coord location;
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
		
		shared_ptr<Cyt_Node> cyt = make_shared<Cyt_Node>(location,this_cell);
		cyt->set_Damping(new_damping);
		cyt_nodes.push_back(cyt);
		num_cyt_nodes++;
	}
	
	//update equilibrium angle
	update_Wall_Equi_Angles();
	//update wall angles
	update_Wall_Angles();	
	//is_divided = false;
	return;
}

// Destructor
Cell::~Cell() {
	//not needed using smartpointers
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
void Cell::set_Damping(double new_damping) {
	this->damping = new_damping;
	return;
}
void Cell::update_Life_Length() {
	life_length++;
	return;
}
void Cell::reset_Life_Length(){
	life_length = 0;
	return;
}
int Cell::get_Node_Count() {
	return num_wall_nodes + num_cyt_nodes;
}
void Cell::get_Wall_Nodes_Vec(vector<shared_ptr<Wall_Node>>& walls) {
	walls = this->wall_nodes;
	return;
}
void Cell::add_wall_node_vec(shared_ptr<Wall_Node> curr) {
	this->wall_nodes.push_back(curr);
	this->num_wall_nodes++;
	return;
}
void Cell::get_Cyt_Nodes_Vec(vector<shared_ptr<Cyt_Node>>& cyts) {
	cyts = cyt_nodes;
	return;
}
void Cell::update_cyt_node_vec(shared_ptr<Cyt_Node> new_node){
	this->cyt_nodes.push_back(new_node);
	this->num_cyt_nodes++;
	return;
}
void Cell::reset_Cell_Progress(){
	this->Cell_Progress = 0;
	return;
}
void Cell::update_Cell_Progress() {
	this->Cell_Progress++;
	return;
}
void Cell::calc_WUS(Coord L1_AVG) {
	
	//new data from eric
	//CZ ~5 cells wide
	//layer 1
	/*if(this->rank == 0){
		this->wuschel = 11;
	}
	else if(this->rank == 1){
		this->wuschel = 17; 
	}
	else if(this->rank == 2){
		this->wuschel = 16; 
	}
	else if(this->rank == 5){
		this->wuschel = 12; 
	}
	else if(this->rank == 6){
		this->wuschel = 14;; 
	}
	//layer 2
	else if(this->rank == 9){
		this->wuschel = 9;
	}
	else if(this->rank == 10){
		this->wuschel = 18;
	}
	else if(this->rank == 11){
		this->wuschel = 15; 
	}
	else if(this->rank == 14){
		this->wuschel = 20;
	}
	else if(this->rank == 15){
		this->wuschel = 13; 
	}
	//layer 3
	else if(this->rank == 18){
		this->wuschel = 44; 
	}
	else if(this->rank == 19){
		this->wuschel = 23; 
	}
	else if(this->rank == 20){
		this->wuschel = 19; 
	}
	else if(this->rank == 22){
		this->wuschel = 18; 
	}
	else if(this->rank == 23){
		this->wuschel = 23; 
	}
	//layer 4
	else if(this->rank == 25){
		this->wuschel = 43;
	}
	else if(this->rank == 26){
		this->wuschel = 21; 
	}
	else if(this->rank == 27){
		this->wuschel = 50; 
	}
	else if(this->rank == 29){
		this->wuschel = 14; 
	}
	else if(this->rank == 30){
		this->wuschel = 11;
	}
	//layer 5
	else if(this->rank == 32){
		this->wuschel = 44;
	}
	else if(this->rank == 33){
		this->wuschel = 21;
	}
	else if(this->rank == 35){
		this->wuschel = 27;
	}
	//layer 6
	else if(this->rank == 37){
		this->wuschel = 36;
	}
	else if(this->rank == 38){
		this->wuschel = 35;
	}
	else if(this->rank == 40){
		this->wuschel = 23;
	}
	//layer 7
	else if(this->rank == 42){
		this->wuschel = 23;
	}
	else if(this->rank == 43){
		this->wuschel = 27;
	}
	else if(this->rank == 45){
		this->wuschel = 22;
	}
	else{
	//pz the rest of the cells
	//function for pz cells is:
	this->wuschel = 0;
	//this->wuschel = -0.01624*5*pow((cell_center-Coord(0,-24)).length(),2) + 1.519*5*(cell_center-Coord(0,-24)).length() + 9.253; 
	}*/
	//from 2018 paper
	double distance = (cell_center-(L1_AVG-Coord(0,21))).length();
	//if(distance < 140*.15){
		this->wuschel = 84.6*exp(-0.01573*(distance));
	//}
	//else {
	//	this->wuschel = 9.36*exp(0.01573*(-distance/.15+280));
	//}
	return;
}
void Cell::calc_CK(Coord L1_AVG) {
	double distance = (cell_center-(L1_AVG-Coord(0,21))).length();
	if((this->get_Layer() == 1)||(this->get_Layer() ==2)){
		this->cytokinin = 0;
	}
	else{// if((this->get_Layer() >2) && (this->get_Layer() < 6)){
		this->cytokinin = 110*exp(-0.01637*distance);
	}
	//else {
	//	this->cytokinin = 70;
	//}
	/*if(this->layer==1){
		this->cytokinin = 10;
	}
	else if(this->layer==2){
		this->cytokinin = 10;
	}
	else if(this->layer==3){
		this->cytokinin = 40;
	}
	else if(this->layer==4){
		this->cytokinin = 50;
	}
	else if(this->layer==5){
		this->cytokinin = 40;
	}
	else if(this->layer==6){
		this->cytokinin = 20;
	}
	else{
		this->layer = 10;
	}*/
	return;
}
void Cell::set_growth_rate() {
	//this->growth_rate = unifRandInt(5000,30000);
	if(this->wuschel < 46){
		this->growth_rate = unifRandInt(2000,10000);
	}
	else if((this->wuschel >= 46)&&(this->wuschel < 50)){
		this->growth_rate = unifRandInt(10000,12510);
	}
	else if((this->wuschel >=53.8) && (this->wuschel < 57)){
		this->growth_rate = unifRandInt(12510,15012);
	}
	else if((this->wuschel >= 57)&&(this->wuschel <60.8)){
		this->growth_rate = unifRandInt(15012,17514);
	}
	else if((this->wuschel >=60.8) &&(this->wuschel <64)){
		this->growth_rate = unifRandInt(17514,20016);
	}
	else if((this->wuschel >=64) &&(this->wuschel <67.8)){
		this->growth_rate = unifRandInt(20016,22518);
	}
	else if((this->wuschel >=67.8) &&(this->wuschel <71)){
		this->growth_rate = unifRandInt(22518,25020);
	}
	else if((this->wuschel >=71) &&(this->wuschel <74.8)){
		this->growth_rate = unifRandInt(25020,27522);
	}
	else if((this->wuschel >=74.8) &&(this->wuschel <78)){
		this->growth_rate = unifRandInt(27522,30024);
	}
	else if((this->wuschel >=78) && (this->wuschel <81.8)){
		this->growth_rate = unifRandInt(30024,35526);
	}
	else{
		this->growth_rate = unifRandInt(37530,40032);
	}
	if((this->cytokinin >= 80)){

		this->growth_rate = growth_rate*.1;
	}

	//this->growth_rate = 5000;
	//2018 paper	
	/*if(this->wuschel < 11){
		this->growth_rate = unifRandInt(5000,8000);
	}
	else if((this->wuschel >= 12) &&(this->wuschel <24)) {
		this->growth_rate = unifRandInt(8000,10000);
	}
	else if((this->wuschel >= 24) && (this->wuschel <36)){
		this->growth_rate = unifRandInt(10000,12000);
	}
	else if ((this->wuschel >= 36) && (this->wuschel <48)){
		this->growth_rate = unifRandInt(12000,14000);
	}
	else if ((this->wuschel >= 48) && (this->wuschel < 60)){
		this->growth_rate = unifRandInt(14000,16000);
	}	
	else if ((this->wuschel >= 60) && (this->wuschel <72)){
		this->growth_rate = unifRandInt(16000,17000);
	}
	else if ((this->wuschel >= 72) && (this->wuschel < 84)){
		this->growth_rate = unifRandInt(17000,18000);
	}
	else if ((this->wuschel >= 84) && (this->wuschel < 96)){
		this->growth_rate = unifRandInt(18000,20000);
	}
	else if((this->wuschel >=96)&&(this->wuschel < 108)) {
		this->growth_rate = unifRandInt(200000,22000);
	}
	else if((this->wuschel >=108)&&(this->wuschel < 120)) {
		this->growth_rate = unifRandInt(22000,24000);
	}
	else if((this->wuschel >=120)&&(this->wuschel < 132)) {
		this->growth_rate = unifRandInt(24000,27000);
	}
	else if(this->wuschel>= 132) {
		this->growth_rate = unifRandInt(27000,30000);
	}*/
	//if(this->cytokinin > 1200){
	//	this->growth_rate = unifRandInt(2000,5000);
	//}

	return;
}
void Cell::update_growth_direction(){
	//signaling stuff
	if((this->layer == 1)||(this->layer ==2)){
		this->growth_direction = Coord(1,0);
	}
	else if(this->wuschel > this->cytokinin){
                 this->growth_direction = Coord(0,0);
        }
        else{
	 	this->growth_direction = Coord(0,1);
	}

	this->update_node_parameters_for_growth_direction();
	return;
}
void Cell::update_node_parameters_for_growth_direction(){
	vector<shared_ptr<Wall_Node>> walls;
	double k_bend;
	this->get_Wall_Nodes_Vec(walls);
	for(unsigned int i = 0; i < walls.size();i++) {	
		k_bend = compute_k_bend(walls.at(i));
		walls.at(i)->set_K_BEND(k_bend);
	}
	this->update_Wall_Equi_Angles();
	return;
}
void Cell::set_growth_direction(Coord gd){
	this->growth_direction = gd;
	return;
}
void Cell::get_Neighbor_Cells(vector<shared_ptr<Cell>>& cells) {
	cells = this->neigh_cells;
	return;
}
void Cell::set_Left_Corner(shared_ptr<Wall_Node> new_left_corner) {
	this->left_Corner = new_left_corner;
	return;
}
void Cell::set_Wall_Count(int number_nodes) {
	this->num_wall_nodes = number_nodes;
	return;
}
double Cell::compute_membr_thresh(shared_ptr<Wall_Node> current) {
	//ran simulations changing equilibrium length to 
	//introduce a growth bias and it did not have a 
	//big effect 
	//could do more sensitivity analysis with this but 
	//for now all nodes will have same equilibrium length
	/*double l_thresh = 0;
        double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	curr_len = curr_vec.length();	
	growth_len = 1;
	costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	theta = acos( min( max(costheta,-1.0), 1.0) );
	if(this->growth_direction == Coord(0,1)){
		if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
			l_thresh = Membr_Equi_Len_Long;
		}
		else { 
			l_thresh = Membr_Equi_Len_Short;
		}
	}*/
	double l_thresh = Membr_Equi_Len_Short;

	return l_thresh;
}

double Cell::compute_k_lin(shared_ptr<Wall_Node> current) {
	//had the idea that changin linear springs would create a
	//growth direction bias
	//but linear springs did not have the larges effect
	//for now all spring have the same linear spring constant
	/*double k_lin = 0;
        double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	curr_len = curr_vec.length();	
	growth_len = 1;
	costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	theta = acos( min( max(costheta,-1.0), 1.0) );
	k_lin = K_LINEAR_LOOSE + K_LINEAR_STIFF*(1-pow(costheta,2));*/
	
	double k_lin = K_LINEAR_LOOSE;
	
	return k_lin;
}
double Cell::compute_k_bend(shared_ptr<Wall_Node> current) {
	//coefficient of bending spring is very important
	//nodes that are parrallel to growth direction have 
	//high bending coefficient
	//nodes that are perpendicular to growth direction have 
	//low bending coefficient
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))||(growth_direction == Coord(0,0))) {
		//fine
	}
	else{
		cout << "No growth direction assigned" << endl;
		exit(1);
	}
	double k_bend = 0;

       	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))){
		double theta = 0;
        	double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		growth_len = 1;
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "Theta: " << theta << endl;
		if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
			k_bend = K_BEND_STIFF;
		}
		else { 
			k_bend = K_BEND_LOOSE;
		}
	}
	else{
		//cout << "elsesssseses" << endl;
		k_bend = K_BEND_UNIFORM;
	}
	//if((layer == 1)||(layer == 2)){

	//	k_bend = K_BEND_UNIFORM;
	//}
	//cout << "K bend: " << k_bend << endl;
	return k_bend;
}
double Cell::compute_k_bend_div(shared_ptr<Wall_Node> current) {
	//coefficient of bending spring is very important
	//nodes that are parrallel to growth direction have 
	//high bending coefficient
	//nodes that are perpendicular to growth direction have 
	//low bending coefficient
	//cout << "k bend div" << endl;
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0)) || (growth_direction == Coord(0,0))) {
		//fine
	}
	else{
		cout << "No growth direction assigned" << endl;
		exit(1);
	}
	double k_bend = 0;
        if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))){
		double theta = 0;
        	double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		growth_len = 1;
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "Theta: " << theta << endl;
		if((theta < ANGLE_FIRST_QUAD_Div) || (theta > ANGLE_SECOND_QUAD_Div)){
			k_bend = K_BEND_STIFF;
		}
		else { 
			k_bend = K_BEND_LOOSE;
		}
	}
	else{
		k_bend = K_BEND_UNIFORM;
	}
	//cout << "K bend: " << k_bend << endl;
	return k_bend;
}

void Cell::update_Wall_Angles() {
	//cout << "wall angles" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i< walls.size();i++) {
		//cout<< "updating" <<endl;
		walls.at(i)->update_Angle();
	}
	//cout << "Success" << endl;
	return;
}
void Cell::update_Wall_Equi_Angles() {
	//cout << "equi angles" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel 
	{
        	double theta = 0;
        	double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		int counter = 0;
		double new_equi_angle = 0; 
		double circle_angle  = (this->num_wall_nodes-2)*pi/(this->num_wall_nodes);
		#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < walls.size();i++) {	
			if(this->growth_direction != Coord(0,0)){
				curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
				curr_len = curr_vec.length();	
				growth_len = this->growth_direction.length();
				costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
				theta = acos( min( max(costheta,-1.0), 1.0) );
				if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
					new_equi_angle = pi;
				}
				else{
					counter++;
					new_equi_angle = circle_angle;
	
				}
			}
			else{
				//cout << "elseeeeee" << endl;
				new_equi_angle = circle_angle;
			}
			//if((layer == 1)||(layer == 2)){
			//	new_equi_angle = circle_angle;
			//}
					
			
		//this was an idea to make the round part of the cell
		//smaller in radius but not necessary
		//if(new_equi_angle != pi) {
		//	new_equi_angle =  (counter*2-2)*pi/(counter*2);
		//}
	
			walls.at(i)->update_Equi_Angle(new_equi_angle);
		}
	}
	return;
}
void Cell::update_Wall_Equi_Angles_Div() {
	//cout << "equi angles div" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel 
	{
        	double theta = 0;
        	double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		//int counter = 0;
		double new_equi_angle = 0; 
		double circle_angle  = (this->num_wall_nodes-2)*pi/(this->num_wall_nodes);
		#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < walls.size();i++) {	
			if(this->growth_direction != Coord(0,0)){
				curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
				curr_len = curr_vec.length();	
				growth_len = this->growth_direction.length();
				costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
				theta = acos( min( max(costheta,-1.0), 1.0) );
				if((theta < ANGLE_FIRST_QUAD_Div) || (theta > ANGLE_SECOND_QUAD_Div)){
					new_equi_angle = circle_angle;
				}
				else{
					//counter++;
					new_equi_angle = circle_angle;
				}
			}
			else {
				new_equi_angle = circle_angle;
			}
		//this was an idea to make the round part of the cell
		//smaller in radius but not necessary
		//if(new_equi_angle != pi) {
		//	new_equi_angle =  (counter*2-2)*pi/(counter*2);
		//}
	
			walls.at(i)->update_Equi_Angle(new_equi_angle);
		}
	}
	return;
}
void Cell::update_Cell_Center() {
	vector<shared_ptr<Wall_Node>> walls;
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
void Cell::update_Linear_Bending_Springs(){
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//double new_damping = this->get_Damping();
        //double k_lin = 0;
	double k_bend = 0;
	//double l_thresh = 0;

        for(unsigned int i = 0; i < walls.size();i++) {	
		//walls.at(i)->set_Damping(new_damping);
		//walls.at(i)->set_membr_len(MembrEquLen);
		//k_lin = compute_k_lin(walls.at(i));
		k_bend = compute_k_bend(walls.at(i));
		//l_thresh = compute_k_bend(walls.at(i));
		//walls.at(i)->set_K_LINEAR(k_lin);
		walls.at(i)->set_K_BEND(k_bend);
		//walls.at(i)->set_membr_len(l_thresh);
	}
	update_Wall_Equi_Angles();
	return;
}
//=============================================================
//=========================================
// Keep Track of neighbor cells and Adhesion springs
//=========================================
//=============================================================
void Cell::update_Neighbor_Cells() {
	//clear prev vector of neigh cells
	//cout << "clear" << endl;
	neigh_cells.clear();
	//grab all cells from tissue
	//cout << "cleared" << endl;
	vector<shared_ptr<Cell>> all_Cells;
	my_tissue->get_Cells(all_Cells);

	// Empty variables for holding info about other cells
	double prelim_threshold = 20;
	shared_ptr<Cell> sp_this = shared_from_this();
	//cout << "made pointer to cell" << endl;
	// iterate through all cells
	#pragma omp parallel
	{
		Coord curr_Cent;
		Coord distance;
		#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < all_Cells.size(); i++) {
			shared_ptr<Cell> curr = all_Cells.at(i);
			//cout << "made pointer to current neighbor" << endl;
			if (curr != sp_this) {
				curr_Cent = curr->get_Cell_Center();
//				cout << "got center" << endl;
				// Check if cell centers are close enough together
				distance = sp_this->cell_center - curr_Cent;
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
//each cell wall node holds a vector of adhesion 
//connections and this function clears that for
//all cell wall nodes in the cell
void Cell::clear_adhesion_vectors() {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel
	{	

		#pragma omp for schedule(static,1)	
		for(unsigned int i=0; i< walls.size();i++) {
			walls.at(i)->clear_adhesion_vec();
		}
	}
	return;
}
//for each cell wall node on current cell 
//this function searches through all the cell wall nodes on neighboring cells
//if a cell wall node on a neighboring cell is within the ADHthresh
//this function updates adhesion vector which is a private
//member variable for each cell wall node on the current cell
//this function pushes the current wall node on the neighboring cell
//onto adhesion vector
void Cell::update_adhesion_springs() {
	//get wall nodes for this cell
	vector<shared_ptr<Wall_Node>> current_cell_walls;
	this->get_Wall_Nodes_Vec(current_cell_walls);
	vector<shared_ptr<Wall_Node>> nghbr_walls_total;
	vector<shared_ptr<Wall_Node>> nghbr_walls_current;
	//int counter;
	//get all neighboring cells to this cell
	vector<shared_ptr<Cell>> neighbors;
	this->get_Neighbor_Cells(neighbors);
	for(unsigned int i = 0; i < neighbors.size(); i++) {
		neighbors.at(i)->get_Wall_Nodes_Vec(nghbr_walls_current);
		nghbr_walls_total.insert(nghbr_walls_total.end(), nghbr_walls_current.begin(), nghbr_walls_current.end());
	}	
//	#pragma omp parallel 
//	{
//		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < current_cell_walls.size(); i++) {
			//counter++;
			//cout<< counter << endl;
		        //cout << "Wall node" << current_cell_walls.at(i) << endl;
			current_cell_walls.at(i)->make_connection(nghbr_walls_total);
			//cout << "connection made" << endl;
		}
//	}
	//for all cell wall nodes
	//look at adh vector, is a nodes is in the curr cell wall
	//nodes adh vector make sure the current cell wall node
	//is in that nodes adh vector
	for(unsigned int i = 0; i < current_cell_walls.size(); i++) {
		current_cell_walls.at(i)->one_to_one_check();

	}
	return;
}
//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
//calls calc_forces
void Cell::calc_New_Forces(int Ti) {
	//cout << "cyts forces" << endl;	
	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);
	
	#pragma omp parallel for
	for (unsigned int i = 0; i < cyts.size(); i++) {
		cyts.at(i)->calc_Forces(Ti);
	}
	//cout << "cyts done" << endl;
	//calc forces on wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//cout<< "walls  forces" << endl;
	#pragma omp parallel
	{
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			walls.at(i)->calc_Forces(Ti);
		}	
	}

	return;
}
//calls update location
//calls update cell center
//calls update wall angles
void Cell::update_Node_Locations() {
	//update cyt nodes
	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);
	#pragma omp parallel 
	{
		#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < cyts.size(); i++) {
			  cyts.at(i)->update_Location();
		}	
	}

	//update wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	#pragma omp parallel 
	{	
		#pragma omp for schedule(static,1)
		for(unsigned int i=0; i< walls.size();i++) {
			//cout << "update locaation" << endl;
			walls.at(i)->update_Location();
		}
	}
	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	if((this->life_length == 2000)) {
	update_Wall_Equi_Angles();
	}
	update_Wall_Angles();
	//cout << "done" << endl;
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
	//stem and boundary?	
	if((Ti%growth_rate == (growth_rate -1))){
		//cout << "cyt node added "<< endl;
		this->add_Cyt_Node();
	  	this->Cell_Progress++;
	}
	return;
}
void Cell::division_check(){
	vector<shared_ptr<Cell>> neighbor_cells;
	//cout <<"Before div progress" << Cell_Progress << endl;	
	if(this->Cell_Progress >= 30){


		cout << "dividing cell" << this->rank <<  endl;
		//orientation of division should be 
		//fed to the division  function here
		shared_ptr<Cell> new_Cell= this->division();
		cout << "division success" << endl;
		//cout << "Parent cell prog" << Cell_Progress<< endl;
		//cout << "Sister cell prog" << new_Cell->get_Cell_Progress()<< endl;
		this->my_tissue->update_Num_Cells(new_Cell);
		//setting info about new cell
		//cout << "Num cells" << this->my_tissue->get_num_cells() << endl;
		new_Cell->set_Rank(this->my_tissue->get_num_cells()-1);
		//cout << "set rank" << endl;
		cout << "Parent rank: " << this->rank << endl;
		cout << "sister rank: " << new_Cell->get_Rank() << endl;
		cout << new_Cell->get_wall_count() << endl;
		cout << new_Cell->get_cyt_count() << endl;
		cout << this->get_wall_count() << endl;
		cout << this->get_cyt_count() << endl;
		cout << "parent" << this << endl;
		cout << "Parent progress: " << this->get_Cell_Progress() << endl;
		cout << "new cell" << new_Cell << endl;
		cout << "New progress: " << new_Cell->get_Cell_Progress() << endl;

		//layer in division function		
		//damping in division function
		//boundary needs to be figured out
		//life length set to 0 in constructor
		//life length of parent cell reset in 
		//division function
		//cyt nodes in divison function
		//wall nodes in division function
		//all cell progress set to 0 in division 
		//cell center in division function
		//cyt and wus in division function
		//growth rate in div function
		//growth direction inherited in div function
		//left corner in divison function  
		//cout << "adhesion division" << endl;
		new_Cell->update_Neighbor_Cells();
		new_Cell->update_adhesion_springs();
		new_Cell->get_Neighbor_Cells(neighbor_cells);
		for(unsigned int i =0; i < neighbor_cells.size(); i++) {
			neighbor_cells.at(i)->update_Neighbor_Cells();
			neighbor_cells.at(i)->clear_adhesion_vectors();
			neighbor_cells.at(i)->update_adhesion_springs();

		}
	}
	return;
}

double Cell::calc_Area() {
	shared_ptr<Wall_Node> curr = left_Corner;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;
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
double Cell::get_curr_perimeter() {
	//measure perimeter
	vector<shared_ptr<Wall_Node>> walls; 
	this->get_Wall_Nodes_Vec(walls);
	shared_ptr<Wall_Node> current;
	Coord curr_vec;
	double curr_len;
	double curr_perimeter = 0;
	shared_ptr<Wall_Node> left_neighbor;
	//#pragma omp parallel for reduction(+:curr_perimeter)
	current = walls.at(0);
	shared_ptr<Wall_Node> start = current;
	do{
		left_neighbor = current->get_Left_Neighbor();
		curr_vec = left_neighbor->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		curr_perimeter += curr_len;
		current = left_neighbor;
	}while(current != start);

	return curr_perimeter;
}
void Cell::set_perimeter(double new_perimeter){
	this->perimeter = new_perimeter;
	return;
}
void Cell::add_wall_Node_Check(int Ti) {
	//cout << "adding a wall node" << endl;
	//#pragma omp for schedule(static,1)
	if(this->life_length < 1000){
		//do nothing
	}
	else {
	double curr_perim = this->get_curr_perimeter();
	double increase = curr_perim - this->get_perimeter();
	//cout << "curr perim " << curr_perim << endl;
	//cout << "old perim" << this->get_perimeter()<< endl;
	//cout << "increase" << increase << endl;
	this->set_perimeter(curr_perim);
	if(increase > PERIM_INCREASE){
		add_Wall_Node(Ti);
	}
	}
	return;
}
void Cell::delete_wall_Node_Check(int Ti){
	delete_Wall_Node(Ti);
	return;
}
void Cell::add_Wall_Node(int Ti) {

//find node to the right of largest spring
	shared_ptr<Cell> this_cell= shared_from_this();
	shared_ptr<Wall_Node>right = NULL;
	//vector<pair<shared_ptr<Wall_Node>,double>> nodes;
	//cout  << "Find largest length" << endl;
	find_Largest_Length(right);
	//cout << "Largest found" << endl;
	shared_ptr<Wall_Node> left;
	Coord location;
	double l_thresh;
	double k_lin;
	double k_bend;
	//if(nodes.size() >0) {
	if(right != NULL){
//find location and set neighbors for new node
		//for(int i = 0; i< nodes.size(); i++){
		//right = nodes[i].first;
		//cout << "adding node" << endl;
		left = right->get_Left_Neighbor();
		location  = (right->get_Location() + left->get_Location())*0.5;
		shared_ptr<Wall_Node> added_node = make_shared<Wall_Node>(location, this_cell, left, right);
		this->add_wall_node_vec(added_node);
		double new_damping = this->get_Damping();
		right->set_Left_Neighbor(added_node);
		left->set_Right_Neighbor(added_node);
		//set the variables for the new node
		l_thresh = compute_membr_thresh(added_node);
		k_lin = compute_k_lin(added_node);
		k_bend = compute_k_bend(added_node);
		added_node->set_Damping(new_damping);
		added_node->set_K_LINEAR(k_lin);
		added_node->set_K_BEND(k_bend);
		added_node->set_membr_len(l_thresh);	
		added_node->set_added(1);
		//adhesion for the new node
		vector<shared_ptr<Cell>> neighbors;
		this->get_Neighbor_Cells(neighbors);
		//cout << "adh added node find closest" << endl;
		vector<shared_ptr<Wall_Node>> nghbr_walls_total;
		vector<shared_ptr<Wall_Node>> nghbr_walls_current;
		for(unsigned int i = 0; i < neighbors.size(); i++) {
			neighbors.at(i)->get_Wall_Nodes_Vec(nghbr_walls_current);
			nghbr_walls_total.insert(nghbr_walls_total.end(), nghbr_walls_current.begin(), nghbr_walls_current.end());
		}
		added_node->make_connection(nghbr_walls_total);
		//cout << "adh added node success" << endl;
		//update angles
		//should i update k_bennnnndddddd?????
		//update_Linear_Bending_Springs();
		update_Wall_Equi_Angles();
		update_Wall_Angles();
		//}
	}
	return;
}
void Cell::delete_Wall_Node(int Ti) {
	shared_ptr<Wall_Node> left = NULL;
	shared_ptr<Wall_Node> right = NULL;
	shared_ptr<Wall_Node> small = NULL;
	//vector<Cell*>neighbors;

	this->find_Smallest_Length(small);
	if(small !=NULL) {
		cout << "delete initiated" << endl;
		left = small->get_Left_Neighbor();
		right = small->get_Right_Neighbor();
	
		//if small is the left corner cell reassign
		if(this->left_Corner == small) {
			//cout << " set left corner" << endl;
			this->set_Left_Corner(left);
		}
		//need to make sure all nodes connected to small
		//via adhesion are erased
		//small->remove_from_adh_vecs();
		//small->clear_adh_vec();

		//set new neighbors so nothing points at small
		left->set_Right_Neighbor(right);
		right->set_Left_Neighbor(left);
		//reindex
		this->wall_nodes.clear();
		this->num_wall_nodes = 0;
	
		shared_ptr<Wall_Node> curr = this->left_Corner;
		shared_ptr<Wall_Node> next = NULL;
		shared_ptr<Wall_Node> orig = curr;
	
		do {
			this->wall_nodes.push_back(curr);
			next = curr->get_Left_Neighbor();
			num_wall_nodes++;
			curr = next;
		}while(next != orig);
		//cout << "update equi angles" << endl;
	
		update_Wall_Equi_Angles();
		
		//cout << "update angles" << endl;
		update_Wall_Angles();
	}
	return;
}
//finds right neighbor node of smallest length on membrane
void Cell::find_Smallest_Length(shared_ptr<Wall_Node>& right) {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 100;
	//#pragma omp parallel
	//{
		shared_ptr<Wall_Node> left_neighbor;
		double curr_len = 0;
		//#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < walls.size();i++) {
			left_neighbor = walls.at(i)->get_Left_Neighbor();
			curr_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
			if(curr_len < .05){
				if(curr_len < max_len) {
					//#pragma omp critical
					max_len = curr_len;
					right = walls.at(i);
				}
			}
		}
	//}
	return;
}
//finds right neighbor node of largest length on membrane
void Cell::find_Largest_Length(shared_ptr<Wall_Node>& node) {
	vector<shared_ptr<Wall_Node>> walls; 
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 0;
	//double second_max_len = 0;
	shared_ptr<Wall_Node> biggest;
	//shared_ptr<Wall_Node> second_biggest;
	//double temp_len;
	//shared_ptr<Wall_Node> temp_pointer;
	// UNUSED double theta = 0;
        // UNUSED double costheta = 0;
	double curr_len = 0;
	// UNUSED double growth_len = 0;
	Coord curr_vec;	
		
	//#pragma omp parallel 
	//{

	int start = unifRandInt(0,num_wall_nodes-1); 
	shared_ptr<Wall_Node> starter = walls.at(start);
	shared_ptr<Wall_Node> left_neighbor;
	shared_ptr<Wall_Node> current = starter;
	//double average_length;
	//#pragma omp for schedule(static,1)
	do {
		left_neighbor = current->get_Left_Neighbor();
		curr_vec = left_neighbor->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		// UNUSED growth_len = 1;
		// UNUSED costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		// UNUSED theta = acos( min( max(costheta,-1.0), 1.0) );
		if(curr_len > max_len){
			max_len = curr_len;
			biggest = current;
		}
		//other stuff for vector of nodes etc
		//average_length = average_length + curr_len;
		//if((theta < ADD_WALL_NODE_ANGLE_FIRST_QUAD) ||(theta > ADD_WALL_NODE_ANGLE_SECOND_QUAD)){
			//if(curr_len > MEMBR_THRESH_LENGTH) {			
				//cout << "got in here" << endl;
				//if(curr_len > second_max_len){
				//	if(curr_len > max_len){
						//temp_len = max_len;
						//temp_pointer = biggest;
						//max_len = curr_len;
						//biggest = current;
						//second_max_len = temp_len;
						//second_biggest = temp_pointer;
					//}
					//else {
						//second_max_len = curr_len;
						//second_biggest = current;
					//}
					//cout << "max_len" << max_len << endl;
					//cout << "second_max_len"<< second_max_len << endl;
				//}
			//}	

		//}
		current = left_neighbor;
	}while (left_neighbor != starter);
	//cout << "have max lens" << endl;
	//nodes.push_back(make_pair(biggest,max_len));
	//nodes.push_back(make_pair(second_biggest, second_max_len));
	//cout << "average length" << average_length/num_wall_nodes<< endl;
	node = biggest;

return;

}
Coord Cell::compute_direction_of_highest_tensile_stress(){
	//average position of all cell wall nodes
	vector<shared_ptr<Wall_Node>> wall_nodes;
	this->get_Wall_Nodes_Vec(wall_nodes);
	Coord next_coord;
	Coord curr_coord;
	Coord direction_vec;
	shared_ptr<Wall_Node> curr = wall_nodes.at(0);
	shared_ptr<Wall_Node> orig = curr;
	shared_ptr<Wall_Node> next;
	double delta_x = 0;
	double delta_y = 0;
	double x = 0;
	double y = 0;
	double average_x;
	double average_y;
	int counter = 0;
	double curr_length;
	double strain;
	do{
		next = curr->get_Left_Neighbor();
		curr_coord = curr->get_Location();
		next_coord = next->get_Location();
		curr_length = (next->get_Location() - curr->get_Location()).length();
		delta_x = (next_coord.get_X() - curr_coord.get_X())/curr_length;
		delta_y = (next_coord.get_Y() - curr_coord.get_Y())/curr_length;
		if(delta_x < 0) {
			delta_x = delta_x*-1;
		}
		if(delta_y < 0) {
			delta_y = delta_y*-1;
		}
		
		strain = (curr_length - Membr_Equi_Len_Long)/Membr_Equi_Len_Long; 
		cout << "strain" << strain << endl;
		x = x+ strain*delta_x;
		y = y+ strain*delta_y;
		counter++;
		curr = next;
	} while(next != orig);
	average_x = x/counter;
	average_y = y/counter;
	
	direction_vec = Coord(-average_y,average_x);
	return direction_vec;
}
/*Coord Cell::compute_point_on_line(double t){
	Coord r_0 = this->get_Cell_Center();
	Coord v = this->compute_longest_axis();
	Coord p = Coord(r_0.get_X() + v.get_X()*t, r_0.get_Y() + v.get_Y()*t);
	return p;
}*/
void Cell::add_Cyt_Node() {
	//cout << "cyt" << endl;
	double new_damping = this->get_Damping();
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Cyt_Node> cyt = make_shared<Cyt_Node>(cell_center, this_cell);
	cyt_nodes.push_back(cyt);
	cyt->set_Damping(new_damping);

	num_cyt_nodes++;
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
	//cout << "ID before: " << id << endl;
	int rel_cnt = 0;

	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do { 
		curr_wall->update_VTK_Id(id);
		id++;
		for(unsigned int i = 0; i < curr_wall->get_adh_vec().size(); i++){
			rel_cnt++;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);
	
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->update_VTK_Id(id);
		id++;
	}
	//cout << "ID after: " << id << endl;
	return rel_cnt;
}
void Cell::print_VTK_Adh(ofstream& ofs) {

	int my_id, nei_id;
	shared_ptr<Wall_Node> neighbor = NULL;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	vector<shared_ptr<Wall_Node>> nodes;
	
	do {
		for(unsigned int i = 0; i < curr_wall->get_adh_vec().size(); i++){
			nodes = curr_wall->get_adh_vec();
			neighbor = nodes.at(i);
			if(neighbor != NULL) {
				my_id = curr_wall->get_VTK_Id();
				nei_id = neighbor->get_VTK_Id();
				ofs << 2 << ' ' << my_id << ' ' << nei_id << endl;
			}
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while(curr_wall != left_Corner);
	return;
}
Coord Cell::average_coordinates(){
	Coord direction = Coord(0,0);
	for(unsigned int i = 0; i <wall_nodes.size() ;i++){
		direction = direction + wall_nodes.at(i)->get_Location();
	}
	double num_nodes = (double) wall_nodes.size();
	direction = direction/num_nodes;
	return direction;
}
void Cell::print_direction_vec(ofstream& ofs){
	//average coordinates of all nodes
	Coord sum = this->average_coordinates();
	//point 1 is center + 1*direction vector
	Coord point1 = this->cell_center + sum*.05;
	Coord point2 = this->cell_center - sum*.05;
	ofs << point1.get_X() << ' ' << point1.get_Y() << ' ' << 0 << endl;
	ofs << point2.get_X() << ' ' << point2.get_Y() << ' ' << 0 << endl;
	
	return;
}
void Cell::print_locations(ofstream& ofs) {
	ofs << this->get_Rank() << endl;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	shared_ptr<Wall_Node> orig = curr_wall;
	//	cout << "knows left corner" << endl;
	do {


		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 <<' '<< 1 << endl;
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();
	
	//cout << "did it  " << count << endl;
	} while (curr_wall != orig);
	
	//cout << "walls worked" << endl;
	/*for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
	
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << 0 << endl;
	}*/
	return;
}

void Cell::print_VTK_Points(ofstream& ofs, int& count) {

	shared_ptr<Wall_Node> curr_wall = left_Corner;
	shared_ptr<Wall_Node> orig = curr_wall;
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
	};
	
	//cout << "points worked" << endl;
	return;
}
/*void Cell::print_VTK_Scalars_Wall_Pressure(ofstream& ofs){
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
}*/
void Cell::print_VTK_Scalars_Average_Pressure(ofstream& ofs) {
	//float pressure = this->average_Pressure();
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		//concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
	//	ofs << pressure << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	//for(unsigned int i=0;i<wall_nodes.size();i++){
	//	ofs << pressure << endl;
	//}
	for(unsigned int i=0;i<cyt_nodes.size();i++){
	//	ofs<< pressure << endl;
	}
	return;
}
/*void Cell::print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs) {
	double pressure = this->average_Pressure();
	for(unsigned int i=0;i<1;i++){
		ofs << pressure << endl;
	}
	return;
}*/
	
void Cell::print_VTK_Scalars_WUS(ofstream& ofs) {

	double concentration = 0;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;
	}
	return;
}
void Cell::print_VTK_Scalars_CK(ofstream& ofs) {

	double concentration = 0;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;
	}
	return;
}
/*void Cell::print_VTK_Scalars_WUS_cell(ofstream& ofs) {
	for(unsigned int i = 0; i < 1; i++) {
		
		double concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();
		
		ofs << concentration << endl;
	}
	return;
}
	
void Cell::print_VTK_Scalars_CYT(ofstream& ofs) {

//	double concentration = 0;

	shared_ptr<Wall_Node> curr_wall = left_Corner;
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
}*/
void Cell::print_VTK_Scalars_Node(ofstream& ofs) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		if(currW->get_added()==1){
			color = 30.0;
		} else {
			color = 0.0;
		}
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while (currW != left_Corner);
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		color = 0.0;
		ofs << color << endl;
	}
	return;
}

void Cell::print_VTK_Tensile_Stress(ofstream& ofs) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		color = currW->calc_Tensile_Stress();
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		color = CYT_COLOR;
		ofs << color << endl;
	}
	return;
}

void Cell::print_VTK_Shear_Stress(ofstream& ofs) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		color = currW->calc_Shear_Stress();
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		color = CYT_COLOR;
		ofs << color << endl;
	}
	return;
}
/*		   
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
}*/





//////////////////////////////////
