//cell_div.cpp
//========================
//Forward Declarations
//


//========================
//Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//==========================

Cell* Cell::divide() {
	Cell* sister = NULL;
	//calculate area
//	cout << "in division function" << endl;
//	double area = this->calc_Area();
//	cout << "area calculated " << area << endl;
//	if(area > AREA_DOUBLED) {
		//add in iff statement for amount of cytokinin and wuschel and logintudinal vs. radial pressure
	//	if(this->rank != 1) {
			cout << "Cell " << this->rank << "  passed area threshold for division lengthwise" << endl;
			sister = this->division();
			cout << "divided" << endl;
	//	}
	//	else if(layer == 3) {
	//		cout << "Cell " << this->rank << " passed area threshold for division widthwise" << endl;
	//		sister = this->divide_width_wise();
	//		cout << "divided" << endl;
	//	}
	return sister;
}

void Cell::find_Largest_Length_Div(Wall_Node*& right_one, Wall_Node*& right_two) {
	Wall_Node* curr = left_Corner;
	Wall_Node* biggest = NULL;
	Wall_Node* second_biggest = NULL;
	Wall_Node* orig = curr;
	Coord left_Neighb_loc;
	Coord curr_Loc;
	Coord diff_vect;
	double max_len = 0;
	double second_len = 0;
	double len;
	double counter = 0;
//	double temp_len;
//	Wall_Node* temp_node;
	//int big_gaps = 0;
	//loop through all possible Cell Wall 'links' to find biggest
	do {
	//	cout << "finding current lengths and comparing" << endl;
		left_Neighb_loc = curr->get_Left_Neighbor()->get_Location();
		curr_Loc = curr->get_Location();
		diff_vect = left_Neighb_loc - curr_Loc;
		len = diff_vect.length();
		//if (len > MEMBR_THRESH_LENGTH) {
		//	big_gaps++;
		if(len > max_len) {
				
//				second_len = max_len;
				max_len = len;
//				right_two = right_one;
				right_one = curr;
		}
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);
	
	Wall_Node* new_start = right_one;
	Wall_Node* new_end = right_one;
	for(unsigned int i = 0; i < 80; i++) {
			new_start = new_start->get_Left_Neighbor();
			new_end = new_end->get_Right_Neighbor();
	}
	max_len = 0;
	curr = new_start;
	do {
	//	cout << "finding current lengths and comparing" << endl;
		left_Neighb_loc = curr->get_Left_Neighbor()->get_Location();
		curr_Loc = curr->get_Location();
		diff_vect = left_Neighb_loc - curr_Loc;
		len = diff_vect.length();
		//if (len > MEMBR_THRESH_LENGTH) {
		//	big_gaps++;
		if(len > max_len) {
				
//				second_len = max_len;
				max_len = len;
//				right_two = right_one;
				right_two = curr;
		}
		curr = curr->get_Left_Neighbor();
	} while (curr != new_end);
	//cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;
	return;
}

Cell* Cell::division() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
	Cell* sister = new Cell(my_tissue);
	double center_x = this->cell_center.get_X();
	double center_y = this->cell_center.get_Y();

	cout << "Made new cell pointer" << endl;
	
	Wall_Node* closest = NULL;
	Wall_Node* next_closest = NULL;
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	

	Wall_Node* start = NULL;
	Wall_Node* end = NULL;
	
	cout << "Find gaps " << endl;
	find_Largest_Length_Div(start,end);
	double counter = 0;
	curr = start;
	do { 
		next = curr->get_Left_Neighbor();
		counter++;
		curr = next;
	} while(next!=end);
	cout << "counter: " << counter << endl;

	if(start == NULL) {
		cout << "start bad" << endl;
	}
	if(end == NULL) {
		cout << end << endl;
	}
	cout << "Delete for goodness" << endl;
	Wall_Node* left_start = start->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node*  right_start  = start->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	
	delete start->get_Right_Neighbor()->get_Right_Neighbor();
	delete start->get_Right_Neighbor();
	delete start->get_Left_Neighbor()->get_Left_Neighbor();
	delete start->get_Left_Neighbor();
	delete start;

	cout << "deleted" << endl;

	Wall_Node* left_end = end->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	cout << "deleted 2" << endl;
	Wall_Node* right_end = end->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();

	delete end->get_Left_Neighbor()->get_Left_Neighbor();
	delete end->get_Right_Neighbor()->get_Right_Neighbor();
	delete end->get_Left_Neighbor();
	delete end->get_Right_Neighbor();
	delete end;

	cout << "Deleted" << endl;
	cout << "lengths" << endl;
	double left_length = (left_end->get_Location() - left_start->get_Location()).length();
	double right_length = (right_end->get_Location() - right_start->get_Location()).length();

	double total_num_left_dbl = (left_length/(MembrEquLen*4));
	int total_num_left = (int) total_num_left_dbl -1;
	double total_num_right_dbl = (right_length/(MembrEquLen*4));
	int total_num_right = (int) total_num_right_dbl -1;
	cout << "Total num right" << total_num_right << endl;	
	cout << "Total num left" << total_num_left << endl;	
	
	//bool left = false;
	//bool right = false;

	double x_left = left_end->get_Location().get_X() - left_start->get_Location().get_X();
	double x_right = right_end->get_Location().get_X() - right_start->get_Location().get_X();
	//cout << x_left << endl;
	double y_left = left_end->get_Location().get_Y() - left_start->get_Location().get_Y();
	double y_right = right_end->get_Location().get_Y() - right_start->get_Location().get_Y();
	 
	double delta_x_left = (double)x_left/(double)total_num_left;
	double delta_x_right = (double)x_right/(double)total_num_right;
	//cout <<delta_x_left<< endl;
	double delta_y_left = (double)y_left/(double)total_num_left;
	double delta_y_right = (double)y_right/(double)total_num_left;

	double start_x_left = left_start->get_Location().get_X();
	double start_x_right = right_start->get_Location().get_X();
	double start_y_left = left_start->get_Location().get_Y();
	double start_y_right = right_start->get_Location().get_Y();

//	cout << "make left side" << endl;
	curr = NULL;
	double curr_x = start_x_left;
	double curr_y = start_y_left;
	Coord curr_coord = Coord(curr_x+delta_x_left,curr_y+delta_y_left);
	Wall_Node* prev = left_start;
	
	for(unsigned int i = 0; i < total_num_left-2; i++) {
		curr_coord = curr_coord + Coord(delta_x_left, delta_y_left);
		cout << curr_coord << endl;
		curr = new Wall_Node(curr_coord, this);
//		cout << "Setting neighbors: " << i << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(curr);
	left_Corner = left_start;

//	cout << "make right side" << endl;;
	prev = right_start;
	curr_y = start_y_right;
	curr_x = start_x_right;
	curr_coord = Coord(curr_x+delta_x_right, curr_y+delta_y_right);
//	cout << total_num_right << endl;
	for(unsigned int j = 0; j< total_num_right-2; j++) {
		curr_coord = curr_coord + Coord(delta_x_right,delta_y_right);;
		cout << curr_coord << endl;
		curr = new Wall_Node(curr_coord, sister);
//		cout << "Setting neighbors: " << j << endl;
		prev->set_Right_Neighbor(curr);
		//cout << "prev" << endl;
		curr->set_Left_Neighbor(prev);
		//cout<< "curr" << endl;
		prev = curr;
		//cout << "prev = curr" << endl;
	}
	curr->set_Right_Neighbor(right_end);
	right_end->set_Left_Neighbor(curr);
	sister->set_Left_Corner(right_start);

	//count wall nodes
//	cout << "begin count wall nodes" << endl;
	curr = this->left_Corner;
	next = NULL;
	orig = curr;
	int number_nodes_A = 0;
//	cout << "cell a counting" << endl;
	do {
		number_nodes_A++;
		curr->update_Cell(this);
		next = curr->get_Left_Neighbor();
	//	cout << number_nodes_A << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while(next != orig);
//	cout << "done cell a counting" << endl;	
	this->set_Wall_Count(number_nodes_A);
	cout << "Parent nodes: " << number_nodes_A << endl;	
	curr = sister->get_Left_Corner();
	next = NULL;
	orig = curr;
	int number_nodes_B = 0;
//	cout << "Cell b counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		next = curr->get_Left_Neighbor();
	//	cout << number_nodes_B << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while (next != orig);
//	cout << "done counting cell b" << endl;
	sister->set_Wall_Count(number_nodes_B);
	cout << "Sister: " << number_nodes_B << endl;
//	cout << "updating angles" << endl;
	//update wall angles
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();
//	cout << "updating center" << endl;
	//update cell center
	this->update_Cell_Center();
	sister->update_Cell_Center();
	this->calc_WUS();
	sister->calc_WUS();
	this->calc_CYT();
	sister->calc_CYT();
	this->calc_Total_Signal();
	sister->calc_Total_Signal();
	sister->set_Layer(this->layer);

	double K_LINEAR = -3.3673*(this->cytokinin) + 5.7335*(this->wuschel) + 269.4673;
	double K_LINEAR_X;
	double K_LINEAR_Y;
		
	if(this->layer == 1) {
		K_LINEAR_Y = 600;
		K_LINEAR_X = 150;
	}	
	else {
		K_LINEAR_X = 600;
		K_LINEAR_Y = 150;
	}	

	this->K_LINEAR = Coord(K_LINEAR_X, K_LINEAR_Y);
	
	sister->set_K_LINEAR(K_LINEAR_X,K_LINEAR_Y);
	double new_damping = this->get_Damping();
	sister->set_Damping(new_damping);

	//delete all old cyt nodes
	int num_cyts = cyt_nodes.size();
	int new_cyt_cnt = num_cyts/2;
	Cyt_Node* cyt = NULL;
	while(!cyt_nodes.empty()) {
		cyt = cyt_nodes.at(cyt_nodes.size() -1);
		delete cyt;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
//	cout << "Finished deleting old cyt nodes" << endl;
	
//	cout << "get most up/down and left/right for radius" << endl;
/*	Wall_Node* up_sis = NULL;
	Wall_Node* down_sis = NULL;
	Wall_Node* left_sis = NULL;
	Wall_Node* right_sis = NULL;
	this->closest_node_top(up);
	this->closest_node_bottom(down);
	if(up == NULL) {

	cout << "this up null" << endl;
	}
	if(down == NULL) {
		cout << "this down null" << endl;
	}
	sister->closest_node_top(up_sis);
	sister->closest_node_bottom(down_sis);
	if(up_sis == NULL) {
		cout << "up sis null" << endl;
	}
	if(down_sis == NULL) {
		cout << "down sis null" << endl;
	}
	this->closest_node_left(left);
	this->closest_node_right(right);
	if(left == NULL) {
		cout << "left null" << endl;
	}
	if(right==NULL) {
		cout << "right null" << endl;
	}
	sister->closest_node_left(left_sis);
	sister->closest_node_right(right_sis);
	if(left_sis == NULL) {

	cout << "left sis null" << endl;
	}
	if(right_sis == NULL) {
		cout << "right sis null" << endl;
	}


	double radius_x = ((left_end->get_Location() - right->get_Location()).length())*0.5;
	double radius_y = ((up->get_Location() - down->get_Location()).length())*0.5;
//	cout << "first radius" << endl;
	double radius_x_s = ((left_sis->get_Location() - right_sis->get_Location()).length())*0.5;
	double radius_y_s = ((up_sis->get_Location() - down_sis->get_Location()).length())*0.5;
//	cout << "second radius" << endl;
	
*/
	//old cell get radius
	//sister get radius
	double old_cell_radius = this->find_radius();
	double sister_cell_radius = sister->find_radius();
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		//cout << "adding cytoplasm" << i << endl;
		this->add_Cyt_Node_Div(old_cell_radius);
		sister->add_Cyt_Node_Div(sister_cell_radius);
	}
	return sister;
}

/*Cell* Cell::divide_width_wise() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
	bool islength = false;
//	cout << "Made new cell pointer" << endl;
	Cell* sister = new Cell(my_tissue);
	Wall_Node* up = NULL;
	Wall_Node* down = NULL;
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	this->closest_node_left(left);
	if(up == NULL) {
		cout << "Top NULL" << endl;
		exit(1);
	}
	this->closest_node_right(right);*	if(down == NULL) {
		cout << "Bottom NULL" << endl;
		exit(1);
	}

	Wall_Node* top_start = left->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node* bottom_start = left->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	
	delete left->get_Right_Neighbor()->get_Right_Neighbor();
	delete left->get_Right_Neighbor();
	delete left->get_Left_Neighbor()->get_Left_Neighbor();
	delete left->get_Left_Neighbor();
	delete left;

	Wall_Node* top_end = right->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	Wall_Node* bottom_end = right->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();

	delete right->get_Left_Neighbor()->get_Left_Neighbor();
	delete right->get_Right_Neighbor()->get_Right_Neighbor();
	delete right->get_Left_Neighbor();
	delete right->get_Right_Neighbor();
	delete right;

	double top_slope = (top_end->get_Location().get_Y() - top_start->get_Location().get_Y())/(top_end->get_Location().get_X() - top_start->get_Location().get_X());
	double bottom_slope = (bottom_end->get_Location().get_Y() - bottom_start->get_Location().get_Y())/(bottom_end->get_Location().get_X() - bottom_start->get_Location().get_X());
	
	double b_top = top_end->get_Location().get_Y()-top_slope*top_end->get_Location().get_X();
	double b_bottom = bottom_end->get_Location().get_Y()-bottom_slope*bottom_end->get_Location().get_X();

	double top_length = (top_end->get_Location() - top_start->get_Location()).length();
	double bottom_length = (bottom_end->get_Location() - bottom_start->get_Location()).length();

	int total_num_top = (int) (top_length/(MembrEquLen*2));
	int total_num_bottom = (int) (bottom_length/(MembrEquLen*2));

	double x_top = sqrt(pow(top_end->get_Location().get_X() - top_start->get_Location().get_X(),2));
	double x_bottom = sqrt(pow(bottom_end->get_Location().get_X() - bottom_start->get_Location().get_X(),2));

	double delta_x_top = x_top/total_num_top;
	double delta_x_bottom = x_bottom/total_num_bottom;

	double start_x_top = top_start->get_Location().get_X();
	double start_x_bottom = bottom_start->get_Location().get_X();
	
//	cout << "make left side" << endl;
	Wall_Node* curr = NULL;
	Coord curr_coord;
	Wall_Node* prev = top_start;
	double curr_x = start_x_top;
	for(unsigned int i = 0; i < total_num_top; i++) {
		curr_x = curr_x + delta_x_top;
		curr_coord = Coord(curr_x + 0.04, top_slope*curr_x + b_top);
		curr = new Wall_Node(curr_coord, this);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(top_end);
	top_end->set_Right_Neighbor(curr);
	left_Corner = top_start;

//	cout << "make right side" << endl;;
	prev = bottom_start;
	curr_x = start_x_bottom;

	for(int i = 0; i< total_num_bottom; i++) {
		curr_x = curr_x + delta_x_bottom;
		curr_coord = Coord(curr_x - .04, bottom_slope*curr_x + b_bottom);
		curr = new Wall_Node(curr_coord, sister);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Right_Neighbor(curr);
		curr->set_Left_Neighbor(prev);
		prev = curr;
	}
	curr->set_Right_Neighbor(bottom_end);
	bottom_end->set_Left_Neighbor(curr);
	sister->set_Left_Corner(bottom_start);

	//count wall nodes
	//cout << "begin count wall nodes" << endl;
	curr = this->left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	int number_nodes_A = 0;
	//cout << "cell a counting" << endl;
	do {
		number_nodes_A++;
		curr->update_Cell(this);
		next = curr->get_Left_Neighbor();
		//cout << number_nodes_A << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while(next != orig);
	//cout << "done cell a counting" << endl;	
	this->set_Wall_Count(number_nodes_A);
	
	curr = sister->get_Left_Corner();
	next = NULL;
	orig = curr;
	int number_nodes_B = 0;
	//cout << "Cell b counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		next = curr->get_Left_Neighbor();
		//cout << number_nodes_B << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while (next != orig);
	//cout << "done counting cell b" << endl;
	sister->set_Wall_Count(number_nodes_B);
	
	//cout << "updating angles" << endl;
	//update wall angles
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();
	//cout << "updating center" << endl;
	//update cell center
	this->update_Cell_Center();
	sister->update_Cell_Center();
	this->calc_WUS();
	sister->calc_WUS();
	this->calc_CYT();
	sister->calc_CYT();
	double K_LINEAR_Y = .1540*pow(wuschel,3) + -4.8350*pow(wuschel,2) + 54.2901*wuschel + -50.7651;
	double K_LINEAR_X = -13.2177*wuschel + 473.7440;
	this->K_LINEAR = Coord(K_LINEAR_X, K_LINEAR_Y);
	//cout << "update layer" << endl;
	//update layer information
	sister->set_Layer(this->layer);
	//cout << "update growth rate" << endl;
//	sister->set_growth_rate(this->growth_rate);
	
//	this->area = calc_Area();
//	double new_area = sister->calc_Area();
//	sister->set_Area(new_area);	
	
//	this->reset_Cell_Progress();
//	sister->reset_Cell_Progress();
//	this->Cell_Progress_add_node = area;
//	sister->set_Cell_Progress_add_node(new_area);
	//cout << "Updated Angles and cell centers"<< endl;
	
	double new_damping = this->get_Damping();
	sister->set_Damping(new_damping);
	//distribute cyt nodes between sister cells
	//cout << "deleting cyt nodes" << endl;
	
	//delete all old cyt nodes
	int num_cyts = cyt_nodes.size();
	int new_cyt_cnt = num_cyts/2;
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
	//cout << "Finished deleting old cyt nodes" << endl;
	
//	cout << "get most up/down and left/right for radius" << endl;
	Wall_Node* up_sis = NULL;
	Wall_Node* down_sis = NULL;
	Wall_Node* left_sis = NULL;
	Wall_Node* right_sis = NULL;
	this->closest_node_top(up);
	this->closest_node_bottom(down);
	sister->closest_node_top(up_sis);
	sister->closest_node_bottom(down_sis);
	this->closest_node_left(left);
	this->closest_node_right(right);
	sister->closest_node_left(left_sis);
	sister->closest_node_right(right_sis);

	double radius_x = ((left->get_Location() - right->get_Location()).length())*0.5;
	double radius_y = ((up->get_Location() - down->get_Location()).length())*0.5;
	double radius_x_s = ((left_sis->get_Location() - right_sis->get_Location()).length())*0.5;
	double radius_y_s = ((up_sis->get_Location() - down_sis->get_Location()).length())*0.5;
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div(radius_x, radius_y, islength);
		sister->add_Cyt_Node_Div(radius_x_s, radius_y_s, islength);
	}
	return sister;
}
*/
