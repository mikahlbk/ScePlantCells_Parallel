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
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//==========================
//find the nodes to set up the division plane
shared_ptr<Wall_Node> Cell::find_nodes_for_div_plane(Coord& orientation) {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	double curr_theta = 0;
        double curr_diff = 0;
	double target_theta = 1.57;
	double smallest_diff = 100;
	double costheta = 0;
	double curr_len = 0;
	double orientation_len = 0;
	Coord curr_vec;	
	shared_ptr<Wall_Node> node;

	for(unsigned int i = 0; i< walls.size(); i++) {
		curr_vec = (walls.at(i)->get_Left_Neighbor()->get_Location()-walls.at(i)->get_Location());
		curr_len = curr_vec.length();
		orientation_len = 1;
		costheta = orientation.dot(curr_vec)/(curr_len*orientation_len);
		curr_theta = acos( min( max(costheta,-1.0), 1.0) );
	 	curr_diff = abs(curr_theta - target_theta);
		if(curr_diff < smallest_diff){
			node = walls.at(i);
			smallest_diff = curr_diff;
		} 
	}
	if(node == NULL){
		cout << "Did not pick up div plane nodes" << endl;
		exit(1);
	}
	return node;
}

shared_ptr<Cell> Cell::division() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as parent cell
	//	this function will create a sister cell
	//	and return it to the tissue
	//cout << "New Cell" << endl;
	shared_ptr<Cell> sister = make_shared<Cell> (this->my_tissue);
	sister->set_Layer(this->layer);
	sister->set_growth_direction(this->get_growth_direction());
	
	//this is the vector for the desired division orientation
	Coord orientation = Coord(0,1);
	//finds node on one side of cell
	shared_ptr<Wall_Node> first = find_nodes_for_div_plane(orientation);
	//give opposite vector to find node on other side of the cell
	orientation = orientation*-1;
	shared_ptr<Wall_Node> second = find_nodes_for_div_plane(orientation);
	
	//***Warning***make sure the nodes you are getting rid of dont have 
	//any other nodes pointing to them
	//left and right neighbors will be reset below
	//other pointers come from adhesion!!!!
	//cout << "adhesion pointers" << endl;
	first->remove_from_adh_vec();
	first->clear_closest_in_adh_vec();
	second->remove_from_adh_vec();
	second->clear_closest_in_adh_vec();
	//cout << "adhesion success" << endl;
	
	//set the starting point of daughter cell one
	shared_ptr<Wall_Node> start_daughter_one = first->get_Left_Neighbor();
	//move over one so that cells have enough room between them
	//but first need to get rid of pointers in adhesion	
	start_daughter_one->remove_from_adh_vec();
	start_daughter_one->clear_closest_in_adh_vec();	
	start_daughter_one = start_daughter_one->get_Left_Neighbor();
	//move over one so that cells have enough room between them
	//but first need to get rid of pointers in adhesion	
	start_daughter_one->remove_from_adh_vec();
	start_daughter_one->clear_closest_in_adh_vec();	
	//move happens here
	start_daughter_one = start_daughter_one->get_Left_Neighbor();
	
	//set ending point for daughter cell one
	shared_ptr<Wall_Node>end_daughter_one = second->get_Right_Neighbor();
	//move over one so that cells have enough room between them
	//but first need to get rid of pointers in adhesion
	end_daughter_one->remove_from_adh_vec();
	end_daughter_one->clear_closest_in_adh_vec();	
	end_daughter_one = end_daughter_one ->get_Right_Neighbor();
	//move over one so that cells have enough room between them
	//but first get rid of adhesion pointers
	end_daughter_one->remove_from_adh_vec();
	end_daughter_one->clear_closest_in_adh_vec();	
	//move happens here
	end_daughter_one = end_daughter_one ->get_Right_Neighbor();
	
	//set starting point for daughter cell two	
	shared_ptr<Wall_Node>start_daughter_two = second->get_Left_Neighbor();
	start_daughter_two->remove_from_adh_vec();
	start_daughter_two->clear_closest_in_adh_vec();	
	start_daughter_two = start_daughter_two->get_Left_Neighbor();
	start_daughter_two->remove_from_adh_vec();
	start_daughter_two->clear_closest_in_adh_vec();	
	start_daughter_two = start_daughter_two->get_Left_Neighbor();
	
	//set ending point for daughter cell two
	shared_ptr<Wall_Node>end_daughter_two = first->get_Right_Neighbor();
	end_daughter_two->remove_from_adh_vec();
	end_daughter_two->clear_closest_in_adh_vec();	
	end_daughter_two = end_daughter_two ->get_Right_Neighbor();
	end_daughter_two->remove_from_adh_vec();
	end_daughter_two->clear_closest_in_adh_vec();	
	end_daughter_two = end_daughter_two ->get_Right_Neighbor();


	//this is the length of the division plane
	double daughter_one_length = (start_daughter_one->get_Location() - end_daughter_one->get_Location()).length();
	double daughter_two_length = (start_daughter_two->get_Location() - end_daughter_two->get_Location()).length();

	//this is how many new nodes will be made for each cell
	int total_num_one = static_cast<int>(daughter_one_length/(Membr_Equi_Len_Short*4));
	int total_num_two = static_cast<int>(daughter_two_length/(Membr_Equi_Len_Short*4));
		
	//this is how much change there is in the x direction of divison plane for each cell
	double x_one = end_daughter_one->get_Location().get_X() - start_daughter_one->get_Location().get_X();
	double x_two = end_daughter_two->get_Location().get_X() - start_daughter_two->get_Location().get_X();
	
	//this is how much change there is in the y direction of division plane for each cell
	double y_one = end_daughter_one->get_Location().get_Y() - start_daughter_one->get_Location().get_Y();
	double y_two = end_daughter_two->get_Location().get_Y() - start_daughter_two->get_Location().get_Y();
	 
	//this is how much x should increment by when
	//making new cell wall nodes for 
	//each cell
	double delta_x_one = (double)x_one/(double)total_num_one;
	double delta_x_two = (double)x_two/(double)total_num_two;

	
	//this is how much y should increment by when 
	//making new cell wall nodes for 
	//each cell
	double delta_y_one = (double)y_one/(double)total_num_one;
	double delta_y_two = (double)y_two/(double)total_num_two;

	//some variables to be used for making new nodes
	shared_ptr<Wall_Node> currW = start_daughter_one;
	Coord location = start_daughter_one->get_Location();
	shared_ptr<Wall_Node> prevW = start_daughter_one;
	shared_ptr<Cell> this_cell = shared_from_this();
	double curr_X = start_daughter_one->get_Location().get_X()+2*delta_x_one;
	double curr_Y = start_daughter_one->get_Location().get_Y()+2*delta_y_one;
	//this counter was for debugging
	//int counter = 0;
	
	cout << "Make duaghter one new cell wall nodes" << endl;
	for(unsigned int i = 0; i< total_num_one-4; i++){
		curr_X = curr_X + delta_x_one;
		curr_Y = curr_Y + delta_y_one;
		location = Coord(curr_X,curr_Y);
		//node made 
		cout << "making new wall node here" << endl;
		currW = make_shared<Wall_Node>(location,this_cell);
		//set neighbors
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;
		//debugging
		//counter++;
		//cout << counter << endl;	
	}
	currW->set_Right_Neighbor(end_daughter_one);
	end_daughter_one->set_Left_Neighbor(currW);
	this->left_Corner = start_daughter_one;
	
	cout << "Make daughter two new cell wall nodes" << endl;
	//counter = 0;
	currW = start_daughter_two;
	location = start_daughter_two->get_Location();;
	prevW = start_daughter_two;
	curr_X = location.get_X()+2*delta_x_two;
	curr_Y = location.get_Y()+2*delta_y_two;
	for(unsigned int i = 0; i< total_num_two-4; i++){
		curr_X = curr_X + delta_x_two;
		curr_Y = curr_Y + delta_y_two;
		location = Coord(curr_X,curr_Y);
		//node made
		cout << "make new node here" << endl;
		currW =make_shared<Wall_Node>(location,sister);
		//set neighbors		
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;	
		//debugging
		//counter++;
		//cout<< counter << endl;

	}
	currW->set_Right_Neighbor(end_daughter_two);
	end_daughter_two->set_Left_Neighbor(currW);
	sister->set_Left_Corner(start_daughter_two);
	
	cout << "Made wall nodes" << endl;
	//here i am refilling
	//the wall node vectors for each cell 
	//and setting private member variables for 
	//each new wall node
	double new_damping = this->get_Damping();	
	sister->set_Damping(new_damping);
	shared_ptr<Wall_Node> curr = sister->get_Left_Corner();
	shared_ptr<Wall_Node> orig = curr;
	double l_thresh = 0;
	double k_lin = 0;
	double k_bend = 0;
	int number_nodes_B = 0;
	//cout << "Cell two counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		curr->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(curr);
		k_lin = compute_k_lin(curr);
		k_bend = compute_k_bend(curr);
		curr->set_membr_len(l_thresh);
		curr->set_K_LINEAR(k_lin);
		curr->set_K_BEND(k_bend);

		sister->add_wall_node_vec(curr);
		curr = curr->get_Left_Neighbor();
		//cout << number_nodes_B << endl;
		//cout << curr->get_Location() << endl;
		if(curr == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
	} while (curr != orig);
	//cout << "Done counting cell two" << endl;
	sister->set_Wall_Count(number_nodes_B);
	//cout << "Sister: " << number_nodes_B << endl;
	
	this->wall_nodes.clear();

	curr = this->left_Corner;
	orig = this->left_Corner;
	int number_nodes_A = 0;
	//cout << "Celll one counting" << endl;
	do {
		number_nodes_A++;
		curr->update_Cell(this_cell);
		curr->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(curr);
		k_lin = compute_k_lin(curr);
		k_bend = compute_k_bend(curr);
		curr->set_membr_len(l_thresh);
		curr->set_K_LINEAR(k_lin);
		curr->set_K_BEND(k_bend);


		this->add_wall_node_vec(curr);
		curr = curr->get_Left_Neighbor();
		//cout << number_nodes_A << endl;
		//cout << curr->get_Location() << endl;
		if(curr == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
	} while(curr != orig);
	//cout << "Done counting cell one" << endl;	
	this->set_Wall_Count(number_nodes_A);
	//cout << "Parent nodes: " << number_nodes_A << endl;	
	
	//update cell level variables
	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();;
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	this->update_Cell_Center();
	sister->update_Cell_Center();
	//this->calc_WUS();
	//sister->calc_WUS();
	

	this->set_growth_rate();
	sister->set_growth_rate();
	this->life_length = 0;	
	this->Cell_Progress =0;		
	sister->reset_Cell_Progress();
	sister->reset_Life_Length();
	

	//cout << "Reassign cyt nodes" << endl;
	vector<shared_ptr<Cyt_Node>> temp_cyts;
	this->get_Cyt_Nodes_Vec(temp_cyts);
	this->cyt_nodes.clear();
	this->num_cyt_nodes = 0;
	
	//this is the center of the new cell wall for each cell
	Coord center_1 = (start_daughter_one->get_Location() + end_daughter_one->get_Location())*.5;
	Coord center_2 = (start_daughter_two->get_Location() + end_daughter_two->get_Location())*.5;	
	double length_1;
	double length_2;
	//for debugging
	//int counter = 0;
	for(unsigned int i = 0; i < temp_cyts.size(); i++) {
		length_1 = (temp_cyts.at(i)->get_Location()-center_1).length();
		length_2 = (temp_cyts.at(i)->get_Location()-center_2).length();	
		if(length_1 < length_2) {
			temp_cyts.at(i)->update_Cell(this_cell);
			this->update_cyt_node_vec(temp_cyts.at(i));
			//cout << temp_cyts.at(i)->get_Location() << endl;
			this->Cell_Progress++;
		}
		else{
			temp_cyts.at(i)->update_Cell(sister);
			sister->update_cyt_node_vec(temp_cyts.at(i));
			//cout << temp_cyts.at(i)->get_Location() << endl;
			sister->update_Cell_Progress();
		}
		//counter++;
		//cout << counter << endl;
	}
	//move cyt nodes that are too close 
	//to the new cell wall
	//cout << "Move cyts here" << endl;
	this->move_cyt_nodes(center_1);
	sister->move_cyt_nodes(center_2);
	
	//cout << "Finished deleting old cyt nodes" << endl;
	
	return sister;
}
void Cell::move_cyt_nodes(Coord center_point){
	Coord vector_from_center;
	double length_from_center_pt;
	Coord location;
	for(unsigned int i=0; i<cyt_nodes.size(); i++){
		length_from_center_pt = (cyt_nodes.at(i)->get_Location()-center_point).length();
		vector_from_center = cyt_nodes.at(i)->get_Location() - this->cell_center;
		//if(length_from_center_pt< 3) {
			location = cell_center + vector_from_center*.6;
			cyt_nodes.at(i)->new_location(location);
			//cout << cyt_nodes.at(i)->get_Location() << endl;
		//
	}
	return;
}
//////=======================================================================================================
//end of cell_div.cpp
