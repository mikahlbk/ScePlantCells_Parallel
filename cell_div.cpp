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
void Cell::find_nodes_for_div_plane(Coord& orientation, vector<shared_ptr<Wall_Node>>& nodes) {
	//we will find the longest axis of the cell and feed in
	//the vector perpendicular to that axis to this function
	//this input  vector is the orientation vector above
	//this function finds the two cell wall nodes
	//on the current cell where the line between them is
	//closest representation to the orientation vector
	//need to know all the cell wall nodes of curr cell
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//this vector will hold the distance of each node
	//the the division line and will be sorted on distance
	vector<pair<double,shared_ptr<Wall_Node>>> pairs;
	//need to find two points on line that goes
	//through cell center and parallel to orientation vector
	//so we get all distances and put in the above vector
	//then we will take the smallest as first node and more
	//later about how to find the second node
	Coord center = this->get_Cell_Center();
	Coord q = Coord(center.get_X(), center.get_Y());
	Coord r = Coord(center.get_X() + orientation.get_X(), center.get_Y() + orientation.get_Y());
	Coord q_to_r = Coord(r.get_X() - q.get_X(), r.get_Y() - q.get_Y());
	Coord q_to_curr_wall;
	double q_to_r_cross_q_to_curr_wall;
	double curr_distance=0;
	shared_ptr<Wall_Node> first;
	shared_ptr<Wall_Node> second;
	for(unsigned int i = 0; i< walls.size(); i++) {
	
		q_to_curr_wall = Coord(walls.at(i)->get_Location().get_X() - q.get_X(), walls.at(i)->get_Location().get_Y()-q.get_Y());
		q_to_r_cross_q_to_curr_wall = q_to_r.cross(q_to_curr_wall);
		curr_distance = abs(q_to_r_cross_q_to_curr_wall)/q_to_r.length();
		pairs.push_back(make_pair(curr_distance,walls.at(i)));
	} 
	sort(pairs.begin(), pairs.end()); 
	first = pairs[0].second;
	shared_ptr<Wall_Node> curr = first;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = NULL;
	Coord a_i;
	Coord a_j;
	double area_1 = 0;
	double area_2 = 0;
	double curr_area = 0;
	double diff = 100;
	double curr_diff;
	//cout << "Second" << endl;
	for(unsigned int i = 1; i<11; i++){
		orig = pairs[i].second;
		do {
			next = curr->get_Left_Neighbor();
			a_i = curr->get_Location() - cell_center;
			a_j = next->get_Location() - cell_center;
			curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
			area_1 += curr_area;
			curr = next;
		} while(next != orig);
		cout << "area 1" << endl;
		cout << area_1 << endl;
		curr = orig;
		orig = first;
		do {
		
		
			next = curr->get_Left_Neighbor();
			a_i = curr->get_Location() - cell_center;
			a_j = next->get_Location() - cell_center;
			curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
			area_2 += curr_area;
			curr = next;
		} while(next != orig);
		//cout << "area two" << endl;
		//cout << area_2 << endl;
		curr_diff = abs(area_1 - area_2);
		//cout << "curr diff " << curr_diff << endl;
		if(curr_diff < diff){
			second = pairs[i].second;
			diff = curr_diff;
		}

	
	}
	//cout << "diff" << diff << endl;
	//cout << "end second" << endl;
	//cout << "first" << first << endl;
	//cout << "Second"  << second << endl;
	if((first == NULL)||(second == NULL)){
	
		cout << "Did not pick up div plane nodes" << endl;
		
		exit(1);


	}
	nodes.push_back(first);
	nodes.push_back(second);
	//cout << "nodes size " << nodes.size() << endl;
	return;
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
	this->set_growth_rate();
	sister->set_growth_rate();
	this->life_length = 0;	
	this->Cell_Progress =0;		
	sister->reset_Cell_Progress();
	sister->reset_Life_Length();
	
	//this is the vector for the desired division orientation
	Coord orientation;
	if((this->layer == 1) || (this->layer == 2)){
		orientation = Coord(0,1);
	}
	else {
		orientation = Coord(1,0);
	}
	//finds node on one side of cell
	vector<shared_ptr<Wall_Node>> nodes;
	//cout << "Nodes before" << nodes.size() << endl;
	find_nodes_for_div_plane(orientation, nodes);
	//cout << "Nodes after" << nodes.size() << endl;

	shared_ptr<Wall_Node> first  = nodes.at(0);
	//finds node on other side of cell
	shared_ptr<Wall_Node> second = nodes.at(1); 	
	//***Warning***make sure the nodes you are getting rid of dont have 
	//any other nodes pointing to them
	//left and right neighbors will be reset below
	//other pointers come from adhesion!!!!
	//cout << "adhesion pointers" << endl;
	//to delete node we have to clear everything that points to it
	//neighbors will be reassigned later
	//other pointing is adhesion vec
	//every node in adhesion pairs vec has this node in its
	//adhesions pairs vec
	first->remove_from_adh_vecs();
	first->clear_adh_vec();
	second->remove_from_adh_vecs();
	second->clear_adh_vec();
	//cout << "adhesion success" << endl;
	
	//set the starting point of daughter cell one
	shared_ptr<Wall_Node> start_daughter_one = first->get_Left_Neighbor();
	//move over one so that cells have enough room between them
	//but first need to get rid of pointers in adhesion	
	start_daughter_one->remove_from_adh_vecs();
	start_daughter_one->clear_adh_vec();	
	start_daughter_one = start_daughter_one->get_Left_Neighbor();
	//move over one so that cells have enough room between them
	//but first need to get rid of pointers in adhesion	
	//start_daughter_one->remove_from_adh_vec();
	//start_daughter_one->clear_closest_in_adh_vec();	
	//move happens here
	//start_daughter_one = start_daughter_one->get_Left_Neighbor();
	
	//set ending point for daughter cell one
	shared_ptr<Wall_Node>end_daughter_one = second->get_Right_Neighbor();
	//move over one so that cells have enough room between them
	//but first need to get rid of pointers in adhesion
	end_daughter_one->remove_from_adh_vecs();
	end_daughter_one->clear_adh_vec();	
	end_daughter_one = end_daughter_one ->get_Right_Neighbor();
	//move over one so that cells have enough room between them
	//but first get rid of adhesion pointers
	//end_daughter_one->remove_from_adh_vec();
	//end_daughter_one->clear_closest_in_adh_vec();	
	//move happens here
	//end_daughter_one = end_daughter_one ->get_Right_Neighbor();
	
	//set starting point for daughter cell two	
	shared_ptr<Wall_Node>start_daughter_two = second->get_Left_Neighbor();
	start_daughter_two->remove_from_adh_vecs();
	start_daughter_two->clear_adh_vec();	
	start_daughter_two = start_daughter_two->get_Left_Neighbor();
	//start_daughter_two->remove_from_adh_vec();
	//start_daughter_two->clear_closest_in_adh_vec();	
	//start_daughter_two = start_daughter_two->get_Left_Neighbor();
	
	//set ending point for daughter cell two
	shared_ptr<Wall_Node>end_daughter_two = first->get_Right_Neighbor();
	end_daughter_two->remove_from_adh_vecs();
	end_daughter_two->clear_adh_vec();	
	end_daughter_two = end_daughter_two ->get_Right_Neighbor();
	//end_daughter_two->remove_from_adh_vec();
	//end_daughter_two->clear_closest_in_adh_vec();	
	//end_daughter_two = end_daughter_two ->get_Right_Neighbor();


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
	int counter = 0;
	
	//cout << "Make daughter one new cell wall nodes" << endl;
	for(unsigned int i = 0; i< total_num_one-4; i++){
		curr_X = 
		curr_X + delta_x_one;
		curr_Y = curr_Y + delta_y_one;
		location = Coord(curr_X,curr_Y);
		//node made 
		//cout << "making new wall node here" << endl;
		currW = make_shared<Wall_Node>(location,this_cell);
		//set neighbors
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;
		//debugging
		counter++;
		//cout << counter << endl;	
	}
	currW->set_Right_Neighbor(end_daughter_one);
	end_daughter_one->set_Left_Neighbor(currW);
	this->left_Corner = start_daughter_one;
	//cout << "one new" << counter << endl;
	//cout << "Make daughter two new cell wall nodes" << endl;
	counter = 0;
	currW = start_daughter_two;
	location = start_daughter_two->get_Location();;
	prevW = start_daughter_two;
	curr_X = location.get_X()+2*delta_x_two;
	curr_Y = location.get_Y()+2*delta_y_two;
	for(unsigned int i = 0; i< total_num_two-4; i++){
	
		curr_X = curr_X + delta_x_two;
		curr_Y =
		curr_Y + delta_y_two;
		location = Coord(curr_X,curr_Y);
		//node made
		//cout << "make new node here" << endl;
		currW =make_shared<Wall_Node>(location,sister);
		//set neighbors		
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;	
		//debugging
		counter++;
		//cout<< counter << endl;

	}
	currW->set_Right_Neighbor(end_daughter_two);
	end_daughter_two->set_Left_Neighbor(currW);
	sister->set_Left_Corner(start_daughter_two);
	//cout << "two new" << counter << endl;
	//cout << "Made wall nodes" << endl;
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

	Coord vector_from_center;
	Coord new_location;
	//for debugging
	counter = 0;
	Coord center_1 = (end_daughter_one->get_Location() - start_daughter_one->get_Location())*.5;
	Coord center_2 = (end_daughter_two ->get_Location()-start_daughter_two->get_Location())*.5;
	double length_1;
	double length_2;
	for(unsigned int i = 0; i < temp_cyts.size(); i++) {
		length_1 = (temp_cyts.at(i)->get_Location() - this->get_Cell_Center()).length();
		length_2 = (temp_cyts.at(i)->get_Location() - sister->get_Cell_Center()).length();
		if(length_1 < length_2) {
			temp_cyts.at(i)->update_Cell(this_cell);
			this->update_cyt_node_vec(temp_cyts.at(i));
			//cout << temp_cyts.at(i)->get_Location() << endl;
			this->Cell_Progress++;
			//cout << "this" << endl;
		}
		else{
			temp_cyts.at(i)->update_Cell(sister);
			sister->update_cyt_node_vec(temp_cyts.at(i));
			//cout << temp_cyts.at(i)->get_Location() << endl;
			sister->update_Cell_Progress();
			//cout << "sister" << endl;
		}
		counter++;
		cout << "counter" <<  counter << endl;
	}
	//move cyt nodes that are too close 
	//to the new cell wall
	//cout << "Move cyts here" << endl;
	//Coord center_1 = Coord(0,0);
	//Coord center_2 = Coord(0,0);
	this->move_cyt_nodes(center_1);
	sister->move_cyt_nodes(center_2);
	
	cout << "Finished deleting old cyt nodes" << endl;
	

return sister;
}
void Cell::move_cyt_nodes(Coord center_point){
	Coord vector_from_center;
	double length_from_center_pt;
	Coord location;
	for(unsigned int i=0; i<cyt_nodes.size(); i++){
		length_from_center_pt = (cyt_nodes.at(i)->get_Location()-center_point).length();
		vector_from_center = cyt_nodes.at(i)->get_Location() - this->cell_center;
		//if(length_from_center_pt< 4) {
	
			location = cell_center + vector_from_center*.8;
			cyt_nodes.at(i)->new_location(location);
			//cout << cyt_nodes.at(i)->get_Location() << endl;
		//}
	}
	return;
}
//////=======================================================================================================
//end of cell_div.cpp
