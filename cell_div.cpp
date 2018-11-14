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
void Cell::find_nodes_for_div_plane(Coord& direction, vector<shared_ptr<Wall_Node>>& nodes) {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
        double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	double curr_diff = 100;
	double smallest_1 = 100;
	double smallest_2 = 100;
	double index_first;
	double index_second;
	//find first node that makes 90 degree angle
	for(unsigned int i = 0; i < walls.size();i++) {	
		curr_vec = walls.at(i)->get_Location() - cell_center;
		//cout << walls.at(i)->get_Location() << endl;
		//cout << "center" << cell_center << endl;
		//cout << "curr vec" << curr_vec << endl;
		curr_len = curr_vec.length();	
		//cout << "curr len" << curr_len << endl;
		growth_len = direction.length();
		//cout << "growth len" << growth_len << endl;
		costheta = direction.dot(curr_vec)/(curr_len*growth_len);
		//cout << "cos" << costheta << endl;
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "theta " << theta << endl;
		
		curr_diff = fabs(theta-1.5708);
		//cout << curr_diff << "curr_diff" << endl;
	//	if(curr_diff < .1){
			if(curr_diff<smallest_1){
				smallest_1 = curr_diff;
				index_first = i;
			}
	//	}
	
	}
	Coord flipped_dir = direction*(-1);
	for(unsigned int i = 0; i < walls.size();i++) {	
		curr_vec = walls.at(i)->get_Location() - cell_center;
		//cout << walls.at(i)->get_Location() << endl;
		//cout << "center" << cell_center << endl;
		//cout << "curr vec" << curr_vec << endl;
		curr_len = curr_vec.length();	
		//cout << "curr len" << curr_len << endl;
		growth_len = flipped_dir.length();
		//cout << "growth len" << growth_len << endl;
		costheta = flipped_dir.dot(curr_vec)/(curr_len*growth_len);
		//cout << "cos" << costheta << endl;
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "theta " << theta << endl;
		
		curr_diff = fabs(theta-1.5708);
		//cout << curr_diff << "curr_diff" << endl;
	//	if(curr_diff < .1){
			if(curr_diff<smallest_2){
				if(i != index_first){
				smallest_2 = curr_diff;
				index_second = i;
				}
			}
	//	}
	}
	//cout << "first" << first << endl;
	//cout << "second" << second << endl;
	nodes.push_back(walls.at(index_first));
	nodes.push_back(walls.at(index_second));
	return;	
}

void Cell::division(shared_ptr<Cell>& sister) {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
	//shared_ptr<Cell>  sister = make_shared<Cell>(this->my_tissue);
	//pointers for while loops throughout function
	//cout << num_wall_nodes << endl;
	//shared_ptr<Wall_Node> curr = NULL;
	//Wall_Node* next = NULL;
	//Wall_Node* orig = NULL;
	
	vector<shared_ptr<Wall_Node>> nodes;
	//Wall_Node* start = NULL;
	//Wall_Node* end = NULL;
	Coord horizontal = Coord(1,0);
	Coord vertical = Coord(0,1);
	
	if((this->get_Layer() != 1)&&(this->get_Layer()!=2)) {
		find_nodes_for_div_plane(vertical,nodes);
	}
	else {
		find_nodes_for_div_plane(horizontal,nodes);
	}
	shared_ptr<Wall_Node> start = nodes.at(0);
	shared_ptr<Wall_Node> end = nodes.at(1);
	//cout << "got start and end" << endl;	
	if((start==NULL)||(end==NULL)) {
		cout << "pointer problem in division"<< endl;
	}
	//cout << "Setting new start and end Points" << endl;
	shared_ptr<Wall_Node> left_start = start;
	//cout << "left start broke" << left_start << endl;
	shared_ptr<Wall_Node> left_end = end;
	//cout << "left end broke" << left_end << endl;
	shared_ptr<Wall_Node> right_start = start;
	//cout << "right start broke" << right_start << endl;
	shared_ptr<Wall_Node> right_end = end;
	//cout << "right end broke" << right_end << endl;
	//move start and end pointers to correct location
	//dont woory about deletion
	//node gets deleted when smart pointers dont point anymore
	vector<shared_ptr<Wall_Node>> adhesion_vec;
	for(unsigned int i = 0; i < 4; i++) {
	//cout << "left start" << left_start << endl;
	//shared_ptr<Wall_Node> currD = left_start;
	adhesion_vec = left_start->get_adhesion_vec();
	for(unsigned int j=0; j<adhesion_vec.size(); j++){
		adhesion_vec.at(j)->set_Closest(NULL,100);
	}
	
	left_start = left_start->get_Right_Neighbor();

	//cout << "right start" << right_start << endl;
	adhesion_vec = right_start->get_adhesion_vec();
	for(unsigned int j=0; j<adhesion_vec.size(); j++){
		adhesion_vec.at(j)->set_Closest(NULL,100);
	}
	
	right_start = right_start->get_Left_Neighbor();
	
	//currD = left_end;	
	left_end = left_end->get_Left_Neighbor();
	adhesion_vec = left_end->get_adhesion_vec();
	for(unsigned int j=0; j<adhesion_vec.size(); j++){
		adhesion_vec.at(j)->set_Closest(NULL,100);
	}
	
	right_end = right_end->get_Right_Neighbor();
	adhesion_vec = right_end->get_adhesion_vec();
	for(unsigned int j=0; j<adhesion_vec.size(); j++){
		adhesion_vec.at(j)->set_Closest(NULL,100);
	}
	
	}
	
		
	//cout << "lengths" << endl;
	//Coord left_center = (left_start->get_Location() + left_end->get_Location())*.5;
	//Coord right_center = (right_start->get_Location() + right_end->get_Location())*.5;
	double left_length = (left_end->get_Location() - left_start->get_Location()).length();
	double right_length = (right_end->get_Location() - right_start->get_Location()).length();

	//double total_num_left_dbl = (left_length/(MembrEquLen*2));
	int total_num_left = static_cast<int>(left_length/(MembrEquLen*3));
	//double total_num_right_dbl = (right_length/(MembrEquLen*2));
	int total_num_right = static_cast<int>(right_length/(MembrEquLen*3));
	//cout << "Total num right" << total_num_right << endl;	
	//cout << "Total num left" << total_num_left << endl;	
	
	double x_left = left_end->get_Location().get_X() - left_start->get_Location().get_X();
	double x_right = right_end->get_Location().get_X() - right_start->get_Location().get_X();
	//cout << x_left << endl;
	double y_left = left_end->get_Location().get_Y() - left_start->get_Location().get_Y();
	double y_right = right_end->get_Location().get_Y() - right_start->get_Location().get_Y();
	 
	double delta_x_left = (double)x_left/(double)total_num_left;
	double delta_x_right = (double)x_right/(double)total_num_right;
	//cout <<delta_x_left<< endl;
	//cout << delta_x_right << endl;
	double delta_y_left = (double)y_left/(double)total_num_left;
	double delta_y_right = (double)y_right/(double)total_num_left;
	//cout << delta_y_left << endl;
	//cout << delta_y_right << endl;
	//double start_x_left = left_center.get_X();
	//double start_x_right = right_center.get_X();
	//double start_y_left = left_center.get_Y();
	//double start_y_right = right_center.get_Y();
	
	//here i will make line from center to left start
	//left side
	//Wall_Node* left_center_node = new Wall_Node(left_center,this);
	shared_ptr<Wall_Node> currW = left_start;
	Coord location = left_start->get_Location();
	shared_ptr<Wall_Node> prevW = left_start;
	shared_ptr<Cell> this_cell = shared_from_this();
	double curr_X = left_start->get_Location().get_X()+2*delta_x_left;
	double curr_Y = left_start->get_Location().get_Y()+2*delta_y_left;
	//int counter = 0;
	//int num_nodes = floor(total_num_left*.5)-3;
	for(unsigned int i = 0; i< total_num_left-4; i++){
		curr_X = curr_X + delta_x_left;
		curr_Y = curr_Y + delta_y_left;
		location = Coord(curr_X,curr_Y);
	//	cout << "location" << location << endl;
		currW = make_shared<Wall_Node>(location,this_cell);
		currW->set_Right_Neighbor(prevW);
		prevW->set_Left_Neighbor(currW);
		prevW = currW;
		//counter++;
		//cout << counter << endl;	
	}
	currW->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(currW);
	this->left_Corner = left_start;
	//here i will make line from center to left end*/
	/*keeeeeep
 	currW = left_center_node;
	location = left_center;
	prevW = currW;
	curr_X = start_x_left;
	curr_Y = start_y_left;
	for(unsigned int i = 0; i< num_nodes; i++){
		curr_X = curr_X + delta_x_left;
		curr_Y = curr_Y + delta_y_left;
		location = Coord(curr_X,curr_Y);
		currW = new Wall_Node(location,this);
		currW->set_Right_Neighbor(prevW);
		prevW->set_Left_Neighbor(currW);
		prevW = currW;	
	}
	currW->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(currW);
	this->left_Corner = left_start;KEEEEEEEEP*/

	//right side
	//line from center to right start
	//Wall_Node* right_center_node = new Wall_Node(right_center,sister);
	currW = right_start;
	location = right_start->get_Location();;
	prevW = right_start;
	curr_X = location.get_X()+2*delta_x_right;
	curr_Y = location.get_Y()+2*delta_y_right;
	//num_nodes = floor(total_num_right*.5)-3;
	for(unsigned int i = 0; i< total_num_right-4; i++){
		curr_X = curr_X + delta_x_right;
		curr_Y = curr_Y + delta_y_right;
		location = Coord(curr_X,curr_Y);
	//	cout << "location" << location << endl;
		currW =make_shared<Wall_Node>(location,sister);
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;	
	}
	currW->set_Right_Neighbor(right_end);
	right_end->set_Left_Neighbor(currW);
	sister->set_Left_Corner(right_start);
	//here i will make line from center to left end*/
	/*keeeeeep
 	currW = right_center_node;
	location = right_center;
	prevW = currW;
	curr_X = start_x_right;
	curr_Y = start_y_right;
	for(unsigned int i = 0; i< num_nodes; i++){
		curr_X = curr_X + delta_x_right;
		curr_Y = curr_Y + delta_y_right;
		location = Coord(curr_X,curr_Y);
		currW = new Wall_Node(location,sister);
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;	
	}
	currW->set_Right_Neighbor(right_end);
	right_end->set_Left_Neighbor(currW);
	sister->set_Left_Corner(right_start);KEEEEEEP	*/
	//here i will make the left side
	//find the angle of the starting node
	/*keep 
 	Coord vector_to_start_node = left_start->get_Location()-cell_center;
	Coord vector_to_end_node = left_end->get_Location()-cell_center;
	Coord x_axis = Coord(1,0);
	double start_len = vector_to_start_node.length();	
	double end_len= vector_to_end_node.length();
	double start_costheta = x_axis.dot(vector_to_start_node)/(start_len);
	double end_costheta =2 x_axis.dot(vector_to_end_node)/(end_len);	
	double start_theta = acos( min( max(start_costheta,-1.0), 1.0) );
	double end_theta = acos( min( max(end_costheta,-1.0), 1.0) );	
	double total_angle = start_theta+end_theta;
	double angle_increment = total_angle/total_num_left;
	
	double curr_theta = start_theta;
	double curr_X = 0;
	double curr_Y = 0;
	double center_x = (left_start->get_Location().get_X() + left_end->get_Location().get_X())*.5;
	double center_y = (left_start->get_Location().get_Y() + left_end->get_Location().get_Y())*.5;
	Coord location;
	Wall_Node* currW = NULL;
	Wall_Node* prevW = left_start;
	double width = min(left_start->get_Location().get_X()-center_x,left_end->get_Location().get_X()-center_x)-.5; 
	double length = .01;//(fabs(left_start->get_Location().get_Y()-cell_center.get_Y())+ fabs(left_end->get_Location().get_Y()-cell_center.get_Y()))*.5;
	
	for(unsigned int i = 0; i < total_num_left; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = center_x + width*cos(curr_theta);
		curr_Y = center_y + length*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		currW = new Wall_Node(location, this);
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}
	currW->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(currW);
	left_Corner = left_start;
	
	vector_to_start_node = right_start->get_Location()-cell_center;
	vector_to_end_node = right_end->get_Location()-cell_center;
	x_axis = Coord(1,0);
	start_len = vector_to_start_node.length();	
	end_len= vector_to_end_node.length();
	start_costheta = x_axis.dot(vector_to_start_node)/(start_len);
	end_costheta = x_axis.dot(vector_to_end_node)/(end_len);	
	start_theta = acos( min( max(start_costheta,-1.0), 1.0) );
	end_theta = acos( min( max(end_costheta,-1.0), 1.0) );	
	total_angle = start_theta+end_theta;
	angle_increment = total_angle/total_num_right;
	
	curr_theta = end_theta;
	curr_X = 0;
	curr_Y = 0;
	center_x = (right_start->get_Location().get_X() + right_end->get_Location().get_X())*.5;
	center_y = (right_start->get_Location().get_Y() + right_end->get_Location().get_Y())*.5;
	location;
	currW = NULL;
	prevW = right_end;
	width = min(right_start->get_Location().get_X()-center_x, right_end->get_Location().get_X()-center_x)-.5; 
	length = .01;//(fabs(right_start->get_Location().get_Y()-cell_center.get_Y())+ fabs(right_end->get_Location().get_Y()-cell_center.get_Y()))*.5;
	
	for(unsigned int i = 0; i < total_num_right; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = center_x + width*cos(curr_theta);
		curr_Y = center_y + length*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		currW = new Wall_Node(location, this);
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}
	currW->set_Left_Neighbor(right_start);
	right_start->set_Right_Neighbor(currW);
	sister->set_Left_Corner(right_start);
keep*/
	//count wall nodes
	//cout << "begin delete wall nodes" << endl;
	//Wall_Node* wall = NULL;
//	cout << "made wall nodes" << endl;
	shared_ptr<Wall_Node> curr = sister->get_Left_Corner();
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;
	double k_lin = 0;
        double theta = 0;
        double costheta = 0;
	double curr_len = 0;
	double growth_len = 0;
	Coord curr_vec;	
	sister->set_growth_direction(this->get_growth_direction());
	
	int number_nodes_B = 0;
//	cout << "Cell b counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		//curr_vec = curr->get_Left_Neighbor()->get_Location() - curr->get_Location();
		//curr_len = curr_vec.length();	
		//growth_len = this->growth_direction.length();
		//costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		//theta = acos( min( max(costheta,-1.0), 1.0) );
		//curr->set_membr_len(MembrEquLen);
		curr->set_membr_len(MembrEquLen);
		k_lin = compute_k_lin(curr);//150 + 500*(1-pow(costheta,2));
		curr->set_K_LINEAR(k_lin);

		//k_lin = 150 + 500*pow(costheta,2);//162.1142*pow(theta,2) + -509.2962*theta + 650;
		//curr->set_K_LINEAR(k_lin);
	
		//cout << "add new node sister " << endl;
		sister->add_wall_node_vec(curr);
		//cout << "before vec" << endl;
		next = curr->get_Left_Neighbor();
		//cout << number_nodes_B << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while (next != orig);
//	cout << "done counting cell b" << endl;
	sister->set_Wall_Count(number_nodes_B);
//	cout << "Sister: " << number_nodes_B << endl;
	
	//this->wall_nodes.clear();
	//num_wall_nodes = 0;
	//shared_ptr<Wall_Node> current = NULL;
	while(!wall_nodes.empty()) {
		//current = current.at(wall_nodes.size() -1);
		wall_nodes.pop_back();
		num_wall_nodes--;
	}	
	curr = this->left_Corner;
	//next = NULL;
	orig = this->left_Corner;
	int number_nodes_A = 0;
	//cout << "cell a counting" << endl;
	do {
		number_nodes_A++;
		curr->update_Cell(this_cell);
		//cout << "add new node " << endl;
		//curr_vec = curr->get_Left_Neighbor()->get_Location() - curr->get_Location();
		//curr_len = curr_vec.length();	
		//growth_len = this->growth_direction.length();
		//costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		//theta = acos( min( max(costheta,-1.0), 1.0) );
		//curr->set_membr_len(MembrEquLen);
		curr->set_membr_len(MembrEquLen);
		k_lin = compute_k_lin(curr);//150 + 500*(1-pow(costheta,2));
		curr->set_K_LINEAR(k_lin);


		this->add_wall_node_vec(curr);
		curr = curr->get_Left_Neighbor();
		//k_lin = 150 + 500*pow(costheta,2);//162.1142*pow(theta,2) + -509.2962*theta + 650;
		//curr->set_K_LINEAR(k_lin);
	
		//cout << number_nodes_A << endl;
		if(curr == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		//curr = next;
	} while(curr != orig);
	//cout << "done cell a counting" << endl;	
	this->set_Wall_Count(number_nodes_A);
	//cout << "Parent nodes: " << number_nodes_A << endl;	
	//cout << "updating angles" << endl;
	//update wall angles
	//cout << "Equi" << endl;

	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();
	//cout << "Parent angles" << endl;
	this->update_Wall_Angles();
	//cout << "Sister angles" << endl;
	sister->update_Wall_Angles();
	//cout << "updating center div" << endl;
	//update cell center
	this->update_Cell_Center();
	sister->update_Cell_Center();
	this->calc_WUS();
	sister->calc_WUS();
//	this->calc_CYT();
//	sister->calc_CYT();
//	this->calc_Total_Signal();
//	sister->calc_Total_Signal();
	sister->set_Layer(this->layer);
	this->set_growth_rate();
	sister->set_growth_rate();
	//sister->set_K_LINEAR(K_LINEAR_X,K_LINEAR_Y);
	double new_damping = this->get_Damping();
	sister->set_Damping(new_damping);
	this->life_length = 0;	
	this->Cell_Progress =0;		
	//this->Cell_Progress_add_node = 0;
	sister->reset_Cell_Progress();
	//double new_rank = this->get_Tissue()->get_num_cells();
	//sister->set_Rank(new_rank);
	//this->get_Tissue()->update_Num_Cells(sister);	
	//delete all old cyt nodes*/
	/*keep 
 	int num_cyts = cyt_nodes.size();
	int new_cyt_cnt = num_cyts/2;
	Cyt_Node* cyt = NULL;
	
	while(!cyt_nodes.empty()) {
		cyt = cyt_nodes.at(cyt_nodes.size() -1);
		delete cyt;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}keeep*/
	vector<shared_ptr<Cyt_Node>> temp_cyts;
	this->get_Cyt_Nodes_Vec(temp_cyts);
	//this->cyt_nodes.clear();
	this->num_cyt_nodes = 0;
	//Cyt_Node* cyt = NULL;
	//temp_cyts.clear();
	while(!cyt_nodes.empty()) {
		cyt_nodes.pop_back();
		//num_cyt_nodes--;
	}
//	cout << "cyt nodes" << endl;
	double length_1;
	double length_2;
	double counter = 0;
	for(unsigned int i = 0; i < temp_cyts.size(); i++) {
		length_1 = (temp_cyts.at(i)->get_Location()-this->cell_center).length();
		length_2 = (temp_cyts.at(i)->get_Location()-sister->get_Cell_Center()).length();	
		if(length_1 < length_2) {
			temp_cyts.at(i)->update_Cell(this_cell);
			this->update_cyt_node_vec(temp_cyts.at(i));
			this->Cell_Progress++;
		}
		else{
			temp_cyts.at(i)->update_Cell(sister);
			sister->update_cyt_node_vec(temp_cyts.at(i));
			//sister->update_Cell_Progress_Var();
		}
		counter++;
		//cout << counter << endl;
	}
	this->move_cyt_nodes();
	sister->move_cyt_nodes();
	//cout << "Finished deleting old cyt nodes" << endl;
	
//	cout << "get most up/down and left/right for radius" << endl;*/
/*keep	Wall_Node* up_sis = NULL;
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
//	double old_cell_radius_x;
//	double old_cell_radius_y;
//	this->find_radius(old_cell_radius_x,old_cell_radius_y);
//	old_cell_radius_y = ((this->find_closest_node_top()->get_Location() - this->find_closest_node_bottom()->get_Location()).length())*0.5;	
//	old_cell_radius_x = ((this->find_closest_node_left()->get_Location() - this->find_closest_node_right()->get_Location()).length())*0.5;
//	double sister_cell_radius_x;
//	double sister_cell_radius_y;

//	sister_cell_radius_y = ((sister->find_closest_node_top()->get_Location() - sister->find_closest_node_bottom()->get_Location()).length())*0.5;	
//	sister_cell_radius_x = ((sister->find_closest_node_left()->get_Location() - sister->find_closest_node_right()->get_Location()).length())*0.5;

//	sister->find_radius(sister_cell_radius_x,sister_cell_radius_y)	
	//create new ones for each cell
//	for(int i = 0; i < new_cyt_cnt; i++) {
		//cout << "adding cytoplasm" << i << endl;
//		this->add_Cyt_Node_Div(old_cell_radius_x,old_cell_radius_y);
//		sister->add_Cyt_Node_Div(sister_cell_radius_x,sister_cell_radius_y);
//	}
	//cout << "division function" << endl;
	//for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
	//		cout << cyt_nodes.at(i)->get_Location() << endl;
	//}
	/*vector<Wall_Node*>walls;
	sister->get_Wall_Nodes_Vec(walls);
	int counter = 0;	
	cout << "sister" << endl;
	for(unsigned int i = 0; i < walls.size(); i++) {
		cout << walls.at(i)->get_Location() << endl;
		counter++;
	}
	cout << counter << endl;
	cout << "this" << endl;
	counter = 0;
	for(unsigned int i =0; i < wall_nodes.size(); i++) {
		cout << walls.at(i)->get_Location() << endl;
		counter++;
	}
	cout << counter << endl;keep*/
//	return sister;
}
void Cell::move_cyt_nodes(){
	Coord vector;
	Coord location;
	for(unsigned int i=0; i<cyt_nodes.size(); i++){
		vector = cyt_nodes.at(i)->get_Location()-cell_center;
	//	if(vector.length() < .4) {
			location = cell_center + vector*.6;
			cyt_nodes.at(i)->new_location(location);
	//	}
	}
	return;
}
/*keep Cell* Cell::divide_width_wise() {
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
