//tissue.cpp
//=========================

//=========================
//Include Dependencies
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>
#include <random>
#include "phys.h"
#include "coord.h"
#include "cell.h"
#include "tissue.h"
//=========================
// Public Member Functions for Tissue.cpp
//tissue constructor makes new tissue from 
//.txt file that is read in
//TEMPLATE CONSTRUCTOR

Tissue::Tissue(string filename) {
	num_cells = 0;
	num_deleted = 0;
	set_up_counts();
	ifstream ifs(filename.c_str());

	if(!ifs) {
		cout << filename << " is not available" << endl;
		return;
	}

	stringstream ss;
	string line;
	string temp;
	char trash;
	int rank;
	int layer;
	int boundary;
	int stem;
	double radius;
	Coord center;
	double x, y;
	Tissue* my_tissue = this;
	while (getline(ifs,line)) {
		ss.str(line);

		getline(ss,temp,':');

		if (temp == "CellRank") {
			ss >> rank;
		}
		else if (temp == "Center") {
			ss >> x >> trash >> y;
			Coord loc(x,y);
			center = loc;
		}
		else if (temp == "Radius") {
			ss >> radius;
		}
		else if (temp == "Layer") {
			ss >> layer;
		}
		else if (temp == "Boundary"){
			ss >> boundary;
		}
		else if(temp == "Stem"){
			ss >> stem;
		}
		else if (temp == "End_Cell") {
			//create new cell with collected data 
			//and push onto vector that holds all cells in tissue 
			//cout<< "making a cell" << endl;
			if (layer >= STEM_LAYER) { 
				stem = 1;
			} else { 
				stem = 0;
			}
			shared_ptr<Cell> curr= make_shared<Cell>(rank, center, radius, my_tissue, layer,boundary, stem);
			//give that cell wall nodes and internal nodes
			curr->make_nodes(radius);
			//curr->make_nodes_experimental("experimental_nodes.txt");
			//cout<< "make nodes" << endl;
			num_cells++;
			cells.push_back(curr);
		}
		ss.clear();
	}
	ifs.close();
}
///EXPERIMENTAL CONSTRUCTOR
/*
Tissue::Tissue(string filename) {
	num_cells = 0;
	set_up_counts();
	ifstream ifs(filename.c_str());

	if(!ifs) {
		cout << filename << " is not available" << endl;
		return;
	}
	shared_ptr<Wall_Node> lc;
	shared_ptr<Wall_Node> curr;
	shared_ptr<Wall_Node> temp;
	Coord loc;

	stringstream ss;
	string line;
	string temp;
	char trash;
	int rank;
	int this_rank = -1;
	int layer;
	int boundary;
	int stem;
	double radius;
	Coord center;
	double x, y;
	Tissue* my_tissue = this;

	do { 
		ss.str(line);
		if (ss.eof()) { 
			//Do stuff because we're done
			break;
		}
		ss << x << y << rank;
		if (rank != this_rank) { 
			//x,y belong to new cell. Finish
			//up old cell (if it exists), 
			//and continue.
			if (this_rank == -1) { 
				//do nothing
			} else { 
				//Finish up old cell. 
			}
		}



	} while();
	
	while (getline(ifs,line)) {
		ss.str(line);

		getline(ss,temp,':');

		if (temp == "CellRank") {
			ss >> rank;
		}
		else if (temp == "Center") {
			ss >> x >> trash >> y;
			Coord loc(x,y);
			center = loc;
		}
		else if (temp == "Radius") {
			ss >> radius;
		}
		else if (temp == "Layer") {
			ss >> layer;
		}
		else if (temp == "Boundary"){
			ss >> boundary;
		}
		else if(temp == "Stem"){
			ss >> stem;
		}
		else if (temp == "End_Cell") {
			//create new cell with collected data 
			//and push onto vector that holds all cells in tissue 
			//cout<< "making a cell" << endl;
			shared_ptr<Cell> curr= make_shared<Cell>(rank, center, radius, my_tissue, layer,boundary, stem);
			//give that cell wall nodes and internal nodes
			curr->make_nodes(radius);
			//curr->make_nodes_experimental("experimental_nodes.txt");
			//cout<< "make nodes" << endl;
			num_cells++;
			cells.push_back(curr);
		}
		ss.clear();
	}
	ifs.close();
}*/

Tissue::~Tissue() {
	//not necessary because 
	//using smartpointers
}

//*********functions for tissue to return information i.e.************//
//*********updating or returning the number of cells in the tissue*******//
void Tissue::get_Cells(vector<shared_ptr<Cell>>& cells) {
	cells = this->cells;
	return;
}
void Tissue::update_Num_Cells(shared_ptr<Cell>& new_Cell) {
	num_cells++;
	//cout << num_cells << endl;
	if (isnan(new_Cell->get_Left_Corner()->get_Location().get_X())) { 
		cout << "Passed a bad cell!" << endl;
		cout << "Num_Cells:" << endl;
		exit(1);
	}
	cells.push_back(new_Cell);
	return;
}
Coord Tissue::Compute_L1_AVG(){
	Coord avg;
	double avgx = 0;
	double avgy = 0;

	//int counter = 1;
	for(unsigned int i= 0; i< cells.size(); i++){
		if(cells.at(i)->get_Layer() == 1){
			//cout << "counter: " << counter <<  endl;
			//cout << cells.at(i)->get_Cell_Center().get_X()<< endl;
			avgx = avgx + cells.at(i)->get_Cell_Center().get_X();
			//cout << cells.at(i)->get_Cell_Center().get_Y() << endl;
			avgy = avgy + cells.at(i)->get_Cell_Center().get_Y();
			//cout << avgx << "avg x" << endl;
			//cout << avgy << "avg y" << endl;
			//counter++;
		}
	}
	//cout << "really?? " << avgx << endl;
	//cout << "really?? " << avgy << endl;
	avgx = avgx/cells.size();
	avgy = avgy/cells.size();
	avg = Coord(avgx,avgy);

	return avg;
}
//**********functions for tissue to perform on cells********//
//updates current neighbors of each cell
void Tissue::assign_dist_vecs(vector<int> dist1, vector<int> dist2, vector<int> dist3, vector<int> dist4){
	this->dist1 = dist1;
	this->dist2 = dist2;
	this->dist3 = dist3;
	this->dist4 = dist4;
	for(unsigned int i = 0; i<4; i++){
		counts.push_back(0);
	}
}
void Tissue::set_up_counts(){
	for(unsigned int i = 0; i < 4; i++){
		counts.push_back(0);
	}
}
int Tissue::return_counts(int index){
	return counts.at(index);
}
void Tissue::set_counts(int index){
	counts.at(index) = counts.at(index) + 1;
}
int Tissue::get_next_random(int dist, int count){
	if(dist == 1){
		return this->dist1.at(count);
	}
	else if (dist == 2){
		return this->dist2.at(count);
	}
	else if (dist == 3){
		return this->dist3.at(count);
	}
	else{
		return this->dist4.at(count);
	}
}
void Tissue::update_Signal(){
	
	Coord L1_AVG = this->Compute_L1_AVG();
	for(int i = 0; i < num_cells; i++){
		//cout<< "WUS" << endl;
		cells.at(i)->calc_WUS(L1_AVG);
		//cout << "CK" << endl;
		cells.at(i)->calc_CK(L1_AVG);
		//cout << "GROWTH RATE" << endl;
		cells.at(i)->set_growth_rate();
		//cout<< "growth rate: " << i << " " << cells.at(i)->get_growth_rate() << endl;
	}
	return;

}
void Tissue::update_growth_direction(){
	
	for(int i = 0; i < num_cells; i++){
		cells.at(i)->update_growth_direction();
	}
	return;
}
void Tissue::update_Neighbor_Cells() {
	//update vectors of neighboring cells
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Neighbor_Cells();	
	}
	return;
}
//adds node to cell wall for each cell
void Tissue::add_Wall(int Ti) {
	//#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->add_Wall_Node_Check(Ti);
	}
	//cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_Wall_Count() << endl;
	return;
}
void Tissue::delete_Wall(int Ti) {
 	// (pragma) pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		//cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_wall_count() << endl;
		cells.at(i)->delete_Wall_Node_Check(Ti);
		//cout<< "Wall Count Cell " << i << ": " << this->cells.at(i)->get_wall_count() << endl;
	}
	return;
}

void Tissue::inc_Num_Deleted() { 
	this->num_deleted++; 
	cout << "Total Number of Wall Nodes Deleted: " << this->num_deleted << endl;
	return;
}

//updates adhesion springs for each cell
void Tissue::update_Adhesion() {
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0;i<cells.size();i++) {
		//cout << "Updating adhesion for cell" << i <<  endl;
		cells.at(i)->clear_adhesion_vectors();
	}
	for(unsigned int i = 0;i<cells.size();i++) {
		cells.at(i)->update_adhesion_springs();
		cells.at(i)->update_Adh_Neighbors();
	}
	return;
}
//this function is not in use
void Tissue::update_Linear_Bending_Springs(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->update_Linear_Bending_Springs();
	}
	return;
}	
//**********updates cell cycle of each cell************//
void Tissue::update_Cell_Cycle(int Ti) {
	//cout << "Current number of cells: " << cells.size() << endl; 
	// UNUSED int number_cells = cells.size();
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {

		//cout << "updating cell" << i << endl;

		cells.at(i)->update_Cell_Progress(Ti);
	}
	//cout << "Number cells is: " << cells.size() << endl;
	return;

}
void Tissue::division_check(){
	// UNUSED int number_cells = cells.size();
	//#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		//cout << "dating cell" << i << endl;
		cells.at(i)->division_check();
	}
	return;
}
	
//calculates the forces for nodes of  each cell 
void Tissue::calc_New_Forces(int Ti) {
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		
		//cout << "Calc forces for cell: " << i << endl;
		cells.at(i)->calc_New_Forces(Ti);
		//cout << "success for cell: " << i << endl;
	
	}
	return;
}
//updates the location of all the nodes of each cell
void Tissue::update_Cell_Locations(int Ti) {
	#pragma omp parallel for schedule(static,1)	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Node_Locations(Ti);
		//if(Ti%5000 == 0){
		//vector<pair<shared_ptr<Wall_Node>,double>> nodes;
		//cells.at(i)->find_Largest_Length(nodes);
		//}
	}

	return;
}

void Tissue::locations_output(ofstream& ofs, bool cytoplasm){
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_locations(ofs,cytoplasm);
	}
return;
}
/*void Tissue::nematic_output(ofstream& ofs){
	Coord average;
	double angle;
	for(unsigned int i = 0; i<cells.size();i++){
		ofs << cells.at(i)->get_Layer()<<endl;
		ofs << cells.at(i)->get_Rank() << endl;
	}*/
	//ofs << "average vec" << endl;
	/*for(unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->nematic(average, angle);
		ofs<< average.get_X() << endl;
	}
	//ofs << "average vec" << endl;
	for(unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->nematic(average, angle);
		ofs << average.get_Y() << endl;
	}
	//ofs << "angles" << endl;
	for(unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->nematic(average, angle);
		ofs<< angle << endl;
	}
	//ofs<< "centers2" << endl;
	//for(unsigned int i = 0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
	//		ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	}
	//}
	//for(unsigned int i = 0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
	//		ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
	//	}
	//}
	//ofs << "centers1" << endl;
	//for(unsigned int i = 0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==1){
	//		ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	}
	//}
	//ofs << "centers1" << endl;
	//for(unsigned int i = 0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==1){
	//		ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
	//	}
	//}
	//ofs << "centers2" << endl;
	for(unsigned int i = 0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
			ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	}
	}
	//ofs << "centers1" << endl;
	for(unsigned int i = 0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
			ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
	//	}
	}
	//for(unsigned int i = 0; i< cells.size();i++) {
	//	ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	ofs <<cells.at(i)->get_Cell_Center().get_Y() << endl;
	//}
	for (unsigned int i = 0; i < cells.size(); i++) {
		ofs << cells.at(i)->get_WUS_concentration()<<endl;
	}
	for (unsigned int i = 0; i < cells.size(); i++) {
	//	ofs << cells.at(i)->average_Pressure() << endl;
	}
		

	
	return;
}*/
//***Functions for VTK output****//
int Tissue::update_VTK_Indices(bool cytoplasm) {

	int id = 0;
	int rel_cnt = 0;

	//iterate through cells to reassign vtk id's - starting at 0
	for (unsigned int i = 0; i < cells.size(); i++) {
		//iterate through
		rel_cnt += cells.at(i)->update_VTK_Indices(id, cytoplasm);
	}

	//cout << "final ID: " << id << endl;
	//cout << "rel_cnt: " << rel_cnt << endl;

	return rel_cnt;
}

void Tissue::print_VTK_Direction_File(ofstream& ofs){
	
	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Points representing direction vector for Cells" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	//Good up to here

	//Need total number of points for all cells
	int num_Points = 0;
	for (unsigned int i = 0; i < cells.size(); i++){
		num_Points = num_Points+2;;
	}

	ofs << "POINTS " << num_Points << " float" << endl;
	
	for(unsigned int i = 0; i < cells.size(); i++){
		cells.at(i)->print_direction_vec(ofs);
	}

	ofs << endl;
	
	ofs << "CELLS " << cells.size() << ' ' << 2*cells.size() << endl;
	
	int k = 0;
	for(unsigned int i = 0; i< cells.size() ; i++) {
		ofs << 1 << ' ' << k << endl;
		k++;
	}

	ofs << endl;

	ofs << "CELL_TYPES " << cells.size() << endl;
	for (unsigned int i = 0; i < cells.size(); i++){
		ofs << 1 << endl;
	}
	return;
}

void Tissue::print_VTK_File(ofstream& ofs, bool cytoplasm) {
	//Argument of update_VTK_indices is cytoplasm, whether or not to index cyt nodes.
	int rel_cnt = update_VTK_Indices(cytoplasm);
	
	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Sub_cellular elem model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	//Good up to here
	
	//Need total number of points for all cells
	int num_Points = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		if (cytoplasm) { 
			num_Points += cells.at(i)->get_Node_Count();
		} else { 
			num_Points += cells.at(i)->get_Wall_count();
		}
	}
	
	ofs << "POINTS " << num_Points << " float64" << endl;
	
	vector<int> start_points;
	vector<int> end_points;
	int count = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		start_points.push_back(count);
		cells.at(i)->print_VTK_Points(ofs,count,cytoplasm);
		end_points.push_back(count - 1);
	}

	ofs << endl;
	
	//to be used without visualizing ADH springs
	//ofs << "CELLS " << cells.size()<< ' ' << (num_Points + start_points.size())  << endl;
	
	//to be used for visualizing adh springs
	ofs << "CELLS " << cells.size() + rel_cnt<< ' ' << 
		(num_Points + start_points.size()) + (rel_cnt*3)  << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		if (cytoplasm) { 
			ofs << cells.at(i)->get_Node_Count();
		} else { 
			ofs << cells.at(i)->get_Wall_count();
		}

		for (int k = start_points.at(i); k <= end_points.at(i); k++) {
			ofs << ' ' << k;
		}
		ofs << endl;
	}
	
	//output pairs of node indices to draw adh line
	for(unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Adh(ofs);
	}

	ofs << endl;

	//no adh visualization
	//ofs << "CELL_TYPES " << start_points.size() << endl;
	
	//adh visualization
	ofs << "CELL_TYPES " << start_points.size()+rel_cnt << endl;
	for (unsigned int i = 0; i < start_points.size(); i++) {
		ofs << 2 << endl;
	}
	
	for(int i = 0; i < rel_cnt; i++) {
		//type for adh relationship
		ofs << 3 << endl;
	}

	ofs << endl;



	ofs << "POINT_DATA " << num_Points << endl;
	ofs << "SCALARS WUS float64 " << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_WUS(ofs, cytoplasm);
	}

	ofs << endl;

	ofs << "SCALARS CK  float64 " << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_CK(ofs, cytoplasm);
	}

	ofs << endl;

	/*ofs << "Scalars average_pressure float" << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		//cells.at(i)->print_VTK_Scalars_Average_Pressure(ofs);
	}
	ofs << endl;
	*/

	ofs << "Scalars wall_pressure float64" << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_Node(ofs, cytoplasm);
	}

	ofs << endl;
	
	ofs << "Scalars tensile_stress float64" << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Tensile_Stress(ofs, cytoplasm);
	}

	ofs << endl;

	ofs << "Scalars shear_stress float64" << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Shear_Stress(ofs, cytoplasm);
	}

	ofs << endl;

	ofs << "Scalars Neighbors float64" << 1 << endl;
	ofs << "LOOKUP_TABLE discrete_colors" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Neighbors(ofs, cytoplasm);
	}

	ofs << endl;

	ofs << "Scalars Growth_Direction float64" << 1 << endl;
	ofs << "LOOKUP_TABLE discrete_colors" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Growth_Dir(ofs, cytoplasm);
	}
	ofs << endl;

	ofs << "Scalars Curved_Walls float64" << 1 << endl;
	ofs << "LOOKUP_TABLE discrete_colors" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Curved(ofs, cytoplasm);
	}
	ofs << endl;

	ofs << "Scalars Corners float64" << 1 << endl;
	ofs << "LOOKUP_TABLE discrete_colors" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Corners(ofs, cytoplasm);
	}
	ofs << endl;

	ofs << "Scalars Mother_Daughter float64" << 1 << endl;
	ofs << "LOOKUP_TABLE discrete_colors" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_MD(ofs, cytoplasm);
	}
	ofs << endl;


	return;
}

void Tissue::NAN_CATCH(int Ti) {
	#pragma omp parallel for schedule(static,1)	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->NAN_CATCH(Ti);
	}
	return;
}
 
void Tissue::BAD_CATCH(int call, int Ti) { 
	unsigned int ui_num_cells = num_cells;
	if (cells.size() != ui_num_cells) { 
		cout << "BAD_CATCH num_cells disparity triggered at location " << call << " at Time = " << Ti << endl;
	}
	if (num_cells < 0) { 
		cout << "Negative num_cells at location " << call << " at Time = " << Ti << endl;
	}
	double x_curr, y_curr, x_neigh, y_neigh;

	for (unsigned int i = 0; i < cells.size(); i++) { 
		x_curr = cells.at(i)->get_Left_Corner()->get_Location().get_X();
		y_curr = cells.at(i)->get_Left_Corner()->get_Location().get_Y();
		x_neigh = cells.at(i)->get_Left_Corner()->get_Right_Neighbor()->get_Location().get_X();
		y_neigh = cells.at(i)->get_Left_Corner()->get_Right_Neighbor()->get_Location().get_Y();

		if (isnan(x_curr) || isnan(y_curr)) { 
			cout << "Bad X-coord in left corner of cell " << cells.at(i)->get_Rank() << " call " << call << " Time = " << Ti;
		}
		if (x_curr == x_neigh && y_curr == y_neigh) { 
			cout << "Overlapping nodes in cell " << cells.at(i)->get_Rank() << " call " << call << " Time = " << Ti;
			cout << "x_curr = " << x_curr << endl << "; y_curr = " << y_curr << endl;
		}
	}
	return;
}
//=========================
//End of tissue.cpp
