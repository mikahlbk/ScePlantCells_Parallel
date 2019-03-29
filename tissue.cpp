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
#include "phys.h"
#include "coord.h"
#include "cell.h"
#include "tissue.h"
//=========================
// Public Member Functions for Tissue.cpp
//tissue constructor makes new tissue from 
//.txt file that is read in
Tissue::Tissue(string filename) {
	num_cells = 0;
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
			shared_ptr<Cell> curr= make_shared<Cell>(rank, center, radius, my_tissue, layer,boundary, stem);
			//give that cell wall nodes and internal nodes
			curr->make_nodes(radius);
			num_cells++;
			cells.push_back(curr);
		}

		ss.clear();
	}

	ifs.close();
}

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
	cells.push_back(new_Cell);
	return;
}
//**********functions for tissue to perform on cells********//
//updates current neighbors of each cell
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
		cells.at(i)->add_Wall_Node(Ti);
	}
	//cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_Wall_Count() << endl;
	return;
}
void Tissue::delete_Wall(int Ti) {
 	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		//cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_wall_count() << endl;
		//cells.at(i)->delete_Wall_Node(Ti);
		//cout<< "Wall Count Cell " << i << ": " << this->cells.at(i)->get_wall_count() << endl;
	}
	return;
}
//updates adhesion springs for each cell
void Tissue::update_Adhesion() {
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0;i<cells.size();i++) {
		//cout << "Updating adhesion for cell" << i <<  endl;
		cells.at(i)->clear_adhesion_vectors();
	}
	for(unsigned int i=0;i<cells.size();i++) {
		cells.at(i)->update_adhesion_springs();
	}
	return;
}
//this function is not in use
void Tissue::update_Linear_Bending_Springs(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i = 0;i<cells.size();i++){
		cells.at(i)->update_Linear_Bending_Springs();
	}
	return;
}	
//**********updates cell cycle of each cell************//
void Tissue::update_Cell_Cycle(int Ti) {
	//cout << "Current number of cells: " << cells.size() << endl; 
	int number_cells = cells.size();
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {

		//cout << "updating cell" << i << endl;

		cells.at(i)->update_Cell_Progress(Ti);
	}
	//cout << "Number cells is: " << cells.size() << endl;
	return;

}
void Tissue::division_check(){
	int number_cells = cells.size();
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
void Tissue::update_Cell_Locations() {
	#pragma omp parallel for schedule(static,1)	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Node_Locations();
	}

	return;
}

void Tissue::locations_output(ofstream& ofs){
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_locations(ofs);
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
	/*for(unsigned int i=0; i < cells.size(); i++) {
		cells.at(i)->nematic(average, angle);
		ofs<< average.get_X() << endl;
	}
	//ofs << "average vec" << endl;
	for(unsigned int i=0; i < cells.size(); i++) {
		cells.at(i)->nematic(average, angle);
		ofs << average.get_Y() << endl;
	}
	//ofs << "angles" << endl;
	for(unsigned int i=0; i < cells.size(); i++) {
		cells.at(i)->nematic(average, angle);
		ofs<< angle << endl;
	}
	//ofs<< "centers2" << endl;
	//for(unsigned int i=0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
	//		ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	}
	//}
	//for(unsigned int i=0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
	//		ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
	//	}
	//}
	//ofs << "centers1" << endl;
	//for(unsigned int i=0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==1){
	//		ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	}
	//}
	//ofs << "centers1" << endl;
	//for(unsigned int i=0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==1){
	//		ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
	//	}
	//}
	//ofs << "centers2" << endl;
	for(unsigned int i=0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
			ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
	//	}
	}
	//ofs << "centers1" << endl;
	for(unsigned int i=0; i< cells.size();i++) {
	//	if(cells.at(i)->get_Layer() ==2){
			ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
	//	}
	}
	//for(unsigned int i=0; i< cells.size();i++) {
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
int Tissue::update_VTK_Indices() {

	int id = 0;
	int rel_cnt = 0;

	//iterate through cells to reassign vtk id's - starting at 0
	for (unsigned int i = 0; i < cells.size(); i++) {
		//iterate through
		rel_cnt += cells.at(i)->update_VTK_Indices(id);
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

void Tissue::print_VTK_File(ofstream& ofs) {
	int rel_cnt = update_VTK_Indices();

	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Sub_cellular elem model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	//Good up to here

	//Need total number of points for all cells
	int num_Points = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		num_Points += cells.at(i)->get_Node_Count();
	}

	ofs << "POINTS " << num_Points << " float" << endl;
	
	vector<int> start_points;
	vector<int> end_points;
	int count = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		start_points.push_back(count);
		cells.at(i)->print_VTK_Points(ofs,count);
		end_points.push_back(count - 1);
	}

	ofs << endl;
	
	//to be used without visualizing ADH springs
	//ofs << "CELLS " << cells.size()<< ' ' << (num_Points + start_points.size())  << endl;
	
	//to be used for visualizing adh springs
	ofs << "CELLS " << cells.size()+rel_cnt<< ' ' << (num_Points + start_points.size())+(rel_cnt*3)  << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		ofs << cells.at(i)->get_Node_Count();

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

	for(unsigned int i = 0; i < rel_cnt; i++) {
		//type for adh relationship
		ofs << 3 << endl;
	}

	ofs << endl;



	ofs << "POINT_DATA " << num_Points << endl;
	ofs << "SCALARS WUS  double " << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_WUS(ofs);
	}

	ofs << endl;

	ofs << "SCALARS CK  double " << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_CK(ofs);
	}

	ofs << endl;
	/*ofs << "Scalars average_pressure float" << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		//cells.at(i)->print_VTK_Scalars_Average_Pressure(ofs);
	}
	ofs << endl;
	*/

	ofs << "Scalars wall_pressure float" << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_Node(ofs);
	}
	return;
}


//=========================
//End of tissue.cpp
