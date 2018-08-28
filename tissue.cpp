//tissue.cpp
//=========================

//=========================
//Include Dependencies
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <omp.h>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "cell.h"
#include "tissue.h"
//=========================
// Public Member Functions for Tissue.cpp

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
	double radius;
	Coord center;
	double x, y;
	//cout << "hey" << endl;
	//shared_ptr<Tissue> sp_this1 = make_shared<Tissue>();
	//shared_ptr<Tissue> sp_this2 = sp_this1->getptr();
	Tissue* sp_this2 = this;
	//cout << "made it" << endl;
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
		else if (temp == "End_Cell") {
			//create new cell with collected data and push onto vector 
			//cout<< "making a cell" << endl;
			shared_ptr<Cell> curr= make_shared<Cell>(rank, center, radius, sp_this2, layer);
			num_cells++;
			cells.push_back(curr);
		}

		ss.clear();
	}

	ifs.close();
}

Tissue::~Tissue() {
	
	//shared_ptr<Cell> curr = NULL;
	//while ( !cells.empty() ) {
	//	curr = cells.at(cells.size() - 1);
	//	delete curr;
	//	cells.pop_back();
	//}
}
//*********functions for tissue to return information i.e.************//
//*********updating or returning the number of cells in the tissue*******//
void Tissue::get_Cells(vector<shared_ptr<Cell>>& cells) {
	cells = this->cells;
	return;
}
//shared_ptr<Tissue> Tissue::getptr(){
//	return shared_from_this();
//	return;
//}
void Tissue::update_Num_Cells(shared_ptr<Cell>& new_Cell) {
	num_cells++;
	cout << num_cells << endl;
	cells.push_back(new_Cell);
	return;
}
//**********function for tissue to perform on cells********//
//**********updates cell cycle of each cell************//
void Tissue::update_Cell_Cycle(int Ti) {
	//cout << "Current number of cells: " << cells.size() << endl; 
	int number_cells = cells.size();
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < number_cells; i++) {
		//cout << "updating cell" << i << endl;
		cells.at(i)->update_Cell_Progress(Ti);
	}
	//cout << "Number cells is: " << cells.size() << endl;
	return;
}
//adds node to cell wall if needed for each cell
void Tissue::add_Wall(int Ti) {
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->add_Wall_Node(Ti);
	}
	//cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_Wall_Count() << endl;
	return;
}
void Tissue::delete_Wall(int Ti) {
 	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_wall_count() << endl;
		
		cells.at(i)->delete_Wall_Node(Ti);

		cout<< "Wall Count Cell " << i << ": " << this->cells.at(i)->get_wall_count() << endl;
	}
	return;
}
//calculates the forces for nodes of  each cell 
void Tissue::calc_New_Forces(int Ti) {
	if(Ti>15000){
	//#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		
		cout << "Calc forces for cell: " << i << endl;
		
		cells.at(i)->calc_New_Forces(Ti);
		
		cout << "success for cell: " << i << endl;
	
	}
	}	
	else{
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		
		cout << "Calc forces for cell: " << i << endl;
		
		cells.at(i)->calc_New_Forces(Ti);
		
		cout << "success for cell: " << i << endl;
	
	}	
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
//updates current neighbors of each cell
void Tissue::update_Neighbor_Cells() {
	//update vectors of neighboring cells
	#pragma omp parallel for schedule(static,1)
	for (unsigned int i = 0; i < cells.size(); i++) {
		cout<< "cell" << i << "neighbors" << endl;
		cells.at(i)->update_Neighbor_Cells();	
		cout << "cell" << i << "done" << endl;
	}
	return;
}
//updates adhesion springs for each cell
void Tissue::update_Adhesion() {
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0;i<cells.size();i++) {
		//cout << "Updating adhesion for cell" << endl;
		cells.at(i)->update_adhesion_springs();
	}
	return;
}
/*void Tissue::update_Microfibrils(int Ti) {
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i < cells.size(); i++) {
		cells.at(i)->update_microfibril_springs();
	}
	return;
}*/
//this is taken care of in cell progress function in cell class
/*void Tissue::add_cyt_node(){
	#pragma omp parallel for schedule(static,1)
	for(unsigned int i =0; i < cells.size(); i++) {
		cells.at(i)->add_Cyt_Node();
	}
	return;
}*/
/*void Tissue::nematic_output(ofstream& ofs){
	Coord average;
	double angle;
	//ofs << "average vec" << endl;
	for(unsigned int i=0; i < cells.size(); i++) {
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
	for(unsigned int i=0; i< cells.size();i++) {
		if(cells.at(i)->get_Layer() ==2){
			ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
		}
	}
	for(unsigned int i=0; i< cells.size();i++) {
		if(cells.at(i)->get_Layer() ==2){
			ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
		}
	}
	//ofs << "centers1" << endl;
	for(unsigned int i=0; i< cells.size();i++) {
		if(cells.at(i)->get_Layer() ==1){
			ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
		}
	}
	//ofs << "centers1" << endl;
	for(unsigned int i=0; i< cells.size();i++) {
		if(cells.at(i)->get_Layer() ==1){
			ofs<< cells.at(i)->get_Cell_Center().get_Y() << endl;
		}
	}
	for(unsigned int i=0; i< cells.size();i++) {
		ofs<< cells.at(i)->get_Cell_Center().get_X() << endl;
		ofs <<cells.at(i)->get_Cell_Center().get_Y() << endl;
	}
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_WUS_cell(ofs);
	}
//	for (unsigned int i = 0; i < cells.size(); i++) {
//		cells.at(i)->print_VTK_Scalars_Average_Pressure_cell(ofs);
//	}
		

	
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

//	cout << "final ID: " << id << endl;
//	cout << "rel_cnt: " << rel_cnt << endl;

	return rel_cnt;
}

void Tissue::print_VTK_File(ofstream& ofs) {
	int rel_cnt = update_VTK_Indices();


	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Sub_cellular elem model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	// Good up to here

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
		cells.at(i)->print_VTK_Points(ofs, count);
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
	
//	//output pairs of node indices to draw adh line
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
	//ofs << "SCALARS magnitude double " << 1 << endl;
	//ofs << "LOOKUP_TABLE default" << endl;
	//for (unsigned int i = 0; i < cells.size(); i++) {
	//	cells.at(i)->print_VTK_Scalars_WUS(ofs);
	//}

	//ofs << endl;

	//ofs << "Scalars average_pressure float" << endl;
	//ofs << "LOOKUP_TABLE default" << endl;
	//for (unsigned int i = 0; i < cells.size(); i++) {
	//	cells.at(i)->print_VTK_Scalars_Average_Pressure(ofs);
//	}
//	ofs << endl;

	ofs << "Scalars wall_pressure float" << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_Node(ofs);
	}
	return;
}


//=========================
//End of tissue.cpp
