//tissue.h
//=================
//Include Guards
#ifndef _TISSUE_H_INCLUDED_
#define _TISSUE_H_INCLUDED_
//=========================
//forward declarations

//=======================
//Include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
//=========================
// Tissue Class Declaration

//class Tissue: public enable_shared_from_this<Tissue> {
class Tissue{
	private:
		// We'll need a better data structure for later
		vector<shared_ptr<Cell>> cells;
		vector<shared_ptr<Cell>>top_cells;
		int num_cells;
	public:
		Tissue(string filename);
	//	shared_ptr<Tissue> getptr();
		void get_Cells(vector<shared_ptr<Cell>>& cells);
		void get_top_cells(vector<shared_ptr<Cell>>& cells);
		int  get_num_cells() {return num_cells;}
		void update_Num_Cells(shared_ptr<Cell>& new_Cell);
		void update_WUS(int Ti);
		void update_Cell_Cycle(int Ti);
		void add_Wall(int Ti);
		void delete_Wall(int Ti);
		void update_Linear_Bending_Springs();
		void calc_New_Forces(int Ti);
		void update_Cell_Locations();
		void update_Neighbor_Cells();
		void update_Adhesion();
		void update_Microfibrils(int Ti);
		void compression_Test();
		void pressure();
		void add_cyt_node();
		void set_Stationary_Points(int Ti);
		void stretching_Test();
		void elastic_mod_measurements();
		void cell_area();
		void nematic_output(ofstream& ofs);
		//void cell_strain();
		void make_Vectors();
		void print_Data_Output(ofstream& ofs);
		void locations_output(ofstream& ofs);
		int update_VTK_Indices();
		void print_VTK_File(ofstream& ofs);
		int get_Num_Cells() {return num_cells;}
		//Destructor
		~Tissue();
};


//===========================
//End of file

#endif
