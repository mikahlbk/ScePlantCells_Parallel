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

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
//=========================
// Tissue Class Declaration

class Tissue {

	private:
		// We'll need a better data structure for later
		vector<Cell*> cells;
		int num_cells;
	public:
		Tissue(string filename);
		void get_Cells(vector<Cell*>& cells);
		void update_Num_Cells(Cell*& new_Cell);
		void update_Cell_Cycle(int Ti);
		void update_Wall();
		void calc_New_Forces(int Ti);
		void update_Cell_Locations();
		void update_Neighbor_Cells();
		void update_Adhesion();
		void compression_Test();
		void pressure();
		void add_cyt_node();
		void set_Stationary_Points(int Ti);
		void stretching_Test();
		void elastic_mod_measurements();
		void cell_area();
		//void cell_strain();
		void make_Vectors();
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices();
		void print_VTK_File(ofstream& ofs);
		int get_Num_Cells() {return num_cells;}
		//Destructor
		~Tissue();
};


//===========================
//End of file

#endif
