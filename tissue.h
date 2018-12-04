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
		int num_cells;
	public:
		Tissue(string filename);
		void get_Cells(vector<shared_ptr<Cell>>& cells);
		//set/get the number of cells in the tissue
		void update_Num_Cells(shared_ptr<Cell>& new_Cell);
		int  get_num_cells() {return num_cells;}
		void update_Neighbor_Cells();
		void add_Wall(int Ti);
		void delete_Wall(int Ti);
		void update_Adhesion();
		//not in use
		void update_Linear_Bending_Springs();
		
		void update_Cell_Cycle(int Ti);
		void division_check();
		void calc_New_Forces(int Ti);
		void update_Cell_Locations();
		
		//stuff for data output
		void nematic_output(ofstream& ofs);
		void print_Data_Output(ofstream& ofs);
		void locations_output(ofstream& ofs);
		int update_VTK_Indices();
		void print_VTK_File(ofstream& ofs);
		
		//Destructor
		~Tissue();
};


//===========================
//End of file

#endif
