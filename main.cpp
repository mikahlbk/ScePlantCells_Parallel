// Main.cpp
//============================

//===========================
// Include Dependencies
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <omp.h>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
#include "rand.h"
//==============================

using namespace std;

//============================

int main(int argc, char* argv[]) {


	string anim_folder = argv[1];

	int start = clock();
	
	string init_tissue = "cell_maker.txt";
	
	//make new cell objects in tissue
	Tissue growing_Tissue(init_tissue);
	cout << "Finished creating Cells" << endl;
	//parameters for time step
	double numSteps = 500;
	seed();
	// Variable for dataoutput
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Plant_Cell_";
	string Filename;
	ofstream ofs;
	int out = 0; //counter for creating output/vtk files

	//loop for time steps
	for(int Ti = 0; Ti*dt< numSteps; Ti++) {
		//loop through all cells
		//for now only one cell
		//cout << "Ti = " << Ti << endl;
		//Print to dataOutput and VTK files
		if (Ti % 1000 == 0) {
			
			digits = ceil(log10(out + 1));
			if (digits == 1 || digits == 0) {
				Number = "0000" + to_string(out);
			}	
			else if (digits == 2) {
				Number = "000" + to_string(out);
			}	
			else if (digits == 3) {
				Number = "00" + to_string(out);
			}
			else if (digits == 4) {
				Number = "0" + to_string(out);
			}

			Filename = anim_folder+ initial + Number + format;

			ofs.open(Filename.c_str());
			growing_Tissue.print_VTK_File(ofs);
			ofs.close();	
			out++;
		}
		
		if (Ti % 1000 == 0) {
			cout << "Simulation still running. Ti: " << Ti << endl;
		}
		
	
		// Tissue Growth
		
		//cout << "update cell cycle of each cell" << endl;
		//this includes a check for division and addition
		//of new internal nodes according to growth rate
	//	if(Ti > 2500){
			growing_Tissue.update_Cell_Cycle(Ti);
	//	}
		//cout << ":wq
		//updated cell cycle" << endl;	
		//cout << "add new cell wall nodes if needed" << endl;
		//adds one new cell wall node in the biggest gap
		
		//if(Ti< 2000) {
		//	growing_Tissue.update_Wall();
		//}
	//	if(Ti%200 ==0) {
		//if(Ti%200==0) {
			growing_Tissue.update_Wall();
		//}	
	//	}
		if (Ti%200 == 0) {
			//cout << "Find Neighbors" << endl;
			growing_Tissue.update_Neighbor_Cells();
		}
		//if(Ti< 30000) {
			//cout << "Finding adhesion points" << endl;
		//if(Ti%200 ==199) {
		//	growing_Tissue.update_Adhesion();
		//}
	//	}
	//	if(Ti>= 30000) {
		if (Ti%200 == 0) {
			growing_Tissue.update_Adhesion();
		}
	//	}
		
		//if(Ti%1000 == 999) {
		//	growing_Tissue.add_cyt_node();
		//}
		//Calculate new forces on cells and nodes
		//cout << "forces" << endl;
		growing_Tissue.calc_New_Forces(Ti);
	
		//cout << "locations" << endl;
		//Update node positions
		growing_Tissue.update_Cell_Locations();
		//cout << "Finished" << endl;
		//growing_Tissue.cell_area();
}
	
	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;

	return 0;
}











