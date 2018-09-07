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
#include <memory>
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

	if (argc != 3) {
	
	}
	
	// reads in name of folder that 
	// stores vtk files, given in run.sh
	string anim_folder = argv[1];
	//reads in name of folder that 
	//stores data output , given in run.sh
	string nematic_folder = argv[2];

	int start = clock();	
	//.txt file that tells how cells start
	string init_tissue = "cell_maker.txt";
	//cout << "Read in cell starter" << endl;	
	//make new cell objects in tissue
	Tissue growing_Tissue(init_tissue);
	//Tissue* growing_Tissue= new Tissue(init_tissue);
	//cout << "Finished creating Cells" << endl;
	
	//parameters for time step
	double numSteps = 500;
	
	// Variable for dataoutput
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Plant_Cell_";
	string Filename;
	ofstream ofs_anim;
	int out = 0;

	ofstream ofs_nem;
	int Number2 =0;
	string nem_Filename;
	string nem_initial = "/Nematic_";
	
	//growing_Tissue.update_Neighbor_Cells();
		
	//loop for time steps
	for(int Ti = 0; Ti*dt< numSteps; Ti++) {
		// print to dataOutput and vtk files
		if(Ti%5000 ==0) {
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

			ofs_anim.open(Filename.c_str());
			growing_Tissue.print_VTK_File(ofs_anim);
			ofs_anim.close();	
			out++;
		}

			
		//loop through all cells
		if (Ti % 10000 == 0) {
			cout << "Simulation still running. Ti: " << Ti << endl;
		}
		// Tissue Growth
		//fills vector of neighbor cells for each cell
		if (Ti%10000 == 0) {
			//cout << "Find Neighbors" << endl;
			growing_Tissue.update_Neighbor_Cells();
		}	
		
		//cout << "add new cell wall nodes if needed" << endl;
		//adds one new cell wall node in the biggest gap
		if(Ti%1000 == 0){	
			//cout << "add wall " << endl;
			growing_Tissue.add_Wall(Ti);
		}
		//deletes a cell wall node if too close together
		//if(Ti%1000 == 0){	
			//cout << "delete wall" << endl;
			growing_Tissue.delete_Wall(Ti);
		//}
		//matches wall nodes with adhesion pairs
		if(Ti < 2500) {
		if(Ti%1000 == 0) {
			//cout << "adhesion"<< endl;
			growing_Tissue.update_Adhesion();
		}
		}
		else {
		if(Ti%10000 == 0) {
			//cout << "adhesion"<< endl;
			growing_Tissue.update_Adhesion();
		}
		}


		//cout << "cell cycle" << endl;
		growing_Tissue.update_Cell_Cycle(Ti);

		
		//Calculate new forces on cells and nodes
		//cout << "forces" << endl;
		growing_Tissue.calc_New_Forces(Ti);
	
		//cout << "locations" << endl;
		//Update node positions
		growing_Tissue.update_Cell_Locations();	
		//cout << "Finished" << endl;
		
		if(Ti%5000 == 1){
			nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
			ofs_nem.open(nem_Filename.c_str());
			growing_Tissue.nematic_output(ofs_nem);
			ofs_nem.close();
			Number2++;
		}

}
	nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
	ofs_nem.open(nem_Filename.c_str());
	growing_Tissue.nematic_output(ofs_nem);
	ofs_nem.close();
	
	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;

	return 0;
}











