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
	//stores data output, given in run.sh
	string nematic_folder = argv[2];
	//this is a folder that holds node
	//locations in the fashion that weitao 
	//asked for to couple the models
	string locations_folder = argv[3];
	//keep track of time
	int start = clock();	
	//.txt file that tells configuration 
	//that cells start in
	string init_tissue = "cell_maker.txt";
	
	//cout << "Read in cell starter" << endl;	
	//make new cell objects in tissue
	Tissue growing_Tissue(init_tissue);
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
	
	ofstream ofs_loc;
	string locations_Filename;
	string locations_initial = "/Locations_";
	int Number3 = 0;
	
	//loop for time steps
	//which matlab file tells you how many
	//seconds each time step reprsents (2.5?)
	for(int Ti = 0; Ti*dt< numSteps; Ti++) {
		// print to dataOutput and vtk files
		if(Ti%1000 ==0) {
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

			
		//keep track of simulation runs
		if (Ti % 1000 == 0) {
			cout << "Simulation still running. Ti: " << Ti << endl;
		}
	
		// Tissue Growth

		//fills vector of neighbor cells for each cell
		//in tissue class this goes through each cell and calls
		//updated neighbor on each cell
	
		if(Ti==0) {
			cout << "Find Neighbors" << endl;
			growing_Tissue.update_Neighbor_Cells();
		}	
		
		//adds one new cell wall node per cell everytime it is called
		//dont call it right away to give cell time to find initial configuration
	
		if(Ti%2000 == 0){	
			cout << "add new cell wall nodes if needed" << endl;
			growing_Tissue.add_Wall(Ti);
		}
	
		
		//deletes a cell wall node if too close together
		//if(Ti%1000 == 0){	
			//cout << "delete wall" << endl;
			//growing_Tissue.delete_Wall(Ti);
		//}
	
		if(Ti < 999){
			if(Ti%200 == 0) {
				growing_Tissue.update_Adhesion();
			}
		}
		else{	
			if(Ti%1000 == 0) {
				cout << "adhesion"<< endl;
				growing_Tissue.update_Adhesion();
			}
		}
		
		cout << "cell cycle" << endl;
		//this is where the cell decides if 
		//it will divide
		growing_Tissue.update_Cell_Cycle(Ti);
		//cout << "divide necessary cells" << endl;
		growing_Tissue.division_check();
		
		//Calculate new forces on cells and nodes
		cout << "forces" << endl;
		growing_Tissue.calc_New_Forces(Ti);
	
		cout << "locations" << endl;
		//Update node positions
		growing_Tissue.update_Cell_Locations();	
		
		cout << "Finished" << endl;
		
		//data output from simulations
		//for cell center etc
		/*if(Ti%5000 == 1){
			nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
			ofs_nem.open(nem_Filename.c_str());
			growing_Tissue.nematic_output(ofs_nem);
			ofs_nem.close();
			Number2++;
		}
		if(Ti%1000 == 1){
			locations_Filename = locations_folder + locations_initial + to_string(Number3) + ".txt";
			ofs_loc.open(locations_Filename.c_str());
			growing_Tissue.locations_output(ofs_loc);
			ofs_loc.close();
			Number3++;
		}*/

}
	/*nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
	ofs_nem.open(nem_Filename.c_str());
	growing_Tissue.nematic_output(ofs_nem);
	ofs_nem.close();*/
	
	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;

	return 0;
}











