//cell.h
//===================
// Inlcude Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_
//===================
// forward declarations
class Tissue;
//===================
// include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "phys.h"
#include "coord.h"
#include "node.h"
//===================
// Cell Class Declaration

class Cell {

	private:
		int rank;
		double damping;
		Tissue* my_tissue;
		double Cell_Progress;
		double Cell_Progress_add_node;
		int Cell_Progress_div;
		Coord K_LINEAR;
		int counter;
		int num_cyt_nodes;
		int layer;
		int life_length;
		Coord cell_center;
		double cytokinin;
		double wuschel;
		double total_signal;
		int num_wall_nodes;
		vector<double> strain_vec;
		vector<double> stress_vec;
		double curr_area;
		Wall_Node* top;
		Wall_Node* bottom;
		int counter_left;
		int counter_right;
		vector<Cyt_Node*> cyt_nodes;
		vector<Cell*> neigh_cells;
		Wall_Node* left_Corner;	
	public:
		
		Cell(Tissue* tissue);
		Cell(int rank, Coord center, double radius, Tissue* tiss, int layer);

		// Destructor
		~Cell();

		// Getters and Setters
		int get_Rank() {return rank;}
		double get_Damping() {return damping;}
		Tissue* get_Tissue() {return my_tissue;}
		double get_Cell_Progress() {return Cell_Progress;}
		double get_Cell_Progress_Add_Node() {return Cell_Progress_add_node;}
		Coord get_K_LINEAR() {return K_LINEAR;}
		int get_Cytoplasm_Count() {return num_cyt_nodes;}
		int get_Layer() {return layer;}
		int get_life_length() {return life_length;}
		Coord get_Cell_Center() {return cell_center;}
		double get_WUS_concentration() {return wuschel;}
		double get_CYT_concentration() {return cytokinin;}
		double get_total_concentration() {return total_signal;}
		int get_Wall_Count() {return num_wall_nodes;}
		Wall_Node* get_Wall_Nodes() {return left_Corner;}
		Wall_Node* get_Left_Corner() {return left_Corner;}		
		int get_Node_Count();
		void get_Cyt_Nodes(vector<Cyt_Node*>& cyts);
		void get_Neighbor_Cells(vector<Cell*>& cells);
		void set_div_time(int& Ti);	
		void set_K_LINEAR(double& x, double& y);
		void set_Damping(double& new_damping);
		void set_Rank(const int id);
		void set_Layer(int layer);
		void reset_Cell_Progress();
		void set_Left_Corner(Wall_Node*& new_left_corner);
		void set_Wall_Count(int& number_nodes);
		
		void update_Cell_Progress_add_node(int& time);
		void update_Life_Length();
	
		// Keep track of neighbor cells
		void update_Neighbor_Cells();
		void update_adhesion_springs();
	
		// Forces and Positionsing
		void calc_New_Forces(int Ti);
		void update_Node_Locations();
		void update_Wall_Angles();
		void update_Wall_Equi_Angles();
		void update_Cell_Center();
	
		//Growth of a cell
		void update_Cell_Progress(int Ti);
		double calc_Area();
		void wall_Node_Check();
		void add_Wall_Node();
		void find_Largest_Length(Wall_Node*& right);
		void find_Largest_Length_Div(Wall_Node*& right_one, Wall_Node*& right_two);
		void add_Cyt_Node();
		
		//Functions for Division
		double find_radius();
		void add_Cyt_Node_Div(double radius);
		void stress_Tensor_Eigenvalues(double& a, double& b, double& c, double& d, vector<double>& eigen_Max);
		double compute_Stress_Tensor_XY();
		double compute_Stress_Tensor_Y();
		double compute_Stress_Tensor_X();
		//Functions for calibration
		void get_Strain(vector<double>& strain);
		void get_Stress(vector<double>& stress);
		void set_isStationary();
		void get_stretch_nodes();
		void compress();
	//	void extensional_strain();
	//	void tensile_Stress();	
		//Output Functions
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices(int& id);
		void print_VTK_Adh(ofstream& ofs);
		void print_VTK_Points(ofstream& ofs, int& count);
		void print_VTK_Scalars_Force(ofstream& ofs);
		void print_VTK_Scalars_WUS(ofstream& ofs);
		void print_VTK_Scalars_CYT(ofstream& ofs);
		void print_VTK_Scalars_Total(ofstream& ofs);
		void print_VTK_Vectors(ofstream& ofs);
		
		//Not in use
	//	void most_Left_Right(Wall_Node*& left, Wall_Node*& right);
		Wall_Node* closest_node_top();
		Wall_Node* closest_node_bottom();
	//	void closest_node_left(Wall_Node*& left);
	//	void closest_node_right(Wall_Node*& right);
	//	void closest_node(Wall_Node*& closest);
	
		void set_Stationary_Points(int Ti);
		double compute_pressure();	
		void calc_WUS();
		void calc_CYT();
		void calc_Total_Signal();
	
		//Division 
		Cell* divide();
		Cell* divide_length_wise();
		Cell* divide_width_wise();
	};


// End Cell Class
//===================

#endif

