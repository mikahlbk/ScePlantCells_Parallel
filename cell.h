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
		Tissue* my_tissue;
		int rank;
		int layer;
		double damping;
		int life_length;
		int growth_rate;
		int num_cyt_nodes;
		int num_wall_nodes;
		double Cell_Progress;
		double Cell_Progress_add_node;
		int Cell_Progress_div;
		Coord cell_center;
		double cytokinin;
		double wuschel;
		double total_signal;
		Coord K_LINEAR;
		vector<Wall_Node*> wall_nodes;
		vector<Cyt_Node*> cyt_nodes;
		vector<Cell*> neigh_cells;
		bool is_deleted;
		Wall_Node* left_Corner;	
	public:
		
		Cell(Tissue* tissue);
		Cell(int rank, Coord center, double radius, Tissue* tiss, int layer);

		// Destructor
		~Cell();

		// Getters and Setters
		Tissue* get_Tissue() {return my_tissue;}
		int get_Rank() {return rank;}
		void set_Rank(const int id);
		int get_Layer() {return layer;}
		void set_Layer(int layer);
		double get_Damping() {return damping;}
		void set_Damping(double& new_damping);
		int get_life_length() {return life_length;}
		void update_Life_Length();
		int get_Node_Count();
		int get_wall_count() {return num_wall_nodes;}
		int get_cyt_count() {return num_cyt_nodes;}
		bool return_is_deleted() {return is_deleted;}
		void get_Wall_Nodes_Vec(vector<Wall_Node*>& walls);
		void add_wall_node_vec(Wall_Node* curr);
		void get_Cyt_Nodes_Vec(vector<Cyt_Node*>& cyts);
		double get_Cell_Progress() {return Cell_Progress;}
		void reset_Cell_Progress();
		double get_Cell_Progress_Add_Node() {return Cell_Progress_add_node;}
		void update_Cell_Progress_add_node(double& add_node_prog);
		int get_Cell_Progress_div() {return Cell_Progress_div;}
		void update_Cell_Progress_div(int& div_time);	
		Coord get_Cell_Center() {return cell_center;}
		double get_WUS_concentration() {return wuschel;}
		double get_CYT_concentration() {return cytokinin;}
		double get_total_concentration() {return total_signal;}
		void set_growth_rate();
		Coord get_K_LINEAR() {return K_LINEAR;}
		void set_K_LINEAR(double& x, double& y);
		Wall_Node* get_Wall_Nodes() {return left_Corner;}
		Wall_Node* get_Left_Corner() {return left_Corner;}		
		void get_Neighbor_Cells(vector<Cell*>& cells);
		void set_Left_Corner(Wall_Node*& new_left_corner);
		void set_Wall_Count(int& number_nodes);
		void calc_WUS();
		void calc_CYT();
		void calc_Total_Signal();

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
		void update_Cell_Progress(int& Ti);
		double calc_Area();
		void add_wall_Node_Check();
		void delete_wall_Node_Check();
		void add_Wall_Node();
		void delete_Wall_Node();
		void compute_Main_Strain_Direction(double& x_length, double& y_length); 
		Wall_Node* find_closest_node_top(); 
		Wall_Node* find_closest_node_bottom();
		Wall_Node* find_closest_node_left();
		Wall_Node* find_closest_node_right();
		void find_Smallest_Length(Wall_Node*& right);
		void find_Largest_Length(Wall_Node*& right);
		void find_Largest_Length_Div(Wall_Node*& right, Wall_Node*& second_right);
		void add_Cyt_Node();
		
		//Functions for Division
		double find_radius();
		void add_Cyt_Node_Div(double radius_x,double radius_y);
		void stress_Tensor_Eigenvalues(double& a, double& b, double& c, double& d, vector<double>& eigen_Max);	
		double compute_Stress_Tensor_XY();
		double compute_Stress_Tensor_Y();
		double compute_Stress_Tensor_X();
		
		
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
		
		//Division 
		Cell* divide();
		//void find_Largest_Length_Div(Wall_Node*& right_one, Wall_Node*& right_two);
		Cell* division();
	};


// End Cell Class
//===================

#endif

