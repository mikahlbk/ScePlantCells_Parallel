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
		Coord cell_center;
		double cytokinin;
		double wuschel;
		Coord growth_direction;
		vector<Wall_Node*> wall_nodes;
		vector<Cyt_Node*> cyt_nodes;
		vector<Cell*> neigh_cells;
		//vector<Coord> nematic_vec;
		//vector<double> angle_vec;
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
		void update_Life_Length();
		int get_life_length() {return life_length;} 
		int get_Node_Count();
		void get_Wall_Nodes_Vec(vector<Wall_Node*>& walls);
		void add_wall_node_vec(Wall_Node* curr);
		void get_Cyt_Nodes_Vec(vector<Cyt_Node*>& cyts);
		void reset_Cell_Progress(int cyt_size);
		void update_cyt_node_vec(Cyt_Node* new_node);
		int get_wall_count() {return num_wall_nodes;}
		int get_cyt_count() {return num_cyt_nodes;}
		double get_Cell_Progress() {return Cell_Progress;}
		Coord get_Cell_Center() {return cell_center;}
		double get_WUS_concentration() {return wuschel;}
		double get_CYT_concentration() {return cytokinin;}
		void set_growth_rate();
		void set_growth_direction(Coord gd);
		Coord get_growth_direction(){return growth_direction;}
		Wall_Node* get_Wall_Nodes() {return left_Corner;}
		Wall_Node* get_Left_Corner() {return left_Corner;}		
		void move_cyt_nodes();
		void get_Neighbor_Cells(vector<Cell*>& cells);
		void set_Left_Corner(Wall_Node*& new_left_corner);
		void set_Wall_Count(int& number_nodes);
		void calc_WUS();
		void calc_WUSwildtype();
		void calc_WUSBAP12hr();
		void calc_WUSBAP24hr();
		void calc_CYT();
		void calc_Total_Signal();

		// Keep track of neighbor cells
		void update_Neighbor_Cells();
	//	void update_adhesion_springs();
	//	void update_microfibril_springs();
	//	void make_top_bottom_vectors(vector<Wall_Node*>&top,vector<Wall_Node*>&bottom);
	//	void make_left_right_vectors(vector<Wall_Node*>&left,vector<Wall_Node*>&right);
		// Forces and Positionsing
		void calc_New_Forces(int Ti);
		void update_Node_Locations();
		void update_Wall_Angles();
		void update_Wall_Equi_Angles();
		void update_Cell_Center();
	
		//Growth of a cell
		void update_Cell_Progress(int& Ti);
		double calc_Area();
		void nematic(Coord& avg_vec, double& angle);
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
		void find_nodes_for_div_plane(Coord& direction,vector<Wall_Node*>& nodes);
		//void find_nodes_for_div_plane_anticlinal(Coord& direction,vector<Wall_Node*>& nodes);
		void add_Cyt_Node_Div(double radius_x,double radius_y);
		void stress_Tensor_Eigenvalues(double& a, double& b, double& c, double& d, vector<double>& eigen_Max);	
		double compute_Stress_Tensor_XY();
		double compute_Stress_Tensor_Y();
		double compute_Stress_Tensor_X();
		
		
		//Output Functions
		double average_Pressure();
		void wall_Pressure();
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices(int& id);
		void print_VTK_Adh(ofstream& ofs);
		void print_VTK_Points(ofstream& ofs, int& count);
		void print_VTK_Scalars_Wall_Pressure(ofstream& ofs);
		void print_VTK_Scalars_Average_Pressure(ofstream& ofs);
		void print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs);
		void print_VTK_Scalars_WUS(ofstream& ofs);
		void print_VTK_Scalars_WUS_cell(ofstream& ofs);
		void print_VTK_Scalars_CYT(ofstream& ofs);
		void print_VTK_Scalars_Total(ofstream& ofs);
		void print_VTK_Vectors(ofstream& ofs);
		void print_VTK_Scalars_Node(ofstream& ofs);	
		//Division 
		Cell* divide();
		//void find_Largest_Length_Div(Wall_Node*& right_one, Wall_Node*& right_two);
		Cell* division();
	};


// End Cell Class
//===================

#endif

