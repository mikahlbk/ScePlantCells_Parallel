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
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
//===================
// Cell Class Declaration

class Cell: public enable_shared_from_this<Cell> {

	private:
		//shared_ptr<Tissue> my_tissue;
		Tissue* my_tissue;
		int rank;
		int layer;
		//int boundary;
		double damping;
		int life_length;
		int num_cyt_nodes;
		int num_wall_nodes;
		double Cell_Progress;
		Coord cell_center;
		//double cytokinin;
		double wuschel;
		int growth_rate;
		Coord growth_direction;
		vector<shared_ptr<Wall_Node>> wall_nodes;
		vector<shared_ptr<Cyt_Node>> cyt_nodes;
		vector<shared_ptr<Cell>> neigh_cells;
		shared_ptr<Wall_Node> left_Corner;	
	public:
		
		Cell(Tissue* tissue);
		Cell(int rank, Coord center, double radius, Tissue* tiss, int layer, int boundary);
		void make_nodes(double radius);
		
		// Destructor
		~Cell();

		// Getters and Setters
		//set tissue
		Tissue* get_Tissue() {return my_tissue;}
		//set/get rank
		void set_Rank(const int id);
		int get_Rank() {return rank;}
		//set/get layer
		void set_Layer(int layer);
		int get_Layer() {return layer;}
		//set/get damping
		void set_Damping(double new_damping);
		double get_Damping() {return damping;}
		//set/get life length
		void update_Life_Length();
		int get_life_length() {return life_length;} 
		//get total node count
		int get_Node_Count();
		//get counts of individual nodes
		int get_wall_count() {return num_wall_nodes;}
		int get_cyt_count() {return num_cyt_nodes;}
		//get wall nodes
		void get_Wall_Nodes_Vec(vector<shared_ptr<Wall_Node>>& walls);
		//add new wall node
		void add_wall_node_vec(shared_ptr<Wall_Node> curr);
		//get cyt nodes
		void get_Cyt_Nodes_Vec(vector<shared_ptr<Cyt_Node>>& cyts);
		//add new cyt node
		void update_cyt_node_vec(shared_ptr<Cyt_Node> new_node);
		//reset cell_Progress
		void reset_Cell_Progress();
		//set/get cell progress
		void update_Cell_Progress();
		double get_Cell_Progress() {return Cell_Progress;}
		//get cell center
		Coord get_Cell_Center() {return cell_center;}
		//set/get WUS conc
		void calc_WUS(int Ti);
		double get_WUS_concentration() {return wuschel;}
		//double get_CYT_concentration() {return cytokinin;}
		//set growth rate based on WUS
		void set_growth_rate();
		//set/get growth direction
		void set_growth_direction(Coord gd);
		Coord get_growth_direction(){return growth_direction;}
		//get current neighbor cells		
		void get_Neighbor_Cells(vector<shared_ptr<Cell>> cells);
		//set/get left corner
		void set_Left_Corner(shared_ptr<Wall_Node> new_left_corner);
		shared_ptr<Wall_Node> get_Left_Corner() {return left_Corner;}			      //is this necessary?
		void set_Wall_Count(int number_nodes);
		double compute_k_lin(shared_ptr<Wall_Node> current);
		//void calc_CYT();
	
		
		//Keep track of neighbor cells
		void update_Neighbor_Cells();
		void update_adhesion_springs();
		
		//Forces and Positionsing
		void calc_New_Forces(int Ti);
		void update_Node_Locations();
		void update_Wall_Angles();
		void update_Wall_Equi_Angles();
		void update_Cell_Center();
	
		//Growth of a cell
		void update_Cell_Progress(int& Ti);
		double calc_Area();
		void add_wall_Node_Check(int Ti);
		void delete_wall_Node_Check(int Ti);
		void add_Wall_Node(int Ti);
		void delete_Wall_Node(int Ti);
		
		void find_Smallest_Length(shared_ptr<Wall_Node>& right);
		void find_Largest_Length(shared_ptr<Wall_Node>& right);
		void find_Largest_Length_Div(shared_ptr<Wall_Node>& right, shared_ptr<Wall_Node>& second_right);
		void add_Cyt_Node();
		
		//Functions for Division
		double find_radius();
		void find_nodes_for_div_plane(Coord& direction,vector<shared_ptr<Wall_Node>>& nodes);
		void find_nodes_for_div_plane_anticlinal(Coord& direction,vector<shared_ptr<Wall_Node>>& nodes);
		void add_Cyt_Node_Div(double radius_x,double radius_y);
		void move_cyt_nodes();
	
		//Output Functions
		double average_Pressure();
		void wall_Pressure();
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices(int& id);
		void print_VTK_Adh(ofstream& ofs);
		void print_locations(ofstream& ofs);
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
		shared_ptr<Cell> divide();
		//void find_Largest_Length_Div(Wall_Node*& right_one, Wall_Node*& right_two);
		void division(shared_ptr<Cell>& new_cell);
	};


// End Cell Class
//===================

#endif

