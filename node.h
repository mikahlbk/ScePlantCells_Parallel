//node.h
//=====================
// Include Guards
#ifndef _NODE_H_INCLUDED_  //if node.h hasn't been included yet
#define _NODE_H_INCLUDED_  //    define it so the compiler knows 
//=====================
// Forward Declarations 
class Wall_Node;
class Cell;
//=====================
// Include Declarations
#include <iostream>
#include <vector>
#include "phys.h"
#include "coord.h"
//=====================

class Node {
    protected:
    //variables that will be shared by all nodes
		Coord my_loc;
		Coord new_force;
		int vtk_id;
   public:
    //functions that you will want performed on all nodes
        //Constructor
        Node(Coord loc);
        //some functions you can define in base class because 
        //    all nodes will use the exact same function
        	virtual Coord get_Location();
		virtual Coord get_Force();
		virtual void update_Location(double& new_damping);
       		virtual void update_VTK_Id(int id);
		virtual int get_VTK_Id() {return vtk_id;}
	//other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        	//virtual void calc_Forces();
		
		virtual ~Node();
};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank
		Cell* my_cell;
    public:
        	//constructor
		Cyt_Node(Coord loc, Cell* my_cell);
		void update_Cell(Cell* cell);
		void new_location(Coord location); 
		Cell* get_My_Cell(){return my_cell;}
		void calc_Forces(int Ti);
		Coord calc_Morse_II(int Ti);
		Coord calc_Morse_MI(Wall_Node* orig, int Ti);
		Coord morse_Equation(Cyt_Node* cyt, int Ti);
		Coord morse_Equation(Wall_Node* wall, int Ti);
		~Cyt_Node();
};

class Wall_Node: public Node {
    protected:
    //variables that will be shared by all wall nodes
       		Wall_Node* left;
        	Wall_Node* right;
		Cell* my_cell;
		double membr_equ_len;
		//double pressure;
		//bool is_new;
		double my_angle;
		double K_LINEAR;
		double equi_angle;
		double cross_Prod;
		Coord cyt_force;
		//Wall_Node* closest;
		//double closest_len;
		//Wall_Node* microfibril_pair;
		//double curr_slope;
		bool is_added;
    public:
    //function that you want performed on all wall nodes
		// Constructors
       		Wall_Node(Coord loc, Cell* my_cell);
        	Wall_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right);
        // Getters and Setters
		void set_is_new();
		bool get_status() {return is_added;}
		Wall_Node* get_Left_Neighbor() {return left;}
		Wall_Node* get_Right_Neighbor() {return right;}
		double get_Angle() {return my_angle;}
		double get_Equi_Angle() {return equi_angle;}
		double get_k_lin(){return K_LINEAR;}
		double get_membr_len(){return membr_equ_len;}
		Cell* get_My_Cell() {return my_cell;}
		Coord get_CytForce() {return cyt_force;}
      				
	//	bool get_color() {return is_new;}
		void set_K_LINEAR(double& k_lin);
		void set_Left_Neighbor(Wall_Node* new_Left);
		void set_Right_Neighbor(Wall_Node* new_Right);
		//void set_Delete(int y);
		//void set_is_new(bool yes);
		void set_membr_len(double length);
		void update_Angle();
		void update_Equi_Angle(double new_theta);
		void update_Cell(Cell* new_cell);
		//Wall_Node* get_Closest(){return closest;}
		//double get_Closest_Len() {return closest_len;}
		//Wall_Node* find_Closest_Node(vector<Cell*>& neighbors);
		//void make_Connection(Wall_Node* curr_Closest);
		//Coord get_Ext_Force() {return f_EXT;}
		//void set_Closest(Wall_Node* closest, double closest_len);
		//void set_Closest_Vec(Wall_Node* closest);
		//void clear_Closest_Vec();
		//void find_microfibril_pair_horiz(vector<Wall_Node*>side2);
		//void find_microfibril_pair_vert(vector<Wall_Node*>side2);
		//void set_microfibril_pair(Wall_Node* pair,double curr_slope);
		//Wall_Node* get_micro_pair(){return microfibril_pair;}
	//	double get_pressure(){return pressure;}
	//	void set_pressure(double& new_press);
		//Force Calculations
		void calc_Forces(int Ti);
		Coord calc_Morse_SC(int Ti);
		Coord calc_Morse_DC(int Ti);
		Coord neighbor_nodes(Cell* neighbor,int Ti);
		Coord calc_Bending();
		Coord calc_Linear();
			
		// Mathematical Force Equations
		Coord morse_Equation(Cyt_Node* cyt, int Ti);
		Coord morse_Equation(Wall_Node* wall, int Ti);
		Coord bending_Equation_Center();
		Coord bending_Equation_Left();
		Coord bending_Equation_Right();
		Coord linear_Equation(Wall_Node* wall);
//		Coord linear_Equation_ADH(Wall_Node*& wall);
//		Coord linear_Equation_microfibril(Wall_Node*& wall);
		~Wall_Node();

};
//===========================
#endif  
