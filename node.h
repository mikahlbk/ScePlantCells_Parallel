//node.h
//=====================
// Include Guards
#ifndef _NODE_H_INCLUDED_  //if node.h hasn't been included yet
#define _NODE_H_INCLUDED_  //    define it so the compiler knows 
//=====================
// Forward Declarations 
class Wall_Node;
class Cell;
class Cyt_Node;
//=====================
// Include Declarations
#include <iostream>
#include <vector>
#include <memory>
#include "phys.h"
#include "coord.h"
//=====================

class Node{
    protected:
    //variables that will be shared by all nodes
		Coord my_loc;
		Coord new_force;
		double damping;
		int vtk_id;
   public:
    //functions that you will want performed on all nodes
        //Constructor
        Node(Coord loc);
        //some functions you can define in base class because 
        //    all nodes will use the exact same function
        	virtual Coord get_Location() {return my_loc;}
		virtual Coord get_Force() {return new_force;}
		virtual void set_Damping(double new_damping);
		virtual void update_Location();
       		virtual void update_VTK_Id(int id);
		virtual int get_VTK_Id() {return vtk_id;}
	//other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        	virtual ~Node();
};

class Cyt_Node: public Node, public enable_shared_from_this<Cyt_Node>{

   protected: 
    //if don't need to keep any more information then just leave blank
		shared_ptr<Cell> my_cell;
    public:
        	//constructor
		Cyt_Node(Coord loc, shared_ptr<Cell> my_cell);
		void update_Cell(shared_ptr<Cell> cell);
		void new_location(Coord location); 
		shared_ptr<Cell> get_My_Cell(){return my_cell;}
		void calc_Forces(int Ti);
		Coord calc_Morse_II(int Ti);
		Coord calc_Morse_MI(shared_ptr<Wall_Node> orig, int Ti);
		Coord morse_Equation(shared_ptr<Cyt_Node> cyt, int Ti);
		Coord morse_Equation(shared_ptr<Wall_Node> wall, int Ti);
		~Cyt_Node();
};

class Wall_Node: public Node, public enable_shared_from_this<Wall_Node> {
    protected:
    //variables that will be shared by all wall nodes
       		shared_ptr<Wall_Node> left;
        	shared_ptr<Wall_Node> right;
		shared_ptr<Cell> my_cell;
		double membr_equ_len;
		double my_angle;
		double K_LINEAR;
		double equi_angle;
		double cross_Prod;
		Coord cyt_force;
		
		//variables used for adhesion
		shared_ptr<Wall_Node> curr_closest;
		shared_ptr<Wall_Node> closest;
		double closest_len;
		vector<shared_ptr<Wall_Node>> adhesion_pairs;

    public:
    //function that you want performed on all wall nodes
		
    		// Constructors
       		Wall_Node(Coord loc, shared_ptr<Cell> my_cell);
        	Wall_Node(Coord loc, shared_ptr<Cell> my_cell, shared_ptr<Wall_Node> left, shared_ptr<Wall_Node> right);
        	
		// Getters and Setters
		//set left/right neighbor
		void set_Left_Neighbor(shared_ptr<Wall_Node> new_Left);
		void set_Right_Neighbor(shared_ptr<Wall_Node> new_Right);
		//get left/right neighbor
		shared_ptr<Wall_Node> get_Left_Neighbor() {return left;}
		shared_ptr<Wall_Node> get_Right_Neighbor() {return right;}
		//set/get cell
		void update_Cell(shared_ptr<Cell> new_cell);
		shared_ptr<Cell> get_My_Cell() {return my_cell;}
		//set/get membrane length
		void set_membr_len(double length);
		double get_membr_len(){return membr_equ_len;}
		//set/get angle
		void update_Angle();
		double get_Angle() {return my_angle;}
		//set/get K_linear
		void set_K_LINEAR(double k_lin);
		double get_k_lin(){return K_LINEAR;}
		//set/get equi angle
		void update_Equi_Angle(double new_theta);
		double get_Equi_Angle() {return equi_angle;}
		//crossprod necessary?
		Coord get_CytForce() {return cyt_force;}
	
		//Adhesion functions
		shared_ptr<Wall_Node> find_Closest_Node(vector<shared_ptr<Cell>> neighbors);
		void set_curr_Closest(shared_ptr<Wall_Node> curr_closest);
		shared_ptr<Wall_Node> get_curr_Closest(){return curr_closest;}
		void make_Connection(shared_ptr<Wall_Node> curr_Closest);
		void set_Closest(shared_ptr<Wall_Node> closest, double closest_len);
		shared_ptr<Wall_Node> get_Closest(){return closest;}
		double get_Closest_Len() {return closest_len;}
		void add_adh_pair(shared_ptr<Wall_Node> pair);
		vector<shared_ptr<Wall_Node>>get_adhesion_vec(){return adhesion_pairs;}
      		void clear_adh_vec();
		
		//functions for calculating forces
		void calc_Forces(int Ti);
		Coord calc_Morse_SC(int Ti);
		Coord calc_Linear();
		Coord calc_Bending();
		Coord calc_Morse_DC(int Ti);
		Coord neighbor_nodes(shared_ptr<Cell> neighbor,int Ti);
		
		// Mathematical Force Equations
		Coord morse_Equation(shared_ptr<Cyt_Node> cyt, int Ti);
		Coord morse_Equation(shared_ptr<Wall_Node> wall, int Ti);
		Coord bending_Equation_Center();
		Coord bending_Equation_Left();
		Coord bending_Equation_Right();
		Coord linear_Equation(shared_ptr<Wall_Node> wall);
		Coord linear_Equation_ADH(shared_ptr<Wall_Node>& wall);
		~Wall_Node();

};
//===========================
#endif  
