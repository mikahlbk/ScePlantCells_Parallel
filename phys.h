//phys.h
#ifndef _PHYS_H_INCLUDED_
#define _PHYS_H_INCLUDED_
//===================
//// Forward declarations
//
////====================
//// Include dependencies
#include <math.h>
#include "coord.h"
////=====================
//// Simulation Constants
const double dt = .0003;
const int Init_Num_Cyt_Nodes = 20;
const int Init_Wall_Nodes = 100;
const double pi = acos(-1.0);
//these are used to control equi angles
//and bending spring constants
const double ANGLE_FIRST_QUAD = 0.785398;
const double ANGLE_SECOND_QUAD = 2.35619;
////// Cell wall mechanical parameters
const double K_BEND_STIFF = 60;
const double K_BEND_LOOSE =12;
//K_LINEAR_STIFF not used
const double K_LINEAR_STIFF = 0;
const double K_LINEAR_LOOSE = 50;
////Adhesion spring mechanical params
const double K_ADH = 35;
const double K_ADH_L1 = 35;
const double MembrEquLen_ADH = .8;
const double ADHThresh = 2.5;
//equilibrium length of linear springs
//add node threshold length not being used
const double MEMBR_THRESH_LENGTH = .7; 
const double Membr_Equi_Len_Long = .07;
const double Membr_Equi_Len_Short = .07;
//
/////// Subcellular element parameters for membrane - membrane interactions
const double U_MM = 15;
const double W_MM =  0;
const double xsi_MM = 0.5;
const double gamma_MM = 1.5625;	
//
/////// Subcellular element parameters for membrane  - internal interactions
const double U_MI = 75;
const double W_MI = 0;
const double xsi_MI = .3;
const double gamma_MI = 1.34;
const double xsi_MI_div = .2;
//
/////// Subcellular element parameters for internal - internal interactions
const double U_II = 95;
const double W_II = 6.71;
const double xsi_II = .8;
const double gamma_II = 1.34;
const double xsi_II_div = .4;
////=====================
#endif
