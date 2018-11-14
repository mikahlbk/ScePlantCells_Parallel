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
//
//// Simulation Constants
const double dt = .0003;
const int Init_Num_Cyt_Nodes = 15;
const int Init_Wall_Nodes = 150;
const double pi = acos(-1.0);
const double pi_fourths = acos(sqrt(2)/2.0);
const double three_pi_fourths = acos(-sqrt(2)/2.0);
////// Cell wall mechanical parameters
const double K_BEND_STIFF = 12;
const double K_BEND_LOOSE =12;
const double K_LINEAR_STIFF = 500;
const double K_LINEAR_LOOSE = 150;
////Adhesion spring mechanical params
const double K_ADH = 0;
const double K_ADH_L1 = 0;
const double MembrEquLen_ADH = .8;
const double ADHThresh = 2;
/////Microfibril spring mechanical params
const double K_microfibril = 15;
const double MembrEquLen_microfibril = 8;

//
////linear spring equilibrium length
const double MembrEquLen = .07;
const double MEMBR_THRESH_LENGTH = 0.15; 
//
/////// Subcellular element parameters for membrane - membrane interactions
const double U_MM = 100;
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
const double U_II = 120;
const double W_II = 6.71;
const double xsi_II = .8;
const double gamma_II = 1.34;
const double xsi_II_div = .4;
////=====================
#endif
