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
//
////// Cell wall mechanical parameters
const double K_BEND = 10;
//const double K_BEND_L1 =2;
//
////Adhesion spring mechanical params
const double K_ADH = 20;
const double K_ADH_L1 = 20;
const double MembrEquLen_ADH = .8;
const double ADHThresh = 2.5;
//
////linear spring equilibrium length
const double MembrEquLen = .07;
const double MEMBR_THRESH_LENGTH = 0.25; 
//
/////// Subcellular element parameters for membrane - membrane interactions
const double U_MM =  5;
const double W_MM =  0;
const double xsi_MM = .8;
const double gamma_MM = 1.5625;	
//
/////// Subcellular element parameters for membrane  - internal interactions
const double U_MI = 25;
const double W_MI = 0;
const double xsi_MI = 1.1;
const double gamma_MI = 1.34;
//
/////// Subcellular element parameters for internal - internal interactions
const double U_II = 50;
const double W_II = 0;
const double xsi_II = 1.1;
const double gamma_II = 1.34;
////=====================
#endif
