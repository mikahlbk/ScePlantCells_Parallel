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
const double dt = .00075;
const int Init_Num_Cyt_Nodes = 20;
const int Init_Wall_Nodes = 300;
//
const double AREA_THRESH = 64;
const double GROWTH_RATE = .000001;
//// Global Physics Constants
const double pi = acos(-1.0);
const int calibStart = 50000;
/////// Cell wall mechanical parameters
const double K_BEND = 2;
//
////Adhesion spring mechanical params
const double K_ADH = 20;
const double MembrEquLen_ADH = 0.6;
const double ADHThresh = .78125;
//
////Calibration parameters
const double EXTERNAL_FORCE = .001;
const double INIT_AREA = .018;
const double HIGH_GROWTH = 0.0000012;
const double LOW_GROWTH = 0.00000093;
////linear spring equilibrium length
const double MembrEquLen = .07; 
const double MEMBR_THRESH_LENGTH = 0.2; 
//
//
/////// Subcellular element parameters for membrane - membrane interactions
const double U_MM =  3.9;
const double W_MM =  0;
const double xsi_MM = 0.125;
const double gamma_MM = 1.5625;	
//
/////// Subcellular element parameters for membrane  - internal interactions
const double U_MI = 45;
const double W_MI = 0;
const double xsi_MI = .15;
const double gamma_MI = 1.55;
//
/////// Subcellular element parameters for internal - internal interactions
const double U_II = 75;
const double W_II = 6.71;
const double xsi_II = .4;
const double gamma_II = 1.34;
////const double ALPHA = 2.5;
////const double DELTA = .35;
////const double MORSE_EQ = 1.65;
////const double ALPHA_MI = 1.9;
////const double DELTA_MI = .15;
////const double MORSE_EQ_MI = 1.25;
////=====================
#endif
