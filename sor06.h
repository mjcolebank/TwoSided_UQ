/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.7 2005/07/07 22:19:51 heine Exp $
// Last updated on October 23, 2014 by M. Umar Qureshi

#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves, N_aorta;               // Number of vessels in the tree.

int    tmstps = 8192,                 // The number of timesteps per period.
	   plts   = 512;                 // Number of plots per period.

// const char*  CO_filename = "Qin.dat";   // Input flow file at the heart.
// //const char*  CO_filename = "p0_in4096.dat"; //Changeed to "p0_in4096.dat" if using pressure as an input condition.
// const char*  PL_filename = "LA_pressure.dat";       // Outflow (Static pressure) file.

const char*  CO_filename;   // Input flow file at the heart.
const char*  PL_filename;       // Outflow (Static pressure) file.
     double vel_power = 9.0; // Exponent for the power law profile (2 poiseuille, 9 plug flow)

double conv   = 1333.220,             // Conversion from mmHg to SI-units.
       rho    = 1.055,                // Density of blood [g/cm^3].
       mu     = 0.032,                // Viscosity of blood [g/cm/s]. WITH DIMENSIONS
       nu     = mu/rho,             // Kinematic viscosity [cm^2/s] WITH DIMENSIONS
       Tper   = 0.85,//0.9,                 // The period of one heart beat [s].

       
       Fcst   = 2*0.07*sqrt(2*M_PI/nu/Tper),//50, //17.7778,        // Determines the damping coeff. Increases the peak pressure only
                                      // for the friction.
       Lr     = 1.0,                  // characteristic radius of the
                                      // vessels in the tree [cm].
       Lr2    = sq(Lr),               // The squared radius [cm2].
       Lr3    = cu(Lr),               // The radius to the third power [cm^3].
       g      = 981.0,                // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,             // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),    // The squared Froudes number.
       Re     = q*rho/mu/Lr,          // Reynolds number.
       Period = Tper*q/Lr3,           // The dimension-less period.
       nu_pl  = nu,//nu*q/Lr2,               // Dimension-less kinematic viscosity
       mu_pl  = mu,//*rho*g*Lr,       // Dimension-less viscosity
       k      = Period/tmstps,        // Length of a timestep.
       Deltat = Period/plts,         // Interval between each point plottet.
       bound_thick = sqrt(Tper*nu/(2.0*M_PI))/Lr, // Boundary layer thickness


       p0     = 0.0;//-10/rho/g/Lr*conv;      // Ensures a certain diastolic pressure. (Unstressed pressure)

double       *fjac[18],   // Work space used by bound_bif.
             *fj[8];      // Work space used by bound_match.

#endif
