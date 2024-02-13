/* The sor06.C main program */

// $Id: sor06.C,v 1.10 2005/10/14 18:05:59 heine Exp $
// Last updated on February 14, 2020, MJC

#include "sor06.h"
#include "tools.h"
#include "arteries.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

// The vessel-network is specified and initialized, and the flow and
// pressures are to be determed. All global constants must be defined
// in the header file. That is the number of dots per cm, and the number
// of vessels.
int main(int argc, char *argv[])
{

    /*
     MJC: Set up the code to feed in the parameters describing vascular stiffness.
     THERE ARE NOT OPTIONAL ARGUMENTS: All parameters are needed to run the model. Default
     values will be enforced in future versions. First input argument is the executable file ./sor06
     */

    if (argc<6) {
        printf("ERROR: Not enough input arguments. Exiting.\n");
        return 0;
    }

    

    double Af1,Af2,Af3;                // Large artery stiffness
    double fs1,fs2,fs3;                // Small vessel stiffness
    double Vf1,Vf2,Vf3;                // Large venous stiffness
    double alpha_b,beta_b,lrrA,lrrV,rm;
                                    //artery/vein length to radius ratio, and minimum radius.
    int num_pts, cycles, fileID;            // Number of spatial points and number of output cycles

    Af1       = 0.0;
    Af2       = 0.0;
    Af3       = atof(argv[1]);
    fs1       = 0.0;
    fs2       = 0.0;
    fs3       = atof(argv[2]);
    Vf1       = 0.0;
    Vf2       = 0.0;
    Vf3       = atof(argv[3]);
    alpha_b   = atof(argv[4]);
    beta_b    = atof(argv[5]);
    lrrA      = atof(argv[6]);
    lrrV      = atof(argv[7]);
    rm        = atof(argv[8]);
    fileID    = atoi(argv[9]);
    num_pts   = 10;
    cycles    = 1;

  double tstart, tend, finaltime;
    char namepu1[24],nameinflow[24],nameoutflow[24];
    sprintf(namepu1, "output_%d.2d", fileID);
    sprintf(nameinflow, "Qin_%d.dat", fileID);
    sprintf(nameoutflow, "LAout_%d.dat", fileID);

	FILE *fpu1 = fopen (namepu1, "w");

   CO_filename  = nameinflow;   // Input flow file at the heart.
    PL_filename = nameoutflow;       // Outflow (Static pressure) file.


  // Workspace used by bound_match
  for(int i=0; i<8; i++) fj[i] = new double[8];

  // Workspace used by bound_bif
  for(int i=0; i<18; i++) fjac[i] = new double[18];

//   clock_t c1 = clock();        // Only used when timing the program.
  nbrves    = 27; //11;              // Total number of large arteries and veins in the network
//   nbrves    = 11;
    tstart    = 0.0;           // Starting time.

  // The number of vessels in the network is given when the governing array of
  // vessels is declared.

  Tube   *Arteries[nbrves];                    // Array of blood vessels.

  // *---------Initialisation of the vessel network.---------*
  // init == 1 implies initial (inflow) artery.
  // init == 2 implies initial (outflow) vein.
  // init == 3 implies terminal artery to be matched to vein,
  //  (with LD == 0, RD == pointer to vein).
  // init == 4 implies monof,
  // r_min == 0 implies vessel has two daughter vessels.
  // r_min != 0 implies 'terminal' vessel,
  //  (terminal arteries - one daughter,
  //   terminal veins - no daughters).

// Parameters required to initiate class Tube (Length,topradius,botradius,LeftDaughter,RightDaughter,rmin, points,
                                               //init,K,f1,f2,f3,fs1,fs2,fs3,Vf1,Vf2,Vf3,alpha_b,beta_b,lrrA,lrrV);


  // Initialization of the Veins. Approximated dimensions: Data 1 in PV_dimensionV2.xlsx

    // 27 vessels
    // Veins
    Arteries[26]  = new Tube( 1.17, 0.610,0.610, 0, Arteries[14], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[25]  = new Tube( 1.73, 0.460,0.460, 0, Arteries[13], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[24]  = new Tube( 1.43, 0.562,0.562, 0, Arteries[12], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[23]  = new Tube( 2.11, 0.829,0.829, 0, Arteries[11], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[22]  = new Tube( 0.95, 0.293,0.293, 0, Arteries[10], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[21]  = new Tube( 1.10, 0.433,0.433, 0, Arteries[ 9], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[20]  = new Tube( 1.31, 0.514,0.514, 0, Arteries[ 8], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[19]  = new Tube( 1.93, 0.757,0.757, 0, Arteries[ 7], rm, num_pts,0,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[18]  = new Tube( 1.92, 0.824,0.824, Arteries[25], Arteries[26], 0, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[17]  = new Tube( 2.35, 0.864,0.864, Arteries[23], Arteries[24], 0, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[16]  = new Tube( 1.23,  0.716,0.716, Arteries[21], Arteries[22], 0, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[15]  = new Tube( 2.15, 0.641,0.641, Arteries[19], Arteries[20], 0, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    // Arteries
    Arteries[14]  = new Tube( 1.55, 0.610,0.610, 0,  Arteries[26], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[13]  = new Tube( 1.17, 0.460,0.460, 0,  Arteries[25], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[12]  = new Tube( 1.43, 0.562,0.562, 0,  Arteries[24], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[11]  = new Tube( 2.11, 0.829,0.829, 0,  Arteries[23], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[10]  = new Tube( 0.95, 0.293,0.293, 0,  Arteries[22], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 9]  = new Tube( 1.10, 0.433,0.433, 0,  Arteries[21], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 8]  = new Tube( 1.31, 0.514,0.514, 0,  Arteries[20], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 7]  = new Tube( 1.93, 0.757,0.757, 0,  Arteries[19], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 6]  = new Tube( 1.92, 0.755,0.755, Arteries[13],  Arteries[14], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 5]  = new Tube( 2.35, 0.922,0.922, Arteries[11],  Arteries[12], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 4]  = new Tube( 1.23, 0.481,0.481, Arteries[ 9],  Arteries[ 10], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 3]  = new Tube( 2.15, 0.842,0.842, Arteries[ 7],  Arteries[ 8], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 2]  = new Tube( 5.75, 1.10,1.10,  Arteries[ 5],  Arteries[ 6], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[ 1]  = new Tube( 2.50, 0.90,0.90,  Arteries[ 3],  Arteries[ 4], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
    Arteries[0] = new Tube( 4.30, 1.35,1.35,  Arteries[ 1], Arteries[ 2], 0, num_pts,1,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3,  alpha_b, beta_b, lrrA, lrrV,1);

  // 11 Vessels
    // Veins
//      Arteries[10]  = new Tube( 1.90, 0.755,0.755, 0, 0, rm, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[9]  = new Tube( 2.40,  0.922,0.922, 0, 0, rm, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[8]  = new Tube( 1.20, 0.481,0.481, 0, 0, rm, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[7]  = new Tube( 2.10, 0.842,0.842, 0, 0, rm, num_pts,2,0,Vf1,Vf2,Vf3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);

    // Arteries
//     Arteries[6]  = new Tube( 1.92, 0.755,0.755, 0,  Arteries[10], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[5]  = new Tube( 2.35, 0.922,0.922, 0,  Arteries[ 9], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[4]  = new Tube( 1.23, 0.481,0.481, 0,  Arteries[ 8], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[3]  = new Tube( 2.15, 0.842,0.842, 0,  Arteries[ 7], rm, num_pts,3,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[2]  = new Tube( 5.75, 1.10,1.10,  Arteries[ 5],  Arteries[ 6], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[1]  = new Tube( 2.50, 0.90,0.90,  Arteries[ 3],  Arteries[ 4], 0, num_pts,0,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3, alpha_b, beta_b, lrrA, lrrV,2);
//     Arteries[0] = new Tube( 4.30, 1.35,1.35,  Arteries[ 1], Arteries[ 2], 0, num_pts,1,0,Af1,Af2,Af3,fs1,fs2,fs3,fs1,fs2,fs3,  alpha_b, beta_b, lrrA, lrrV,1);

  int period_counter = 1; // Count the number of periods you have solved for
  double norm_sol = 1e+6;
  double sol_tol  = 1e+0;
//   printf("NORM_SOL: %f\n",norm_sol);
  double sol_p1[tmstps],sol_p2[tmstps];
  tend      = Deltat;
    fflush(stdout);


  // SOLVE THE MODEL ONCE
  // Note: Only want to test the pressure at the inlet
  int sol_ID = 0;
  while (tend<=period_counter*Period)
  {
  solver (Arteries, tstart, tend, k);
  sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0]); // for printing
  sol_p1[sol_ID] *= rho*g*Lr/conv;
  tstart = tend;
  tend   = tend + Deltat; // The current ending time is increased by Deltat.
  sol_ID++;
  }


  // LOOP FOR CONVERGENCE
  double sse;
  while (norm_sol>=sol_tol)
  {
      sol_ID = 0;
      sse    = 0;
      period_counter++;
      while (tend<=period_counter*Period)
      {
          solver (Arteries, tstart, tend, k);
          sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0]); // for printing
          sol_p2[sol_ID] *= rho*g*Lr/conv;
          sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
          tstart = tend;
          tend   = tend + Deltat; // The current ending time is increased by Deltat.
          sol_ID++;
      }
      norm_sol = sse;
      memcpy (sol_p1, sol_p2, sizeof(sol_p2));

  }

  // The loop is continued until the final time
  // is reached. If one wants to make a plot of
  // the solution versus x, tend is set to final-
  // time in the above declaration
    period_counter++;
    finaltime = (period_counter+(cycles-1))*Period;
      while (tend <= finaltime)
        {
          for (int j=0; j<nbrves; j++)
          {
            int ArtjN = Arteries[j]->N;
            for (int i=0; i<ArtjN; i++)
            {
              Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
              Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
            }
          }

          // Solves the equations until time equals tend.
          solver (Arteries, tstart, tend, k);
            // Single artery and vein
          for (int save_id=0; save_id<nbrves; save_id++)
         {
             Arteries[ save_id] -> printALLt (fpu1, tend, 0);
            
         }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printALLt(fpu1, tend, Arteries[save_id]->N/2);
        }
        for (int save_id=0; save_id<nbrves; save_id++)
        {
            Arteries[ save_id] -> printALLt (fpu1, tend, Arteries[save_id]->N);
        }

            fflush(fpu1);

         // The time within each print is increased.
         tstart = tend;
         tend   = tend + Deltat; // The current ending time is increased by Deltat.
        }


  // In order to termate the program correctly the vessel network and hence
  // all the vessels and their workspace are deleted.
//   for (int i=0; i<nbrves; i++) delete Arteries[i];

  // Matrices and arrays are deleted
//   for (int i=0; i<18; i++) delete[] fjac[i];
//   for (int i=0; i<8;  i++) delete[] fj[i];


 return 2;
}
