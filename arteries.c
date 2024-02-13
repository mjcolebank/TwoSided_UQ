
/***************************************************************************/
/*                                                                         */
/*  Program: arteries.C                                                    */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This module can predict the flow and pressure in an tree of elastic    */
/*  vessels as described in IMFUFATEKST NR 297, and D2.1-4. The dependen-  */
/*  cies of the vessels in the tree must be specified in the main module   */
/*  according to the tree in question (for further details see documenta-  */
/*  tion in problem.default.cxx).                                          */
/*  This module includes all the functions needed to solve the system      */
/*  of equations. That is the description of all functions in the class    */
/*  containing the vessel (for further details see arteries.h), and in     */
/*  particular the functions needed to solve the system of equations nu-   */
/*  merically.                                                             */
/*                                                                         */
/*  The module is dependent on the utilities in tools.C, and               */
/*  their corresponding h-files, and also arteries.h that includes the     */
/*  declaration of the vessel-object.                                      */
/*                                                                         */
/***************************************************************************/

// Last updated on February 14, 2020, MJC

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>

#include "tools.h"
#include "arteries.h"

using namespace std;

extern int nbrves;
extern int tmstps;
extern char* CO_filename;
extern char* PL_filename;

// A f90 subroutine that determines the impedance at the root of a
// structured tree. Used as BC for the tree consisting of the larger
// arteries.

extern "C" void impedance_driver_(int *tmstps, double *Period,
                                  double *rho, double *mu,
                                  double *r_root, double *rmin,
                                  double *y11, double *y12, double *y21, double *y22,
                                  double *Lr, double *q, double *g,
                                  double *fa1, double *fa2, double *fa3,
                                  double *fv1, double *fv2, double *fv3,
                                  double *alpha_b, double *beta_b,
                                  double *lrrA, double *lrrV, int *term_ID);


/* Methods of class Tube, see arteries.h for description of this. */

// The constructor. When an object is made this function will initialize
// all the attributes of the specific tube. The parameters for the length
// of the specific vessel, the top and bottom radii, and if applicable
// the pointers to the daughter arteries will be initialized according
// to the actual parameters passed in the call.
// If the tube is terminal then the peripheral resistance must be set,
// and the daughter vessels should be NIL. Otherwise the pointers to
// the daughter vessels must be given.
// Further all the work arrays are declared and initialized, and the
// initial condition for the system equations is applied.
Tube :: Tube (double Length,
              double topradius, double botradius,
              Tube *LeftDaughter, Tube *RightDaughter,
              double rmin, double points, int init, double K,
              double f1, double f2, double f3,
              double fa1, double fa2, double fa3,
              double fv1, double fv2, double fv3,
              double alpha_b, double beta_b,
              double lrrA, double lrrV,
              int term_ID):
	L(Length),
	rtop(topradius),
	rbot(botradius),
	LD(LeftDaughter),
	RD(RightDaughter),
  rm(rmin),
	pts(points),
  init(init),
	K_loss(K),
	f1(f1),
	f2(f2),
	f3(f3),
  fa1(fa1),
  fa2(fa2),
  fa3(fa3),
  fv1(fv1),
  fv2(fv2),
  fv3(fv3)
{
  // Initialization of the basic parameters
  N	  = int(pts*L);
  h	  = 1.0/pts/Lr;

  //fprintf(stdout, "K_loss=%8.4f\n",K_loss);

  // Declaration and Initialization of the needed intermediate arrays.
  Qnew	  = new double[N+1];
  Anew	  = new double[N+1];
  Qold	  = new double[N+1];
  Aold	  = new double[N+1];
  Qprv	  = new double[N+1];
  Aprv	  = new double[N+1];
  R1	    = new double[N+1];
  R2	    = new double[N+1];
  S1	    = new double[N+1];
  S2	    = new double[N+1];
  r0	    = new double[N+1];
  r0h	    = new double[N+2];
  dr0dx   = new double[N+1];
  dr0dxh  = new double[N+2];
  wom     = new double[N+1];
  A0      = new double[N+1];
  A0h     = new double[N+2];
  fr      = new double[N+1];
  frh     = new double[N+2];
  dfrdr0  = new double[N+1];
  dfrdr0h = new double[N+2];
  p1      = new double[N+1];
  p1h     = new double[N+2];
  dp1dr0  = new double[N+1];
  dp1dr0h = new double[N+2];
  Ah	  = new double[N];
  Qh	  = new double[N];
  R1h	  = new double[N];
  R2h	  = new double[N];
  S1h	  = new double[N];
  S2h	  = new double[N];
  pL	  = new double[tmstps];
  y11	  = new double[tmstps];
  y12	  = new double[tmstps];
  y21	  = new double[tmstps];
  y22	  = new double[tmstps];

  double rgLr  = 1.0/rho/g/Lr;

    // MJC: Define Nondimensional stiffnesses
    double f1ND = f1*rgLr;
    double f2ND = f2*Lr;
    double f3ND = f3*rgLr;
    double fa1ND = fa1;//*rgLr;
    double fa2ND = fa2;//*Lr;
    double fa3ND = fa3;//*rgLr;
    double fv1ND = fv1;//*rgLr;
    double fv2ND = fv2;//*Lr;
    double fv3ND = fv3;//*rgLr;

  // Vessel geometry is tabulated and initial conditions are applied
  for (int i=0; i<=N; i++)
  {
    r0 [i]     = rtop*exp(i*log(rbot/rtop)/N)/Lr;
    r0h[i]     = rtop*exp((i-0.5)*log(rbot/rtop)/N)/Lr;
    dr0dx [i]  = (log(rbot/rtop)/h/N)*r0 [i];
    dr0dxh[i]  = (log(rbot/rtop)/h/N)*r0h[i];
    wom[i]     = r0[i]*sqrt(2.0*M_PI/Period/nu);
    A0 [i]     = M_PI*sq(r0 [i]);
    A0h[i]     = M_PI*sq(r0h[i]);
    fr [i]     = 4.0/3.0*(f1ND*exp(f2ND*r0 [i])+f3ND);
    frh[i]     = 4.0/3.0*(f1ND*exp(f2ND*r0h[i])+f3ND);
    dfrdr0 [i] = 4.0/3.0*f1ND*f2ND*exp(f2*r0 [i]);// /Lr;
    dfrdr0h[i] = 4.0/3.0*f1ND*f2ND*exp(f2*r0h[i]);///Lr;
    p1 [i]     = fr [i]/M_PI;
    p1h[i]     = frh[i]/M_PI;
    dp1dr0 [i] = dfrdr0 [i]/M_PI;
    dp1dr0h[i] = dfrdr0h[i]/M_PI;
    Qnew[i]    = 1.0;
    Anew[i]    = A0[i];
  }
  r0h[N+1]     = rtop*exp((N+0.5)*log(rbot/rtop)/N)/Lr;
  dr0dxh[N+1]  = log(rbot/rtop)/h/N*r0h[N+1];
  A0h[N+1]     = M_PI*sq(r0h[N+1]);

//   fr[N+1]    = 4.0/3.0*(f1ND*exp(f2ND*r0 [N+1])+f3ND);
  frh[N+1]   = 4.0/3.0*(f1ND*exp(f2ND*r0h[N+1])+f3ND);
  p1h[N+1]     = frh[N+1]/M_PI;
  dp1dr0h[N+1] = dfrdr0h[N+1]/M_PI;
  Ah05         = A0[N];
  Qh05         = 1.0;

  // Read from file data for the inflow profile.
  if (init == 1)
  {
    Q0 = new double[tmstps+1]; // Reading inflow from "gvPA_4096.dat".
                               //Comment if useing pressure as an input condition
    //Ps = new double[tmstps+1]; //Reading input pressure from "p0_in4096.dat".
    FILE *fi = fopen (CO_filename, "r");
//     if (fi) fprintf(stdout, "Input file opened\n");
//     else error ("arteries.C"," INPUT FILE NOT OK");

    for (int i=0; i<=tmstps; i++)
    {
      fscanf(fi,"%lf",&Q0[i]);   //Comment out if using pressure input
      Q0[i] = Q0[i]/q; // If the indata have dimensions they should be made non-dimensional.

      //fscanf(fi,"%lf",&Ps[i]); //Uncomment for pressure input
      //Ps[i] = Ps[i]/(rho*g*Lr/conv); //Uncomment for pressure input
    }
      fclose(fi);
  }

  // Static pressure condition at the outlet of four large veins
   if (init == 2)
   {
    Pout = new double[tmstps+1];

     FILE *fi = fopen (PL_filename, "r" );
//      if (fi) fprintf(stdout, "PL opened\n");
//      else error ("arteries.C"," PL NOT OK");

        for (int i=0; i<=tmstps; i++)
     {
        fscanf(fi,"%lf",&Pout[i]);
        Pout[i] =( Pout[i])/(rho*g*Lr/conv); // Pout[i]+a, where a is the non-zero impost outlet pressure
         // If the indata have dimensions they should be made non-dimensional.
     }
   }

  // In case of an end-tube evaluate the impedances for the boundary condition.
  // This is done by calling the f90 routine root_imp which calculates the
  // impedance at the root of a structured tree. The underscores is sensitive
  // to the compiler but can be seen at the bottom of the file root_imp.o.
  if (init == 3)
  {
//      fprintf(stdout,"Calling f90 subroutines\n");
//      impedance_driver_(&tmstps,&Period,&rho,&mu_pl,&rbot,&rmin,y11,y12,y21,y22,&Lr,&q,&g,&fa1,&
//                        fa2,&fa3,&fv1,&fv2,&fv3,&alpha_b,&beta_b,&lrrA,&lrrV, &term_ID);
//      printf("Finished with f90 subroutines.\n");
      // MJC
//       fprintf(stdout,"Calling f90 subroutines\n");
      impedance_driver_(&tmstps,&Period,&rho,&mu_pl,&rbot,&rmin,y11,y12,y21,y22,&Lr,&q,&g,&fa1ND,
                        &fa2ND,&fa3ND,&fv1ND,&fv2ND,&fv3ND,&alpha_b,&beta_b,&lrrA,&lrrV, &term_ID);
//       printf("Finished with f90 subroutines.\n");

      // Initialize the array pL used when determining the convolution
      // in the right boundary condition (see the subroutine bound_right).
     for (int j=0; j<tmstps; j++)
      {
        pL[j] = p1[N]; //Runs from 0 to N-1
      }
  }
    /*
     ADDED BY MJC: In the case we have a single inlet, we need to call both bound-left and bound_match
     */
     if (init == 1 && LD==0)
     {
//         fprintf(stdout,"Calling f90 subroutines\n");
//         impedance_driver_(&tmstps,&Period,&rho,&mu_pl,&rbot,&rmin,y11,y12,y21,y22,&Lr,&q,&g,&fa1,&fa2,&fa3,&fv1,&fv2,&fv3,&alpha_b,&beta_b,&lrrA,&lrrV, &term_ID);
//         printf("Finished with f90 subroutines.\n");
         // MJC
//          fprintf(stdout,"Calling f90 subroutines\n");
               impedance_driver_(&tmstps,&Period,&rho,&mu_pl,&rbot,&rmin,y11,y12,y21,y22,&Lr,&q,&g,&fa1ND,
                           &fa2ND,&fa3ND,&fv1ND,&fv2ND,&fv3ND,&alpha_b,&beta_b,&lrrA,&lrrV, &term_ID);
//          printf("Finished with f90 subroutines.\n");

         // Initialize the array pL used when determining the convolution
         // in the right boundary condition (see the subroutine bound_right).
        for (int j=0; j<tmstps; j++)
         {
           pL[j] = p1[N];//Runs from 0 to N-1
         }
     }
}

// The destructor. When the tube-objects terminates, all arrays are deleted,
// in order to free the memory occupied by the object.
Tube :: ~Tube ()
{
  delete[] Anew;
  delete[] Qnew;
  delete[] Aold;
  delete[] Qold;
  delete[] Aprv;
  delete[] Qprv;
  delete[] Ah;
  delete[] Qh;
  delete[] y11;
  delete[] y12;
  delete[] y21;
  delete[] y22;
  delete[] pL;
  delete[] R1h;
  delete[] R2h;
  delete[] S1h;
  delete[] S2h;
  delete[] R1;
  delete[] R2;
  delete[] S1;
  delete[] S2;
  delete[] r0;
  delete[] r0h;
  delete[] dr0dx;
  delete[] dr0dxh;
  delete[] A0;
  delete[] A0h;
  delete[] fr;
  delete[] frh;
  delete[] dfrdr0;
  delete[] dfrdr0h;
  delete[] p1;
  delete[] p1h;
  delete[] dp1dr0;
  delete[] dp1dr0h;
}

// ----------------------PLOTTING ROUTINES WITH DIMENSIONS ------------

void Tube :: printQ0 (FILE *fd)
{
  for (int i=0; i<=tmstps; i++)
  {
    fprintf (fd, "%15.10f\n", Q0[i]*q);
  }
}


// The following functions prints p, q(x,t) in terms of the re-dimensionalized
// variables. The parameters for the function are the  position (x),
// and the time (t).
void Tube :: printPt (FILE *fd, double t, int i)
{
 fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, (P(i,Anew[i])+p0)*rho*g*Lr/conv);
}

void Tube :: printQt (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, Qnew[i]*q);
//fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, Qnew[N-1]*q);
}

void Tube :: printAt (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, Anew[i]*Lr2);
}

void Tube :: printPUt (FILE *fd, double t, int i)// printing pressure, velocity and the wavespeed at a given position
{
	fprintf (fd, "%13.10f %15.10f %15.10f %17.10f\n", t*Lr3/q, (P(i,Anew[i])+p0)*rho*g*Lr/conv, (Qnew[i]*q)/(Anew[i]*Lr2), c(i, Anew[i])*Fr2);
}

void Tube :: printFt (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t*Lr3/q, F(Qnew[i],Anew[i])*sq(q)/Lr3);
}

// The following functions prints P, Q, A, and F as functions of
// (x, t). This is done in terms of the re-dimensionalized variables.
// In this case the functions is plotted for a
// fixed time, but for all x along the vessel in question. Since the
// doesn't have to be the first vessel in the tree, it would have
// some offset from the heart. Which determines the position for x.
// Therefore there are two arguments passed to this function the time
// and the offset.
// Umar modified the print routines to write importnat parameters in one file. Convenienrt files for WIA.
void Tube :: printALLt (FILE *fd,double t, int i)
{
    fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f\n", t*Lr3/q, i*h*Lr, (P(i,Anew[i])+p0)*rho*g*Lr/conv, Qnew[i]*q, Anew[i]*Lr2, c(i, Anew[i]));// NEED TO REDIMENSIONALIZE WAVESPEED
}



void Tube :: printPxt (FILE *fd, double t, int offset)
{
  if (offset == 0) fprintf (fd, "\n");
  for (int i=0; i<N; i++)
  {

      // //Original Pxt files
fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f\n",
	 t*Lr3/q, (i+offset)*h*Lr, (P(i,Anew[i])+p0)*rho*g*Lr/conv, Qnew[i]*q, Anew[i]*Lr2, c(i, Anew[i])*Fr2); // modifeid to print more variables
    // Time, Vessel Length, Pressure [mmHg], Flow rate [cm^3/sec], Cross-sectional Area [cm^2], Pulse wave velocity [m/sec]
      // commentout the above line to print out velocity instead of flow
  }
}

void Tube :: printQxt (FILE *fd, double t, int offset)
{
  if (offset == 0) fprintf (fd, "\n");
  for (int i=0; i<N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, Qnew[i]*q);
  }
}

void Tube :: printAxt (FILE *fd, double t, int offset)
{
  for (int i=0; i<=N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, Anew[i]*Lr2);
  }
}

void Tube :: printFxt (FILE *fd, double t, int offset)
{
  for (int i=0; i<=N; i++)
  {
    fprintf (fd, "%13.10f %13.10f %15.10f\n",
             t*Lr3/q, (i+offset)*h*Lr, F(Qnew[i],Anew[i])*sq(q)/Lr3);
  }
}

// A function that prints p(Q) for all t. This is done in terms of
// the re-dimensionalized variables. In this case the plot is made for a
// fixed point in space, but for all t along the vessel in question.
void Tube :: printPQ (FILE *fd, int i)
{
  fprintf (fd,"%15.10f %15.10f\n",
           Qnew[i]*q, (P(i,Anew[i])+p0)*rho*g*Lr/conv);
}

void Tube :: printPA (FILE *fd, int i)
{
  fprintf (fd,"%15.10f %15.10f\n",
           (P(i,Anew[i])+p0)*rho*g*Lr/conv, Anew[i]*Lr2);
}

// Plotting the terms in the continuity equation on dimension-less form.
void Tube :: printdQdx (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Qnew[i+1]-Qnew[i-1])/2.0/h);

}

void Tube :: printdAdt (FILE *fd, double t, int i, double Aprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Anew[i]-Aprev)/tmst);
}

void Tube :: printTotConRes (FILE *fd, double t, int i, double Aprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
          (Qnew[i+1] - Qnew[i-1])/2.0/h +
          (Anew[i]-Aprev)/tmst);
}

// Plotting the terms in the momentum equation on dimension-less form.
void Tube :: printdQdt (FILE *fd, double t, int i, double Qprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (Qnew[i]-Qprev)/tmst);
}

void Tube :: printddxQ2divA (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           (sq(Qnew[i+1])/Anew[i+1] - sq(Qnew[i-1])/Anew[i-1])/2.0/h);

}

void Tube :: printdPdx (FILE *fd, double t, int i)
{
    fprintf (fd, "%13.10f %15.10f\n", t,
             Anew[i]*(P(i+1,Anew[i+1])-P(i-1,Anew[i-1]))/Fr2/2.0/h);
}

void Tube :: printFric (FILE *fd, double t, int i)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           F(Qnew[i],Anew[i]));
}

void Tube :: printTotMomRes (FILE *fd, double t, int i, double Qprev, double tmst)
{
  fprintf (fd, "%13.10f %15.10f\n", t,
           Anew[i]*(P(i+1,Anew[i+1])-P(i-1,Anew[i-1]))/Fr2/2.0/h +
           (sq(Qnew[i+1])/Anew[i+1] - sq(Qnew[i-1])/Anew[i-1])/2.0/h +
           (Qnew[i]-Qprev)/tmst + F(Qnew[i],Anew[i]));
}

// Further print functions can be added, and they would look similar to
// the two functions above!

// The next function returns the pressure p as a function of a fixed x,
// and the corresponding cross-sectional area A. The pressure is defined
// according to the  mathematical model, described in IMFUFATEKST no 297,
// and D2.1-4.
double Tube :: P (int i, double A)
{
//   double pold = fr[i]*(1-sqrt(A0[i]/A));
  double pold = fr[i]*(sqrt(A/A0[i])-1.0);  // Linear B

  return pold;
}

double Tube :: dPdA (int i, double A)
{
//    double pold = 0.5*fr[i]*sqrt(A0[i]/cu(A));
    double pold = 0.5*fr[i]/sqrt(A0[i]*A);  // Linear B

  return pold;
}

double Tube :: dPdx1(int i, double A)
{
//   double pold = (dfrdr0[i]*(1.0-sqrt(A0[i]/A))-fr[i]*sqrt(M_PI/A))*dr0dx[i];
    double pold = (dfrdr0[i]*(sqrt(A/A0[i])-1)-fr[i]*sqrt(M_PI*A)/A0[i])*dr0dx[i];  // Linear B

  return pold;
}

double Tube :: B (int i, double A)
{
//   double pold = fr[i]*(sqrt(A0[i]*A)-A0[i])/Fr2;
    double pold = fr[i]*(sqrt(cu(A)/A0[i])-A0[i])/Fr2/3.0;  // Linear B

  return pold;
}

double Tube :: Bh (int i, double A)
{
   int ip1 = i+1;
//    double pold =  frh[ip1]*(sqrt(A0h[ip1]*A)-A0h[ip1])/Fr2;
    double pold = frh[ip1]*(sqrt(cu(A)/A0h[ip1])-A0h[ip1])/Fr2/3.0;  // Linear B

   return pold;
}

double Tube :: dBdx1 (int i, double A)
{
   double dfr = dfrdr0[i];
//   double pold = dr0dx[i]/Fr2*(
//          sqrt(A)*(2.0*sqrt(M_PI)*fr[i]+2.0*sqrt(A0[i])*dfr)-
//          A*dfr-2.0*M_PI*r0[i]*fr[i]-A0[i]*dfr);
    double inteval = 2.0*M_PI*r0[i]*fr[i]+A0[i]*dfr;
           double pold = dr0dx[i]*(2.0*sqrt(cu(A))*(sqrt(M_PI)*fr[i]-sqrt(A0[i])*dfr)/A0[i]/3.0
                                   + A*dfr - inteval/3.0)/Fr2; // Linear B

  return pold;
}

double Tube :: dBdx1h (int i, double A)
{
  int ip1 = i+1;

  double dfr = dfrdr0h[ip1];

//   double pold = dr0dxh[ip1]/Fr2*(
//          sqrt(A)*(2.0*sqrt(M_PI)*frh[ip1]+2.0*sqrt(A0h[ip1])*dfr)-
//          A*dfr-2.0*M_PI*r0h[ip1]*frh[ip1]-A0h[ip1]*dfr);

    double inteval = 2.0*M_PI*r0h[ip1]*frh[ip1]+A0h[ip1]*dfr;
       double pold = dr0dxh[ip1]*(2.0*sqrt(cu(A))*(sqrt(M_PI)*frh[ip1]-sqrt(A0h[ip1])*dfr)/A0h[ip1]/3.0
                                   + A*dfr - inteval/3.0)/Fr2; // Linear B

  return pold;
}

double Tube :: dBdAh (int i, double A)
{
  int ip1      = i+1;
//   double pold = 0.5*frh[ip1]*sqrt(A0h[ip1]/A)/Fr2;
    double pold = 0.5*frh[ip1]*sqrt(A/A0h[ip1])/Fr2;  // Linear B

  return pold;
}

double Tube :: d2BdAdxh (int i, double A)
{
   int ip1      = i+1;
   double dfr = dfrdr0h[ip1];
//    double pold = (-dfr+1.0/sqrt(A)*(sqrt(M_PI)*frh[ip1]+
//            sqrt(A0h[ip1])*dfr))*dr0dxh[ip1]/Fr2;

    double pold = (dfr+sqrt(A)*(sqrt(M_PI)*frh[ip1] - sqrt(A0h[ip1])*dfr)/A0h[ip1])*dr0dxh[ip1]/Fr2;  // Linear B

   return pold;
}

// When determining or checking the step-size (k) the CFL-condition is applied.
// This is determined according to the result reached from the analysis
// made using the method of characteristics (See IMFUFATEKST no 297).
// In this function the minimal step-size fulfilling this condition for this
// tube is returned.

double Tube :: CFL () // The CFL-condition
{
  double minimum = 64000000.0;
  for (int i=0; i<=N; i++)
  {
    double c_tmp = c(i, Anew[i]);
    double Vnew  = Qnew[i]/Anew[i];
    double temp = min (h / fabs (Vnew - c_tmp),
                h / fabs (Vnew + c_tmp));
    if (temp < minimum) minimum = temp;
  }
  return (minimum);
}

// When taking a Lax-Wendroff step, the flux of the system must be determined.
// This is evaluated at i + j/2, and the prediction is given as described
// in IMFUFATEKST no 297 and D2.1-4. The integer k determines whether we deal
// with the first or the second component of the vector.
double Tube :: Rvec (int k, int i, int j, double Q, double A)
{
      double qA2_term=((vel_power+2.0)/(vel_power+1.0)); // Power law

  if(k==1) return(Q); else
  if(k==2) return(qA2_term*sq(Q)/A + ((j==0)?B(i,A):Bh(i,A)));
  else error ("arteries.cxx","Call of non-existing vector-component of R");
  return(0);
}

// Similarly the right hand side of the system of equations must be determined
// at i + j/2. Also in this case the function is given as stated in
// the mathematical model, and also in this case k states the needed component
// of the vector.
double Tube :: Svec (int k, int i, int j, double Q, double A)
{
  if(k==1) return(0.0); else
  if(k==2) return(F(Q,A) + ((j==0)?dBdx1(i,A):dBdx1h(i,A)));
  else error ("arteries.cxx","Call of non-existing vector-component of S");
  return(0);
}

// The solutions of Anew and Qnew are found for all interior points
// of the vessel at (t+k), where k is the length of the current
// time-step. This function saves the results in the arrays Anew and
// Qnew, and the function is made according to Lax-Wendroff's method
// as described in IMFUFATEKST no 297 and D2.1-4.
void Tube :: step (double k)
{
  double theta = k/h;  // Theta is determined.
  double gamma = 0.5*k;  // Gamma is determined.

  for (int i=0; i<=N; i++)  // Remember the values at this time level.
  {
    Qold[i] = Qnew[i];
    Aold[i] = Anew[i];
  }

  // Anew and Qnew are predicted at the new time level (t+k).
  for (int i=0; i<=N; i++)
  {
    R1[i] = Rvec(1,i,0,Qold[i],Aold[i]);
    R2[i] = Rvec(2,i,0,Qold[i],Aold[i]);
    S1[i] = Svec(1,i,0,Qold[i],Aold[i]);
    S2[i] = Svec(2,i,0,Qold[i],Aold[i]);
  }

  for (int i=0; i<N; i++)
  {
    Ah[i]  = 0.5*(Aold[i+1]+Aold[i]) - 0.5*theta*(R1[i+1]-R1[i]) + 0.5*gamma*(S1[i+1]+S1[i]);
    //Qh[i]  = 0.5*(Qold[i+1]+Qold[i]) - 0.5*theta/Cu[i]*(R2[i+1]-R2[i]) + 0.5*gamma/Cu[i]*(S2[i+1]+S2[i]);
    Qh[i]  = 0.5*(Qold[i+1]+Qold[i]) - 0.5*theta*(R2[i+1]-R2[i]) + 0.5*gamma*(S2[i+1]+S2[i]);
    R1h[i] = Rvec(1,i,1,Qh[i],Ah[i]);
    R2h[i] = Rvec(2,i,1,Qh[i],Ah[i]);
    S1h[i] = Svec(1,i,1,Qh[i],Ah[i]);
    S2h[i] = Svec(2,i,1,Qh[i],Ah[i]);
  }
  for (int i=1; i<N; i++)
  {
    Anew[i] = Aold[i] - theta*(R1h[i]-R1h[i-1]) + gamma*(S1h[i]+S1h[i-1]);
    //Qnew[i] = Qold[i] - theta/Cu[i]*(R2h[i]-R2h[i-1]) + gamma/Cu[i]*(S2h[i]+S2h[i-1]);
    Qnew[i] = Qold[i] - theta*(R2h[i]-R2h[i-1]) + gamma*(S2h[i]+S2h[i-1]);
  }
}

//==========INPUT CONDITION (Adjust routines accordingly to use Q or P as an INPUT)====================

// The left boundary (x=0) uses this function to model an inflow into
// the system. The actual parameter given to the function is the model time.
// As stated in the mathematical model the constants of the function are
// chosen in order to ensure a certain CO (specified in main.h). Hence we have
// the specified value of b. Further the period (dimension-less) is assumed
// to be Period.
double Tube :: Q0_init (double t, double k, double Period)  //Coment for pressure input
{
  if (t <= Period) return (Q0[int(t/k)]); else
  if (t >  Period) return (Q0_init((t-Period),k,Period));
  else return (0);
}

/*double Tube :: Ps_init (double t, double k, double Period) //Uncoment for pressure input
{
	if (t <= Period) return (Ps[int(t/k)]); else
    if (t >  Period) return (Ps_init((t-Period),k,Period));
    else return (0);
}*/


// Update of the left boundary at time t. This function uses Q0 to determine
// the flow rate at the next time-step. From this the value of A is predicted
// using Lax-Wendroff's numerical scheme. This function is only relevant
// when the tube is an inlet vessel.
void Tube :: bound_left (double t, double k, double Period)
{
    double Q05 = 0.0;
    Qnew[0] = Q0_init(t,k,Period);
    Q05     = Q0_init(t-k,k,Period);

    if (int(t/k) < 0)
    printf("t/k negative in bound_left\n");

  double Qhm05 = Qnew[0]+Q05 - Qh[0];
  double R1hm05    = Qhm05;
  Anew[0]   = Aold[0] - k*(R1h[0] - R1hm05)/h;
}

//=====Uncomment the follwoing routine to use pressure as an input condition=====

/*void Tube :: bound_left (double t, double k, double Period)
{
	double Pnew   = Ps_init(t,k,Period);
	double P05    = Ps_init(t-k,k,Period);

	if (int(t/k) < 0)
    printf("t/k negative in bound_right\n");

	Anew[0] = A0[0]/(sq(1 - (Pnew/fr[0])));
	double A05 = A0[0]/(sq(1 - (P05/fr[0])));
	double Ahm05 = 2*A05 - Ah[0];
	double Qhm05 = R1h[0] + h*(Anew[0] - Aold[0])/k;
	double R2hm05 = Rvec(2,-1,1,Qhm05,Ahm05);
	double S2hm05 = Svec(2,-1,1,Qhm05,Ahm05);
	Qnew[0] = Qold[0] - k*(R2h[0] - R2hm05)/h + k*(S2h[0] + S2hm05)/2;
}*/
//==============================================================================

// The value at the right boundary at time t is predicted. NB: This should
// only be used with terminal vessels, i.e. for vessels that don't bifurcate
// into further branches.
// In that situation the bifurcation boundary function should be called
// instead. Again the procedure specified is given according to the mathemati-
// cal theory presented in IMFUFATEKST no 297 and D2.1-4.

double Tube :: c (int i, double A) // The wave speed through aorta.
{
//   double cnst =  0.5*fr[i]*sqrt(A0[i]/A)/Fr2;
  double  cnst =  0.5*fr[i]*sqrt(A/A0[i])/Fr2;  //Linear B

  return sqrt (cnst);
}

// The matching boundary.

void Tube :: bound_match (int qLnb, double t, double k, double theta, double gamma)
{
    double PA;
    double PV;
    double PAh;
    double PVh;
    int j = 1;
    int ok = false;
    int ntrial = 1000000;

    int qLnb_1 = qLnb + 1;

    // Make sure that qLnb_1 runs in the interval [0:tmstps-1].
    if (qLnb_1 == (int) tmstps)
    {
      qLnb_1 = 0;
    }

    // In order to make a relation between P(x_L, t+dt) and Q(x_L, t+dt), and
    // P(x_L, t+dt/2) and Q(x_L,t+dt/2) we need to extract the term involving
    // y[0] (see mathematical derivation). This term corresponds to the peripheral
    // The remaining terms in the convolution present at the boundary,
    // see mathematical derivation.

  double paterms = 0.0;
  double pvterms = 0.0;
  int M = RD->N;
  double RDtheta = k/RD->h;

  if (t > Period)
  {
    for (int m=1; m<tmstps; m++)
    {
      int pindex  = (qLnb_1 + tmstps - m) % tmstps;
      paterms  = paterms  + (pL[pindex])*y11[m]  + (RD->pL[pindex])*y12[m];
      pvterms  = pvterms  + (pL[pindex])*y21[m]  + (RD->pL[pindex])*y22[m];
    }
    paterms  = k*paterms;
    pvterms  = k*pvterms;
  }

  double k1   = Aold[N] + theta*R1h[N-1];
  double k2   = Qold[N] + theta*R2h[N-1] + gamma*S2h[N-1];
  double k3   = paterms;
  double k4   = k*y11[0];
  double k5   = k*y12[0];
  double k6   = RD->Aold[M] + RDtheta*(RD->R1h[M-1]);
  double k7   = RD->Qold[M] + RDtheta*(RD->R2h[M-1]) + gamma*(RD->S2h[M-1]);
  double k8   = pvterms;
  double k9   = k*y21[0];
  double k10  = k*y22[0];
  double k11  = k3 - 0.5*Qh[N-1];
  double k12  = Ah[N-1];
  double k13  = RD->Ah[M-1];
  double k14  = k8 - 0.5*(RD->Qh[M-1]);

  //Unknowns declared, and initial guesses applied
  double xb[8];

 if (t <= Period)
 {
   xb[0] = Ah[N-1];
   xb[1] = Qh[N-1];
   xb[2] = Aold[N];
   xb[3] = Qold[N];
   xb[4] = RD->Ah[M-1];
   xb[5] = RD->Qh[M-1];
   xb[6] = RD->Aold[M];
   xb[7] = RD->Qold[M];
 }

 if (t > Period)
 {
    xb[0] = Aold[N];
    xb[1] = Qold[N];
    xb[2] = Ah05; //Aold[N];
    xb[3] = Qh05; //Qold[N];
    xb[4] = RD->Aold[M];
    xb[5] = RD->Qold[M];
    xb[6] = RD->Ah05; //RD->Aold[M];
    xb[7] = RD->Qh05; //RD->Qold[M];
 }


  while (j <= ntrial && ok==false)
  {
    double fvec[8];

    PA = P(N,xb[0]);
    PV = RD->P(M,xb[4]);
    PAh = P(N,0.5*(k12 + xb[2]));
    PVh = RD->P(M,0.5*(k13 + xb[6]));

    fvec[0] = k1 - xb[0] - theta*xb[3];
    fvec[1] = k2 - xb[1] - theta*(sq(xb[3])/xb[2] + Bh(N,xb[2])) + gamma*(F(xb[3],xb[2])+dBdx1h(N,xb[2]));
    fvec[2] = k3 - xb[1] + k4*PA + k5*PV;
    fvec[3] = k6 - xb[4] - RDtheta*xb[7];
    fvec[4] = k7 - xb[5] - RDtheta*(sq(xb[7])/xb[6] + RD->Bh(M,xb[6])) + gamma*(RD->F(xb[7],xb[6])+RD->dBdx1h(M,xb[6]));
    fvec[5] = k8 - xb[5] + k9*PA + k10*PV;
    fvec[6] = k11 - xb[3]/2.0 + k4*PAh + k5*PVh;
    fvec[7] = k14 - xb[7]/2.0 + k9*PAh + k10*PVh;


    for (int row = 0; row < 8; row++)
    {
      for (int col = 0; col < 8; col++)
      {
        fj[row][col] = 0.0;
      }
    }

    fj[0][0]  = -1.0;
    fj[0][3]  = -theta;

    fj[1][1]  = -1.0;
    fj[1][2]  = theta*(sq(xb[3]/xb[2]) - dBdAh(N,xb[2])) + gamma*(dFdA(xb[3],xb[2]) + d2BdAdxh(N,xb[2]));
    fj[1][3]  = -2.0*theta*xb[3]/xb[2] + gamma*dFdQ(xb[2]);

    fj[2][0]  = k4*dPdA(N,xb[0]);
    fj[2][1]  = -1.0;
    fj[2][4]  = k5*RD->dPdA(M,xb[4]);

    fj[3][4]  = -1.0;
    fj[3][7]  = -RDtheta;

    fj[4][5]  = -1.0;
    fj[4][6]  = RDtheta*(sq(xb[7]/xb[6]) - RD->dBdAh(M,xb[6])) + gamma*(dFdA(xb[7],xb[6]) + RD->d2BdAdxh(M,xb[6]));
    fj[4][7]  = -2.0*RDtheta*xb[7]/xb[6] + gamma*dFdQ(xb[6]);

    fj[5][0]  = k9*dPdA(N,xb[0]);
    fj[5][4]  = k10*RD->dPdA(M,xb[4]);
    fj[5][5]  = -1.0;

    fj[6][2]  = k4*dPdA(N,(k12+xb[2])/2);
    fj[6][3]  = -0.5;
    fj[6][6]  = k5*RD->dPdA(M,(k13+xb[6])/2);

    fj[7][2]  = k9*dPdA(N,(k12+xb[2])/2);
    fj[7][6]  = k10*RD->dPdA(M,(k13+xb[6])/2);
    fj[7][7]  = -0.5;

    int ch = zero (xb, 8, 1.0e-8, 1.0e-8, fvec, fj);
    if (ch == 1) ok = true;

    j = j+1;
  }

  Anew[N]        = xb[0];
  Qnew[N]        = xb[1];
  Ah05           = xb[2];
  Qh05           = xb[3];
  RD->Anew[M]    = xb[4];
  RD->Qnew[M]    = xb[5]; // -1.0;
  RD->Ah05       = xb[6];
  RD->Qh05       = xb[7];
  pL[qLnb_1]     = P(N,Anew[N]);
  RD->pL[qLnb_1] = RD->P(M,xb[4]);

  if (j >=ntrial)
  {
    error ("arteries.C","Root not found in the matching");
    printf("t=%10.15f, t/(k*tmstps)=%10.15f \n",t,t/(k*tmstps));
    printf("AM=%10.15f, QM=%10.15f \n",xb[0],xb[1]);
    printf("AL=%10.15f, QL=%10.15f \n",xb[4],xb[5]);
    printf("rbot=%10.15f \n",rbot);
  }
//    printf("QA: %lf, QV: %lf\n",xb[1],xb[5]);

 }


 // The value at the bifurcation point at time t is predicted. NB: This should
 // only be done for tubes that do bifurcate into further branches. If
 // this is not the case we have a terminal vessel and bound_right should be
 // called instead. The procedure operates according to the specifications
 // in the mathematical model as a link between this tube and its daughters.
 // Therefore there will be three tubes involved in this function.
 // One problem is however, that the rather complicated system of equations does
 // not converge for all choices of parameters (the peripheral resistance, the
 // top radius, and the bottom radius).
 void Tube :: bound_bif (double theta, double gamma)
 {
   double PN;
   int j = 1;
   int ok = false;
     const int ntrial = 4000;//40;

   double g1   = Qold[N]     + theta*R2h[N-1] + gamma*S2h[N-1];
   double g2   = LD->Qold[0] - theta*(LD->R2h[0]) + gamma*(LD->S2h[0]);
   double g2a  = RD->Qold[0] - theta*(RD->R2h[0]) + gamma*(RD->S2h[0]);

   double k1   = Aold[N]     + theta*R1h[N-1];
   double k2   = LD->Aold[0] - theta*(LD->R1h[0]);
   double k2a  = RD->Aold[0] - theta*(RD->R1h[0]);

   double k3   = Qh[N-1]/2.0;
   double k4   = LD->Qh[0]/2.0;
   double k4a  = RD->Qh[0]/2.0;

   double k5   = Ah[N-1]/2.0;
   double k6   = LD->Ah[0]/2.0;
   double k6a  = RD->Ah[0]/2.0;

   double xb[18];

   // The approximative initial guesses are applied.
   xb[ 0] =  Qh[N-1];                      //Initial guess for Q1_xb n+1
   xb[ 1] = (Qold[N-1] + Qold[N])/2.0;       //Initial guess for Q1_xb^n+0.5
   xb[ 2] =  Qold[N];                      //Initial guess for Q1_xb+0.5 n+0.5
   xb[ 3] =  LD->Qh[0];                    //Initial guess for Q2_xb n+1
   xb[ 4] = (LD->Qold[0] + LD->Qold[1])/2.0; //Initial guess for Q2_xb n+0.5
   xb[ 5] =  LD->Qold[0];                  //Initial guess for Q2_xb+0.5 n+0.5
   xb[ 6] =  RD->Qh[0];                    //Initial guess for Q3_xb n+1
   xb[ 7] = (RD->Qold[0] + RD->Qold[1])/2.0; //Initial guess for Q3_xb n+0.5
   xb[ 8] =  RD->Qold[0];                  //Initial guess for Q3_xb+0.5 n+0.5
   xb[ 9] =  Ah[N-1];                      //Initial guess for A1_xb n+1
   xb[10] = (Aold[N-1] + Aold[N])/2.0;       //Initial guess for A1_xb^n+0.5
   xb[11] =  Aold[N];                      //Initial guess for A1_xb+0.5 n+0.5
   xb[12] =  LD->Ah[0];                    //Initial guess for A2_xb n+1
   xb[13] = (LD->Aold[0] + LD->Aold[1])/2.0; //Initial guess for A2_xb n+0.5
   xb[14] =  LD->Aold[0];                  //Initial guess for A2_xb+0.5 n+0.5
   xb[15] =  RD->Ah[0];                    //Initial guess for A3_xb n+1
   xb[16] = (RD->Aold[0] + RD->Aold[1])/2.0; //Initial guess for A3_xb n+0.5
   xb[17] =  RD->Aold[0];                  //Initial guess for A3_xb+0.5 n+0.5

   double k7nh  = LD->K_loss/2.0; //32*mu/(2*LD->rtop*rho*q);
   double k7n   = LD->K_loss/2.0; //32*mu/(2*RD->rtop*rho*q);
   double k7anh = RD->K_loss/2.0; //32*mu/(2*LD->rtop*rho*q);
   double k7an  = RD->K_loss/2.0; //32*mu/(2*RD->rtop*rho*q);

   // The residuals (fvec), and the Jacobian is determined, and if possible
   // the system of equations is solved.
   while (j <= ntrial && ok==false) // Find the zero
   {
     double fvec[18];
     // The residuals.

     //look at Bh/ dBx1h
     fvec[0]  = g1  - xb[0] -
                theta*(sq(xb[2])/xb[11] + Bh(N,xb[11])) +
                gamma*(F(xb[2],xb[11])+dBdx1h(N,xb[11]));

     fvec[1]  = g2  - xb[3] +
            theta*(sq(xb[5])/xb[14] + LD->Bh(-1,xb[14])) +
            gamma*(F(xb[5],xb[14])  + LD->dBdx1h(-1,xb[14]));

     fvec[2]  = g2a - xb[6] +
            theta*(sq(xb[8])/xb[17] + RD->Bh(-1,xb[17])) +
            gamma*(F(xb[8],xb[17])  + RD->dBdx1h(-1,xb[17]));

     fvec[3]  = - theta*xb[2] - xb[9]  + k1;
     fvec[4]  =   theta*xb[5] - xb[12]  + k2;
     fvec[5]  =   theta*xb[8] - xb[15]  + k2a;
     fvec[6]  = - xb[ 1] + xb[ 2]/2.0 + k3;
     fvec[7]  = - xb[ 4] + xb[ 5]/2.0 + k4;
     fvec[8]  = - xb[ 7] + xb[ 8]/2.0 + k4a;
     fvec[9]  = - xb[10] + xb[11]/2.0 + k5;
     fvec[10] = - xb[13] + xb[14]/2.0 + k6;
     fvec[11] = - xb[16] + xb[17]/2.0 + k6a;
     fvec[12] = - xb[ 1] + xb[ 4]   + xb[7];
     fvec[13] = - xb[ 0] + xb[ 3]   + xb[6];

     PN    = P(N,xb[10]);
     double sq211 = sq(xb[1]/xb[10]);

     if (xb[1] > 0)
     {
       fvec[14] =  - PN + LD->P(0,xb[13]) + k7nh*sq211;
 //        - 0.5*(sq211 - sq(xb[4]/xb[13]));
       fvec[15] =  - PN + RD->P(0,xb[16]) + k7anh*sq211;
 //        - 0.5*(sq211 - sq(xb[7]/xb[16]));
     } else
     {
       fvec[14] =  - PN + LD->P(0,xb[13]) - k7nh*sq211;
 //        - 0.5*(sq211 - sq(xb[4]/xb[13]));
       fvec[15] =  - PN + RD->P(0,xb[16]) - k7anh*sq211;
 //        - 0.5*(sq211 - sq(xb[7]/xb[16]));
     };

     PN    = P(N,xb[9]);
     double sq110 = sq(xb[0]/xb[9]);
     if (xb[0] > 0)
     {
       fvec[16] = - PN + LD->P(0,xb[12]) + k7n*sq110;
 //           - 0.5*(sq110 - sq(xb[3]/xb[12]));
       fvec[17] = - PN + RD->P(0,xb[15]) + k7an*sq110;
 //           - 0.5*(sq110 - sq(xb[6]/xb[15]));
     } else
     {
       fvec[16] = - PN + LD->P(0,xb[12]) - k7n*sq110;
 //           - 0.5*(sq110 - sq(xb[3]/xb[12]));
       fvec[17] = - PN + RD->P(0,xb[15]) - k7an*sq110;
 //           - 0.5*(sq110 - sq(xb[6]/xb[15]));
     };

     for (int row = 0; row < 18; row++)
       for (int col = 0; col < 18; col++)
         fjac[row][col] = 0.0;

     // The Jacobian.
     fjac[ 0][0]  = -1.0;
     fjac[13][0]  = -1.0;
     if (xb[0] > 0)
     {
       fjac[16][0] = xb[0]/sq(xb[9])*(2.0*k7n);//-1);
       fjac[17][0] = xb[0]/sq(xb[9])*(2.0*k7an);//-1);
     } else
     {
       fjac[16][0] = xb[0]/sq(xb[9])*(-2*k7n);//-1);
       fjac[17][0] = xb[0]/sq(xb[9])*(-2*k7an);//-1);
     };
     fjac[ 6][1] = -1.0;
     fjac[12][1] = -1.0;
     if (xb[1] > 0)
     {
       fjac[14][1] = xb[1]/sq(xb[10])*(2.0*k7nh);//-1);
       fjac[15][1] = xb[1]/sq(xb[10])*(2.0*k7anh);//-1);
     } else
     {
       fjac[14][1] = xb[1]/sq(xb[10])*(-2.0*k7nh);//-1);
       fjac[15][1] = xb[1]/sq(xb[10])*(-2.0*k7anh);//-1);
     };
     fjac[ 0][2] = -2.0*theta*xb[2]/xb[11] + gamma*dFdQ(xb[11]);
     fjac[ 3][2] = -theta;
     fjac[ 6][2] =  0.5;

     fjac[ 1][3] = -1.0;
     fjac[13][3] =  1.0;
 //  fjac[16][3] = xb[3]/sq(xb[12]);

     fjac[ 7][4] = -1.0;
     fjac[12][4] =  1.0;
 //  fjac[14][4] = xb[4]/sq(xb[13]);

     //fjac[ 1][5] =  2*theta/Cu[0]*xb[5]/xb[14] + gamma/Cu[0]*dFdQ(xb[14]);
     fjac[ 1][5] =  2.0*theta*xb[5]/xb[14] + gamma*dFdQ(xb[14]);
     fjac[ 4][5] =  theta;
     fjac[ 7][5] =  0.5;

     fjac[ 2][6] = -1.0;
     fjac[13][6] =  1.0;
 //  fjac[17][6] = xb[6]/sq(xb[15]);

     fjac[ 8][7] = -1.0;
     fjac[12][7] =  1.0;
 //  fjac[15][7] = xb[7]/sq(xb[16]);

     fjac[ 2][8] = 2.0*theta*xb[8]/xb[17] + gamma*dFdQ(xb[17]);
     fjac[ 5][8] = theta;
     fjac[ 8][8] = 0.5;

     fjac[ 3][9] = -1.0;
     if (xb[0] > 0)
     {
       fjac[16][9] = - dPdA(N,xb[9])
              + sq(xb[0])/cu(xb[9])*(-2.0*k7n);//+1);
       fjac[17][9] = - dPdA(N,xb[9])
              + sq(xb[0])/cu(xb[9])*(-2.0*k7an);//+1);
     } else
     {
       fjac[16][9] = - dPdA(N,xb[9])
              + sq(xb[0])/cu(xb[9])*(2.0*k7n);//+1);
       fjac[17][9] = - dPdA(N,xb[9])
              + sq(xb[0])/cu(xb[9])*(2.0*k7an);//+1);
     };
     fjac[9][10] = -1.0;
     if (xb[1] > 0)
     {
       fjac[14][10] = - dPdA(N,xb[10])
              + sq(xb[1])/cu(xb[10])*(-2.0*k7nh);//+1);
       fjac[15][10] = - dPdA(N,xb[10])
                  + sq(xb[1])/cu(xb[10])*(-2.0*k7anh);//+1);
     } else
     {
       fjac[14][10] = - dPdA(N,xb[10])
                  + sq(xb[1])/cu(xb[10])*(2.0*k7nh);//+1);
       fjac[15][10] = - dPdA(N,xb[10])
                  + sq(xb[1])/cu(xb[10])*(2.0*k7anh);//+1);
     };
     //fjac[ 0][11] = theta/Cu[N]*(  sq(xb[2]/xb[11]) - dBdAh(N,xb[11])) +
     //               gamma/Cu[N]*(dFdA(xb[2],xb[11]) + d2BdAdxh(N,xb[11]));
     fjac[ 0][11] = theta*(  sq(xb[2]/xb[11]) - dBdAh(N,xb[11])) +
                    gamma*(dFdA(xb[2],xb[11]) + d2BdAdxh(N,xb[11]));
     fjac[ 9][11] = 0.5;

     fjac[ 4][12] = -1.0;
     fjac[16][12] = LD->dPdA(0,xb[12]);// - sq(xb[3])/cu(xb[12]);

     fjac[10][13] = -1.0;
     fjac[14][13] = LD->dPdA(0,xb[13]);// - sq(xb[4])/cu(xb[13]);

     fjac[ 1][14] = theta*( -sq(xb[5]/xb[14]) + LD->dBdAh(-1,xb[14])) +
                    gamma*(dFdA(xb[5],xb[14]) + LD->d2BdAdxh(-1,xb[14]));
     fjac[10][14] = 0.5;

     fjac[ 5][15] = -1.0;
     fjac[17][15] = RD->dPdA(0,xb[15]);// - sq(xb[6])/cu(xb[15]);

     fjac[11][16] = -1.0;
     fjac[15][16] = RD->dPdA(0,xb[16]);// - sq(xb[7])/cu(xb[16]);

     fjac[ 2][17] = theta*( -sq(xb[8]/xb[17]) + RD->dBdAh(-1,xb[17])) +
                    gamma*(dFdA(xb[8],xb[17]) + RD->d2BdAdxh(-1,xb[17]));
     fjac[11][17] = 0.5;

     // Check whether solution is close enough. If not run the loop again.
     int ch = zero (xb, 18, 1.0e-12, 1.0e-12, fvec, fjac);

     if (ch == 1) ok = true;

     j = j+1;
   }

   // Solutions is applied, and right boundary is updated.
   Anew[N]     = xb[ 9];
   Qnew[N]     = xb[ 0];
   LD->Anew[0] = xb[12];
   LD->Qnew[0] = xb[ 3];
   RD->Anew[0] = xb[15];
   RD->Qnew[0] = xb[ 6];

   if (j >=ntrial) error ("arteries.C","Root not found in the bifurcation");
 }

 // The right boundary (x=L) uses this function to model an inflow into
 // the system. The actual parameter given to the function is the model time.
 // As stated in the mathematical model the constants of the function are
 // chosen in order to ensure a certain CO (specified in main.h). Hence we have
 // the specified value of b. Further the period (dimension-less) is assumed
 // to be Period.
 double Tube :: PL_init (double t, double k, double Period)
 {
   if (t <= Period) return (Pout[int(t/k)]); else
   if (t >  Period) return (PL_init((t-Period),k,Period));
   else return (0);
 }





 // Update of the right boundary at time t. This function uses PL to determine
 // A at the next time-step. From this the value of Q is predicted
 // using Lax-Wendroff's numerical scheme. This function is only relevant
 // when the tube is an outlet vessel.
 void Tube :: bound_right (double t, double k, double Period)
 {
   double Pnew   = PL_init(t,k,Period);
   //double P05    = (Pnew + PL_init(t-k,k,Period))/2;
   double P05    = PL_init(t-k,k,Period);

   if (int(t/k) < 0)
     printf("t/k negative in bound_right\n");

   Anew[0] = A0[0]/(sq(1 - (Pnew/fr[0])));

   double A05 = A0[0]/(sq(1 - (P05/fr[0])));
   double Ahm05 = 2.0*A05 - Ah[0];
   double Qhm05 = R1h[0] + h*(Anew[0] - Aold[0])/k;
   double R2hm05 = Rvec(2,-1,1,Qhm05,Ahm05);
   double S2hm05 = Svec(2,-1,1,Qhm05,Ahm05);

   Qnew[0] = Qold[0] - k*(R2h[0] - R2hm05)/h + k*(S2h[0] + S2hm05)/2.0;

 }
 // Solves the non-linear PDE's (momentum and continuity eqn's.
 // from t = tstart to t= tend.
 //
 // This function checks the maximal possible size of the next time-step,
 // reduces it to make sure that we don't walk to far, and takes the
 // next step. This is done by executing the step routine, then updating
 // the left boundary and finally updating bifurcation points and the
 // right boundaries. This is carried out as long as the time hasn't passed
 // the desired ending time (tend) which is passed to the function as a
 // parameter.
 void solver (Tube *Arteries[], double tstart, double tend, double k)
 {
   // The following definitions only used when a variable time-stepping is
   // used.

     

   double t    = tstart;
   int qLnb = (int) fmod(t/k,tmstps);

   // As long as we haven't passed the desired ending time do:
   while (t < tend)
   {
     // Check that the step we take is valid. If this error occurs when forcing
     // a constant step-size the program must be terminated.
     if (t+k > tend)
     {
       double kold = k;
       k = tend - t;
       printf("ERROR (arteries.C): Step-size changed, t+k=%10.15f, tend=%10.15f k=%10.15f kold=%10.15f\n",t+kold,tend,k,kold);
     }
     
     // Check that the CFL-condition applies.
     for (int i=0; i<nbrves; i++)
     {
//          fprintf(stdout,"CFL - vessel %d, k = %lf\n",i,k);
//            fflush(stdout);
       if (k > Arteries[i] -> CFL())
       {
         error("arteries.C","Step-size too large CFL-condition violated\n");
       }
     }
      

     // solve for interior points, by calling step.
     for (int i=0; i<nbrves; i++)
     {
		Arteries[i] -> step (k);
     }
//        fprintf(stdout,"Done with step\n");
//        fflush(stdout);
     // Update left and right boundaries, and the bifurcation points.
     for (int i=0; i<nbrves; i++)
     {
       if (Arteries[i] -> init == 1)
       {
         Arteries[i] -> bound_left(t+k, k, Period);
//         fprintf(stdout,"Done with inlet %d\n",i);
//        fflush(stdout);
         if (Arteries[i]->LD==0)
         {
            double theta = k/Arteries[i]->h;
            double gamma = k/2;
            Arteries[i] -> bound_match (qLnb, t, k, theta, gamma);
         }

       };

         // MJC 11/10/2022: want to get junction points before we impose pressure?
       if (Arteries[i] -> init == 2)
       {
         Arteries[i] -> bound_right(t+k, k, Period);
//            fprintf(stdout,"Done with pressure BC %d\n",i);
//            fflush(stdout);
       };
       
      
      
       if (Arteries[i] -> rm == 0)
       {
         double theta = k/Arteries[i]->h;
         double gamma = k/2;
         Arteries[i] -> bound_bif (theta, gamma);
//         fprintf(stdout,"Done with bif %d\n",i);
//         fflush(stdout);
         //If init == 4, then mono...
       };

         if (Arteries[i] -> init == 3)
       {
         double theta = k/Arteries[i]->h;
         double gamma = k/2;
         Arteries[i] -> bound_match (qLnb, t, k, theta, gamma);
//        fprintf(stdout,"Done with Admittance %d\n",i);
//        fflush(stdout);
       }
        
     }
     // Update the time and position within one period.
     t = t + k;
     qLnb = (qLnb + 1) % tmstps;
   }

}
