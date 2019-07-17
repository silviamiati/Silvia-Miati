/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000,m_bins=100;
int n_props, ie, ih=0, igofr, igofp, length, nbins;
double bin_size;
double walker[m_props],norm[m_props];

bool restart;

// averages
double blk_av[m_props],histo[m_bins],blk_norm,accepted,attempted,d;
double glob_av[m_props],glob_av2[m_props], glob_is[m_bins], glob_is2[m_bins];
double stima_e,err_e,err_gdir,stima_h;

//configuration
double xold,mu,sigma2;

// thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk;
double delta;

//pigreco
const double pi=3.1415927;
const double sigma=0.34*1E-9;
const double kb=1.38054*1E-23;
const double epsilon_kb=120;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(void);
void Measure(void);
double Wave2(double, double, double);
double Wave(double, double, double);
double norma(double, double, double);
double Pbc(double);
double Error(double,double,int);
double Acceptance();
void Histo();
void EUltimo();
void Eblock(int);



#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
