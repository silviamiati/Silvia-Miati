/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"
#include <string>

//parameters, observables
Random rnd;
const int m_props=1000;
bool restart,Realsimulation;
int n_props, L;
double bin_size,nbins;
int iv,ik,it,ie,ip, F, E;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_p,stima_pr,stima_k,stima_t,stima_et,err_k,err_t,err_et,err_p,err_pr,err_gdir;


//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part],x_n[m_part],y_n[m_part],z_n[m_part];
double vx[m_part],vy[m_part],vz[m_part],vx_n[m_part],vy_n[m_part],vz_n[m_part] ;

// thermodynamical state
int npart;
int nblk,igofr;
int block=1;
double energy,temp,vol,rho,box,rcut; 
double ave_pot,ave_kin,ave_etot,ave_temp,ave_pot2,ave_kin2,ave_temp2,ave_etot2,ave_pres,ave_pres2;
 
//characteristic type of gas
double massa, fattoret;
 
// simulation
int nstep, iprint, seed;
double delta;

const double pi=3.1415927;
double sigma; 
double kb;
double epsilon_kb;


//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Accumulate(void);
void Block(int numero, int i);
double Error(double, double, int);
void Averages(int);
void Reset(int);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
