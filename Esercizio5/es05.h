#include "random.h"


int seed[4];
Random rnd;

double xin, yin, zin;
double xpre, ypre, zpre;
double delta;

bool restart;

const int m_props=10000;
int block=1;
int n_props, ir;
int nblock, nstep;

double ave_raggio,media;
double blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props],err_r;

bool passo, gr;

const double pi=3.1415927;

void Input(void);
double Funz(double, double, double);
void Block(int);
double Error(double,double,int);
void Reset();
void Move(void);
void Measure(void);
double Acceptance();
void ConfFinal(void);
void ConfXYZ(void);

