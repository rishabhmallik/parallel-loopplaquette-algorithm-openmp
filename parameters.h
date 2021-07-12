

float eta=4.01,lmd=1.0;   //theory parameters
const int d=2, Ns=10, Nt=10;    //lattice size and dimensions
float dr=0.01,inf=10; //simulation parameters
const int configs=400,gaps=10,equil=200;    //simulation parameters
const int mu_n=100, int_val=100;
long double rho, I_val[int_val], mu_min=0.0, mu_max=1.0, mu;

float n_avg;
int threads;
int del;
int k[Nt][Ns][d]={0}, a[Nt][Ns][d]={0}
    , a_[Nt][Ns][d]={0}
    , x[d], x1[d], x2[d], x12[d], x_[d];
