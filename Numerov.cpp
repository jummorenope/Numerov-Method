#include <iostream>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>
#include <iomanip>

void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, int N, double l2);
void Potential(std::vector <double>& k2,int N, double ep);
int g2 = 200;

int main (void)
{
    std::setprecision(10);
    int N=1000;
    double l2 = (1.0/(N-1))*(1.0/(N-1));
    double dep=0.1; //Variacion del ep
    double tol=1e-9;
    std::vector <double> k2(N);
    std::vector <double> Psi(N);
    std::vector <double> eps(10);
    double ep = -0.99;
    Psi[0]=0;
    Psi[1]=0.0001;
    double Psi1=0;
    //int ii=1;
    int write=1;
    std::string name="data";
    std::string ext=".txt";
    
    for( int ii=1;ii<=5;ii++)
    {

        Potential(k2, N, ep);
        WaveFunc(Psi,k2,N,l2);
        Psi1=Psi[N-1];
        
        while(std::fabs(Psi1)>tol)
        {
            ep+=dep;
            Potential(k2, N, ep);
            WaveFunc(Psi,k2,N,l2);
            if(Psi1*Psi[N-1]<=0)
            {
                dep=-dep/2.0;
            }
            Psi1=Psi[N-1];
        }
        std::ofstream data;
        data.open(name+std::to_string(ii)+ext);
        for(double jj=0; jj<N; jj++)
        {
            data <<jj/N<<"\t"<<Psi[jj]<<"\n";
        }
        data.close();
        std::cout<<ep<<std::endl;
        dep=0.1;
        ep+=dep;
    }
    
    

    return 0;
}

void Potential(std::vector <double>& k2, int N, double ep)
{
    for(int ii=0; ii<N; ii++)
    {
        k2[ii]=g2*(ep+1);
    }
}

void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, int N, double l2)
{
    for(int ii=2; ii<N; ii++)
    {
        Psi[ii] = (2*(1-(5.0/12)*l2*k2[ii-1])*Psi[ii-1]-(1+(1.0/12)*l2*k2[ii-2])*Psi[ii-2])/(1+(1.0/12)*l2*k2[ii]);
    }
}
