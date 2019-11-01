#include <iostream>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>
#include <iomanip>

void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, double l2);
void WaveNumber (std::vector <double>& k2, std::vector <double> V, double ep);
void Potential (std::vector <double>& V);
int g2 = 200;
int N=1000;

int main (void)
{
    std::cout<<std::setprecision(10);
    double l = (1.0/(N-1));
    double l2 = (1.0/(N-1))*(1.0/(N-1));
    double dep=0.1; //Variacion del ep
    double tol=1e-9;
    std::vector <double> V(N);
    Potential(V);
    std::vector <double> k2(N);
    std::vector <double> Psi(N);
    double ep = -0.99;
    Psi[0]=0;
    Psi[1]=0.0001;
    double ept=0;
    double terror=0;
    double Psi1=0;
    std::string name="data";
    std::string ext=".txt";
    for( int ii=1;ii<=10;ii++)
    {

        double norm=0.0;
        WaveNumber(k2, V, ep);
        WaveFunc(Psi,k2,l2);
        Psi1=Psi[N-1];
        while(std::fabs(Psi1)>tol)
        {
            ep+=dep;
            WaveNumber(k2, V, ep);
            WaveFunc(Psi,k2,l2);
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
            
            data <<jj/N<<"\t"<<Psi[jj]<<"\t"<<Psi[jj]*Psi[jj]<<"\n";
        }
        data.close();
        
        ept=ii*ii*M_PI*M_PI/(g2)-1.0;
        terror=std::fabs(100.0*(ept-ep)/ept);
        std::cout<<ep<<std::endl;
        //std::cout<<ep<<"\t"<<ept<<"\t"<<terror<<'%'<<std::endl;
        dep=0.1;
        ep+=dep;
    }
    return 0;
}

void WaveNumber(std::vector <double>& k2, std::vector <double> V, double ep)
{
    for(int ii=0; ii<N; ii++)
    {
        k2[ii]=g2*(ep-V[ii]);
    }
}

void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, double l2)
{
    for(int ii=2; ii<N; ii++)
    {
        Psi[ii] = (2*(1-(5.0/12)*l2*k2[ii-1])*Psi[ii-1]-(1+(1.0/12)*l2*k2[ii-2])*Psi[ii-2])/(1+(1.0/12)*l2*k2[ii]);
    }
}

void Potential (std::vector <double>& V)
{
    for(double ii=0;ii<N;ii++)
    {
        V[ii]=-1.0;
        //V[ii]=30*(ii/N-0.5)*(ii/N-0.5)-1.0;
        //V[ii]=8*(ii/N-0.5)*(ii/N-0.5)*(ii/N-0.5)*(ii/N-0.5)-1.0;
    }
}
