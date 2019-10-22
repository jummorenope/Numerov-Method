#include <iostream>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cstring>

void PrintVector( const std::vector <int>& a);
void recibirVector( std::vector <int>& a);
void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, int N, double l2);
void Potential(std::vector <double>& k2,int N, double ep);
int g2 = 200;

int main (void)
{
    int N=1000;
    double l2 = (1.0/(N-1))*(1.0/(N-1));
    double dep=0.001; //Variacion del ep
    double tol=1e-05;
    std::vector <double> k2(N);
    std::vector <double> Psi(N);
    double ep = -0.99;
    Psi[0]=0;
    Psi[1]=0.01;
    int ii=1;
    std::string name="data";
    std::string ext=".txt";
    while(ii<=10)
    {
        Potential(k2, N, ep);
        WaveFunc(Psi,k2,N,l2);
        if(Psi[N-1]<tol)
        {
            std::string num=std::to_string(ii);
            std::string file=name+num+ext;
            char cstr[file.size() + 1];
            strcpy(cstr, file.c_str());
            freopen(cstr,"w",stdout);
            for(double jj=0; jj<N; jj++)
            {

                std::cout<<jj/N<<"\t"<<Psi[jj]<<std::endl;
            }
            ii+=1;
        }
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
