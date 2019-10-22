#include <iostream>
#include <vector>
#include <cstdio>

void PrintVector( const std::vector <int>& a);
void recibirVector( std::vector <int>& a);
void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, int N, double l2);
void Potential(std::vector <double>& k2,int N, double ep);
int g2 = 200;

int main (void)
{
    int N=1000;
    std::vector <double> Psi(N);
    double ep = -0.9;
    std::vector <double> k2(N);
    Potential(k2, N, ep);
    double l2 = (1.0/(N-1))*(1.0/(N-1));
    Psi[0]=0;
    Psi[1]=0.0001;
    WaveFunc(Psi,k2,N,l2);
    freopen("data.txt","w",stdout);
    for(double ii=0.0; ii<N; ii++)
    {
        std::cout<<ii/N<<"\t"<<Psi[ii]<<std::endl;
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
