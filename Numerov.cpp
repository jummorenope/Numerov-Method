#include <iostream>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>

void WaveFunc (std::vector <double>& Psi, std::vector <double> k2, int N, double l2);
void Potential(std::vector <double>& k2,int N, double ep);
int g2 = 200;

int main (void)
{
    int N=1000;
    double l2 = (1.0/(N-1))*(1.0/(N-1));
    double dep=0.01; //Variacion del ep
    double tol=1e-9;
    std::vector <double> k2(N);
    std::vector <double> Psi(N);
    std::vector <double> eps(10);
    double ep = -0.99;
    Psi[0]=0;
    Psi[1]=0.0001;
    double Psi1=0;
    int ii=1;
    int write=0;
    std::string name="data";
    std::string ext=".txt";
    while(ii<=10)
    {
        Potential(k2, N, ep);
        WaveFunc(Psi,k2,N,l2);
        //std::cout<<ii<<std::endl;
        
        if(std::fabs(Psi[N-1])>tol&&write)
        {
            //std::cout<<"hola"<<std::endl;
            if(std::signbit(Psi1)!=std::signbit(Psi[N-1])){
                ep-=dep;
                dep=dep/10.0;
            }
            ep+=dep;
            Psi1=Psi[N-1];
        }
        
        if(write!=1)
        {
            Psi1=Psi[N-1];
            write=1;
        }
        
        if(std::fabs(Psi[N-1])<tol)
        {
            std::string num=std::to_string(ii);
            std::string file=name+num+ext;
            char cstr[file.size() + 1];
            strcpy(cstr, file.c_str());
            std::ofstream data;
            data.open(file);
            for(double jj=0; jj<N; jj++)
            {
                data <<jj/N<<"\t"<<Psi[jj]<<"\n";
            }
            data.close();
            std::cout<<ep<<std::endl;
            ii+=1;
            dep=0.001;
            eps[ii-1]=ep;
            ep+=dep;
            write=0;
            
            
        }
        
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
