#include<iostream>
#include<cmath>
#include <iomanip>
#include <fstream>
using namespace std;
const float pi = std::acos(-1);
#define N 80

double matrixsolver(float A[N][N], float B[N], float x[N],int n)
{

    int i,j,flag;
    float sum,xn[N];

    for(i=0;i<n;i++)
    {
        x[i]=0;
    }
    do{
        for(i=0;i<n;i++){
            sum=B[i];
            for(j=0;j<n;j++){
                if(j!=i){
                    sum-=(A[i][j])*x[j];
                }
            }
            xn[i]=sum/(A[i][i]);
        }
        flag=0;
        for(i=0;i<n;i++){
            if(abs(x[i]-xn[i])<0.00001){
                flag++;
            }
        }
        for(i=0;i<n;i++){
            x[i]=xn[i];
        }


    }while(flag<n);
    for(i=0;i<n;i++)
    {
        return x[i];
    }
}
int main()
{
    float dx,dy,s[N], theta[N],mid_x[N],mid_y[N],phi[N],delta[N],beta[N],I[N][N],J[N][N],A[N][N],Z[N],lamda[N],Vt[N],Cp[N];
    float dtheta,r,alhpa,term1,term2,D,E,Q,B,C;
    alhpa=0;
    int i=0,j=0;
    int n=60;
    float x[N],y[N];
    ifstream file;
    file.open("x.txt");

    if(file.is_open())
    {


        for(int i = 0; i < 80; i++)
        {
            file >> x[i];
        }
    }
    file.close();
    ifstream file2;
    file2.open("y.txt");

    if(file2.is_open())
    {


        for(int i = 0; i < 80; i++)
        {
            file2 >> y[i];
        }
    }
    file2.close();
    cout<<"X Y coordinates of control points are:\n";
    for(i=0;i<n;i++)
    {
        mid_x[i]=0.5*(x[i]+x[i+1]);
        mid_y[i]=0.5*(y[i]+y[i+1]);

        dx=x[i+1]-x[i];
        dy=y[i+1]-y[i];
        phi[i]  =atan2(dy,dx);
        s[i]=abs(sqrt((dx*dx)+(dy*dy)));
    }
    for(i=0;i<n;i++)
    {
        if(abs(mid_x[i])<0.000001)
        {
            mid_x[i]=0;

        }
        if(abs(mid_y[i])<0.000001)
        {
            mid_y[i]=0;

        }
        cout<<mid_x[i]<<","<<mid_y[i]<<" "<<"      "<<"\n";
    }


    for(i=0;i<n;i++)
    {
        delta[i]=phi[i]+pi/2;
        beta[i]=delta[i]-alhpa;
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if (j != i)
            {
                Q  = -(mid_x[i]-x[j])*cos(phi[j])-(mid_y[i]-y[j])*sin(phi[j]);
                B  = ((mid_x[i]-x[j])*(mid_x[i]-x[j])) + ((mid_y[i]-y[j])*(mid_y[i]-y[j]));
                C = sin(phi[i]-phi[j]);
                D = -(mid_x[i]-x[j])*sin(phi[i])+(mid_y[i]-y[j])*cos(phi[i]);
                E  = sqrt(abs(B-(Q*Q)));
                if(isnan(E)==true || E==0)
                {
                    I[i][j]=0;
                }
                else
                {
                term1  = 0.5*C*log(((s[j]*s[j] )+ 2*Q*s[j] + B)/B);
                term2  = ((D-Q*C)/E)*(atan2((s[j]+Q),E)-atan2(Q,E)) ;
                I[i][j] = -1*(term1 + term2);
                }

            }
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(j==i)
            {
                A[i][j]=pi;
            }
            else
            {
                A[i][j]=I[i][j];
            }
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
                Q  = -(mid_x[i]-x[j])*cos(phi[j])-(mid_y[i]-y[j])*sin(phi[j]);
                B  = ((mid_x[i]-x[j])*(mid_x[i]-x[j])) + ((mid_y[i]-y[j])*(mid_y[i]-y[j]));
                C = sin(phi[i]-phi[j]);
                D = -(mid_x[i]-x[j])*sin(phi[i])+(mid_y[i]-y[j])*cos(phi[i]);
                E  =abs(sqrt(B-(Q*Q)));
            if(isnan(E)==true)
                {

                    J[i][j]=0;
                }
                else
                {
                term1  = 0.5*((D-Q*C)/E)*log(((s[j]*s[j] )+ 2*Q*s[j] + B)/B);
                term2  = C*(atan2((s[j]+Q),E)-atan2(Q,E));
                J[i][j] = term1 - term2;
                }

        }
    }
    float Vinf=1;
    float adj[N][N];
     for (i = 0;i < n; i++)
     {
         Z[i]= -Vinf*2*pi*cos(beta[i]);
     }
    matrixsolver(A,Z,lamda,n);



    cout<<"\n Cp value are:\n";
    for (i = 0;i < n; i++)
    {
        float val = 0;
        for (j = 0;j < n; j++)
        {
        val = val + (lamda[j]/(2*pi))*J[i][j];
        }
        Vt[i] = Vinf*sin(beta[i]) + val;


        }
    for(i=0;i<n;i++)
    {
        Cp[i] = 1 - (Vt[i]/Vinf)*(Vt[i]/Vinf);
        cout<<Cp[i]<<","<<x[i]<<"\n";
    }
return 0;
}
