#include <stdio.h>
#include <math.h>
float determinant(float [][10000], float);
void cofactor(float [][10000], float);
void transpose(float [][10000], float [][10000], float);


/*Determinant of the Matrix */
float determinant(float a[10000][10000], float k)
{
  float s = 1, det = 0, b[10000][10000];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 
    return (det);
}
 
void cofactor(float num[10000][10000], float f)
{
 float b[10000][10000], fac[10000][10000];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f);
}
/*transpose of matrix*/ 
void transpose(float num[10000][10000], float fac[10000][10000], float r)
{
  int i, j;
  float b[10000][10000], inverse[10000][10000], d;
 
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of matrix is : \n");
  
   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        isnan(inverse[i][j])? printf("\t%f",0);
       printf("\t%f", inverse[i][j]);
        }
    printf("\n");
     }
}
int main()
{
    float Vinf,R,dtheta,alpha,A,B,Cn,Ct,Dn,Dt,E,term1,term2;
    float n;
    float beta[10000],phi[10000],midpoint_x[10000],
    midpoint_y[10000],S[10000], X[10000],Y[10000],theta[10000];
    float I[10000][10000],J[10000][10000];
    const double pi = 4.0 * atan(1.0);
    int index=0,i,j;
    Vinf=30;
    R=1;
    n=8;
    dtheta=2*pi/n;
    alpha=0;
    do
    {
        theta[index]=index*dtheta;
        index++;
    } while (index<n);
    index=0;
    for (index=0;index<n;index++)
    {
        X[index]=R*cos(theta[index]);
        Y[index]=R*sin(theta[index]);
        printf("%f,%f\n",X[index],Y[index]);
    }
    printf("\n\n");
    for (index=0;index<n;index++)
    {
        phi[index]=-alpha+atan2((Y[index+1]-Y[index]),(X[index+1]-Y[index]));
        beta[index]=phi[index]+pi/2;
        midpoint_x[index]=(X[index+1]+X[index])/2;
        midpoint_y[index]=(Y[index+1]+Y[index])/2;
        S[index]=pow((((Y[index+1]-Y[index])*(Y[index+1]-Y[index]))+((X[index+1]-X[index])*(X[index+1]-X[index]))),(1/2));
        printf("%lf,%lf\n",midpoint_x[index],midpoint_y[index]);
    }
    printf("\n\n");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if (j != i)
            {
                A  = -(midpoint_x[i]-X[j])*cos(phi[j])-(midpoint_y[i]-Y[j])*sin(phi[j]);
                B  = ((midpoint_x[i]-X[j])*(midpoint_x[i]-X[j])) + ((midpoint_y[i]-Y[j])*(midpoint_y[i]-Y[j])) ;
                Cn = sin(phi[i]-phi[j]) ;
                Dn = (midpoint_x[i]-X[j])*sin(phi[i])+(midpoint_y[i]-Y[j])*cos(phi[i])  ;
                Ct = -cos(phi[i]-phi[j]);
                Dt = (midpoint_x[i]-X[j])*cos(phi[i])+(midpoint_y[i]-Y[j])*sin(phi[i]) ;
                E  = sqrt(B-(A*A)) ;

                term1  = 0.5*Cn*log(((S[j]*S[j] )+ 2*A*S[j] + B)/B);
                term2  = ((Dn-A*Cn)/E)*(atan2((S[j]+A),E)-atan2(A,E)) ;
                I[i][j] = term1 + term2;
 
                term1  = 0.5*Ct*log(((S[j]*S[j] )+ 2*A*S[j] + B)/B);
                term2  = ((Dt-A*Ct)/E)*(atan2((S[j]+A),E)-atan2(A,E)) ;
                J[i][j] = term1 + term2;
            }
 
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("\t%f,%f",I[i][j],J[i][j]);
        }
        printf("\n");
    }

     float a [10000][10000];                                                  

    for (i=0;i<n;i++)
    {    
        for (j=0;j<n;j++)                                                 
        {
            if (i == j)                                                            
            {
            a[i][j] =pi;
            }                                                      
            else
            {                                                                  
            a[i][j] = I[i][j];
            }
        }

    }
    for (i = 0;i < n; i++)
    {
     for (j = 0;j < n; j++)
       {
         printf("\t%f", a[i][j]);
        }
    printf("\n");
    }
    float d = determinant(a, n);
  if (d == 0)
  { printf("\nInverse of Entered Matrix is not possible\n");}
  else
   {cofactor(a, n);}
}

