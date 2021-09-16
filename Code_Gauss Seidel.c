#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define H 2.0					  // Height of flat plate
#define L 1.0    				  // Length of flat plate
#define I_max 21				  // Maximum vertical gridlines
#define J_max 41				  // Maximum horizontal gridlines
#define Error_mx 0.01
int main() 
{
	int i, j, k=0;
	double dx = L/(I_max-1), dy = H/(J_max-1), beta, Rndm,err;
	beta = dx/dy;
	double x[210], y[210];
	double T1 = 100;          // Bottom wall temp.
	double T2 = 0;          // Left side wall temp.
	double T3 = 0;          // Top wall temp.
	double T4 = 0;           // Right side wall temp.
	double T[210][210], Tn[210][210];
	// initial conditions
	for(i=0; i<I_max; i++)
	{
		for(j=0; j<J_max; j++)
		{
			if(j==0)
			{
				T[i][j]=T1;
			}
			else if(j==J_max-1)
			{
				T[i][j]=T3;
			}
			else if(i==0)
			{
				if(j==0)
				{
					T[i][j]=(T1+T2)/2;
				}
				if(j==J_max-1)
				{
					T[i][j]=(T2+T3)/2;
				}
				else
				{
					T[i][j]=T2;
				}
			}
			else if(i==I_max-1)
			{
				if(j==0)
				{
					T[i][j]=(T1+T4)/2;
				}
				if(j==J_max-1)
				{
					T[i][j]=(T4+T3)/2;
				}
				else
				{
					T[i][j]=T4;
				}
			}
			else
			{
				T[i][j]=0;
			}
			Tn[i][j]=T[i][j];
		}
	}
	// Point Gauss-Seidel iteration
	Rndm=0;
	err=1;
	while(err > Error_mx)
	{
		k = k+1;
		for(i=1; i < I_max-1; i++)
		{
			for(j=1; j < J_max-1; j++)
			{
				Tn[i][j]=((Tn[i+1][j]+Tn[i-1][j])+(pow(beta,2)*(Tn[i][j+1]+Tn[i][j-1])))/(2*(1+pow(beta,2)));
				Rndm = Rndm + fabs(Tn[i][j]-T[i][j]);
				T[i][j]=Tn[i][j];
			}	
		}
		err = Rndm;
		Rndm = 0;
	}
	FILE *fp;
	fp = fopen("gauss_seidel.txt","w");
	printf("Node-x\t\t Node-y\t\tTemp.\n");
	for(j=0; j < J_max; j++)
	{
		for(i=0; i < I_max; i++)
		{
			x[i]=i*dx;
			y[j]=j*dy;
			printf("%0.3lf\t\t%0.3lf\t\t%lf\n", x[i], y[j], T[i][j]);
			fprintf(fp,"%0.3lf\t\t%0.3lf\t\t%lf\n", x[i], y[j], T[i][j]);
		}
	}
	printf("\nTotal iteration = %d", k);
	return 0;
}
