#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "Tangent.h"
#include <math.h>
#include <iomanip>
#include <fstream>


template <class T>
T tersoff_fc(T R0[][3],int total_atoms)
{

	//Lindsay and Broido, Phys. Rev. B 81, 205441 – Published 27 May 2010
	double A(1.3936e3);
	double B(430);
	double lam1(3.4879);
	double lam2(2.2119);
	double beta(1.5724e-7);
	double n(7.2751e-1);
	double c(38049);
	double d(4.3484);
	double h(-0.930);;
	double R(1.95);
	double D(0.15); //Cutoff = R+D

	T rij[3], r, fc, rik[3], r2, rijk, fcik, gijk, psi, bij, aij, fR, fA, Vi1;

	T Etot(0);

	int num_neigh = total_atoms;

	for(int i1=0; i1<total_atoms; i1++)
	{
		Vi1=0;
		for (int i2=0; i2<total_atoms; i2++)
		{
			if(i2!=i1)
			{
				//Distance
				rij[0]=R0[i1][0]-R0[i2][0]; rij[1]=R0[i1][1]-R0[i2][1]; rij[2]=R0[i1][2]-R0[i2][2];
				r = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                
				//Cutoff
				fc=0.5*(1-sin(pi/2*(r-R)/D))*(fabs(r-(R-D))/(r-(R-D))+fabs(R+D-r)/(R+D-r))/2+(fabs(R-D-r)/(R-D-r)+1)/2;

				aij=1;
                
				//Calc attractive and repulsive parts
				fR=A*exp(-lam1*r);
				fA=-B*exp(-lam2*r);

				if (fc==0)
					continue; 
				//Cal psi and eta for all neighbors
				psi=0;

				if (fA == 0)
					bij = 0;
				else
				{
					//Loop over all atoms to compute bij term which depends on coordination between atoms i,j,k
					for(int i3=0; i3<total_atoms; i3++)
					{
						if(i3!=i1 && i3!=i2)
						{
							rik[0]=R0[i1][0]-R0[i3][0]; rik[1]=R0[i1][1]-R0[i3][1]; rik[2]=R0[i1][2]-R0[i3][2];
							r2=sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);
							fcik=0.5*(1-sin(pi/2*(r2-R)/D))*(fabs(r2-(R-D))/(r2-(R-D))+fabs(R+D-r2)/(R+D-r2))/2+(fabs(R-D-r2)/(R-D-r2)+1)/2;
							rijk=rij[0]*rik[0]+rij[1]*rik[1]+rij[2]*rik[2];
							gijk=1+pow(c,2)/pow(d,2)-pow(c,2)/(pow(d,2)+pow(h-rijk/(r*r2),2));
							psi=psi+fcik*gijk;
						}
					}
					bij=pow((1+pow(beta,n)*pow(psi,n)),(-1/(2*n)));
				}

				//Potential contribution
				Vi1=Vi1+fc*(aij*fR+bij*fA);
			}
		}
		Etot=Etot+Vi1/2.0;
	}
	return Etot; 
}
/*
   %Zigzag gnr
   %      o  o      o
   %    o      o  o    3
   %      o  o      o
   %    o      o  o    2
   %      o  o      o
   %    o      o  o    1
   %      o  o      o
   %    o      o  o    4
   %      o  o      o
   %    o      o  o    5
*/
int main(void)
{
	double a=1.475; //e-10;
	double theta = 60*pi/180; //Angle in the hexagon
	int num_uc=5;
	int num_at_per_uc=4;

	int total_atoms=num_uc*num_at_per_uc;
	double coord[4][3][num_uc];
	int dof=3; //x,y,z

	//Initialize Tangent variables
	//Force constant is the second derivative of the total lattice energy with respect to atomic displacements
	//Both the lattice energy (defined in the user defined function tersoff_fc) and atomic displacements should be Tangent variables
	//The output variable fc, where the force constants are stored, is also a Tangent variable
    
	Tangent R0[total_atoms][3]; /* Coordinate data type declared as Tangent instead of traiditional double */
	Tangent fc[num_at_per_uc*dof][num_at_per_uc*dof][num_uc]; /* Force forstant data type declared as Tangent instead of traditional double */

	ofstream forceConstants("force_constant_data_gnr_4_atoms.csv");

	//Generate GNR Coordinates for the center unit cell
	//Variable coord of size coord[num of atoms per unit cell][degrees of freedom][number of unit cells]
    //coordinate generation code omitted for clarity

	//Dynamical Matrix
	int deriv_index1,deriv_index1_coord,deriv_index2,deriv_index2_coord;

	for (int j=0; j<num_at_per_uc*dof; j++)
	{
			int index=0;
			for (int n=0; n<num_uc; n++)
			{
				for (int n1=0; n1<num_at_per_uc; n1++)
				{
                    /* First field of Tangent variable R0 assigned the value of coord*/
					R0[index][0]._v=coord[n1][0][n];
					R0[index][1]._v=coord[n1][1][n];
					R0[index][2]._v=coord[n1][2][n];
					index++;
				}
			} 
			for (int k=0; k<num_at_per_uc*dof; k++)
			{
				for (int inner=0; inner<num_uc; inner++)
				{
					deriv_index1 = floor(j/3);
					deriv_index1_coord = j%3;
					deriv_index2 = floor(k/3);
					deriv_index2_coord = k%3;
					R0[deriv_index1][deriv_index1_coord]._dv1 = 1; /* Define the flag as 1 for the coordinate with which the first partial derivative is needed */
					R0[deriv_index2+num_at_per_uc*inner][deriv_index2_coord]._dv2 = 1; /* Define the flag as 1 for the coordinate with which the second partial derivative is needed */
					fc[j][k][inner] = tersoff_fc<Tangent>(R0, total_atoms); /* Call the tersoff_fc function */
				}
			}
	} 


	return 0;

}


