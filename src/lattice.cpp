#include "../include/lattice.hpp"
#include <cstdlib>

void lattice::init(int Number, int occupancy)
{
	N=Number;
	N_occ=occupancy;

	site Null;
	Null.occ=0;
	Null.angle=0.0;

	spins.resize(Number);
	for( auto &it : spins )
	{
		it.resize(Number, Null);
	}

	int n=0;
	while (n!=occupancy)
	{
		int i=rand() % Number;
		int j=rand() % Number;

		if (spins[i][j].occ==0)
		{	
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			//std::cout<<i<<" "<<j<<" "<<theta<<std::endl;

			spins[i][j].occ=1;
			spins[i][j].angle=theta;

			n++;
		}
	}
}

void lattice::set_const(double j, double k, double frustration)
{
	J=j;
	K=k;
	f=frustration*3.1415962;
}

void lattice::flip(int i, int j)
{
	if (spins[i][j].occ==0) {spins[i][j].occ=1;}
	else {spins[i][j].occ=0;}
}

void lattice::rotate(int i, int j, double theta)
{
	spins[i][j].angle=theta;
}

int lattice::occ(int i, int j)
{
	return spins[i][j].occ;
}

int lattice::how_many()
{
	return N;
}

double lattice::angle(int i, int j)
{
	return spins[i][j].angle;
}

double lattice::H()
{
	double H=0.0;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (i!=N-1)
			{
				int weight=spins[i][j].occ*spins[i+1][j].occ;
				H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; // Still have to add external field...
				H+=K*weight;
			}

			if (j!=N-1)
			{
				int weight=spins[i][j].occ*spins[i][j+1].occ;
				H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; // Still have to add external field...
				H+=K*weight;
			}
		}
	}

	return H;
}

double lattice::H_periodic()
{
	double H=0.0;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (i!=N-1)
			{
				int weight=spins[i][j].occ*spins[i+1][j].occ;
				H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; // Still have to add external field...
				H+=K*weight;
			}
			else
			{
				int weight=spins[i][j].occ*spins[0][j].occ;
				H+=J*cos(spins[i][j].angle-spins[0][j].angle+f*j)*weight; // Still have to add external field...
				H+=K*weight;
			}

			if (j!=N-1)
			{
				int weight=spins[i][j].occ*spins[i][j+1].occ;
				H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; // Still have to add external field...
				H+=K*weight;
			}
			else
			{
				int weight=spins[i][j].occ*spins[i][0].occ;
				H+=J*cos(spins[i][j].angle-spins[i][0].angle-f*i)*weight; // Still have to add external field...
				H+=K*weight;
			}
		}
	}

	return H;
}

double lattice::Bond_Energy(int i, int j, int n, int m)
{
	double B=0.0;

	if (i==n and j==m+1)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle+f*i)*weight;
		B+=K*weight;
	}
	else if (i==n and j==m-1)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle-f*i)*weight;
		B+=K*weight;
	}
	else if (i==n+1 and j==m)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle-f*j)*weight;
		B+=K*weight;
	}
	else if (i==n-1 and j==m)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle+f*j)*weight;
		B+=K*weight;
	}

	return B;
}
