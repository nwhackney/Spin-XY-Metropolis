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

void lattice::restart(int Number, int occupancy, std::string infile)
{

	N=Number;
	N_occ=occupancy;
	
	std::ifstream in;
	in.open(infile);

	if (!in)
	{
     	std::cout << "Unable to open file";
     	exit(1); // terminate with error
     }

     std::vector<double> raw;
     std::string x;

     while (in>>x)
     {
     	std::stringstream item;
     	item<<x;
     	double i=0;
     	item >> i;
     	raw.push_back(i);
     }

     in.close();

     site Null;
	Null.occ=0;
	Null.angle=0.0;

	spins.resize(Number);
	for( auto &it : spins )
	{
		it.resize(Number, Null);
	}

     for (int i=0; i<=raw.size()-4; i+=4)
     {
          if (raw[i+2]!=-1.0)
          {
          	int n=int(raw[i]);
          	int m=int(raw[i+1]);
          	double theta=double(raw[i+2]);
               spins[n][m].occ=1;
               spins[n][m].angle=theta;
          }
     }
}

void lattice::rand_square_init(int Number, int occupancy)
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

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+1][j+1].occ=1;
			spins[i+1][j+1].angle=theta;
		}
	}

	for (int i=0; i<13; i++)
	{
		for (int j=0; j<13; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+22][j+1].occ=1;
			spins[i+22][j+1].angle=theta;
		}
	}

	for (int i=0; i<15; i++)
	{
		for (int j=0; j<15; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+43][j+1].occ=1;
			spins[i+43][j+1].angle=theta;
		}
	}

	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+64][j+1].occ=1;
			spins[i+64][j+1].angle=theta;
		}
	}

	for (int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+1][j+22].occ=1;
			spins[i+1][j+22].angle=theta;
		}
	}

	for (int i=0; i<18; i++)
	{
		for (int j=0; j<18; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+22][j+22].occ=1;
			spins[i+22][j+22].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+43][j+22].occ=1;
			spins[i+43][j+22].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+64][j+22].occ=1;
			spins[i+64][j+22].angle=theta;
		}
	}	

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+1][j+43].occ=1;
			spins[i+1][j+43].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+22][j+43].occ=1;
			spins[i+22][j+43].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+43][j+43].occ=1;
			spins[i+43][j+43].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+64][j+43].occ=1;
			spins[i+64][j+43].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+1][j+64].occ=1;
			spins[i+1][j+64].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+22][j+64].occ=1;
			spins[i+22][j+64].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+43][j+64].occ=1;
			spins[i+43][j+64].angle=theta;
		}
	}

	for (int i=0; i<20; i++)
	{
		for (int j=0; j<20; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+64][j+64].occ=1;
			spins[i+64][j+64].angle=theta;
		}
	}

}

void lattice::square_init(int Number, int length)
{
	N=Number;
	N_occ=length*length;

	site Null;
	Null.occ=0;
	Null.angle=0.0;

	spins.resize(Number);
	for( auto &it : spins )
	{
		it.resize(Number, Null);
	}

	for (int i=1; i<length; i++)
	{
		for (int j=1; j<length; j++)
		{
			double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
			spins[i+20][j+20].occ=1;
			spins[i+20][j+20].angle=theta;
		}
	}
}

void lattice::circle(int Number, int occupancy, double R)
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

	double Nd=(double) N;

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (((double) i-Nd/2.0)*((double) i-Nd/2.0)+((double) j-Nd/2.0)*((double) j-Nd/2.0) <= R*R)
			{
				spins[i][j].angle=3.14;
				spins[i][j].occ=1;
			}
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

int lattice::spin_num()
{
	return N_occ;
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
				 H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; // Original Gauge
				//H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*j)*weight; // New -trial- Gauge
				H+=K*weight;
			}

			if (j!=N-1)
			{
				int weight=spins[i][j].occ*spins[i][j+1].occ;
				H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; //Original Gauge
				//H+=J*cos(spins[i][j].angle-spins[i][j+1].angle)*weight; // New -trial- Gauge
				H+=K*weight;
			}
		}
	}

	return H;
}

double lattice::H_local(int i, int j)
{
	double H=0.0;

	if (i!=0)
	{
		int weight=spins[i][j].occ*spins[i-1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[i-1][j].angle-f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i-1][j].angle-2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	if (i!=N-1)
	{
		int weight=spins[i][j].occ*spins[i+1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	if (j!=0)
	{
		int weight=spins[i][j].occ*spins[i][j-1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][j-1].angle+f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][j-1].angle)*weight; //New Gauge
		H+=K*weight;
	}
	if (j!=N-1)
	{
		int weight=spins[i][j].occ*spins[i][j+1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][j+1].angle)*weight; //New Gauge
		H+=K*weight;
	}

	return H;
}

double lattice::H_local_periodic(int i, int j)
{
	double H=0.0;

	if (i!=0)
	{
		int weight=spins[i][j].occ*spins[i-1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[i-1][j].angle-f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i-1][j].angle-2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	else if (i==0)
	{
		int weight=spins[i][j].occ*spins[N-1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[N-1][j].angle-f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[N-1][j].angle-2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	if (i!=N-1)
	{
		int weight=spins[i][j].occ*spins[i+1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	else if(i==N-1)
	{
		int weight=spins[i][j].occ*spins[0][j].occ;
		H+=J*cos(spins[i][j].angle-spins[0][j].angle+f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[0][j].angle+2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	if (j!=0)
	{
		int weight=spins[i][j].occ*spins[i][j-1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][j-1].angle+f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][j-1].angle)*weight; //New Gauge
		H+=K*weight;
	}
	else if (j==0)
	{
		int weight=spins[i][j].occ*spins[i][N-1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][N-1].angle+f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][N-1].angle)*weight; //New Gauge
		H+=K*weight;
	}
	if (j!=N-1)
	{
		int weight=spins[i][j].occ*spins[i][j+1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][j+1].angle)*weight; //New Gauge
		H+=K*weight;
	}
	else if (j==N-1)
	{
		int weight=spins[i][j].occ*spins[i][0].occ;
		H+=J*cos(spins[i][j].angle-spins[i][0].angle-f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][0].angle)*weight; //New Gauge
		H+=K*weight;	
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
				H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; //OG Gauge
				//H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*j)*weight; //New Gauge
				H+=K*weight;
			}
			else if (i==N-1)
			{
				int weight=spins[i][j].occ*spins[0][j].occ;
				H+=J*cos(spins[i][j].angle-spins[0][j].angle+f*j)*weight; //OG Gauge
				//H+=J*cos(spins[i][j].angle-spins[0][j].angle+2.0*f*j)*weight; //New Gauge
				H+=K*weight;
			}

			if (j!=N-1)
			{
				int weight=spins[i][j].occ*spins[i][j+1].occ;
				H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; //OG Gauge
				//H+=J*cos(spins[i][j].angle-spins[i][j+1].angle)*weight;  //New Gauge
				H+=K*weight;
			}
			else if (j==N-1)
			{
				int weight=spins[i][j].occ*spins[i][0].occ;
				H+=J*cos(spins[i][j].angle-spins[i][0].angle-f*i)*weight; //OG Gauge
				//H+=J*cos(spins[i][j].angle-spins[i][0].angle)*weight; //New Gauge
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
		B+=J*cos(spins[i][j].angle-spins[n][m].angle+f*i)*weight; //OG Gauge
		//B+=J*cos(spins[i][j].angle-spins[n][m].angle)*weight; //New Gauge
		B+=K*weight;
	}
	else if (i==n and j==m-1)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle-f*i)*weight; //OG Gauge
		//B+=J*cos(spins[i][j].angle-spins[n][m].angle)*weight; //New Gauge
		B+=K*weight;
	}
	else if (i==n+1 and j==m)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle-f*j)*weight; //OG Gauge
		//B+=J*cos(spins[i][j].angle-spins[n][m].angle-2.0*f*j)*weight; //New Gauge
		B+=K*weight;
	}
	else if (i==n-1 and j==m)
	{
		int weight=spins[i][j].occ*spins[n][m].occ;
		B+=J*cos(spins[i][j].angle-spins[n][m].angle+f*j)*weight; //OG Gauge
		//B+=J*cos(spins[i][j].angle-spins[n][m].angle+2.0*f*j)*weight; //New Gauge
		B+=K*weight;
	}

	return B;
}

double lattice::H_Neighbor(int i, int j)
{
	double H=0.0;

	if (i!=N-1)
	{
		int weight=spins[i][j].occ*spins[i+1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}

	if (j!=N-1)
	{
		int weight=spins[i][j].occ*spins[i][j+1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][j+1].angle)*weight; //New Gauge
		H+=K*weight;
	}

	return H;
}

double lattice::H_Neighbor_periodic(int i, int j)
{
	double H=0.0;

	if (i!=N-1)
	{
		int weight=spins[i][j].occ*spins[i+1][j].occ;
		H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}
	else
	{
		int weight=spins[i][j].occ*spins[0][j].occ;
		H+=J*cos(spins[i][j].angle-spins[0][j].angle+f*j)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[0][j].angle+2.0*f*j)*weight; //New Gauge
		H+=K*weight;
	}

	if (j!=N-1)
	{
		int weight=spins[i][j].occ*spins[i][j+1].occ;
		H+=J*cos(spins[i][j].angle-spins[i][j+1].angle-f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][j+1].angle)*weight; //New Gauge
		H+=K*weight;
	}
	else
	{
		int weight=spins[i][j].occ*spins[i][0].occ;
		H+=J*cos(spins[i][j].angle-spins[i][0].angle-f*i)*weight; //OG Gauge
		//H+=J*cos(spins[i][j].angle-spins[i][0].angle)*weight; //New Gauge
		H+=K*weight;
	}

	return H;
}

std::vector<double> lattice::Bond_Gauge(int i, int j)
{

	double BER, BED;
	std::vector<double> bonds;

	if (i!=N-1)
	{
		int weight=spins[i][j].occ*spins[i+1][j].occ;
		BER=cos(spins[i][j].angle-spins[i+1][j].angle+f*(double) j)*weight; //OG Gauge
		//BER=cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*(double) j)*weight; //New Gauge
	}
	if (j!=N-1)
	{
		int weight=spins[i][j].occ*spins[i][j+1].occ;
		BED=cos(spins[i][j].angle-spins[i][j+1].angle-f*(double) i)*weight; //OG Gauge
		//BED=cos(spins[i][j].angle-spins[i][j+1].angle)*weight; //New Gauge
	}

	bonds.push_back(BER); bonds.push_back(BED);
	return bonds;
}

std::vector<double> lattice::Bond_Gauge_periodic(int i, int j)
{

	double BER, BED;
	std::vector<double> bonds;

	if (i!=N-1)
	{
		int weight=spins[i][j].occ*spins[i+1][j].occ;
		BER=cos(spins[i][j].angle-spins[i+1][j].angle+f*(double) j)*weight; //OG Gauge
		//BER=cos(spins[i][j].angle-spins[i+1][j].angle+2.0*f*(double) j)*weight; //New Gauge
	}
	else
	{
		int weight=spins[i][j].occ*spins[0][j].occ;
		BER=cos(spins[i][j].angle-spins[0][j].angle+f*(double) j)*weight; //OG Gauge
		//BER=cos(spins[i][j].angle-spins[0][j].angle+2.0*f*(double) j)*weight; //New Gauge
	}
	if (j!=N-1)
	{
		int weight=spins[i][j].occ*spins[i][j+1].occ;
		BED=cos(spins[i][j].angle-spins[i][j+1].angle-f*(double) i)*weight; //OG Gauge
		//BED=cos(spins[i][j].angle-spins[i][j+1].angle)*weight; //New Gauge
	}
	else
	{
		int weight=spins[i][j].occ*spins[i][0].occ;
		BED=cos(spins[i][j].angle-spins[i][0].angle-f*(double) i)*weight; //OG Gauge
		//BED=cos(spins[i][j].angle-spins[i][0].angle)*weight; //New Gauge
	}

	bonds.push_back(BER); bonds.push_back(BED);
	return bonds;
}
