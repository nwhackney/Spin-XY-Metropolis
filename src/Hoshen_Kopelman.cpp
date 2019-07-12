#include "../include/Hoshen_Kopelman.hpp"
#include <fstream>

using namespace std;

int find(int x, std::vector<int> &labels)
{
	int y = x;
	while (labels[y] != y)
		y = labels[y];

	while (labels[x] != x)
	{
		int z = labels[x];
		labels[x] = y;
		x = z;
	}
	return y;
}

int unionize(int x, int y, std::vector<int> &labels)
{
	return labels[find(x, labels)] = find(y, labels);
}

HK::HK(lattice init)
{
	system=init;
	N=system.how_many();
	N_Spins=system.spin_num();	
}

void HK::Find_Cluster()
{
	for (int i=0; i<N_Spins; i++)
	{
		labels.push_back(i);
	}

	matrix.resize(N);
	for (int d=0; d<N; d++)
	{
		matrix[d].resize(N,0);
	}

	for (int n=0; n<N; n++)
	{
		for (int m=0; m<N; m++)
		{
			if (system.occ(n,m) == 0) {continue;}
			else {matrix[n][m] = 1;}
		}
	}

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (matrix[i][j])
			{                        // if occupied ...
				int up = (i==0 ? 0 : matrix[i-1][j]);    //  look up  
				int left = (j==0 ? 0 : matrix[i][j-1]);  //  look left
				
				switch (!!up + !!left)
				{
		  			case 0:
		  				labels[0] ++;
						assert(labels[0] < N_Spins);
						labels[labels[0]] = labels[0];
		  				matrix[i][j] = labels[0];      // a new cluster
		  				break;
		  
					case 1:                              // part of an existing cluster
						matrix[i][j] = max(up,left);    // whichever is nonzero is labelled
						break;
		  
					case 2:                              // this site binds two clusters
						matrix[i][j] = unionize(up, left, labels);
						break;
				}
			}
		}
	}

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (matrix[i][j]==0) {continue;}
			matrix[i][j]=find(matrix[i][j],labels);
		}
	}
}

void HK::print_cluster()
{
	ofstream file;
	file.open("cluster.p");

	file<<"set terminal png"<<endl;
	file<<"set output 'cluster.png'"<<endl;
	file<<"set key off"<<endl;
	file<<"set xrange [0:203]"<<endl;
	file<<"set yrange [0:203]"<<endl;
	file<<"set style arrow 1 head filled size screen 0.03,15 ls 2"<<endl;

	for (int i=0; i<system.how_many(); i++)
	{
		for (int j=0; j<system.how_many(); j++)
		{
			if (matrix[i][j]==0) {continue;}
			
			double x=(i+1)*2.5;
			double y=(j+1)*2.5;

			double theta=system.angle(i,j);
			double dx=cos(theta);
			double dy=sin(theta);

			file<<"set arrow from "<<x<<","<<y<<" to "<<x+dx<<","<<y+dy<<" as 1 lc "<<matrix[i][j]<<endl;
		}
	}
	file<<"plot NaN"<<endl;
	file.close();
}

void HK::clusters_labelled()
{
	ofstream file;
	file.open("labelled.p");

	file<<"set terminal png"<<endl;
	file<<"set output 'labelled.png'"<<endl;
	file<<"set key off"<<endl;
	file<<"set xrange [0:255]"<<endl;
	file<<"set yrange [0:255]"<<endl;
	file<<"set style arrow 1 head filled size screen 0.03,15 ls 2"<<endl;

	for (int i=0; i<system.how_many(); i++)
	{
		for (int j=0; j<system.how_many(); j++)
		{
			if (matrix[i][j]==0) {continue;}
			
			double x=(i+1)*3;
			double y=(j+1)*3;

			double theta=system.angle(i,j);
			double dx=cos(theta);
			double dy=sin(theta);

			file<<"set label '"<<matrix[i][j]<<"' at "<<x<<","<<y<<endl;
		}
	}
	file<<"plot NaN"<<endl;
	file.close();
}

int HK::cluster_size(int label)
{
	int count=0;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (matrix[i][j]==label) {count++;}
		}
	}
	return count;
}

int HK::cluster_count()
{
	int max_label=0;
	for (int i=1; i<N_Spins; i++)
	{
		if (cluster_size(i)!=0) {max_label++;}
	}
	return max_label;
}

double HK::cluster_energy(int label)
{
	double Energy;

	for (int i=0; i<system.how_many(); i++)
	{
		for (int j=0; j< system.how_many(); j++)
		{
			if (matrix[i][j]==label)
			{
				Energy+=system.H_Neighbor(i,j);
			}
		}
	}
	return Energy;
}

std::vector<double> HK::principle_moments(int label)
{
	std::vector<double> moments;
	std::vector<xy> coord;

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (matrix[i][j]==label)
			{
				xy temp;
				temp.x=(double) i; temp.y=(double) j;
				coord.push_back(temp);
			}
		}
	}

	double NL=cluster_size(label);
	double C=1.0/(2.0*((double) NL)*((double) NL));

	double Sxx=0.0,
		  Sxy=0.0,
		  Syx=0.0,
		  Syy=0.0;

	for (int n=0; n<coord.size(); n++)
	{
		for (int m=0; m<coord.size(); m++)
		{
			Sxx+=C*(coord[n].x-coord[m].x)*(coord[n].x-coord[m].x);
			Sxy+=C*(coord[n].x-coord[m].x)*(coord[n].y-coord[m].y);
			Syx+=C*(coord[n].y-coord[m].y)*(coord[n].x-coord[m].x);
			Syy+=C*(coord[n].y-coord[m].y)*(coord[n].y-coord[m].y);
		}
	}

	double Trace = Sxx + Syy;
	double Det = Sxx*Syy-Sxy*Syx;

	double lambda1=0.5*Trace+sqrt(0.25*Trace*Trace - Det);
	double lambda2=0.5*Trace-sqrt(0.25*Trace*Trace - Det);

	if (lambda1<=lambda2)
	{
		moments.push_back(lambda1);
		moments.push_back(lambda2);
	}
	else if(lambda1>lambda2)
	{
		moments.push_back(lambda2);
		moments.push_back(lambda1);
	}

	return moments;
}
