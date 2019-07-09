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
}

void HK::Find_Cluster()
{
	for (int i=0; i<N; i++)
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
						assert(labels[0] < N);
						labels[labels[0]] = labels[0];
		  				matrix[i][j] = labels[0];      // a new cluster
		  				break;
		  
					case 1:                              // part of an existing cluster
						matrix[i][j] = max(up,left);       // whichever is nonzero is labelled
						break;
		  
					case 2:                              // this site binds two clusters
						matrix[i][j] = unionize(up, left, labels);
						break;
				}
			}
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

int HK::cluster_count()
{
	int max_label=0;
	for (int i=0; i<N; i++)
	{
		if (labels[i]>max_label) {max_label++;}
	}
	return max_label;
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

