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
}

void HK::Find_Cluster()
{
	int N=system.how_many();

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
	file.open("thing.p");

	file<<"set terminal png"<<endl;
	file<<"set output 'thing.png'"<<endl;
	file<<"set key off"<<endl;
	file<<"set xrange [0:53]"<<endl;
	file<<"set yrange [0:53]"<<endl;

	for (int i=0; i<system.how_many(); i++)
	{
		for (int j=0; j<system.how_many(); j++)
		{
			if (matrix[i][j]==0) {continue;}
			file<<"set label '"<<matrix[i][j]<<"' at "<<(i+1)*2.5<<","<<(j+1)*2.5<<endl;
		}
	}
	file<<"plot NaN"<<endl;
	file.close();
}