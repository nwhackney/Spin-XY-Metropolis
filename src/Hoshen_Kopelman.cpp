#include "../include/Hoshen_Kopelman.hpp"
#include "Skeleton.cpp"
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

void HK::Find_Cluster_periodic()
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
				int up = (i==0 ? matrix[N-1][j] : matrix[i-1][j]);    //  look up  
				int left = (j==0 ? matrix[i][N-1] : matrix[i][j-1]);  //  look left
				
				switch (!!up + !!left)
				{
		  			case 0:
		  				labels[0]++;
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

				for (int i=0; i<N; i++)
				{
					if (matrix[i][0]!=0 and matrix[i][N-1]!=0) {matrix[i][N-1]=unionize(matrix[i][0],matrix[i][N-1],labels);}
				}

				for (int j=0; j<N; j++)
				{
					if (matrix[0][j]!=0 and matrix[N-1][j]!=0) {matrix[N-1][j]=unionize(matrix[0][j],matrix[N-1][j],labels);}
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
	file<<"set xrange [0:500]"<<endl;
	file<<"set yrange [0:500]"<<endl;
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

int HK::max_label()
{
	int max_label=0;
	for (int i=1; i<N_Spins; i++)
	{
		if (labels[i]>max_label) {max_label=labels[i];}
	}

	return max_label;
}

std::vector<double> HK::mean_distance_to_surface(int label)
{
	std::vector<int> distance;
	std::vector<xy> coord;
	std::vector<double> distr(2,0.0);

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

	double NC= (double) coord.size();

	for (int n=0; n<coord.size(); n++)
	{
		int i=coord[n].x; int j = coord[n].y;
		int nc=0.0;

		if (i== N-1 or i==0.0)
		{
			distance.push_back(0);
			continue;
		}
		else {nc += system.occ(i+1,j)+system.occ(i-1,j);}

		if (j==N-1 or j==0)
		{
			distance.push_back(0);
			continue;
		}
		else {nc += system.occ(i,j+1)+system.occ(i,j-1);}

		if (nc != 4) {distance.push_back(0); continue;}

		int D=2*N;
		for (int m=0; m<coord.size(); m++)
		{

			int s=coord[m].x; int t=coord[m].y;

			int nc2=0.0;
			if (s==N-1 or s==0) {nc2=4;}
			else {nc2 += system.occ(s+1,t)+system.occ(s-1,t);}

			if (t==N-1 or t==0) {nc2=4;}
			else {nc2 += system.occ(s,t+1)+system.occ(s,t-1);}

			int temp=abs(i-s)+abs(j-t);
			if (temp<D and nc2!=4)
			{
				D=temp;
			}

		}
		
		if (D==N) {std::cout<<"DA FUCK?"<<std::endl;}
		distance.push_back(D);
	}

	int sum=0;
	for (int u=0; u<distance.size(); u++)
	{
		sum+=distance[u];
	}
	
	double mean = (double) sum / NC;

	double sstd;
	for (int e=0; e<distance.size(); e++)
	{
		sstd+=(((double) distance[e])-mean)*(((double) distance[e])-mean);
	}
	
	double std=sqrt(sstd/(NC-1.0));
	distr[0]=mean; distr[1]=std;
	
	return distr;
}

std::vector<double> HK::mean_distance_to_surface_periodic(int label)
{
	std::vector<int> distance;
	std::vector<xy> coord;
	std::vector<double> distr(2,0.0);

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

	double NC= (double) coord.size();

	for (int n=0; n<coord.size(); n++)
	{
		int i=coord[n].x; int j = coord[n].y;
		int nc=0.0;

		if (i==0.0) {nc += system.occ(i+1,j)+system.occ(N-1,j);}
		else if (i==N-1.0) {nc += system.occ(0,j)+system.occ(i-1,j);}
		else {nc += system.occ(i+1,j)+system.occ(i-1,j);}

		if (j==0.0) {nc += system.occ(i,j+1)+system.occ(i,N-1);}
		else if (j==N-1.0) {nc += system.occ(i,0)+system.occ(i,j-1);}
		else {nc += system.occ(i,j+1)+system.occ(i,j-1);}

		if (nc != 4) {distance.push_back(0); continue;}

		int D=2*N;
		for (int m=0; m<coord.size(); m++)
		{

			int s=coord[m].x; int t=coord[m].y;

			int nc2=0.0;

			if (s==0) {nc2 += system.occ(s+1,t)+system.occ(N-1,t);}
			else if (s==N-1) {nc2 += system.occ(0,t)+system.occ(s-1,t);}
			else {nc2 += system.occ(s+1,t)+system.occ(s-1,t);}

			if (t==0) {nc2 += system.occ(s,t+1)+system.occ(s,N-1);}
			else if (t==N-1) {nc2 += system.occ(s,0)+system.occ(s,t-1);}
			else {nc2 += system.occ(s,t+1)+system.occ(s,t-1);}

			int temp=abs(i-s)+abs(j-t);
			if (temp<D and nc2!=4)
			{
				D=temp;
			}

		}
		
		distance.push_back(D);
	}

	int sum=0;
	for (int u=0; u<distance.size(); u++)
	{
		sum+=distance[u];
	}
	
	double mean = (double) sum / NC;

	double sstd;
	for (int e=0; e<distance.size(); e++)
	{
		sstd+=(((double) distance[e])-mean)*(((double) distance[e])-mean);
	}
	
	double std=sqrt(sstd/(NC-1.0));
	distr[0]=mean; distr[1]=std;
	
	return distr;
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

double HK::cluster_energy_periodic(int label)
{
	double Energy;

	for (int i=0; i<system.how_many(); i++)
	{
		for (int j=0; j< system.how_many(); j++)
		{
			if (matrix[i][j]==label)
			{
				Energy+=system.H_Neighbor_periodic(i,j);
			}
		}
	}
	return Energy;
}

double HK::cluster_skeletonize(int label)
{
	lattice agg;
	agg.init(system.how_many(),0);
	for (int i=0; i<system.how_many(); i++)
	{
		for (int j=0; j<system.how_many(); j++)
		{
			if (matrix[i][j]==label)
			{
				agg.flip(i,j);
				agg.rotate(i,j,system.angle(i,j));
			}
		}
	}
	stringstream name;
	name<<"Skeleton_"<<label;

	stringstream trim;
	trim<<"Back_Bone_"<<label;

	Skeleton id_clump(agg);
	id_clump.thin(name.str());
	//id_clump.back_bone(trim.str());

	return id_clump.medial_distance();
}

double HK::circularity(int label)
{
	int N=system.how_many();
	double perimeter=0.0;
	double n=0.0;

	for (int i=0; i<N; i++)
	{
		for (int j=0; j< N; j++)
		{
			if (matrix[i][j]==label)
			{
				int neigh=0;
				if (i!=0) {neigh+=system.occ(i-1,j);}
				if (i!=N-1) {neigh+=system.occ(i+1,j);}
				if (j!=0) {neigh+=system.occ(i,j-1);}
				if (j!=N-1) {neigh+=system.occ(i,j+1);}

				if (neigh==4) {n+=1.0;}
				else {perimeter+=1.0;}
			}
		}
	}
	return perimeter/n;
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
