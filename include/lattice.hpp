#include "site.hpp"
#include <vector>
#include <math.h>
#include <iostream>

class lattice
{
	int N;
	int N_occ;

	  double J=-1.0;
	  double K=-0.5;
	  double f=0.0;

	std::vector< std::vector<site> > spins;

	public:

		void set_const(double j, double k, double frustration);
		void init(int Number, int occupancy);
		void init_cut(int Number, int occupancy,int rise, int run);
		void restart(int Number, int occupancy, std::string infile);
		void rand_square_init(int N, int occupancy);
		void square_init(int Number, int length);
		void circle(int N, int occ, double R);
		void flip(int i, int j);
		void rotate(int i, int j, double theta);
		void add_diff(int N, int diff);

		int occ(int i, int j);
		int how_many();
		int spin_num();

		double angle(int i, int j);
		double H();
		double H_local(int i, int j);
		double strain(int i, int j);
		double H_periodic();
		double H_local_periodic(int i, int j);
		double Bond_Energy(int i, int j, int n, int m);
		double H_Neighbor(int i, int j);
		double H_Neighbor_periodic(int i, int j);
		std::vector<double> Bond_Gauge(int i, int j);
		std::vector<double> Bond_Gauge_periodic(int i, int j);

};