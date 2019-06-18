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
		void flip(int i, int j);
		void rotate(int i, int j, double theta);
		/* void print_sys(std::string file_name); */

		int occ(int i, int j);
		int how_many();

		double angle(int i, int j);
		double H();
		double H_periodic();
		double Bond_Energy(int i, int j, int n, int m);

};