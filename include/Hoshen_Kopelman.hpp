#include "xy.hpp"

class HK
{
	int N;
	lattice system;
	std::vector<int> labels;
	std::vector<std::vector<int> > matrix;

public:

	HK(lattice init);
	void Find_Cluster();
	void print_cluster();
	void clusters_labelled();

	int cluster_count();
	int cluster_size(int label);

	double cluster_energy(int label);

	std::vector<double> principle_moments(int label);
};