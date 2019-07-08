
class HK
{
	lattice system;
	std::vector<int> labels;
	std::vector<std::vector<int> > matrix;

public:

	HK(lattice init);
	void Find_Cluster();
	void print_cluster();
};