class Skeleton
{
	int N;
	int N_Spins;
	lattice system;
	std::vector<int> labels;
	std::vector<std::vector<int> > matrix;

public:

	Skeleton(lattice init);
	void thin();
};