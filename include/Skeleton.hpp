class Skeleton
{
	int N;
	lattice system;
	std::vector<std::vector<int> > Bin;
	std::vector<std::vector<int> > Boundary;

public:

	Skeleton(lattice init);
	void thin(std::string file);
	void back_bone(std::string file);
	int boundary();
	double medial_distance();
};