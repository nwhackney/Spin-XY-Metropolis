#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

int main()
{
	for (double f=0.0; f<=0.01; f+=0.0005)
	{
		ofstream config;
		config.open("config.toml");

		stringstream directory;
		directory<<"frust_"<<f;

		stringstream frust;
		frust<<"f="<<f;

		config<<"N=60"<<endl;
		config<<"N_occ=1600"<<endl;
		config<<"J=-1.0"<<endl;
		config<<"K=-0.5"<<endl;
		config<<frust.str()<<endl;
		config<<"Time=100000"<<endl;
		config<<"PBC = 0"<<endl;
		config<<"output = '"<<directory.str()<<"'"<<endl;
		config<<"Bond_Color_Grid = 0"<<endl;
		config<<"Width = 0.4"<<endl;

		config.close();

		stringstream job_name;
		job_name<<"f_"<<f;

		stringstream ssbsub;
		ssbsub<<"bsub -n 1 -R rusage[mem=2048] -q long -W 24:00 -J"<< job_name.str()<<" ./a.out";
		string sbsub=ssbsub.str();
		const char *bscommand = sbsub.c_str(); 
		system(bscommand);

	}

	return 0;
}