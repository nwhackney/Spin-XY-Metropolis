#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include <limits>

#include "lattice.cpp"
#include "../include/tinytoml-master/include/toml/toml.h"
#include "Hoshen_Kopelman.cpp"

using namespace std;

void print_bonds(lattice &system, string file_name)
{
	stringstream file;
	file<<file_name<<".p";

	ofstream out;
	out.open(file.str());

	out<<"set terminal png"<<endl;
	out<<"set output '"<<file_name<<"'"<<endl;
	out<<"set key off"<<endl;
	out<<"set xrange [0:203]"<<endl;
	out<<"set yrange [0:203]"<<endl;
	out<<"set style arrow 1 head filled size screen 0.03,15 ls 2 lc 'black'"<<endl;
	out<<"set style arrow 2 nohead "<<endl;

	double d=2.5;
	double theta,
	       x,
	       y,
	       dx,
	       dy;

	int N=system.how_many();
	int plaq=1;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (system.occ(i,j)==1)
			{
				x=(i+1)*d;
				y=(j+1)*d;

				double Beni=0.0;
				double Benj=0.0;
				if (i!=N-1) {Beni=system.Bond_Energy(i,j,i+1,j);}
				if (j!=N-1) {Benj=system.Bond_Energy(i,j,i,j+1);}

				if (Beni!=0.0)
				{
					double w=-2.0*Beni;
					double R=255.0-205.0*w;
					double G=0.0;
					double B=50.0+205.0*w;

					stringstream Red;
					stringstream Green;
					stringstream Blue;

					Red<<hex<<(int) R;
					Green<<hex<<(int) G;
					Blue<<hex<<(int) B;

					out<< "set arrow from "<<x<<","<<y<<" to "<<x+d<<","<<y<<" as 2 lc rgb \"#"<<Red.str()<<Green.str()<<Blue.str()<<"\""<<endl;
				}
				if (Benj!=0)
				{
					double w=-2.0*Benj;
					double R=255.0-205.0*w*w*w;
					double G=0.0;
					double B=50.0+205.0*w*w*w;

					stringstream Red;
					stringstream Green;
					stringstream Blue;

					Red<<hex<<(int) R;
					Green<<hex<<(int) G;
					Blue<<hex<<(int) B;

					out<< "set arrow from "<<x<<","<<y<<" to "<<x<<","<<y+d<<" as 2 lc rgb \"#"<<Red.str()<<Green.str()<<Blue.str()<<"\""<<endl;
				}
			}
		}
	}

	out<<"plot NaN"<<endl;
}

void print_sys(lattice &system, string file_name)
{
	stringstream file;
	file<<file_name<<".p";

	ofstream out;
	out.open(file.str());

	stringstream png;
	png<<file_name<<".png";

	out<<"set terminal png"<<endl;
	out<<"set output '"<<png.str()<<"'"<<endl;
	out<<"set key off"<<endl;
	out<<"set xrange [0:203]"<<endl;
	out<<"set yrange [0:203]"<<endl;
	out<<"set style arrow 1 head filled size screen 0.03,15 ls 2"<<endl;

	double d=2.5;
	double theta,
	       x,
	       y,
	       dx,
	       dy;

	int N=system.how_many();
	int plaq=1;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (system.occ(i,j)==1)
			{

				x=(i+1)*d;
				y=(j+1)*d;

				theta=system.angle(i,j);
				dx=cos(theta);
				dy=sin(theta);

				double E=system.H_local(i,j);
				double w=(E+2.0)/-4.0;

				double R=255.0-255.0*w*w;
				double G=0.0;
				double B=0.0+255.0*w*w;

				stringstream Red;
				stringstream Green;
				stringstream Blue;

				Red<<hex<<(int) R;
				Green<<hex<<(int) G;
				Blue<<hex<<(int) B;

				//out<<"set arrow from "<<x<<","<<y<<" to "<<x+dx<<","<<y+dy<<" as 1 lc rgb \"#"<<Red.str()<<Green.str()<<Blue.str()<<"\""<<endl;
				out<<"set arrow from "<<x<<","<<y<<" to "<<x+dx<<","<<y+dy<<" as 1 lc 'black'"<<endl;
			}
		}
	}

	out<<"plot NaN"<<endl;
}

void print_sys_data(lattice &system, string file_name)
{
	stringstream file;
	file<<file_name<<"_data.dat";

	ofstream out;
	out.open(file.str());

	int N=system.how_many();
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (system.occ(i,j)==0)
			{
				out<<i<<" "<<j<<" "<<" -1 "<<" 0"<<endl;
			}
			else
			{
				out << i << " " << j << " " << system.angle(i,j) << " " << system.H_local(i,j) << endl;
			}
		}
	}
	out.close();
}

double Box_Muller(double mu, double sigma)
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.1416;

	thread_local double z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

void Metropolis(lattice &system, double T, ofstream &Efile, int pbc=0)
{
	double (lattice::*Hamiltonian)();
	double (lattice::*local_Hamiltonian)(int, int);

	if (pbc==0) {Hamiltonian = &lattice::H; local_Hamiltonian = &lattice::H_local; }
	else if (pbc==1){Hamiltonian = &lattice::H_periodic; local_Hamiltonian = &lattice::H_local_periodic;}

	double E = (system.*Hamiltonian)();
	Efile<<E<<endl;

	int N=system.how_many();
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N;j ++)
		{
			lattice trial = system;
			double trial_E;
			double E_local;

			if (system.occ(i,j)==0) {continue;}

			int flag = rand() % 3;
			//int flag = 0;
			if (flag==0) // rotation
			{
				double width = 0.2*exp(-0.5*T);
				double theta = Box_Muller(system.angle(i,j),width);

				trial.rotate(i,j,theta);
				trial_E=(trial.*local_Hamiltonian)(i,j);
				E_local=(system.*local_Hamiltonian)(i,j);
			}
			else if (flag==1) // translation
			{
				int n = rand() % N;
				int m = rand() % N;

				if (system.occ(n,m)==0)
				{
					trial.flip(i,j);
					trial.flip(n,m);

					double theta=system.angle(i,j);
					trial.rotate(i,j,0.0);
					trial.rotate(n,m,theta);

					trial_E=(trial.*local_Hamiltonian)(i,j);
					trial_E+=(trial.*local_Hamiltonian)(n,m);
					E_local=(system.*local_Hamiltonian)(i,j);
					E_local+=(system.*local_Hamiltonian)(n,m);
				}
			}
			else if (flag==2) // translation+rotation
			{
				int n = rand() % N;
				int m = rand() % N;

				if (system.occ(n,m)==0)
				{
					double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
					
					trial.flip(i,j);
					trial.flip(n,m);

					trial.rotate(i,j,0.0);
					trial.rotate(n,m,theta);

					trial_E=(trial.*local_Hamiltonian)(i,j);
					trial_E+=(trial.*local_Hamiltonian)(n,m);
					E_local=(system.*local_Hamiltonian)(i,j);
					E_local+=(system.*local_Hamiltonian)(n,m);
				}
			}

			double delE = trial_E - E_local;
			double alpha = ((double) rand()/(double)RAND_MAX);

			double U= exp(-1*delE/T);
			if (alpha < min(1.0,U))
			{
				system=trial;
			}
		}
	}
}

void Metropolis_no_sweep(lattice &system, double T, ofstream &Efile, int pbc=0)
{
	double (lattice::*Hamiltonian)();
	double (lattice::*local_Hamiltonian)(int, int);

	if (pbc==0) {Hamiltonian = &lattice::H; local_Hamiltonian = &lattice::H_local; }
	else if (pbc==1){Hamiltonian = &lattice::H_periodic;}

	double E = (system.*Hamiltonian)();
	Efile<<E<<endl;

	int N=system.how_many();
	int count=0;
	while (count<N*N)
	{

		int i=rand() % N;
		int j=rand() % N;

		lattice trial = system;
		double trial_E;
		double E_local;

		if (system.occ(i,j)==0) {continue;}

		int flag = rand() % 3;
		if (flag==0) // rotation
		{
			double width = 0.2*exp(-0.5*T);
			double theta = Box_Muller(system.angle(i,j),width);

			trial.rotate(i,j,theta);
			trial_E=(trial.*local_Hamiltonian)(i,j);
			E_local=(system.*local_Hamiltonian)(i,j);
		}
		else if (flag==1) // translation
		{
			int n = rand() % N;
			int m = rand() % N;

			if (system.occ(n,m)==0)
			{
				trial.flip(i,j);
				trial.flip(n,m);

				double theta=system.angle(i,j);
				trial.rotate(i,j,0.0);
				trial.rotate(n,m,theta);

				trial_E=(trial.*local_Hamiltonian)(i,j);
				trial_E+=(trial.*local_Hamiltonian)(n,m);
				E_local=(system.*local_Hamiltonian)(i,j);
				E_local+=(system.*local_Hamiltonian)(n,m);
			}
		}
		else if (flag==2) // translation+rotation
		{
			int n = rand() % N;
			int m = rand() % N;

			if (system.occ(n,m)==0)
			{
				double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
				
				trial.flip(i,j);
				trial.flip(n,m);

				trial.rotate(i,j,0.0);
				trial.rotate(n,m,theta);

				trial_E=(trial.*local_Hamiltonian)(i,j);
				trial_E+=(trial.*local_Hamiltonian)(n,m);
				E_local=(system.*local_Hamiltonian)(i,j);
				E_local+=(system.*local_Hamiltonian)(n,m);
			}
		}

		double delE = trial_E - E_local;
		double alpha = ((double) rand()/(double)RAND_MAX);

		double U= exp(-1*delE/T);
		if (alpha < min(1.0,U))
		{
			system=trial;
		}

		count++;
	}
}

void run_config()
{
	std::ifstream ifs("config.toml");
	toml::ParseResult pr = toml::parse(ifs);

	if (!pr.valid())
	{
    		cout << pr.errorReason << endl;
     	return;
	}

	const toml::Value& v = pr.value;

	const toml::Value* Np = v.find("N");
	const toml::Value* Occp = v.find("N_occ");
	const toml::Value* Jp = v.find("J");
	const toml::Value* Kp = v.find("K");
	const toml::Value* fp = v.find("f");
	const toml::Value* Tp = v.find("Time");
	const toml::Value* PBC = v.find("PBC");
	const toml::Value* outp = v.find("output");
	const toml::Value* BCG = v.find("Bond_Color_Grid");
	const toml::Value*  W = v.find("Width");

	int N= Np->as<int>();
	int occ= Occp->as<int>();
	int Time=Tp->as<int>();
	int pbc=PBC->as<int>();
	int bcg=BCG->as<int>();
	double w=W->as<double>();
	double J=Jp->as<double>();
	double K=Kp->as<double>();
	double f=fp->as<double>();
	string out_file=outp->as<string>();

	double (lattice::*Hamiltonian)();

	if (pbc==0) {Hamiltonian = &lattice::H;}
	else if (pbc==1){Hamiltonian = &lattice::H_periodic;}

	lattice crystal;
	crystal.set_const(J,K,f);
	crystal.init(N,occ);
	//crystal.rand_square_init(N);

	//cout<<"Initial Energy: "<<(crystal.*Hamiltonian)()<<endl;
	//print_sys(crystal,"init");

	stringstream Efile;
	Efile<<"Energy_"<<out_file<<".dat";
	stringstream out;
	out<<out_file;

	ofstream Edat;
	Edat.open(Efile.str());

	double duration;
	clock_t start;
		start = clock();
	
	double slope,
	       Temp;
	for (int t=0; t<Time; t++)
	{
		slope=10.0/(130000.0);
		Temp=1.0/cosh(w*slope*((double) t));
		Metropolis(crystal,Temp,Edat,pbc);
	}

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	cout<<"Time: "<<duration<<endl;
	cout<<"Final Energy: "<<(crystal.*Hamiltonian)()<<endl;

	HK clump(crystal);
	clump.Find_Cluster();
	clump.print_cluster();
	int NC=clump.cluster_count();

	stringstream info;
	info<<"info_"<<out_file<<".dat";

	ofstream inf;
	inf.open(info.str());
	inf<<N<<"x"<<N<<" lattice"<<endl;
	inf<<occ<<" spin sites"<<endl;
	inf<<Time<<" sweeps"<<endl;
	inf<<"J="<<J<<" K="<<K<<" f="<<f<<endl;
	inf<<"Final Energy: "<<(crystal.*Hamiltonian)()<<endl;
	inf<<"Time: "<<duration<<endl;
	inf<<"Number of Clusters: "<<NC<<endl<<endl;
	for (int n=1; n<=NC;n++)
	{
		int size = clump.cluster_size(n);
		if (size == 0) {continue;}
		inf<<"Cluster "<<n<<":"<<endl;
		inf<<"	Energy: "<<clump.cluster_energy(n)<<endl;
		inf<<"	spin sites: "<<size<<endl;
		vector<double> pm = clump.principle_moments(n);
		inf<<"	principle moment 1: "<<pm[0]<<endl;
		inf<<"	principle moment 2: "<<pm[1]<<endl;
		inf<<"	acylindricity: "<<pm[1]*pm[1]-pm[0]*pm[0]<<endl<<endl;
	}
	inf.close();

	Edat.close();

	print_bonds(crystal,"bonds.png");
	print_sys(crystal,out.str());
	print_sys_data(crystal,out.str());
	clump.clusters_labelled();

}

int main()
{
	srand(time(0));
	run_config();

	return 0;
}