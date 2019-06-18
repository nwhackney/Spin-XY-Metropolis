#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include <limits>

#include "lattice.cpp"
#include "../include/tinytoml-master/include/toml/toml.h"

using namespace std;

void print_sys(lattice &system, string file_name, int bond_flag=0)
{
	stringstream file;
	file<<file_name<<".p";

	ofstream out;
	out.open(file.str());

	out<<"set terminal png"<<endl;
	out<<"set output '"<<file_name<<"'"<<endl;
	out<<"set key off"<<endl;
	out<<"set xrange [0:53]"<<endl;
	out<<"set yrange [0:53]"<<endl;
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

				theta=system.angle(i,j);
				dx=cos(theta);
				dy=sin(theta);

				out<<"set arrow from "<<x<<","<<y<<" to "<<x+dx<<","<<y+dy<<" as 1"<<endl;

				if (bond_flag==1)
				{
					double Beni=0.0;
					double Benj=0.0;
					if (i!=N-1) {Beni=system.Bond_Energy(i,j,i+1,j);}
					if (j!=N-1) {Benj=system.Bond_Energy(i,j,i,j+1);}

					if (Beni!=0.0 and Benj!=0.0)
					{
						double w1=Beni/-1.0;
						double R1=255.0-205.0*w1*w1*w1;
						double G1=0.0;
						double B1=50.0+205.0*w1*w1*w1;
					
						double w2=Benj/-1.0;
						double R2=255.0-205.0*w2*w2*w2;
						double G2=0.0;
						double B2=50.0+205.0*w2*w2*w2;

						stringstream Red;
						stringstream Green;
						stringstream Blue;

						Red<<hex<<(int) (0.5*(R1+R2));
						Green<<hex<<(int) (0.5*(G1+G2));
						Blue<<hex<<(int) (0.5*(B1+B2));

						out<< "set object "<<plaq<<" rect from "<<x<<","<<y<<" to "<<x+d<<","<<y+d<<" back"<<endl;
   						out<<"set object "<<plaq<<" rect fc rgb \"#"<<Red.str()<<Green.str()<<Blue.str()<<"\""<<" fillstyle solid 1.0 noborder"<<endl;
   						plaq++;						
					}
				}
			}
		}
	}

	out<<"plot NaN"<<endl;
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

void Metropolis(lattice &system, double T, ofstream &file, int pbc=0)
{
	double (lattice::*Hamiltonian)();

	if (pbc==0) {Hamiltonian = &lattice::H;}
	else if (pbc==1){Hamiltonian = &lattice::H_periodic;}

	double E = (system.*Hamiltonian)();
	file<<E<<endl;

	int N=system.how_many();
	for (int i=0; i<N-1; i++)
	{
		for (int j=0; j<N-1;j ++)
		{
			lattice trial = system;

			if (system.occ(i,j)==0)
			{
				int flag = rand() % 2;
				if (flag==0) // translation
				{
					int n = rand() % N;
					int m = rand() % N;

					if (system.occ(n,m)==1)
					{
						trial.flip(i,j);
						trial.flip(n,m);

						double theta=system.angle(n,m);
						trial.rotate(i,j,theta);
						trial.rotate(n,m,0.0);
					}
				}
				else if (flag==1) // translation+rotation
				{
					int n = rand() % N;
					int m = rand() % N;

					if (system.occ(n,m)==1)
					{
						double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
						
						trial.flip(i,j);
						trial.flip(n,m);

						trial.rotate(i,j,theta);
						trial.rotate(n,m,0.0);

					}
				}

				double trial_E = (trial.*Hamiltonian)();
				double delE = trial_E - E;
				double alpha = ((double) rand()/(double)RAND_MAX);

				double U= exp(-1*delE/T);
				if (alpha < min(1.0,U))
				{
					E=trial_E;
					system=trial;
				}
			}
			else
			{
				int flag = rand() % 3;
				if (flag==0) // rotation
				{
					double width = 0.2*exp(-0.5*T);
					double theta = Box_Muller(system.angle(i,j),width);

					trial.rotate(i,j,theta);
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

					}
				}

				double trial_E = (trial.*Hamiltonian)();
				double delE = trial_E - E;
				double alpha = ((double) rand()/(double)RAND_MAX);

				double U= exp(-1*delE/T);
				if (alpha < min(1.0,U))
				{
					E=trial_E;
					system=trial;
				}
			}
		}
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

	int N= Np->as<int>();
	int occ= Occp->as<int>();
	int Time=Tp->as<int>();
	int pbc=PBC->as<int>();
	int bcg=BCG->as<int>();
	double J=Jp->as<double>();
	double K=Kp->as<double>();
	double f=fp->as<double>();
	string out_file=outp->as<string>();\

	double (lattice::*Hamiltonian)();

	if (pbc==0) {Hamiltonian = &lattice::H;}
	else if (pbc==1){Hamiltonian = &lattice::H_periodic;}

	for (int i=0; i<5; i++)
	{
		lattice crystal;
		crystal.set_const(J,K,f);
		crystal.init(N,occ);

		//cout<<"Initial Energy: "<<(crystal.*Hamiltonian)()<<endl;

		stringstream Efile;
		Efile<<"Energy_"<<i<<".dat";
		stringstream out;
		out<<out_file<<"_"<<i;
		stringstream info;
		info<<"info_"<<i<<".dat";

		ofstream inf;
		inf.open(info.str());

		ofstream Edat;
		Edat.open(Efile.str());

		//print_sys(crystal,"init");

		inf<<N<<"x"<<N<<" lattice"<<endl;
		inf<<occ<<" spin sites"<<endl;
		inf<<Time<<" sweeps"<<endl;
		
		double slope,
		       Temp;
		for (int t=0; t<Time; t++)
		{
			slope=10.0/((double) Time);
			Temp=1.0/cosh(0.4*slope*((double) t));
			Metropolis(crystal,Temp,Edat,pbc);
		}
		for (int t=0; t<10000; t++)
		{
			Metropolis(crystal,0.0,Edat,pbc);
		}

		cout<<"Final Energy: "<<(crystal.*Hamiltonian)()<<endl;
		inf<<"Final Energy: "<<(crystal.*Hamiltonian)()<<endl;

		Edat.close();
		inf.close();

		print_sys(crystal,out,bcg);

	}
}

int main()
{
	srand(time(0));

	run_config();



	return 0;
}