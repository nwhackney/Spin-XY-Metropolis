#include "../include/Skeleton.hpp"

Skeleton::Skeleton(lattice init)
{
	system=init;
	N=system.how_many();

	std::vector<std::vector<int> > Bin;
	std::vector<std::vector<int> > Boundary;
}

void Skeleton::thin(std::string file)
{
    	Bin.resize(N);
	for( auto &it : Bin )
	{
		it.resize(N, 0);
	}

     for (int i=0; i<N; i++)
     {
     	for (int j=0; j<N; j++)
     	{
     		Bin[i][j]=system.occ(i,j);
     	}
     }

    	Boundary.resize(N);
	for( auto &it : Boundary )
	{
		it.resize(N, 0);
	}

		std::vector<std::vector<int> > copy;
		copy=Bin;

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i-1][j-1]==0 and copy[i][j-1]==0 and copy[i+1][j-1]==0 and copy[i-1][j+1]==1 and copy[i][j+1]==1 and copy[i+1][j+1]==1) // a1
				{
					Boundary[i][j]=1;
				}
			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i][j-1]==0 and copy[i+1][j-1]==0 and copy[i+1][j]==0 and copy[i-1][j]==1 and copy[i][j+1]==1) // b1
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i+1][j-1]==0 and copy[i+1][j]==0 and copy[i+1][j+1]==0 and copy[i-1][j-1]==1 and copy[i-1][j]==1 and copy[i-1][j+1]==1) // a2
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i+1][j]==0 and copy[i][j+1]==0 and copy[i+1][j+1]==0 and copy[i-1][j]==1 and copy[i][j-1]==1) // b2
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i-1][j+1]==0 and copy[i][j+1]==0 and copy[i+1][j+1]==0 and copy[i-1][j-1]==1 and copy[i][j-1]==1 and copy[i+1][j-1]==1) // a3
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i-1][j]==0 and copy[i-1][j+1]==0 and copy[i][j+1]==0 and copy[i][j-1]==1 and copy[i+1][j]==1) // b3
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i-1][j-1]==0 and copy[i-1][j]==0 and copy[i-1][j+1]==0 and copy[i+1][j-1]==1 and copy[i+1][j]==1 and copy[i+1][j+1]==1) // a4
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==0) {continue;}
				if (copy[i-1][j]==0 and copy[i-1][j-1]==0 and copy[i][j-1]==0 and copy[i+1][j]==1 and copy[i][j+1]==1) // b4
				{
					Boundary[i][j]=1;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (copy[i][j]==1)
				{
					copy[i][j]=copy[i][j]-Boundary[i][j];
				}
			}
		}

	int count=0;
	do
	{
		count = 0;
		std::vector<std::vector<int> > temp;
		temp.resize(N);
		for( auto &it : temp )
		{
			it.resize(N, 0);
		}
			
		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i-1][j-1]==0 and Bin[i][j-1]==0 and Bin[i+1][j-1]==0 and Bin[i-1][j+1]==1 and Bin[i][j+1]==1 and Bin[i+1][j+1]==1) // a1
				{
					temp[i][j]=1;
					count++;
				}
			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i][j-1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i-1][j]==1 and Bin[i][j+1]==1) // b1
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0 and Bin[i-1][j-1]==1 and Bin[i-1][j]==1 and Bin[i-1][j+1]==1) // a2
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i+1][j]==0 and Bin[i][j+1]==0 and Bin[i+1][j+1]==0 and Bin[i-1][j]==1 and Bin[i][j-1]==1) // b2
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i-1][j+1]==0 and Bin[i][j+1]==0 and Bin[i+1][j+1]==0 and Bin[i-1][j-1]==1 and Bin[i][j-1]==1 and Bin[i+1][j-1]==1) // a3
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j+1]==0 and Bin[i][j-1]==1 and Bin[i+1][j]==1) // b3
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i+1][j-1]==1 and Bin[i+1][j]==1 and Bin[i+1][j+1]==1) // a4
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}

		///////

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==0) {continue;}
				if (Bin[i-1][j]==0 and Bin[i-1][j-1]==0 and Bin[i][j-1]==0 and Bin[i+1][j]==1 and Bin[i][j+1]==1) // b4
				{
					temp[i][j]=1;
					count++;
				}

			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-temp[i][j];
				}
			}
		}
	}while (count!=0);

	std::stringstream name;
	name<<file<<".p";

	std::ofstream Skull;
	Skull.open(name.str());

	name<<"ng";

	Skull<<"set terminal png"<<std::endl;
	Skull<<"set output '"<<name.str()<<"'"<<std::endl;
	Skull<<"set key off"<<std::endl;
	Skull<<"set xrange [0:86]"<<std::endl;
	Skull<<"set yrange [0:86]"<<std::endl;
	Skull<<"set style arrow 2 nohead ls 10 "<<std::endl;

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N;j++)
		{
			if (Boundary[i][j]==0) {continue;}

			if (Boundary[i-1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i-1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+1<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i-1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+2<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j+2<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i+1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i+1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+1<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i+1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+2<<" as 2 lc rgb 'violet'"<<std::endl;
			}
		}
	}

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (Bin[i][j]==0) {continue;}
			//Skull<<"set label '"<<Bin[i][j]<<"' at "<<i+1<<","<<j+1<<endl;

			if (Bin[i-1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i-1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+1<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i-1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+2<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j+2<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i+1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i+1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+1<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i+1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+2<<" as 2 lc rgb 'black'"<<std::endl;
			}
		}
	}

	Skull<<"plot NaN"<<std::endl;
	Skull.close();
}

void Skeleton::back_bone(std::string file)
{
	std::vector<std::vector<int> > dup;
	dup.resize(N);
	for( auto &it : dup )
	{
		it.resize(N, 0);
	}

	for (int i=1; i<N-1; i++)
	{
		for (int j=1; j<N-1; j++)
		{
			int neigh=Bin[i-1][j-1]+Bin[i-1][j]+Bin[i-1][j+1]+Bin[i][j-1]+Bin[i][j+1]+Bin[i+1][j-1]+Bin[i+1][j]+Bin[i+1][j+1];
			if (neigh>=3)
			{
				dup[i][j]=2;
			}
		}
	}

	std::vector<std::vector<int> > cleaned;
	cleaned.resize(N);
	for( auto &it : cleaned )
	{
		it.resize(N, 0);
	}

	int count;
	do
	{
		count=0.0;
		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]!=1) {continue;}
				if (Bin[i-1][j-1]==1 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j-1]==0 and Bin[i][j+1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==1 and Bin[i-1][j+1]==0 and Bin[i][j-1]==0 and Bin[i][j+1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==1 and Bin[i][j-1]==0 and Bin[i][j+1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j-1]==1 and Bin[i][j+1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j-1]==0 and Bin[i][j+1]==1 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j-1]==0 and Bin[i][j+1]==0 and Bin[i+1][j-1]==1 and Bin[i+1][j]==0 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j-1]==0 and Bin[i][j+1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==1 and Bin[i+1][j+1]==0)
				{
					cleaned[i][j]=1;
					count++;
				}
				else if (Bin[i-1][j-1]==0 and Bin[i-1][j]==0 and Bin[i-1][j+1]==0 and Bin[i][j-1]==0 and Bin[i][j+1]==0 and Bin[i+1][j-1]==0 and Bin[i+1][j]==0 and Bin[i+1][j+1]==1)
				{
					cleaned[i][j]=1;
					count++;
				}	
			}
		}

		for (int i=1; i<N-1; i++)
		{
			for (int j=1; j<N-1; j++)
			{
				if (Bin[i][j]==1)
				{
					Bin[i][j]=Bin[i][j]-cleaned[i][j];
				}
				if (Bin[i][j]==2)
				{
					Bin[i][j]=1;
				}
			}
		}

	} while (count!=0);

	std::stringstream name;
	name<<file<<"_cleaned"<<".p";

	std::ofstream Skull;
	Skull.open(name.str());

	name<<"ng";

	Skull<<"set terminal png"<<std::endl;
	Skull<<"set output '"<<name.str()<<"'"<<std::endl;
	Skull<<"set key off"<<std::endl;
	Skull<<"set xrange [0:86]"<<std::endl;
	Skull<<"set yrange [0:86]"<<std::endl;
	Skull<<"set style arrow 2 nohead ls 10 "<<std::endl;

		for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (Bin[i][j]==0) {continue;}
			//Skull<<"set label '"<<Bin[i][j]<<"' at "<<i+1<<","<<j+1<<endl;

			if (Bin[i-1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i-1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+1<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i-1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+2<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j+2<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i+1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i+1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+1<<" as 2 lc rgb 'black'"<<std::endl;
			}
			if (Bin[i+1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+2<<" as 2 lc rgb 'black'"<<std::endl;
			}
		}
	}

	Skull<<"plot NaN"<<std::endl;
	Skull.close();

}

double Skeleton::medial_distance()
{

	std::vector<double> distance;

	double min;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			if (Bin[i][j]!=1) {continue;}
			min=double(N);
			for (int n=0; n<N; n++)
			{
				for (int m=0; m<N; m++)
				{
					if (Boundary[n][m]==0) {continue;}
					double D=abs(i-n)+abs(j-m);
					if (D<min and D!=0) {min=D;}
				}
			}
			distance.push_back(min);
		}
	}

	double sum=0.0;
	for (int u=0; u<distance.size(); u++)
	{
		sum+=distance[u];
	}

	return sum/(double(distance.size()));
}