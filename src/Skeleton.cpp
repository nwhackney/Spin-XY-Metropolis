#include <set>
#include <iterator>
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

//// Going Back to Old Boundary Finding thing for a test

	//std::vector<std::vector<int> > Boundary;
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

/// End of Old Boundary Finding Thing

    	// std::vector<std::vector<int> > Boundary;
     // Boundary.resize(N);
     // for( auto &it : Boundary )
     // {
     //      it.resize(N, 0);
     // }

     // std::vector<std::vector<int> > copy;
     // std::vector<std::vector<int> > temp_Bin;
     // copy=Bin;

     // for (int i=1; i<N-1; i++)
     // {
     //      for (int j=1; j<N-1; j++)
     //      {
     //           if (copy[i][j]==0) {continue;}
     //           if (copy[i-1][j-1]==0 or copy[i-1][j]==0 or copy[i-1][j+1]==0 or copy[i+1][j-1]==0 or copy[i+1][j]==0 or copy[i+1][j+1]==0)
     //           {
     //                temp_Bin[i][j]=0;
     //           }

     //      }
     // }

     // for (int i=1; i<N-1; i++)
     // {
     //      for (int j=1; j<N-1; j++)
     //      {
     //           Boundary[i][j]=copy[i][j]-temp_Bin[i][j];
     //      }
     // }

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

			// if (Boundary[i-1][j-1]==1)
			// {
			// 	Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j<<" as 2 lc rgb 'violet'"<<std::endl;
			// }
			if (Boundary[i-1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+1<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			// if (Boundary[i-1][j+1]==1)
			// {
			// 	Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+2<<" as 2 lc rgb 'violet'"<<std::endl;
			// }
			if (Boundary[i][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			if (Boundary[i][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j+2<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			// if (Boundary[i+1][j-1]==1)
			// {
			// 	Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j<<" as 2 lc rgb 'violet'"<<std::endl;
			// }
			if (Boundary[i+1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+1<<" as 2 lc rgb 'violet'"<<std::endl;
			}
			// if (Boundary[i+1][j+1]==1)
			// {
			// 	Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+2<<" as 2 lc rgb 'violet'"<<std::endl;
			// }
		}
	}

	for (int i=1; i<N-1; i++)
	{
		for (int j=1; j<N-1; j++)
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

		for (int i=1; i<N-1; i++)
	{
		for (int j=1; j<N-1; j++)
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

int fin(int x, std::vector<int> &labels)
{
     int y = x;
     while (labels[y] != y)
          y = labels[y];

     while (labels[x] != x)
     {
          int z = labels[x];
          labels[x] = y;
          x = z;
     }
     return y;
}

int unioni(int x, int y, std::vector<int> &labels)
{
     return labels[fin(x, labels)] = fin(y, labels);
}

void Find_Cluster(std::vector<std::vector<int> > &matrix,std::vector<std::vector<int> > &Boundary, std::vector<int> &labels, int N_Spins)
{
     int N=Boundary.size();

     for (int i=0; i<N_Spins; i++)
     {
          labels.push_back(i);
     }

     matrix.resize(N);
     for (int d=0; d<N; d++)
     {
          matrix[d].resize(N,0);
     }

     for (int n=0; n<N; n++)
     {
          for (int m=0; m<N; m++)
          {
               if (Boundary[n][m] == 0) {continue;}
               else {matrix[n][m] = 1;}
          }
     }

     for (int i=0; i<N; i++)
     {
          for (int j=0; j<N; j++)
          {
               if (matrix[i][j])
               {                        // if occupied ...
                    int up = (i==0 ? 0 : matrix[i-1][j]);    //  look up  
                    int left = (j==0 ? 0 : matrix[i][j-1]);  //  look left
                    
                    switch (!!up + !!left)
                    {
                         case 0:
                              labels[0] ++;
                              assert(labels[0] < N_Spins);
                              labels[labels[0]] = labels[0];
                              matrix[i][j] = labels[0];      // a new cluster
                              break;
            
                         case 1:                              // part of an existing cluster
                              matrix[i][j] = std::max(up,left);    // whichever is nonzero is labelled
                              break;
            
                         case 2:                              // this site binds two clusters
                              matrix[i][j] = unioni(up, left, labels);
                              break;
                    }
               }
          }
     }

     for (int i=0; i<N; i++)
     {
          for (int j=0; j<N; j++)
          {
               if (matrix[i][j]==0) {continue;}
               matrix[i][j]=fin(matrix[i][j],labels);
          }
     }
}

int Skeleton::boundary()
{
	int n_holes;

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

	std::vector<std::vector<int> > bndries;
     bndries.resize(N);
     for( auto &it : bndries )
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
               if (copy[i-1][j-1]==1 or copy[i-1][j]==1 or copy[i-1][j+1]==1 or copy[i+1][j-1]==1 or copy[i+1][j]==1 or copy[i+1][j+1]==1)
               {
                    Bin[i][j]=0;
               }

          }
     }

     std::vector<std::vector<int> > tmp;
     tmp.resize(N);
     for( auto &it : tmp )
     {
          it.resize(N, 0);
     }

     for (int i=1; i<N-1; i++)
     {
          for (int j=1; j<N-1; j++)
          {
               tmp[i][j]=copy[i][j]-Bin[i][j];
          }
     }

     for (int i=1; i<N-1; i++)
     {
          for (int j=1; j<N-1; j++)
          {
               if (tmp[i][j]==0){bndries[i][j]=1;}
               if (tmp[i][j]==1){bndries[i][j]=0;}
          }
     }

     std::vector<std::vector<int> > matrix;
     std::vector<int> labels;
     Find_Cluster(matrix,bndries,labels,system.how_many());

     std::set <int, std::greater <int> > bounds;
     for (int n=0; n<N; n++)
     {
          for (int m=0; m<N; m++)
          {
               if (matrix[n][m]==0){continue;}
               bounds.insert(matrix[n][m]);
          }
     }
     int n_bounds=bounds.size();

     n_holes=n_bounds-1;

     return n_holes;
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