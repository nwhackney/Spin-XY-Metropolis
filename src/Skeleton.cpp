#include "../include/Skeleton.hpp"

Skeleton::Skeleton(lattice init)
{
	system=init;
	N=system.how_many();
	N_Spins=system.spin_num();	
}

void Skeleton::thin()
{
	vector<vector<int> > Bin;
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

	vector<vector<int> > Boundary;
    	Boundary.resize(N);
	for( auto &it : Boundary )
	{
		it.resize(N, 0);
	}

		vector<vector<int> > copy;
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
		vector<vector<int> > temp;
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

	ofstream Skull;
	Skull.open("Skeleton.p");

	Skull<<"set terminal png"<<endl;
	Skull<<"set output 'Skeleton.png'"<<endl;
	Skull<<"set key off"<<endl;
	Skull<<"set xrange [0:86]"<<endl;
	Skull<<"set yrange [0:86]"<<endl;
	Skull<<"set style arrow 2 nohead ls 10 "<<endl;

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N;j++)
		{
			if (Boundary[i][j]==0) {continue;}

			if (Boundary[i-1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i-1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+1<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i-1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+2<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j+2<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i+1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i+1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+1<<" as 2 lc rgb 'violet'"<<endl;
			}
			if (Boundary[i+1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+2<<" as 2 lc rgb 'violet'"<<endl;
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
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i-1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+1<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i-1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i<<","<<j+2<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+1<<","<<j+2<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i+1][j-1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i+1][j]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+1<<" as 2 lc rgb 'black'"<<endl;
			}
			if (Bin[i+1][j+1]==1)
			{
				Skull<< "set arrow from "<<i+1<<","<<j+1<<" to "<<i+2<<","<<j+2<<" as 2 lc rgb 'black'"<<endl;
			}
		}
	}

	Skull<<"plot NaN"<<endl;
	Skull.close();
}