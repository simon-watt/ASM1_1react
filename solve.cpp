#include<iostream>
#include<fstream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

#include "jacobian.cpp"
#include "feval.cpp"

int main()
{
	VectorXd y(12);
	MatrixXd J,f,delta,invJ;
	double R=0.45,KLA=108,w=0.1,tau=0.5;
	double eps=1e-6,minVal=1e-10;
	ofstream outCT("metrics.dat");
	double bifurH=-1,bifurA=-1;

	tau=0;
	while (tau<1)
	{
		tau+=0.01;
		y << 1,5,100,100,2500,250,200,2,20,1,5,7;
		f=feval(y,tau,R,KLA,w);

		int iter=0,maxIter=1000;
		f=feval(y,tau,R,KLA,w);
		do
		{
			iter++;
			J=jacobian(y,tau,R,KLA,w);
			delta=J.partialPivLu().solve(-f);
			y=y+delta;
			for (int i=0;i<12;i++)
				y(i)=max(y(i),minVal);
			f=feval(y,tau,R,KLA,w);
		} while (f.norm()>eps && iter<maxIter);

		if (iter<maxIter)
		{
			double COD,TN,TIN,TON;
			double SI=y(0),SS=y(1),XI=y(2),XS=y(3),XBH=y(4),XBA=y(5);
			double XP=y(6),SO=y(7),SNO=y(8),SNH=y(9),SND=y(10),XND=y(11);
			double ixb=0.086,ixp=0.06;
			if (R>0)
			{
				COD=SS+SI;
				TN=SNO+SNH+SND;
			}
			else
			{
				COD=SS+SI+(XS+XBA+XBH+XP+XI);
				TN=SNO+SNH+SND+(XND+ixb*(XBH+XBA)+ixp*(XP+XI));
			}
			TIN=SNH+SNO;TON=SND;
			outCT << tau << " " << COD << " " << TN << " " << SNO << " " << SNH << " " << SO << endl;
			if (bifurH<0 && XBH>1e-3)
				bifurH=tau;
			if (bifurA<0 && XBA>1e-3)
				bifurA=tau;
		}
	}
	outCT.close();

	ofstream outB("bifur.dat");outB << bifurH << " " << bifurA << endl;outB.close();
}

