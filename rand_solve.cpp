#include<iostream>
#include<fstream>
#include<eigen3/Eigen/Dense>
#include<ctime>

using namespace Eigen;
using namespace std;
struct paramStruct
{
	double KS,KOA,KOH,KNO,KNH,KX,KA,KH;
	double muh,mua,YH,YA,etag,etah,fp;
	double bh,ba,ixb,ixp,SO_max;
	double SI_in,SS_in,XI_in,XS_in,XBH_in,XBA_in;
	double XP_in,SO_in,SNO_in,SNH_in,SND_in,XND_in;
} param;

#include "jacobian_param.cpp"
#include "feval_param.cpp"

void output(ofstream &out,double x)
{
	out << param.KS << "," << param.KOA << "," << param.KOH << "," << param.KNO << ",";
	out << param.KNH << "," << param.KX << "," << param.KA << "," << param.KH << ",";
	out << param.muh << "," << param.mua << "," << param.YH << "," << param.YA << ",";
	out << param.etag << "," << param.etah << "," << param.fp << "," << param.bh << ",";
	out << param.ba << "," << param.ixb << "," << param.ixp << "," << param.SO_max << ",";
	out << x << endl;
}

int main(int argc, char ** argv)
{
	VectorXd y(12);
	MatrixXd J,f,delta,invJ;
	double R=0.45,KLA=108,w=0.1,tau=0.5;
	double eps=1e-6,minVal=1e-10;
	ofstream outCT("metrics.dat");
	double safety=0.9;
	int count=0;
	ofstream outB("bifur.dat");
	ofstream outSNO4,outSNO6,outSNO8;
	ofstream outSNH4,outSNH6,outSNH8;
	ofstream outSO4,outSO6,outSO8;
	ofstream outBH,outBA;
	char fname[300];
	int id;

	srand(time(0));

	if (argc>1)
		id=atoi(argv[1]);
	else
		id=0;

	sprintf(fname,"data_%i_SNO4.csv",id);outSNO4.open(fname);
	sprintf(fname,"data_%i_SNO6.csv",id);outSNO6.open(fname);
	sprintf(fname,"data_%i_SNO8.csv",id);outSNO8.open(fname);
	sprintf(fname,"data_%i_SNH4.csv",id);outSNH4.open(fname);
	sprintf(fname,"data_%i_SNH6.csv",id);outSNH6.open(fname);
	sprintf(fname,"data_%i_SNH8.csv",id);outSNH8.open(fname);
	sprintf(fname,"data_%i_SO4.csv",id);outSO4.open(fname);
	sprintf(fname,"data_%i_SO6.csv",id);outSO6.open(fname);
	sprintf(fname,"data_%i_SO8.csv",id);outSO8.open(fname);
	sprintf(fname,"data_%i_BH.csv",id);outBH.open(fname);
	sprintf(fname,"data_%i_BA.csv",id);outBA.open(fname);

	outSNO4 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSNO6 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSNO8 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSNH4 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSNH6 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSNH8 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSO4 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSO6 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outSO8 << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outBH << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;
	outBA << "X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,Y" << endl;

	while (clock()<60*CLOCKS_PER_SEC*safety && count<1000)
	{
		cout << count++ << " " << 1.0*clock()/CLOCKS_PER_SEC << endl;

		double bifurH=-1,bifurA=-1;
		double SNH4=-1,SNH6=-1,SNH8=-1;
		double SNO4=-1,SNO6=-1,SNO8=-1;
		double SO4=-1,SO6=-1,SO8=-1;

		double fac;

		param.KS=20;param.KOA=0.4;param.KOH=0.2;param.KNO=0.5;param.KNH=1;param.KX=0.03;param.KA=0.08;param.KH=3;
		param.muh=6;param.mua=0.8;param.YH=0.67;param.YA=0.24;param.etag=0.8;param.etah=0.4;param.fp=0.08;
		param.bh=0.22;param.ba=0.05;param.ixb=0.086;param.ixp=0.06;param.SO_max=10;

		param.SI_in=2;param.SS_in=200;param.XI_in=3;param.XS_in=100;param.XBH_in=0;param.XBA_in=0;
		param.XP_in=0;param.SO_in=2;param.SNO_in=1;param.SNH_in=15;param.SND_in=15;param.XND_in=9;

		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KS*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KOA*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KOH*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KNO*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KNH*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KX*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KA*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.KH*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.muh*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.mua*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.YH*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.YA*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.etag*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.etah*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.fp*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.bh*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.ba*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.ixb*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.ixp*=fac;
		fac=1.0*rand()/RAND_MAX*0.6+0.7;param.SO_max*=fac;

		tau=0;
		while (tau<1)
		{
			tau+=0.001;
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
				if (SNH4<0 && tau>=4.0/24.0)
				{
					SNH4=SNH;SNO4=SNO;SO4=SO;
				}
				if (SNH6<0 && tau>=6.0/24.0)
				{
					SNH6=SNH;SNO6=SNO,SO6=SO;
				}
				if (SNH8<0 & tau>=8.0/24.0)
				{
					SNH8=SNH;SNO8=SNO;SO8=SO;
				}
				if (bifurH<0 && XBH>1e-3)
					bifurH=tau;
				if (bifurA<0 && XBA>1e-3)
					bifurA=tau;
			}
		}
		outB << bifurH << " " << bifurA << endl;
		output(outSNH4,SNH4);output(outSNH6,SNH6);output(outSNH8,SNH8);
		output(outSNO4,SNO4);output(outSNO6,SNO6);output(outSNO8,SNO8);
		output(outSO4,SO4);output(outSO6,SO6);output(outSO8,SO8);
		output(outBH,bifurH);output(outBA,bifurA);
	}
	outCT.close();
	outB.close();
	outSNH4.close();outSNH6.close();outSNH8.close();
	outSNO4.close();outSNO6.close();outSNO8.close();
	outSO4.close();outSO6.close();outSO8.close();
	outBH.close();outBA.close();

}

