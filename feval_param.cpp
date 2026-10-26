MatrixXd feval(VectorXd y,double tau,double R,double KLA,double w)
{
	MatrixXd dy(12,1);

	double SI=y(0),SS=y(1),XI=y(2),XS=y(3),XBH=y(4),XBA=y(5);
	double XP=y(6),SO=y(7),SNO=y(8),SNH=y(9),SND=y(10),XND=y(11);
	/*
	double KS=20,KOA=0.4,KOH=0.2,KNO=0.5,KNH=1,KX=0.03,KA=0.08,KH=3;
	double muh=6,mua=0.8,YH=0.67,YA=0.24,etag=0.8,etah=0.4,fp=0.08;
	double bh=0.22,ba=0.05,ixb=0.086,ixp=0.06,SO_max=10;

	double SI_in=2,SS_in=200,XI_in=3,XS_in=100,XBH_in=0,XBA_in=0;
	double XP_in=0,SO_in=2,SNO_in=1,SNH_in=15,SND_in=15,XND_in=9;
	*/
	double KS=param.KS,KOA=param.KOA,KOH=param.KOH,KNO=param.KNO,KNH=param.KNH,KX=param.KX,KA=param.KA,KH=param.KH;
        double muh=param.muh,mua=param.mua,YH=param.YH,YA=param.YA,etag=param.etag,etah=param.etah,fp=param.fp;
        double bh=param.bh,ba=param.ba,ixb=param.ixb,ixp=param.ixp,SO_max=param.SO_max;

        double SI_in=param.SI_in,SS_in=param.SS_in,XI_in=param.XI_in,XS_in=param.XS_in,XBH_in=param.XBH_in;
	double XBA_in=param.XBA_in,XP_in=param.XP_in,SO_in=param.SO_in,SNO_in=param.SNO_in,SNH_in=param.SNH_in;
	double SND_in=param.SND_in,XND_in=param.XND_in;

	double b1=(1+R)/(R+w);

	double M2=SS/(KS+SS);
	double M8a=SO/(KOA+SO);
	double M8h=SO/(KOH+SO);
	double M9=SNO/(KNO+SNO);
	double M10=SNH/(KNH+SNH);
	double I8=KOH/(KOH+SO);
	double Ksat=XS/(KX*XBH+XS);
	double R_si=0;
	double R_ss=(-muh/YH*M2*(M8h+etag*I8*M9)+KH*Ksat*(M8h+etah*I8*M9))*XBH;
	double R_xi=0;
	double R_xs=(1-fp)*(bh*XBH+ba*XBA)-KH*Ksat*(M8h+etah*I8*M9)*XBH;
	double R_xbh=(muh*M2*(M8h+etag*I8*M9)-bh)*XBH;
	double R_xba=(mua*M10*M8a-ba)*XBA;
	double R_xp=fp*(bh*XBH+ba*XBA);
	double R_so=-muh*(1-YH)/YH*M2*M8h*XBH-mua*(4.57-YA)/YA*M10*M8a*XBA;
	double R_sno=-muh*etag*(1-YH)/(2.86*YH)*M2*I8*M9*XBH+mua/YA*M10*M8a*XBA;
	double R_snh=(-ixb*muh*M2*(M8h+etag*I8*M9)+KA*SND)*XBH-mua*(ixb+1/YA)*M10*M8a*XBA;
	double R_snd=(-KA*SND+KH*Ksat*XND/XS*(M8h+etah*I8*M9))*XBH;
	double R_xnd=(ixb-fp*ixp)*(bh*XBH+ba*XBA)-KH*Ksat*(M8h+etah*I8*M9)*XBH*XND/XS;
	double dSI=1/tau*(SI_in-SI)+R_si;
	double dSS=1/tau*(SS_in-SS)+R_ss;
	double dXI=1/tau*(XI_in-XI)+R/tau*(b1-1)*XI+R_xi;
	double dXS=1/tau*(XS_in-XS)+R/tau*(b1-1)*XS+R_xs;
	double dXBH=1/tau*(XBH_in-XBH)+R/tau*(b1-1)*XBH+R_xbh;
	double dXBA=1/tau*(XBA_in-XBA)+R/tau*(b1-1)*XBA+R_xba;
	double dXP=1/tau*(XP_in-XP)+R/tau*(b1-1)*XP+R_xp;
	double dSO=1/tau*(SO_in-SO)+KLA*(SO_max-SO)+R_so;
	double dSNO=1/tau*(SNO_in-SNO)+R_sno;
	double dSNH=1/tau*(SNH_in-SNH)+R_snh;
	double dSND=1/tau*(SND_in-SND)+R_snd;
	double dXND=1/tau*(XND_in-XND)+R/tau*(b1-1)*XND+R_xnd;

	dy << dSI,dSS,dXI,dXS,dXBH,dXBA,dXP,dSO,dSNO,dSNH,dSND,dXND;

	return dy;
}



