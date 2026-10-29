import sys
from shapleyx import rshdmr
#fname='ishigami_data.csv'
fname1='data_0_SNH4.csv'
fname2='data_0_SNH6.csv'
fname3='data_0_SNH8.csv'
fname4='data_0_SNO4.csv'
fname5='data_0_SNO6.csv'
fname6='data_0_SNO8.csv'
fname7='data_0_SO4.csv'
fname8='data_0_SO6.csv'
fname9='data_0_SO8.csv'
fname10='data_0_BH.csv'
fname11='data_0_BA.csv'
fname=[fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,fname9,fname10,fname11]
index=int(sys.argv[1])
print(fname[index])
model=rshdmr(fname[index],polys=[10,5],method='ard_cv')
sobol,shapley,total=model.run_all()
shapley['rounded']=round(shapley['effect']*100,1)
print(fname[index])
print(shapley[shapley['effect']>0.01][['label','rounded']])
shapley[['label','effect']].to_csv(fname[index]+".dat",sep=" ",index=False,header=False,float_format="%.15e")
