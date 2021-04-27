# from xspec import *
from numpy import *
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special

m_max=19
nu_obs=1

def integ(x):
	i=integrate.quad(lambda y: special.kv(5/3,y),x,np.inf)
	return i[0]

def F(x):
	z=x*integ(x)
	return z

def P(nu,gamma,B):
	nu_ch=B*gamma**2
	x=nu/nu_ch
	return B*F(x)


def f_nu(nu,B,gamma_arr,N):
	tot=0
	for n in range(0,len(gamma_arr)):
		tot+=P(nu,gamma_arr[n],B)*N[n]*gamma_arr[n]

	return tot

slope_loc=[]
f_loc=[]
for m in range(0,m_max):
	gamma_arr=[]
	N=[]
	

	N_evo=open('N_evo_%s.txt' %m,'r')
	k=0
	for line in N_evo:
		row=line.strip()
		row=row.split()
		if float(row[1])!=0:
			gamma_arr.append(float(row[0]))
			N.append(float(row[1]))

		if k==0:
			B_arr=(float(row[2]))
		k+=1
	B=B_arr
	print(B)
	nu_arr=np.logspace(5,11.5,50)

	plt.yscale('log')
	plt.xscale('log')



	F_nu=[]

	for i in range(0, len(nu_arr)):
		F_nu.append(f_nu(nu_arr[i],B,gamma_arr,N))
														
	plt.plot(nu_arr,F_nu)														
	plt.yscale('log')														
	plt.xscale('log')														
														
	count=0														
	while nu_arr[count]<nu_obs:
		count+=1
	slope_loc.append(1-(np.log10(F_nu[count+1])-np.log10(F_nu[count]))/(np.log10(nu_arr[count+1])-np.log10(nu_arr[count])))
	f_loc.append(F_nu[count])
ax1.set_xlabel('$\\nu$ (arbitrary units)', fontsize=28)
ax1.set_ylabel('$F_{\\nu}$', fontsize=28)
plt.show()






