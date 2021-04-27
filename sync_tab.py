from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import scipy.integrate as integrate
import scipy.special as special


R=100.0
p=2.0
p_low=p

index=5.0/3.0

def integ(x):
	i=integrate.quad(lambda y: special.kv(index,y),x,np.inf)
	return i[0]

def F(x):
	z=x*integ(x)
	return z


def h1(nu,R,p):
	p_c=2.0
	i=integrate.quad(lambda y: F(y)*y**((p_c-3)/2),nu/R,nu)
	z1=nu**(-(p_c-1)/2)*i[0]
	return z1



def h2(nu,R,p):
	p_in=p+1
	j=integrate.quad(lambda y: F(y)*y**((p_in-3)/2),0,nu/R)
	z2=nu**(-(p_in-1)/2)*(R)**(p_in/2-1)*j[0]

	return z2

def h_tot(nu,R,p):
	return h1(nu,R,p)+h2(nu,R,p)

h_tot_arr=[]
nu_arr=np.logspace(-3,+6,int(1e2))

for n in range(0,len(nu_arr)):
	h_tot_arr.append(h_tot(nu_arr[n],R,p))
	print(n)

sync_tab= open("sync_tab.txt", 'w')

for n in range(0,len(nu_arr)):
	sync_tab.write('%f\t%f\n' %(nu_arr[n],h_tot_arr[n]))


fig, ax1 = plt.subplots(figsize=(12,9))

plt.plot(nu_arr,h_tot_arr)
plt.xscale('log')
plt.yscale('log')

ax1.set_xlabel('$\\nu/min(\\nu_m,\\nu_c)$', fontsize=28)
ax1.set_ylabel('$F_{\\nu}$', fontsize=28)

plt.show()

