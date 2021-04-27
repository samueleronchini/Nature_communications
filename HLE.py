from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

nu_p=10
theta_max=5*3.14/180


gamma=100.0
beta=1-1/(2*gamma**2)


alpha=-0.67
beta_s=-2.5 

#nu_0_p=1/(2*gamma)*nu_p/(2+alpha) # comment if is a SBPL is used
nu_0_p=1/(2*gamma)*(-(alpha+2)/(beta_s+2))**(-1/(alpha-beta_s))*nu_p  # comment if is a Band is used

def band(nu,nu_0):
	x=nu/nu_0
	if x<(alpha-beta_s):
		y=(x)**(alpha+1)*exp(-x)
	else:
		y=(alpha-beta_s)**(alpha-beta_s)*exp(-alpha+beta_s)*x**(beta_s+1)

	return y

def broken(nu,nu_0):
	x=nu/nu_0
	y=((x)**(-alpha-1)+x**(-beta_s-1))**(-1)
	return y


def D(tau):
	return 2*gamma/tau


def F_ratio(tau,nu):
	tau_in=tau
	nu_in=nu/D(tau_in)
	return broken(nu_in,nu_0_p)  # comment if is a SBPL is used
	# return band(nu_in,nu_0_p) # comment if is a Band is used



tau_max=(1-beta*np.cos(theta_max))/(1-beta)
tau_arr=np.logspace(0,np.log10(tau_max))
gamma_arr=tau_arr*0.0

nu_1=0.5
nu_2=10.0

for n in range(0,len(tau_arr)):
	gamma_arr[n]=1-np.log(F_ratio(tau_arr[n],nu_2)/F_ratio(tau_arr[n],nu_1))/np.log(10/0.5)


def int_f(N0,gamma,tau):

	return N0/(2-gamma)*(20**(2-gamma)-1)*(D(tau))**2*(1-tau*(1-beta))/beta


int_f_arr=tau_arr*0.0

for n in range(0,len(tau_arr)):
	int_f_arr[n]=(int_f(F_ratio(tau_arr[n],nu_1),gamma_arr[n],tau_arr[n])/int_f(F_ratio(tau_arr[0],nu_1),gamma_arr[0],tau_arr[0]))**(-1.)

fig, ax1 = plt.subplots(figsize=(12,9))

plt.plot(int_f_arr,gamma_arr,color='black',ls='-.', linewidth=3.0)

plt.xlim((0.7,10000))
plt.ylim((0.3,3.2))
plt.xscale('log')
ax1.set_xlabel('$F_{max}/F$', fontsize=28)
ax1.set_ylabel('Photon Index', fontsize=28)


plt.show()