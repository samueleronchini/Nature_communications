from numpy import *
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from numpy import NaN, Inf, arange, isscalar, asarray, array,sqrt
import matplotlib.ticker
from matplotlib import gridspec
import random
from scipy.integrate import simps

from scipy import optimize

from matplotlib.axes import Axes

import scipy.integrate as integrate
import scipy.special as special


nu_p=100.0
nu_p_plot=nu_p

f_r=4/3


gamma=100.0
beta=sqrt(1-1/gamma**2)

theta_max=5.*np.pi/180.0

alpha=-0.67
beta_s=-2.5 

nu_p=(-(alpha+2)/(beta_s+2))**(-1/(alpha-beta_s))*nu_p #comment if a Band function is used
#nu_p=nu_p/(2+alpha) #comment if a SBPL is used
nu_p_p=nu_p/(2*gamma)

exp_par=0.0

def broken(nu,nu_0):
	x=nu/nu_0
	y=((x)**(-alpha-1)+x**(-beta_s-1))**(-1.0)*exp(-exp_par*x)
	return y

def band(nu,nu_0):
	x=nu/nu_0
	if x<(alpha-beta_s):
		y=(x)**(alpha+1)*exp(-x)
	else:
		y=(alpha-beta_s)**(alpha-beta_s)*exp(-alpha+beta_s)*x**(beta_s+1)

	return y


def R_R0(theta,tau):
	return tau/(2*gamma**2*(1-beta*cos(theta)))


def D(theta):
	return 1/(gamma*(1-beta*cos(theta)))


def F_nu(tau,nu):

	theta_arr=np.linspace(0.0,theta_max,num=50)
	integrand_arr=theta_arr*0.0

	
	for n in range(0,len(integrand_arr)):
		if R_R0(theta_arr[n],tau)<1.0:

			integrand_arr[n]=(1.0)**(-lam)*broken(nu/(D(theta_arr[n])),nu_p_p*(1.0)**(-f_r-lam))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n]) #comment if a Band function is used
			#integrand_arr[n]=(1.0)**(-lam)*band(nu/(D(theta_arr[n])),nu_p_p*(1.0)**(-f_r-lam))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n]) #comment if a SBPL function is used
		else:
			integrand_arr[n]=(R_R0(theta_arr[n],tau))**(-lam)*broken(nu/(D(theta_arr[n])),nu_p_p*(R_R0(theta_arr[n],tau))**(-f_r-lam))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n]) #comment if a Band function is used
			#integrand_arr[n]=(R_R0(theta_arr[n],tau))**(-lam)*band(nu/(D(theta_arr[n])),nu_p_p*(R_R0(theta_arr[n],tau))**(-f_r-lam))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n]) #comment if a SBPL function is used
	y=simps(integrand_arr,theta_arr)
	return y

fig, ax1 = plt.subplots(figsize=(12,9))

for m in range(0,6):
	lam=float(m)/2.5

	tau_max=10.0*(1-beta*np.cos(theta_max))/(1-beta)

	tau_arr=np.logspace(0,np.log10(tau_max),num=100)
	gamma_arr=tau_arr*0.0

	nu_1=0.5
	nu_2=10.0

	for n in range(0,len(tau_arr)):
		gamma_arr[n]=1-np.log(F_nu(tau_arr[n],nu_2)/F_nu(tau_arr[n],nu_1))/np.log(10/0.5)

	def int_f(N0,gamma_int):
		return N0/(2-gamma_int)*(20**(2-gamma_int)-1)


	F_0=int_f(F_nu(tau_arr[0],nu_1),gamma_arr[0])


	f_ratio_arr=tau_arr*0.0


	for n in range(0,len(tau_arr)):
		f_ratio_arr[n]=int_f(F_nu(tau_arr[n],nu_1),gamma_arr[n])

	int_f_max=max(f_ratio_arr)
	for n in range(0,len(tau_arr)):
		f_ratio_arr[n]=int_f_max/int_f(F_nu(tau_arr[n],nu_1),gamma_arr[n])


	plt.xscale('log')
	plt.xlabel('$F_{max}/F$')
	plt.ylabel('Gamma')
	plt.plot(f_ratio_arr,gamma_arr, linewidth=3.0, label='$\\lambda=$%s' %lam)
plt.title('$\\Delta R\\propto$ R, $\\nu_p=$%s keV' %nu_p_plot, fontsize=20)




plt.xlim((0.7,10000))
plt.ylim((0.3,3.2))
ax1.set_xlabel('$F_{max}/F$', fontsize=28)
ax1.set_ylabel('Photon Index', fontsize=28)

plt.legend()


plt.show()
