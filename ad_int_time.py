from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps


nu_p=100.0
lam=0.4

R0=2e15
gamma=100.0

c=3e10 

nu_p_plot=nu_p


f_r=4/3



beta=sqrt(1-1/gamma**2)

theta_max=5*np.pi/180.0

alpha=-0.667
beta_s=-2.4

nu_p=(-(alpha+2)/(beta_s+2))**(-1/(alpha-beta_s))*nu_p

nu_p_p=nu_p/(2*gamma)


def broken(nu,nu_0):
	x=nu/nu_0
	y=((x)**(-alpha-1)+x**(-beta_s-1))**(-1.0)
	return y


def R_R0(theta,tau):
	return tau/(2*gamma**2*(1-beta*cos(theta)))


def D(theta):
	return 1/(gamma*(1-beta*cos(theta)))


def F_nu(tau,nu):

	theta_arr=np.linspace(0.0, theta_max,num=100)
	integrand_arr=theta_arr*0.0

	
	for n in range(0,len(integrand_arr)):
		if R_R0(theta_arr[n],tau)<1.0:
			integrand_arr[n]=(1.0)**(-lam)*broken(nu/(D(theta_arr[n])),nu_p_p*(1.0)**(-f_r-lam))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n])
		else:
			integrand_arr[n]=(R_R0(theta_arr[n],tau))**(-lam)*broken(nu/(D(theta_arr[n])),nu_p_p*(R_R0(theta_arr[n],tau))**(-f_r-lam))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n])

	y=simps(integrand_arr,theta_arr)
	return y

fig, ax1 = plt.subplots(figsize=(12,9))



plt.xlim((95.0,550))
plt.ylim((1e-3,2.0))

plt.xscale('log')
plt.ylabel('$F/F_{max}$', fontsize=28)
plt.xlabel('$\\delta t_{obs}(s)+100$ s', fontsize=28)
plt.yscale('log')

for m in range(0,6):
	lam=0.4*m

	tau_max=10*(1-beta*np.cos(theta_max))/(1-beta)

	tau_arr=np.logspace(0,np.log10(tau_max),num=50)
	gamma_arr=tau_arr*0.0

	nu_1=0.5
	nu_2=10.0

	for n in range(0,len(tau_arr)):
		gamma_arr[n]=1-np.log(F_nu(tau_arr[n],nu_2)/F_nu(tau_arr[n],nu_1))/np.log(10/0.5)

	def int_f(N0,gamma):
		return N0/(2-gamma)*(20**(2-gamma)-1)


	f_ratio_arr=tau_arr*0.0

	for n in range(0,len(tau_arr)):
		f_ratio_arr[n]=int_f(F_nu(tau_arr[n],nu_1),gamma_arr[n])

	int_f_max=max(f_ratio_arr)
	for n in range(0,len(tau_arr)):
		f_ratio_arr[n]=int_f_max/int_f(F_nu(tau_arr[n],nu_1),gamma_arr[n])


	line1,= plt.plot(R0/c*(1-beta)*(tau_arr-1)+100,f_ratio_arr**(-1.),linewidth=3.0, label='$\\lambda=$%.1f' %lam)

plt.legend()


plt.show()

