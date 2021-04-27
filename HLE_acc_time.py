from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps


nu_p=50
c=3e10 

gamma_0=100.0
R_in=1e16
R_off=5e16

alpha=-0.667
beta_s=-2.5

nu_p=(-(alpha+2)/(beta_s+2))**(1/(alpha-beta_s))*nu_p
nu_p_p=nu_p/(2*gamma_0)

theta_max=5*np.pi/180

def broken(nu,nu_0):
	x=nu/nu_0
	y=((x)**(-alpha-1)+x**(-beta_s-1))**(-1.0)
	return y

def gamma(R):
	if gamma_0*(R/R_in)**k<1.0:
		y=1.0
	else:
		y=gamma_0*(R/R_in)**k
	return y 

def beta(R):
	return sqrt(1-1/(gamma(R))**2)

def D(theta,R):
	return 1/(gamma(R)*(1-beta(R)*cos(theta)))

def R_star(t):
	return ((1+2*(1-2*k)*c*t*gamma_0**2/R_in)**(1/(1-2*k)))*R_in

def th(t,R_em):
	y=arccos(1-c*t/R_em+R_in/(2*gamma_0**2*(1-2*k)*R_em)*((R_em/R_in)**(1-2*k)-1))
	return y


def F_nu(t,nu,k,lam):

	R_arr=np.linspace(R_in,min(R_star(t),R_off),50)
	R_arr=np.flip(R_arr)
	theta_arr=[]
	n=0
	while n<len(R_arr) and th(t,R_arr[n])<theta_max:
		theta_arr.append(th(t,R_arr[n]))
		n+=1		

	theta_arr=np.array(theta_arr)
	integrand_arr=theta_arr*0.0

	if len(theta_arr)<=1:
		y=0.0
	else:
		for n in range(0,len(integrand_arr)):
			integrand_arr[n]=(R_arr[n]/R_off)**(-lam)*broken(nu/D(theta_arr[n],R_arr[n]),nu_p_p*(R_arr[n]/R_off)**(-lam))*(D(theta_arr[n],R_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n])

		y=simps(integrand_arr,theta_arr)
	return y


fig, ax1 = plt.subplots(figsize=(12,9))

k=0.0
lam=0.3

t_arr=np.logspace(-3,3,300) #this is the observer time
gamma_arr=t_arr*0.0

nu_1=0.5
nu_2=10.0

for n in range(0,len(t_arr)):
	if F_nu(t_arr[n],nu_1,k,lam)==0.0:
		gamma_arr[n]=10.0
	else:
		gamma_arr[n]=1-np.log(F_nu(t_arr[n],nu_2,k,lam)/F_nu(t_arr[n],nu_1,k,lam))/np.log(10/0.5)

def int_f(N0,gamma):
	return N0/(2-gamma)*(20**(2-gamma)-1)

int_f_arr=t_arr*0.0

for n in range(0,len(t_arr)):
	int_f_arr[n]=int_f(F_nu(t_arr[n],nu_1,k,lam),gamma_arr[n])


int_f_max=max(int_f_arr)
l=np.argmax(int_f_arr)

for n in range(0,len(t_arr)):
	int_f_arr[n]=int_f_max/int_f(F_nu(t_arr[n],nu_1,k,lam),gamma_arr[n])

plt.plot(t_arr+100.0-t_arr[l],int_f_arr**(-1),label='$\\lambda$=%s, k=%.1f' %(round(lam,1),k),linewidth=3.0)

plt.legend(framealpha=1.0, fontsize=14)
plt.xlim(80,500)
plt.ylim(1e-3,2)
ax1.set_xlabel('$time+100s$', fontsize=28)
ax1.set_ylabel('$F/F_{max}$', fontsize=28)
plt.xscale('log')
plt.yscale('log')

plt.show()
