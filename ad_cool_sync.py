from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

fig, ax1 = plt.subplots(figsize=(12,9))

nu_ch=100

f_r=4/3


theta_max=5.0*np.pi/180.0

lam=0.7

nu_1=0.5
nu_2=10.0


nu_arr=[]
h_tot_arr=[]


gamma=100.0
beta=sqrt(1-1/gamma**2)

sync_tab_R1= open("sync_tab.txt", 'r')

for line in sync_tab_R1:
	app=line.split()
	nu_arr.append(float(app[0]))
	h_tot_arr.append(float(app[1]))

def S_nu(nu,nu_0):
	nu_in=nu/nu_0
	k=0
	while k<len(nu_arr)-2 and nu_in>nu_arr[k]:
		k+=1

	return (h_tot_arr[k]+h_tot_arr[k+1])/2



def R_R0(theta,tau):
	return tau/(2*gamma**2*(1-beta*cos(theta)))


def D(theta):
	return 1/(gamma*(1-beta*cos(theta)))

def F_nu(tau,nu_f,lamb,nu_c_p):
	theta_arr=np.linspace(0.0,theta_max,num=100)
	integrand_arr=theta_arr*0.0

		
	for n in range(0,len(integrand_arr)):
		if R_R0(theta_arr[n],tau)<1.0:
			integrand_arr[n]=(1.0)**(-lamb)*S_nu(nu_f/(D(theta_arr[n])),nu_c_p*(1.0))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n])
		else:
			integrand_arr[n]=(R_R0(theta_arr[n],tau))**(-lamb)*S_nu(nu_f/(D(theta_arr[n])),nu_c_p*(R_R0(theta_arr[n],tau))**(-f_r-lamb))*(D(theta_arr[n]))**3*cos(theta_arr[n])*sin(theta_arr[n])
	y=simps(integrand_arr,theta_arr)
	return y




def int_f(N0,gamma_s):
	return N0/(2-gamma_s)*(20**(2-gamma_s)-1)



nu_c_p= 1/(2*gamma)*nu_ch

tau_max=10*(1-beta*np.cos(theta_max))/(1-beta)
tau_arr=np.logspace(0,np.log10(tau_max),50)
gamma_arr=tau_arr*0.0

nu_1=0.5
nu_2=10.0

for n in range(0,len(tau_arr)):
	gamma_arr[n]=1-np.log(F_nu(tau_arr[n],nu_2,lam,nu_c_p)/F_nu(tau_arr[n],nu_1,lam,nu_c_p))/np.log(10/0.5)


int_f_arr=tau_arr*0.0

for n in range(0,len(tau_arr)):
	int_f_arr[n]=(int_f(F_nu(tau_arr[n],nu_1,lam,nu_c_p),gamma_arr[n])/int_f(F_nu(tau_arr[0],nu_1,lam,nu_c_p),gamma_arr[0]))**(-1.)

plt.plot(int_f_arr,gamma_arr)
plt.xscale('log')
ax1.set_xlabel('$F_{max}/F$', fontsize=28)
ax1.set_ylabel('Photon Index', fontsize=28)

plt.show()
