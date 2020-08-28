import scipy as sp
from scipy.integrate import odeint
import numpy as np
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

#Parametros

#Unidades
g= 9.81 #m/s/s
cm = 0.01 #m
inch = 2.54*cm
km = 1000. #m
mins= 60. #s
h = 60.*mins
#Datos importantes
G=6.67408e-11 #m3/ (kg s2)
U=2*np.pi/(24.*h)
mT=5.972e24 #kg
m=2170. #kg

#Funcion a integrar
#z es el vector de estado
#z=[x,y,vx,vy]
#	dz/dt = satelite(z,t)
#         [    z2    ]
# dz/dt = [          ] (modelo)
#         [ FD/m - g ]

# z[0] -> x
# z[1] -> y
# z[2] -> z
# z[3] -> vx
# z[4] -> vy
# z[5] -> vz
#Matrices de rotacion
def R(t):
	return np.matrix([
	[np.cos(U*t) ,-np.sin(U*t),0],
	[np.sin(U*t) , np.cos(U*t),0],
	[     0      ,      0     ,1]
	])
def Rp(t):
	return U*np.matrix([
	[-np.sin(U*t),-np.cos(U*t),0],
	[ np.cos(U*t), np.sin(U*t),0],
	[     0      ,      0     ,0]
	])
def Rpp(t):
	return (U**2)*np.matrix([
	[-np.cos(U*t), np.sin(U*t),0],
	[-np.sin(U*t),-np.cos(U*t),0],
	[     0      ,      0     ,0]
	])
def satelite(z,t):
	zp=np.zeros(6)
	zp[0:3] = z[3:6]
	r=np.sqrt(np.dot(z[0:3],z[0:3]))

	factor=np.array(Rpp(t)@z[0:3]+2*Rp(t)@z[3:6]).T
	factor1=np.array(-G*mT*np.array(z[0:3])/(r**3)).reshape(3,1)

	igualdad=factor1 - (R(t).T)@factor
	zp[3:6] = (factor1 - (R(t).T)@factor).reshape(1,3)
	return zp

import matplotlib.pylab as plt
#vector de tiempo
t = np.linspace(0,98.6*mins,1001)
#Parte en y = rTierra + 700 km, vx=vx, vy=0
V=np.linspace(0,10000,1001)
#V=[10000.]
#Radio de la tierra=6371km
rTierra=6371.*km
rH=700.*km #altura relativa a la superficie de la tierra(relative height)

for vt in V:
	z0 = np.array([rH+rTierra, 0.,0., 0.,vt,0.])
	sol=odeint(satelite,z0,t)
	(x,y)=(sol[:,0],sol[:,1])
	rho=np.sqrt(x**2 + y**2)
	theta=np.arcsin(y/rho)
	#plt.polar(np.degrees(theta),(rho-rTierra)/km,label=f"Vx = {int(vt)} m/s")
	plt.plot(x/km,y/km,label=f"Vx = {int(vt)} m/s")
	#plt.plot(theta,(rho-rTierra)/km,label=f"Vx = {int(vt)} m/s")
#plt.axis([0,2*np.pi,0,(rH+100.)/km])
plt.figure(1)
plt.legend()
plt.grid(True)
plt.title("Trayectoria para distintas velocidades tangenciales")
#plt.ylabel("Y (km)")
#plt.xlabel("X (km)")
plt.savefig("integracion_edm_basica.png")
plt.show()
