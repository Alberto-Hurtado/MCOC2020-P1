import scipy as sp
from scipy.integrate import odeint
import numpy as np
from time_reader import *
#import warnings
#warnings.filterwarnings("ignore")
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

#Parametros

#Unidades
g= 9.80665 #m/s/s
cm = 0.01 #m
inch = 2.54*cm
km = 1000. #m
mins= 60. #s
h = 60.*mins
#Datos importantes
G=6.67408e-11 #m3/ (kg s2)
U=-2*np.pi/(24.*h)	# radians/s
#Ω = -7.2921150e-5  

#U=0
#Radio de la tierra=6371km
rTierra=6371.*km
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
	return np.array([
	[np.cos(U*t) ,-np.sin(U*t),0],
	[np.sin(U*t) , np.cos(U*t),0],
	[     0      ,      0     ,1]
	],dtype=np.double)
def Rp(t):
	return U*np.array([
	[-np.sin(U*t),-np.cos(U*t),0],
	[ np.cos(U*t),-np.sin(U*t),0],
	[     0      ,      0     ,0]
	],dtype=np.double)
def Rpp(t):
	return (U**2)*np.array([
	[-np.cos(U*t), np.sin(U*t),0],
	[-np.sin(U*t),-np.cos(U*t),0],
	[     0      ,      0     ,0]
	],dtype=np.double)
def satelite(z,t):
	zp=np.zeros(6)
	zp[0:3] = z[3:6]
	r=np.sqrt(np.dot(z[0:3],z[0:3]))


	#Fg=(-G*mT/r**2)*(R(t)@(z[0:3]/r))
	
	#zp[3:6]=R(t).T@(Fg - (2*(Rp(t)@z[3:6]) + (Rpp(t)@z[0:3])))

	factor1=np.array(-G*mT*np.array(z[0:3])/(r**3))
	factor=np.array(Rpp(t)@z[0:3]+2*Rp(t)@z[3:6]).T	
	#igualdad=factor1 - (R.T)@factor
	zp[3:6] = (factor1 - (R(t).T)@factor)
	return zp

#vector de tiempo

t = np.linspace(0,intervalo_en_segundos,1001)
#Parte en y = rTierra + 700 km, vx=vx, vy=0

rH=700.*km #altura relativa a la superficie de la tierra(relative height)

z0 = np.array([-1300354.318621,1000042.793928,6872663.396330,-1692.331483,7262.981572,-1374.216961])
sol=odeint(satelite,z0,t)
zf = np.array([-2034623.889145,-5920783.955308,3289503.706424,-449.399307,3809.729180,6556.395258])


#Presentacion resultados
lector=False
if lector:
	indices=["X ","Y ","Z ","Vx","Vy","Vz"]
	[print(f"delta {indices[i]} = {sol[-1][i] - zf[i]}") for i in range(6)]
else:
	[print(sol[-1][i] - zf[i]) for i in range(6)]


#Graficar
graficar=False
if graficar:
	import plotter_1A
	plotter_1A.graficar_detalle(sol=sol,t=t)
