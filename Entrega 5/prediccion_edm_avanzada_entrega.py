import scipy as sp
from scipy.integrate import odeint
import numpy as np
from leer_eof import leer_eof
from time import perf_counter
import sys
# Como recibir un parametro desde la terminal

nombre = sys.argv[1]
import warnings
warnings.filterwarnings("ignore")
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

#Parametros

#Unidades
g= 9.81998 #m/s/s
cm = 0.01 #m
inch = 2.54*cm
km = 1000. #m
mins= 60. #s
h = 60.*mins
#Datos importantes
G=6.67408e-11 #m3/ (kg s2)
U= 2*np.pi/(24.*h)	# radians/s

#Radio de la tierra=6371km
rTierra=6371.*km
mT=5.972e24 #kg
m=2170. #kg

#Ajuste gravitacional
def P(n,m,r,the): # Polinomio de Legendre
	if   n==2 and m==0:
		pnm = (-1. + 3.*np.sin(the)**2)/(2.*r**3)
	elif n==2 and m==1:
		pnm = (3.*np.sin(the)*np.cos(the))/(r**3)
	elif n==2 and m==2:
		pnm = (3.*np.cos(the)**2.)/(r**3)


	elif n==3 and m == 0:
		pnm = np.sin(the)*(-3. + 5.*np.sin(the)**2)/(2.*r**4)
	elif n==3 and m == 1:
		pnm = np.cos(the)*(-1. + 5.*np.sin(the)**2)*3./(2.*r**4)
	elif n==3 and m == 2:
		pnm = (15.*np.sin(the)*(np.cos(the)**2))/(r**4)
	elif n==3 and m == 3:
		pnm = (15. * np.cos(the)**3)/(r**4)

	return pnm

def C(n,m):
	if n==2:
		if m==1:
			return -0.3504890360e-09
		elif m==2:
			return 0.1574536043e-05
	elif n==3:
		if m==1:
			return 0.2192798802e-05
		elif m==2:
			return 0.3090160446e-06
		elif m==3:
			return 0.1005588574e-06
	elif n==4:
		if   m==1:
			return -0.5087253036e-06
		elif m==2:
			return 0.7841223074e-07
		elif m==3:
			return 0.5921574319e-07
		elif m==4:
			return -0.3982395740e-08	

def S(n,m):
	if n==2:
		if m==1:
			return 0.1635406077e-08
		elif m==2:
			return -0.9038680729e-06
	elif n==3:
		if m==1:
			return 0.2680118938e-06
		elif m==2:
			return -0.2114023978e-06
		elif m==3:
			return 0.1972013239e-06
	elif n==4:
		if   m==1:
			return -0.4494599352e-06
		elif m==2:
			return 0.1481554569e-06
		elif m==3:
			return -0.1201129183e-07
		elif m==4:
			return 0.6525605810e-08

def J(n): # km**n+3  /s**2
	R=6378.1363#km
	if   n==2:
		return -(R**n)*398600.4415*-0.1082635854e-02
	elif n==3:
		return -(R**n)*398600.4415*0.2532435346e-05
	elif n==4:
		return -(R**n)*398600.4415*0.1619331205e-05
	elif n==5:
		return -(R**n)*398600.4415*0.2277161016e-06
	elif n==6:
		return -(R**n)*398600.4415*-0.5396484906e-06
	elif n==7:
		return -(R**n)*398600.4415*0.3513684422e-06
	elif n==8:
		return -(R**n)*398600.4415*0.2025187152e-06


def sumatoria_1(rho,phi,the,Nz=2):
	total=0
	for n in range(Nz+1)[2:]:
		total += J(n)*P(n,0,rho,the)
	return total	

def sumatoria_2(rho,phi,the,Nz=2):
	total=0
	for n in range(Nz+1)[2:]:
		for m in range(n+1)[1:]:
			total += P(n,m,rho,the)*(C(n,m)*np.cos(m*phi) + S(n,m)*np.sin(m*phi))
	return total		

def F2(loc,x,y,z,r):#km/s2
	if loc=='x':
		return J(2)*x*(6.*z**2 - 3.*(x**2 + y**2)/2.)/r**7
	elif loc=='y':
		return J(2)*y*(6.*z**2 - 3.*(x**2 + y**2)/2.)/r**7
	elif loc=='z':
		return J(2)*z*(3.*z**2 - 9.*(x**2 + y**2)/2.)/r**7

def F3(loc,x,y,z,r):#km/s2

	if loc=='x':
		return J(3)*x*z*(10.*z**2 - 15.*(x**2 + y**2)/2.)/r**9
	elif loc=='y':
		return J(3)*y*z*(10.*z**2 - 15.*(x**2 + y**2)/2.)/r**9
	elif loc=='z':
		return J(3)*(4.*z**2 * (z**2 - 3.*(x**2 + y**2)) + 1.5*(x**2 + y**2)**2)/r**9

def eulint(fun,z0,t,Nsubdivisiones = 1.):

	sol=[np.array(z0).copy()]

	for i in range(1,len(t)):
		dt = (t[i] - t[i-1])/Nsubdivisiones

		zi1=sol[-1].copy()
		for n in range(int(Nsubdivisiones)):
			zi1=zi1 +dt*fun(zi1,t[i])
		sol.append(zi1)

	return np.array(sol)

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


def convertir_a_esfericas(z0):
	x=z0[0]
	y=z0[1]
	z=z0[2]
	rho = np.sqrt(x**2 + y**2 + z**2)
	the = np.arcsin(z/rho)
	phi = np.arcsin(y/(rho*np.cos(the)))

	return np.array([rho,phi,the])

def satelite(z,t):
	zp=np.zeros(6)
	zp[0:3] = z[3:6]
	rho,phi,the=convertir_a_esfericas(z)
	
	r=rho

	factor1=np.array(-G*mT*np.array(z[0:3])/(r**3))
	
	factor=np.array(Rpp(t)@z[0:3]+2*Rp(t)@z[3:6]).T	
	zp[3:6] = (factor1 - (R(t).T)@factor)
	zp[3] += F2('x',z[0]/km,z[1]/km,z[2]/km,r/km)*km
	zp[4] += F2('y',z[0]/km,z[1]/km,z[2]/km,r/km)*km
	zp[5] += F2('z',z[0]/km,z[1]/km,z[2]/km,r/km)*km
	zp[3] += F3('x',z[0]/km,z[1]/km,z[2]/km,r/km)*km
	zp[4] += F3('y',z[0]/km,z[1]/km,z[2]/km,r/km)*km
	zp[5] += F3('z',z[0]/km,z[1]/km,z[2]/km,r/km)*km
	return zp

t, x, y, z, vx, vy, vz = leer_eof(nombre)


z0 = np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]])



t1 = perf_counter()
sol=odeint(satelite,z0,t)
#sol=eulint(satelite,z0,t,Nsubdivisiones=10001)
t2 = perf_counter()

#print(f"Elapsed time = {t2-t1} s")
sol_real=sol.copy()
sol_real[:,0]=x
sol_real[:,1]=y
sol_real[:,2]=z
sol_real[:,3]=vx
sol_real[:,4]=vy
sol_real[:,5]=vz

zf = np.array([x[-1],y[-1],z[-1],vx[-1],vy[-1],vz[-1]])

Vsol = np.sqrt(sol[-1][3]**2 + sol[-1][4]**2 +sol[-1][5]**2)
Vf =np.sqrt(zf[3]**2 + zf[4]**2 +zf[5]**2)

#Presentacion resultados
lector=False
if lector:
	indices=["X ","Y ","Z ","Vx","Vy","Vz"]
	[print(f"delta {indices[i]} = {sol[-1][i] - zf[i]}") for i in range(6)]
	print(f"delta V= {Vsol-Vf}")
else:
	[print(sol[-1][i] - zf[i]) for i in range(6)]

#Graficar
graficar=True
if graficar:
	import matplotlib.pylab as plt
	import numpy as np

	vuelta=np.linspace(0,2*np.pi,361)
	mediaVuelta=np.linspace(0,np.pi,101)

	def crearCarpeta(name="New Folder"):
		import os
		try:
			os.mkdir(name)
		except OSError:
			print (f"Creation of the directory {name} failed")

	def graficar_detalle(sol,t,savefig=False,name="prediccion_orbita",mostrarTierra=False,nameFold="Graficos",show=True,apellido=""):
		if savefig:
			crearCarpeta(name=nameFold)
		(x,y,z)=(sol[:,0],sol[:,1],sol[:,2])
		oneT=np.ones(len(t))
		vuelta=np.linspace(0,2*np.pi,361)
		mediaVuelta=np.linspace(0,np.pi,101)
		rho=np.sqrt(x**2 + y**2 + z**2)
		theta=np.arcsin(y/rho)
		phi  =np.arcsin(z/rho)
		V=np.sqrt(sol[0,3]**2+sol[0,4]**2+sol[0,5]**2)
		#plt.polar(np.degrees(theta),(rho-rTierra)/km,label=f"Vx = {sol[0,3]} m/s")

		plt.figure(0)
		plt.plot(x/km,y/km,label=f"Vo {apellido}= {int(V)} m/s")
		if show:
			plt.plot(rTierra*np.cos(vuelta)/km,rTierra*np.sin(vuelta)/km,label="Superficie terrestre")
		plt.legend()
		plt.grid(True)
		plt.title("Trayectoria satelital plano XY")
		plt.ylabel("Y (km)")
		plt.xlabel("X (km)")
		#plt.axis([0,2*np.pi,0,(rH+100.)/km])
		if savefig:
			plt.savefig(f"{nameFold}/{name}_plano_XY.png")

		plt.figure(1)
		#plt.plot(theta,(rho-rTierra)/km)
		plt.plot(t/h,theta)
		plt.grid(True)
		plt.title("Angulo theta versus tiempo")
		plt.ylabel("Theta (radianes)")
		plt.xlabel("Tiempo (horas)")
		if savefig:
			plt.savefig(f"{nameFold}/{name}_theta_vs_t.png")

		plt.figure(2)
		plt.plot(t/h,(rho-rTierra)/km,label=f"Altura satelite {apellido}")
		if show:
			plt.plot(t/h,oneT*700.,"--",label="Altura = 700 km")
			plt.plot(t/h,80.*oneT,label="Atmosfera")
		plt.grid(True)
		plt.legend(loc='lower center')
		plt.title("Altura del satelite versus tiempo")
		plt.ylabel("Altura (km)")
		plt.xlabel("Tiempo (horas)")
		if savefig:
			plt.savefig(f"{nameFold}/{name}_rho_vs_t.png")

		plt.figure(3)
		plt.subplot(3,1,1)
		plt.plot(t/h,x/km,label=f"Eje X {apellido}")
		plt.legend(loc='lower left')
		plt.title("Altura del satelite versus tiempo")
		plt.grid(True)
		plt.subplot(3,1,2)
		plt.plot(t/h,y/km,label=f"Eje Y {apellido}")
		plt.legend(loc='lower left')
		plt.ylabel("Distancia al origen (km)")
		plt.grid(True)
		plt.subplot(3,1,3)
		plt.plot(t/h,z/km,label=f"Eje Z {apellido}")
		plt.legend(loc='lower left')
		plt.grid(True)
		plt.xlabel("Tiempo (horas)")
		if savefig:
			plt.savefig(f"{nameFold}/{name}_h_vs_t.png")

		plt.figure(4)
		axis=plt.axes(projection="3d")
		if mostrarTierra and show:
			for Phi in mediaVuelta:
				axis.plot((rTierra/km)*np.cos(vuelta)*np.sin(Phi),(rTierra/km)*np.sin(vuelta)*np.sin(Phi),np.ones(361)*(rTierra/km)*np.cos(Phi),"g")
		#axis.plot_surface((rTierra/km)*np.cos(vuelta)*np.sin(mediaVuelta),(rTierra/km)*np.sin(vuelta)*np.sin(mediaVuelta),np.ones(361)*(rTierra/km)*np.cos(mediaVuelta))
		axis.plot3D(x/km,y/km,z/km)
		if savefig:
			plt.savefig(f"{nameFold}/{name}_3D.png")
		
		if show:
			plt.show()

	def graficar_deriva(sol,sol_real,t):
		x=sol[:,0]
		y=sol[:,1]
		z=sol[:,2]
		x_real=sol_real[:,0]
		y_real=sol_real[:,1]
		z_real=sol_real[:,2]
		deriva= np.sqrt((x-x_real)**2 + (y-y_real)**2 + (z-z_real)**2)
		plt.plot(t/h,deriva/km)
		plt.grid(True)
		plt.title(f"Deriva en funcion del tiempo, Dmax = {int(max(deriva/km))} km")
		plt.ylabel("Distancia de deriva")
		plt.xlabel("Tiempo (horas)")
		plt.show()

	graficar_deriva(sol,sol_real,t)
	#graficar_detalle(sol=sol,t=t,savefig=False,mostrarTierra=True,show=False)
	#graficar_detalle(sol=sol_real,t=t,savefig=True,mostrarTierra=True,nameFold="Graficos eulint",name="orbita_real",apellido="real")

