import matplotlib.pylab as plt
import numpy as np
from unidades import *

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

def graficar_deriva(sol,sol_real,t,show=True):
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
	if show:
		plt.show()
