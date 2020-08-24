import scipy as sp
from scipy.integrate import odeint
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

#Parametros

#Unidades
g= 9.81 #m/s/s
cm = 0.01 #m
inch = 2.54*cm
#Coeficiente de arrastre
p = 1.225 #kg/m3
cd = 0.47
D = 8.5*inch
r= D/2
A= sp.pi * r**2
CD=0.5*p*cd*A
m=15. #kg
#V = [0,10.,20.] #m/s 
V = 20. #m/s 

#Funcion a integrar
#z es el vector de estado
#z=[x,y,vx,vy]
#	dz/dt = bala(z,t)
#         [    z2    ]
# dz/dt = [          ] (modelo)
#         [ FD/m - g ]

# z[0] -> x
# z[1] -> y
# z[2] -> vx
# z[3] -> vy
def bala(z,t):
	zp=sp.zeros(4)
	zp[0] = z[2]
	zp[1] = z[3]
	v = z[2:4]
	v[0]=v[0]- V
	v2 = sp.dot(v,v)
	vnorm = sp.sqrt(v2)
	FD = -CD*v2*(v/vnorm)
	zp[2] = FD[0]/m
	zp[3] = FD[1]/m  -g
	return zp

#vector de tiempo
t = sp.linspace(0,10,1001)

#Parte en el origen con vx=vy=2 m/s
vi=100.*1000./3600.
z0 = sp.array([0., 0., vi,vi])
soluciones=[]
for V in [0,10.,20.]:
	soluciones.append(odeint(bala,z0,t))
#print(soluciones[0])
import matplotlib.pylab as plt
labels=["V = 0 m/s","V = 10 m/s","V = 20 m/s"]

for i in range(3):
	x=soluciones[i][:,0]
	y= soluciones[i][:,1]
	plt.figure(1)
	plt.plot(x,y,label=labels[i])
	plt.legend()
	plt.grid(True)
	plt.title("Trayectoria para distintos vientos")
	plt.ylabel("Y (m)")
	plt.xlabel("X (m)")
	plt.axis([0, 160, 0, 50])

plt.savefig("balistica.png")
plt.show()
