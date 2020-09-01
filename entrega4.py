import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

t = np.linspace(0,7.5,100)

def funcion(x,t):
	xp=np.zeros(2)
	xp[0]=x[1]
	xp[1]= np.exp(-1.25664*t)*(0.999991*np.cos(1.53906*t) - 3.3816*np.sin(1.53906*t))
	return xp

def eulint(fun,z0,t,Nsubdivisiones = 1.):

	sol=[np.array(z0).copy()]

	for i in range(1,len(t)):
		dt = (t[i] - t[i-1])/Nsubdivisiones

		zi1=sol[-1].copy()
		for n in range(int(Nsubdivisiones)):
			zi1=zi1 +dt*fun(zi1,t[i])
		sol.append(zi1)

	return np.array(sol)

x0=np.array([1,1])

sol_ana=np.exp(-1.25664*t)*(1.46624*np.sin(1.53906*t)+np.cos(1.53906*t))
sol_ode=odeint(funcion,x0,t)
sol_eul  =eulint(funcion,x0,t)
sol_eul10=eulint(funcion,x0,t,Nsubdivisiones=10.)
sol_eul100=eulint(funcion,x0,t,Nsubdivisiones=100.)



plt.plot(t,sol_ana,'k',linewidth=2,label="Sol. Analitica")
plt.plot(t,sol_ode   [:,1],'b',label="Sol. ODE Integrate")
plt.plot(t,sol_eul   [:,1],'g--',label="Sol. Euler Int N=1")
plt.plot(t,sol_eul10 [:,1],'r--',label="Sol. Euler Int N=10")
plt.plot(t,sol_eul100[:,1],color = '#ff7f0e',linestyle='--',label="Sol. Euler Int N=100")
plt.legend()
plt.show()

