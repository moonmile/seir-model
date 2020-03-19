#include package
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#define differencial equation of seir model
def seir_eq(v,t,beta,lp,ip):
    return [
        -beta*v[0]*v[2],
        beta*v[0]*v[2]-(1/lp)*v[1],
        (1/lp)*v[1]-(1/ip)*v[2],
        (1/ip)*v[2]]

#solve seir model

ini_state=[3000,0,5,0]
t_max=100
dt=1
t=np.arange(0,t_max,dt)
plt.plot(t,odeint(seir_eq,ini_state,t,args=(0.001,14,7)))


plt.legend(['Susceptible','Exposed','Infected','Recovered'])
plt.show()
