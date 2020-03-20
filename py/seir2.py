#include package
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#define differencial equation of seir model
def seir_eq(v,t,alpha,beta,gamma,N):

#    ds = - beta * I / N * S             // dS/dt = -βI/N*S
#    de = beta * I / N * S - alpha * E   // dE/dt = βI/N*S-αE
#    di = alpha * E - gamma * I          // dI/dt = αE - γI
#    dr = gamma * I                      // dR/dt = γI 

    return [
        -beta*v[0]*v[2]/N,
        beta*v[0]*v[2]/N-alpha*v[1],
        alpha*v[1]-gamma*v[2],
        gamma*v[2]]

#solve seir model

ini_state=[3000,0,5,0]
t_max=100
dt=1
t=np.arange(0,t_max,dt)

N = ini_state[0] + ini_state[1] + ini_state[2] + ini_state[3]

lp = 14
ip = 7
R0 = 1.0
D  = R0*ip
alpha = 1/lp
gamma = 1/ip
R0 = 10.0
beta = R0*gamma   # 0.001*N
# R0 = beta/gamma 

print('R0 = ', R0 )
plt.plot(t,odeint(seir_eq,ini_state,t,args=(alpha, beta, gamma, N)))


plt.legend(['Susceptible','Exposed','Infected','Recovered'])
plt.show()
