import numpy as np 
import matplotlib.pyplot as plt 

#1.1

G = 9.81 #m/s^2


def diffeq(g,l):
    return 
def euler(theta_0,n,h,l):
    w0 = np.sqrt(G/l)
    theta = np.zeros(n)
    theta[0]=theta_0 
    theta_hjelp = np.zeros(n)

    for i in range(1,len(theta)):
        theta[i]=theta[i-1] + h*theta_hjelp[i-1]
        theta_hjelp[i]=theta_hjelp[i-1]-h*(w0**2)*np.sin(theta[i-1])
    t = np.arange(0, n*h, h)

    return theta,t

ting = euler(.1,1000,.001,.1)
theta = ting[0]
t = ting[1]
plt.plot(t,theta)
plt.show()