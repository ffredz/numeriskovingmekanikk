import numpy as np 
import matplotlib.pyplot as plt 

#1.3 a

G = 9.81 #m/s^2
init_values = [.1,.5,1]


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

def analytisk(theta_0,n,l,h):
    w0 = np.sqrt(G/l)
    t = np.arange(0,n*h,h)
    theta = np.zeros(n)
    theta[0]=theta_0
    for i in range(1,len(theta)):
        theta[i]=theta_0*np.cos(w0*t[i])
    return theta,t


for k in range(len(init_values)):
    stuff = euler(init_values[k], 1000,.001,.1)
    theta_values = stuff[0]/stuff[0][0]
    t = stuff[1]
    plt.plot(t,theta_values,label= f'euler med theta0 = {init_values[k]}')
    
anal = analytisk(1,1000,.1,.001)
plt.plot(anal[1],anal[0],label = "normalisert og linearisert analytisk")
plt.title("Normaliserte numeriske løsninger vs linearisert analytisk")
plt.xlabel("tid")
plt.ylabel("Vinkelutslag i radianer")
plt.grid()
plt.legend()
plt.show()



#1.3 b
T = 10
def energi(theta_0, n,h,l,m):
    w0 = np.sqrt(G/l)
    theta = np.zeros(n)
    theta[0]=theta_0 
    theta_hjelp = np.zeros(n)

    for i in range(1,len(theta)):
        theta_hjelp[i]=theta_hjelp[i-1]-h*(w0**2)*np.sin(theta[i-1])
        theta[i]=theta[i-1] + h*theta_hjelp[i] #symplektisk euler
        
    t = np.arange(0, n*h, h)

    U=m*G*l*(1-np.cos(theta))
    K=.5*m*l*l*theta_hjelp**2
    return U,K,t 


tidssteg = [10**(-3),10**(-2),10**(-1)]

for k in range(len(tidssteg)):
    n = int(T/tidssteg[k])
    energiinfo = energi(.5,n,tidssteg[k],.1,.1)
    potensiell = energiinfo[0]
    kinetisk = energiinfo[1]
    t = energiinfo[2]
    plt.plot(t,potensiell, label=f'Potensiell energi')
    plt.plot(t,kinetisk, label=f'kinetisk energi')
    totalenergi = potensiell+kinetisk
    plt.plot(t,totalenergi, label=f'Total energi')
    plt.grid()
    plt.title(f'Energi med tidssteg {tidssteg[k]}')
    plt.legend()
    plt.show()




