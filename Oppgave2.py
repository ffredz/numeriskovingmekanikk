import numpy as np 
from scipy.integrate import solve_ivp 
import matplotlib.pyplot as plt 

M = 1 #bjelkemasse (kg)
m = .1 #pendelmasse (kg)
L = 1 #pendellengde (m)
G = 9.81 #tyngdeakselerasjon (kgm/s^2 (eller newton da for faen))
M_tot = M+2*m #massen til systemen

def deriverte(t, U):
    theta1, omega1, theta2, omega2 = U 

    #definerer matrisen A
    a11 = L-m*L*(np.cos(theta1)**2)/M_tot
    a12 = -np.cos(theta1)*np.cos(theta2)*m*L/M_tot
    a21 = -np.cos(theta1)*np.cos(theta2)*m*L/M_tot 
    a22 = L-m*L*(np.cos(theta2)**2)/M_tot
    
    A = np.array([[a11,a12],[a21,a22]])

    fellesfaktor = m*L*(np.sin(theta1)*omega1**2+np.sin(theta2)*omega2**2)/M_tot
    
    #så var den den andre vektoren ja 
    b1 = -G*np.sin(theta1)-np.cos(theta1)*fellesfaktor 
    b2 = -G*np.sin(theta2)-np.cos(theta2)*fellesfaktor
    B = np.array([b1,b2])

    akselerasjoner = np.linalg.solve(A,B) #løse for akselerasjonene 
    return [omega1, akselerasjoner[0],omega2, akselerasjoner[1]]

#initialbetingelser
U0 = [.5,0,0,0]
t_span = (0,30)
t_eval = np.linspace(0,30,1000)

resultat = solve_ivp(deriverte,t_span, U0,t_eval=t_eval, method = 'RK45')
print(resultat.message)

plt.figure(figsize=(10,5))
plt.plot(resultat.t, resultat.y[0], label = "theta_1(init .5rad")
plt.plot(resultat.t, resultat.y[2], label = "theta_2(init 0rad)")
plt.title("Oscillasjoner av pendlene")
plt.xlabel("tid")
plt.ylabel("vinkel i rad")
plt.legend()
plt.grid()
#plt.show()


V = -m*L*G*(np.cos(resultat.y[0])+np.cos(resultat.y[2])) #potensiell energi..

vel_bjelke = -m*L*(resultat.y[1]*np.cos(resultat.y[0]) + resultat.y[3]*np.cos(resultat.y[2]))/M_tot #tidsderivert av x_b..

vel_y_1 = L*resultat.y[1]*np.sin(resultat.y[0])
vel_y_2 = L*resultat.y[3]*np.sin(resultat.y[2])

vel_x_1 = vel_bjelke + L*resultat.y[1]*np.cos(resultat.y[0])
vel_x_2 = vel_bjelke + L*resultat.y[3]*np.cos(resultat.y[2])

K = .5*M*vel_bjelke**2 +.5*m*(vel_x_1**2+vel_y_1**2+vel_x_2**2+vel_y_2**2) #kinetisk energi...


plt.plot(resultat.t, V, label = "Potensiell energi")
plt.plot(resultat.t, K, label= "kinetisk energi")
plt.plot(resultat.t, V+K, label="total energi")
plt.title("energiplot")
plt.grid()
plt.legend()
plt.show()



#vi plotter også massesenteret for å bekrefte løsningen

pos_bjelke = -(m*L/M_tot)*np.sin(resultat.y[0]+np.sin(resultat.y[2])) #løst fra bevaringen...

pos_x1 = pos_bjelke +L*np.sin(resultat.y[0])
pos_x2 = pos_bjelke +L*np.sin(resultat.y[2]) #driter i delta siden den er konstant og blir kansellert senere uansett

x_cm = (M*pos_bjelke + m*(pos_x1+pos_x2))/M_tot 

plt.figure(figsize=(10,5))
plt.plot(resultat.t, x_cm, label="Massesenter $x_{cm}$")
plt.axhline(0, color='black', linestyle='--', alpha=0.5) # Referanselinje ved 0
plt.title("Massesenterets posisjon over tid")
plt.xlabel("Tid (s)")
plt.ylabel("Posisjon (m)")
plt.ylim(-0.1, 0.1) # Zoom inn for å se at det er flatt
plt.grid()
plt.legend()
plt.show()