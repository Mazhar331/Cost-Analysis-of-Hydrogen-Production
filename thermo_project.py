#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy.polynomial as poly
import matplotlib.pyplot as plt
import math

R=8.314 # J/(K*mol)
T0=298.15 # K
P0=1 #bar

def mcph(T0,T,A,B,C,D):
    tau=T/T0
    cph_by_R=A+(B/2)*T0*(tau+1)+(C/3)*(T0**2)*(tau**2+tau+1)+D/(tau*(T0**2))
    return(cph_by_R)

def mcps(T0,T,A,B,C,D):
    tau=T/T0
    cps_by_R=A+(B*T0+(C*(T0**2)+D/((tau**2)*(T0**2)))*((tau+1)/2))*((tau-1)/math.log(tau))
    return(cps_by_R)

def solver(P0,P,K):
    # Equilibrium equation involving mole fractions rearranged in the form of a polynomial, feed ratio of H2O/CH4=3 mol/mol
    eq=poly.Polynomial([-48*K*((P0/P)**2),16*K*((P0/P)**2),36*K*((P0/P)**2),0,27-4*K*((P0/P)**2)])
    solns=eq.roots()
    return(solns)

# Reaction: CH4 + H2O (steam) -->  CO + 3H2

CH4_Hfo_T0=-74520 # J/mol
H2Osteam_Hfo_T0=-241818
CO_Hfo_T0=-110525
H2_Hfo_T0=0
deltaHo_T0=3*H2_Hfo_T0+CO_Hfo_T0-CH4_Hfo_T0-H2Osteam_Hfo_T0

CH4_Gfo_T0=-50460 # J/mol
H2Osteam_Gfo_T0=-228572
CO_Gfo_T0=-137169
H2_Gfo_T0=0
deltaGo_T0=3*H2_Gfo_T0+CO_Gfo_T0-CH4_Gfo_T0-H2Osteam_Gfo_T0

for P in [1,3,5,10,15,20]: # bar
    T_array=list()
    epsilon_array=list()
    for T in range(973,1273): # K
        T_array.append(T)
        CH4_cph_by_R=mcph(T0,T,1.702,9.081e-3,-2.164e-6,0)
        H2Osteam_cph_by_R=mcph(T0,T,3.47,1.45e-3,0,0.121e5)
        CO_cph_by_R=mcph(T0,T,3.376,0.557e-3,0,-0.031e5)
        H2_cph_by_R=mcph(T0,T,3.249,0.422e-3,0,0.083e5)
        
        CH4_cps_by_R=mcps(T0,T,1.702,9.081e-3,-2.164e-6,0)
        H2Osteam_cps_by_R=mcps(T0,T,3.47,1.45e-3,0,0.121e5)
        CO_cps_by_R=mcps(T0,T,3.376,0.557e-3,0,-0.031e5)
        H2_cps_by_R=mcps(T0,T,3.249,0.422e-3,0,0.083e5)
        
        # idcph=(summation neu*cph_by_R)*(T-T0)
        idcph=(3*H2_cph_by_R+CO_cph_by_R-CH4_cph_by_R-H2Osteam_cph_by_R)*(T-T0)
        # idcps=(summation neu*cps_by_R)*ln(T/T0)
        idcps=(3*H2_cps_by_R+CO_cps_by_R-CH4_cps_by_R-H2Osteam_cps_by_R)*math.log(T/T0)
        deltaGo_T_by_RT=(deltaGo_T0-deltaHo_T0)/(R*T0)+deltaHo_T0/(R*T)+idcph/T-idcps
        K=math.exp(-deltaGo_T_by_RT)
        
        solns=solver(P0,P,K)
        for i in range(len(solns)):
            if np.isreal(solns[i]) and solns[i]>0 and solns[i]<=1:
                epsilon=solns[i]
        epsilon=epsilon.real
        epsilon_array.append(epsilon)
    T_array=np.array(T_array)
    epsilon_array=np.array(epsilon_array)
    plt.plot(T_array,epsilon_array*100,'-',label='%d bar'%P)

plt.xlabel('Temperature (in K)')
plt.ylabel('Methane Conversion (in %)')
plt.title('Effect of Temperature and Pressure on Equilibrium Conversion when \nFeed Ratio of H2O/CH4 = 3 mol/mol')
plt.legend()
plt.show()


# In[2]:


T=1123.15 # K
P=1 # bar
H2Oliq_Hfo_T0=-285830 # J/mol

CH4_cph_by_R=mcph(T0,T,1.702,9.081e-3,-2.164e-6,0)
H2Osteam_cph_by_R=mcph(T0,T,3.47,1.45e-3,0,0.121e5)
CO_cph_by_R=mcph(T0,T,3.376,0.557e-3,0,-0.031e5)
H2_cph_by_R=mcph(T0,T,3.249,0.422e-3,0,0.083e5)
        
CH4_cps_by_R=mcps(T0,T,1.702,9.081e-3,-2.164e-6,0)
H2Osteam_cps_by_R=mcps(T0,T,3.47,1.45e-3,0,0.121e5)
CO_cps_by_R=mcps(T0,T,3.376,0.557e-3,0,-0.031e5)
H2_cps_by_R=mcps(T0,T,3.249,0.422e-3,0,0.083e5)
        
# idcph=(summation neu*cph_by_R)*(T-T0)
idcph=(3*H2_cph_by_R+CO_cph_by_R-CH4_cph_by_R-H2Osteam_cph_by_R)*(T-T0)
# idcps=(summation neu*cps_by_R)*ln(T/T0)
idcps=(3*H2_cps_by_R+CO_cps_by_R-CH4_cps_by_R-H2Osteam_cps_by_R)*math.log(T/T0)
deltaGo_T_by_RT=(deltaGo_T0-deltaHo_T0)/(R*T0)+deltaHo_T0/(R*T)+idcph/T-idcps
K=math.exp(-deltaGo_T_by_RT)
        
solns=solver(P0,P,K)
for i in range(len(solns)):
    if np.isreal(solns[i]) and solns[i]>0 and solns[i]<=1:
        epsilon=solns[i]
epsilon=epsilon.real  

# Enthalpy Change for feed containing 1 mol CH4 and 3 mol H2O:
deltaH=3*epsilon*(H2_Hfo_T0+R*H2_cph_by_R*(T-T0))+epsilon*(CO_Hfo_T0+R*CO_cph_by_R*(T-T0))+(1-epsilon)*(CH4_Hfo_T0+R*CH4_cph_by_R*(T-T0))+(3-epsilon)*(H2Osteam_Hfo_T0+R*H2Osteam_cph_by_R*(T-T0))-1*CH4_Hfo_T0-3*H2Oliq_Hfo_T0
print('Enthalpy change for Steam Methane Reforming with feed containing 1 mol CH4 and 3 mol H2O at T= %f K and P= %d bar is \n%f J'%(T,P,deltaH))
print('Equilibrium conversion of methane at the above conditions is %f percent'%(100*epsilon))
print('Amount of Hydrogen produced at the above conditions is %f mol'%(3*epsilon))


# In[ ]:




