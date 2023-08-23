# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 12:19:33 2023

@author: Lara Moral
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Cálculo peso molecular

def peso_molecular(X, Y):
    Z = 1-X-Y
    mu = 1/(2*X+(3/4)*Y+(1/2)*Z)
    return mu


# Cálculo radio de cada capa

def radio_capa(R_ini, N):
    h = -R_ini/N
    r = []
    r.append(R_ini)
    for j in range(1, N+1):
        r.append(r[0]+h*j)
    return r


#INTEGRACIÓN DESDE LA SUPERFICIE

# Tres primeras capas del modelo

# Determinación parámetros primeras tres capas
def primeras_capas(X, Y, R_tot, L_tot, M_tot, R_ini, N,n,r):
    mu = peso_molecular(X, Y)
    Z = 1-X-Y
    # Determinamos las temperaturas
    T = []
    A1 = 1.9022*mu*M_tot
    for i in range(n):
        T.append(A1*(1/r[i]-1/R_tot))
    # Determinamos las presiones
    P = []
    A2 = 10.645*np.sqrt(M_tot/(mu*Z*(1+X)*L_tot))
    for i in range(n):
        P.append(A2*(T[i])**(4.25))
    # Ajustamos el número de decimales obtenidos
    r = ['%.5f' % elem for elem in r]
    T = ['%.7f' % elem for elem in T]
    P = ['%.7f' % elem for elem in P]

    L = L_tot*np.ones(n)
    M = M_tot*np.ones(n)
    E=[]
    for i in range(n):
     E.append("--")

    T = [float(t) for t in T]
    P = [float(t) for t in P]
    L = [float(t) for t in L]
    M = [float(t) for t in M]
    return(T, P, L, M, E)




# ENVOLTURA RADIATIVA



# Ritmo de generación de energía
#Se calcula el ritmo de generación de energía tanto para la
#cadena protón-protón como para el ciclo CN. Se escoge el ciclo
#que genere más energía

def ritmo_generacion_energia(X, Y, M, r, T):
    E = []
    t = T*10
    if t < 4:
        E.append("--")
        return (0, 0, 0, 0, E,0)
    if 4 <= t < 12:
        E.append("PP")
        X1 = X2 = X
        if t < 6:
            e1 = 10**(-6.84)
            v = 6
            
        if 6 <= t < 9.5:
            e1 = 10**(-6.04)
            v = 5
        if 9.5 <= t < 12:
            e1 = 10**(-5.56)
            v = 4.5
        e=e1*(X**2)*(t**v)*M/((4/3)*np.pi*(r**3))
    if 12 <= t < 16:
        ePP = (10**-5.02)*(X**2)*M*(t**4)/((4/3)*np.pi*(r**3))
        eCN = (10**-22.2)*X*(1/3)*(1-X-Y)*M*(t**20)/((4/3)*np.pi*(r**3))
        if abs(ePP)-abs(eCN) > 0:
            X1 = X2 = X
            e1 = 10**(-5.02)
            v = 4
            E.append("PP")
            e=ePP
        else:
            E.append("CN")
            X1 = X
            X2 = (1/3)*(1-X-Y)
            e1 = 10**(-22.2)
            v = 20
            e=eCN
    if 16 <= t < 16.5:
        ePP = (10**-5.02)*(X**2)*M*(t**4)/((4/3)*np.pi*(r**3))
        eCN = (10**-19.8)*X*(1/3)*(1-X-Y)*M*(t**18)/((4/3)*np.pi*(r**3))
        if abs(ePP)-abs(eCN) > 0:
            X1 = X2 = X
            e1 = 10**(-5.02)
            v = 4
            E.append("PP")
            e=ePP
        else:
            E.append("CN")
            X1 = X
            X2 = (1/3)*(1-X-Y)
            e1 = 10**(-19.8)
            v = 18
            e=eCN
    if 16.5 <= t < 22.5:
        ePP = (10**-4.40)*(X**2)*M*(t**3.5)/((4/3)*np.pi*(r**3))
        eCN = (10**-19.8)*X*(1/3)*(1-X-Y)*M*(t**18)/((4/3)*np.pi*(r**3))
        if abs(ePP)-abs(eCN) > 0:
            X1 = X2 = X
            e1 = 10**(-4.40)
            v = 3.5
            E.append("PP")
            e=ePP
        else:
            E.append("CN")
            X1 = X
            X2 = (1/3)*(1-X-Y)
            e1 = 10**(-19.8)
            v = 18
            e=eCN
    if 22.5 <= t < 24:
        ePP = (10**-4.40)*(X**2)*M*(t**3.5)/((4/3)*np.pi*(r**3))
        eCN = (10**-17.1)*X*(1/3)*(1-X-Y)*M*(t**16)/((4/3)*np.pi*(r**3))
        if abs(ePP)-abs(eCN) > 0:
            X1 = X2 = X
            e1 = 10**(-4.40)
            v = 3.5
            E.append("PP")
            e=ePP
        else:
            E.append("CN")
            X1 = X
            X2 = (1/3)*(1-X-Y)
            e1 = 10**(-17.1)
            v = 16
            e=eCN
    if 24 <= t < 50:
        E.append("CN")
        X1 = X
        X2 = (1/3)*(1-X-Y)
        if 24 <= t < 27.5:
            e1 = 10**(-17.1)
            v = 16
        if 27.5 <= t < 36:
            e1 = 10**(-15.6)
            v = 15
        if 36 <= t < 50:
            e1 = 10**(-12.5)
            v = 13
        e=e1*X*(1/3)*(1-X-Y)*(t**v)*M/((4/3)*np.pi*(r**3))
    if t >= 50:
        E.append("--")
        return (0, 0, 0, 0, E,0)
    X1 = float(X1)
    X2 = float(X2)
    e1 = float(e1)
    v = float(v)
    return(X1, X2, e1, v, E,e)


# Cálculo primeras derivadas f_i  para los diferentes parámetros 

# factores f_i para P
def f_i_P( P, T, r, L, X, Y, M):
    mu = peso_molecular(X, Y)
    C_p = 8.084*mu
    f_i_P = []
    for i in range(len(P)):
        f_i_P.append(-C_p*P[i]*M[i]/(T[i]*((r[i])**2)))
    return f_i_P

# factores f_i para T

#envoltura radiativa
def f_i_T( P, T, r, L, X, Y, M):
    mu = peso_molecular(X, Y)
    f_i_T = []
    for i in range(len(T)):
      C_t = 0.01679*(1-X-Y)*(1+X)*(mu**2)
      f_i_T.append(-C_t*((P[i])**2)*L[i]/((T[i]**8.5)*(r[i]**2)))
    return f_i_T

#núcleo convectivo
def f_i_T_nucleo( P, T, r, L, X, Y, M,f_T):
    mu = peso_molecular(X, Y)
    f=f_T
    for i in range(len(f_T),len(T)):
         C_t = -3.234*mu
         f.append(C_t*M[i]/(r[i]**2))
    return f


# factores f_i para M
def f_i_M(P, T, r, L, X, Y, M):
    mu = peso_molecular(X, Y)
    C_m = 0.01523*mu
    f_M=[]
    for i in range(len(f_M),len(M)):
      f_M.append(C_m*P[i]*((r[i])**2)/T[i])
    return f_M


# factores f_i para L
def f_i_L(P, T, r, L, X, Y, M):
    mu = peso_molecular(X, Y)
    f_L=[]
    for i in range(len(L)):
        if r[i]==0:
            f_L.append(0)
        else:
         X1 = ritmo_generacion_energia(X, Y, M[i], r[i], T[i])[0]
         X2 = ritmo_generacion_energia(X, Y, M[i], r[i], T[i])[1]
         e1 = ritmo_generacion_energia(X, Y, M[i], r[i], T[i])[2]
         v = ritmo_generacion_energia(X, Y, M[i], r[i], T[i])[3]
         f_L.append(0.01845*e1*X1*X2*(10**v)*(mu**2)*(P[i]**2)*(T[i]**(v-2))*(r[i]**2))
        
    return f_L


# Estimación de la presión

def P_est(i, P, T, r, L, X, Y, M, h):
    f = f_i_P(P, T, r, L, X, Y, M)
    Pest = P[i]+(23/12)*h*f[i]-(4/3)*h*f[i-1]+(5/12)*h*f[i-2]
    return Pest


# Estimación de la temperatura

#Envoltura radiativa
def T_est(i, P, T, r, L, X, Y, M, h):
    f = f_i_T(P, T, r, L, X, Y, M)
    Test = T[i]+(3/2)*h*f[i]-(1/2)*h*f[i-1]
    return Test

#Nucleo convectivo
def T_est_nucleo(i,P, T, r, L, X, Y, M, h,f_T):
    f = f_i_T_nucleo(P, T, r, L, X, Y, M,f_T)
    Test = T[i]+(3/2)*h*f[i]-(1/2)*h*f[i-1]
    return Test


# Cálculo de Mcal
def M_cal(i, P, T, r, L, X, Y, M, h, Pest, Test):
    mu = peso_molecular(X, Y)
    f = f_i_M( P, T, r, L, X, Y, M)
    dM = 0.01523*mu*(Pest/Test)*(r[i+1]**2)
    Mcal = M[i]+(1/2)*h*dM+(1/2)*h*(f[i])
    return Mcal

# Cálculo de Pcal
def P_cal(i, P, T, r, L, X, Y, M, h, Pest, Test, Mcal):
    mu = peso_molecular(X, Y)
    f = f_i_P(P, T, r, L, X, Y, M)
    dP = -8.084*mu*(Pest/Test)*(Mcal/(r[i+1]**2))
    Pcal = P[i]+(1/2)*h*dP+(1/2)*h*(f[i])
    return Pcal

# Cálculo de Lcal

#envoltura radiativa
def L_cal(i, P, T, r, L, X, Y, M, h, Pcal, Test, Mcal):
    mu = peso_molecular(X, Y)
    f = f_i_L( P, T, r, L, X, Y, M)
    X1 = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Test)[0]
    X2 = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Test)[1]
    e1 = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Test)[2]
    v = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Test)[3]
    dL = 0.01845*e1*X1*X2*(10**v)*(mu**2)*(Pcal**2)*(Test**(v-2))*(r[i+1]**2)
    Lcal = L[i]+h*(5/12)*dL+h*((2/3)*f[i]-(1/12)*f[i-1])
    return Lcal

#nucleo convectivo
def L_cal_nucleo(i,P, T, r, L, X, Y, M, h, Pcal, Tcal, Mcal):
    mu = peso_molecular(X, Y)
    f = f_i_L(P, T, r, L, X, Y, M)
    X1 = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[0]
    X2 = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[1]
    e1 = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[2]
    v = ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[3]
    dL = 0.01845*e1*X1*X2*(10**v)*(mu**2)*(Pcal**2)*(Tcal**(v-2))*(r[i+1]**2)
    Lcal = L[i]+h*(5/12)*dL+h*((2/3)*f[i]-(1/12)*f[i-1])
    return Lcal


# Cálculo de Tcal

#envoltura radiativa
def T_cal(i,P, T, r, L, X, Y, M, h, Pcal, Test, Lcal, Mcal):
    mu = peso_molecular(X, Y)
    f = f_i_T( P, T, r, L, X, Y, M)
    dT = -0.01679*(1-X-Y)*(1+X)*(mu**2) *((Pcal**2)/(Test**8.5))*(Lcal/((r[i+1])**2))
    Tcal = T[i]+h*(1/2)*(dT)+(1/2)*h*(f[i])
    Tcal=Tcal
    return (Tcal)

#nucleo convectivo
def T_cal_nucleo(i, P, T, r, L, X, Y, M, h, Mcal,f_T):
   mu = peso_molecular(X, Y)
   f = f_i_T_nucleo( P, T, r, L, X, Y, M,f_T)
   dT = -3.234*mu*Mcal/(r[i+1]**2)
   Tcal = T[i]+h*(1/2)*(dT)+(1/2)*h*(f[i])
   return Tcal


# parámetro n+1
#Calculamos el parámetro n+1 para determinar cuando deja de ser
#válida la hipótesis inicial de transporte radiativo

def parametro_n(i, Pcal, Test, r, Lcal, X, Y, Mcal, h, Pest, Tcal):
    mu = peso_molecular(X, Y)
    dT = -0.01679*(1-X-Y)*(1+X)*(mu**2) *((Pcal**2)/(Test**8.5))*(Lcal/((r[i+1])**2))
    dP = -8.084*mu*(Pest/Test)*(Mcal/(r[i+1]**2))
    # m=n+1
    m = (Tcal/Pcal)*(dP/dT)
    return m



#Cálculo de la P,T,L y M para las capas de la envoltura radiativa
def envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N):
    i = 2  # número de la última capa calculada

# Utilizamos las variables lógicas loop1, loop2 y loop3
# para controlar cada uno de los bucle necesarios para
# programar el algoritmo de la fase A.1. (envoltura radiativa)

    # radios de las capas
    r = radio_capa(R_ini, N)
    # Parámetros de las tres primeras capas
    
    h = -R_ini/N
    T = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,3,r)[0]
    P = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,3,r)[1]
    L = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,3,r)[2]
    M = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,3,r)[3]
    E =primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,3,r)[4]
    n = [0,0,0]

    
    loop1 = True
    while loop1:
        Pest = P_est(i, P, T, r, L, X, Y, M, h)
        Test = T_est(i, P, T, r, L, X, Y, M, h)

        loop2 = True
        while loop2:
            loop3 = True
            while loop3:
                Mcal = M_cal(i, P, T, r, L, X, Y, M, h, Pest, Test)
                Pcal = P_cal(i, P, T, r, L, X, Y, M, h, Pest, Test, Mcal)
                E_relativo_maximo = abs(Pcal-Pest)/Pcal
                if E_relativo_maximo <= 0.0001:
                    decision = "True"
                else:
                    decision = "False"
                if decision == "True":
                    loop3 = False
                else:
                    Pest = Pcal
            Lcal = L_cal(i, P, T, r, L, X, Y, M, h, Pcal, Test, Mcal)
            
            Tcal = T_cal(i, P, T, r, L, X, Y, M, h, Pcal, Test, Lcal,Mcal)
            e_relativo_maximo = abs(Tcal-Test)/Tcal
            if e_relativo_maximo <= 0.0001:
                decision = "True"
            else:
                decision = "False"
            if decision == "True":
                loop2 = False
            else:
                Test = Tcal
        m = parametro_n(i, Pcal, Test, r, Lcal, X, Y, Mcal, h, Pest, Tcal)
        if m <= 2.5:
            decision = "True"
            n.append(m)
        else:
            decision = "False"
        if decision == "True":
            K=Pcal/(Tcal**2.5)    #constante K' del polítropo
            loop1 = False
        else:
            P.append(Pcal)
            T.append(Tcal)
            L.append(Lcal)
            M.append(Mcal)
            n.append(m)
            E.append(ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[4])
            i += 1   # pasamos a la siguiente capa
            
            T = [float(t) for t in T]
            P = [float(t) for t in P]
            L = [float(t) for t in L]
            M = [float(t) for t in M]
            
    return (P, T, L, M, n,i,K,E)




#NÚCLEO CONVECTIVO
#Cálculo de la P,T,L y M para las capas del núcleo convectivo

def nucleo_convectivo( X, Y, R_tot, L_tot, M_tot, R_ini, N):
 #Parámetros envoltura radiativa
 i=envoltura_radiativa(X, Y, R_tot, L_tot, M_tot, R_ini, N)[5]
 P=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[0]
 T=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[1]
 L=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[2]
 M=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[3]
 n=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[4]
 K=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[6]
 E=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[7]
 r = radio_capa(R_ini, N)
 h = -R_ini/N
 f_T=f_i_T(P, T, r, L, X, Y, M)
 
 loop1 = True
 while loop1:
    Test = T_est_nucleo(i, P, T, r, L, X, Y, M, h,f_T)
    loop2 = True
    while loop2:
        Pest=K*(Test**2.5)
        Mcal=M_cal(i, P, T, r, L, X, Y, M, h, Pest, Test)
        if r[i+1]==0:
            Tcal=Test
            loop2=False
        else:
            Tcal=T_cal_nucleo(i, P, T, r, L, X, Y, M, h, Mcal,f_T)
        e_relativo_maximo = abs(Tcal-Test)/Tcal
        if e_relativo_maximo <= 0.0001:
            decision = "True"
        else:
            decision = "False"
        if decision == "True":
            loop2 = False
        else:
            Test=Tcal
    Pcal=K*(Tcal**2.5)
    Lcal=L_cal_nucleo(i, P, T, r, L, X, Y, M, h, Pcal, Tcal, Mcal)
    if r[i+1]<=0:
     decision = "True"
    else:
        decision="False"
    if decision == "True":
        loop1 = False
        P.append(Pcal)
        T.append(Tcal)
        L.append(Lcal)
        M.append(Mcal)
        E.append(ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[4])
        i += 1
    else:
        loop1=False #esta línea detiene el cálculo en la primera capa
                    #que es la única de la que se hará uso más adelante
        P.append(Pcal)
        T.append(Tcal)
        L.append(Lcal)
        M.append(Mcal)
        n.append(0)
        E.append(ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[4])
        i += 1   # pasamos a la siguiente capa
 return (P, T, L, M, E,n)





#Valores en la frontera entre la envoltura radiativa y el núcleo 
#convectivo para la integración desde la superficie

def valores_en_frontera_superficie( X, Y, R_tot, L_tot, M_tot, R_ini, N):
    i=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[5]
    P=nucleo_convectivo( X, Y, R_tot, L_tot, M_tot, R_ini, N)[0]
    T=nucleo_convectivo( X, Y, R_tot, L_tot, M_tot, R_ini, N)[1]
    L=nucleo_convectivo( X, Y, R_tot, L_tot, M_tot, R_ini, N)[2]
    M=nucleo_convectivo( X, Y, R_tot, L_tot, M_tot, R_ini, N)[3]
    n=nucleo_convectivo( X, Y, R_tot, L_tot, M_tot, R_ini, N)[5]
    r=radio_capa(R_ini,N)
    
    P_1=[P[i],P[i+1]]
    T_1=[T[i],T[i+1]]
    L_1=[L[i],L[i+1]]
    M_1=[M[i],M[i+1]]
    n_1=[n[i],n[i+1]]
    r_1=[r[i],r[i+1]]
    
    #Parámetros para n+1=2.5 obtenidos por ajuste lineal

    #y=mx+c
    
    #radio
    [m,c]=np.polyfit(r_1,n_1,deg=1)
    r_frontera=(2.5-c)/m
    
    #presion
    [m,c]=np.polyfit(P_1,n_1,deg=1)
    P_frontera=(2.5-c)/m
    
    #temperatura
    [m,c]=np.polyfit(T_1,n_1,deg=1)
    T_frontera=(2.5-c)/m
    
    #luminosidad
    [m,c]=np.polyfit(L_1,n_1,deg=1)
    L_frontera=(2.5-c)/m
    
    #Masa
    [m,c]=np.polyfit(M_1,n_1,deg=1)
    M_frontera=(2.5-c)/m
    
     
    
    return (r_frontera,P_frontera,T_frontera,L_frontera,M_frontera)





#INTEGRACIÓN DESDE EL CENTRO


# Primeras tres capas
#Cálculo de la P,T,L y M para las tres primeras capas 
#del núcleo convectivo integrando desde el centro

def primeras_tres_capas_centro(X, Y, Tc,R_ini,N,R_tot,L_tot,M_tot):
    mu = peso_molecular(X, Y)
    r =list(np.flip(radio_capa(R_ini, N)))[0:3]
    K=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[6]
    M=[0]
    L=[0]
    T=[Tc]
    P=[]
    
    for i in range(1,3):
     #Calculamos la masa
     M.append(0.005077*mu*K*Tc**(1.5)*(r[i]**3))

     #Calculamos la luminosidad
    for i in range(1,3):
     X1 = ritmo_generacion_energia(X, Y, M[i], r[i], Tc)[0]
     X2 = ritmo_generacion_energia(X, Y, M[i], r[i], Tc)[1]
     e1 = ritmo_generacion_energia(X, Y, M[i], r[i], Tc)[2]
     v = ritmo_generacion_energia(X, Y, M[i], r[i], Tc)[3]
     L.append(0.006150*e1*X1*X2*(10**v)*(mu**2)*(K**2)*(Tc**(3+v))*(r[i]**3))
     #Calculamos la temperatura
     T.append(Tc-0.008207*(mu**2)*K*(Tc**1.5)*(r[i]**2))
     #Calculamos la presión
    for i in range(0,3):
      P.append(K*(T[i]**2.5))

    # Ajustamos el número de decimales obtenidos
    r = ['%.5f' % elem for elem in r]
    T = ['%.7f' % elem for elem in T]
    P = ['%.7f' % elem for elem in P]

    E = ["--", "--", "--"]
    T = [float(t) for t in T]
    P = [float(t) for t in P]
    L = [float(t) for t in L]
    M = [float(t) for t in M]

    return(T, P, L, M, E)






#Capas posteriores
#Cálculo de la P,T,L y M para el resto de capas del núcleo convectivo

def capas_posteriores(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N):
 #Parámetros de las tres primeras capas
 i=2
 P=primeras_tres_capas_centro(X, Y, Tc,R_ini,N,R_tot,L_tot,M_tot)[1]
 T=primeras_tres_capas_centro(X, Y, Tc,R_ini,N,R_tot,L_tot,M_tot)[0]
 L=primeras_tres_capas_centro(X, Y, Tc,R_ini,N,R_tot,L_tot,M_tot)[2]
 M=primeras_tres_capas_centro(X, Y, Tc,R_ini,N,R_tot,L_tot,M_tot)[3]
 K=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[6]
 E=["--","--","--"]
 j=N-envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[5]
 r =list(np.flip( radio_capa(R_ini, N)))
 h = R_ini/N


 loop1 = True
 while loop1:

    Test = T_est_nucleo(i,P, T, r, L, X, Y, M, h,[0])
    
    loop2 = True
    while loop2:
        Pest=K*(Test**2.5)
        Mcal=M_cal(i, P, T, r, L, X, Y, M, h, Pest, Test)
        if i==j:
            Tcal=Test
            loop2=False
        else:
            Tcal=T_cal_nucleo(i, P, T, r, L, X, Y, M, h, Mcal,[0])   
        
        e_relativo_maximo = abs(Tcal-Test)/Tcal
        if e_relativo_maximo <= 0.0001:
            decision = "True"
        else:
            decision = "False"
        if decision == "True":
            loop2 = False
        else:
            Test=Tcal
    Pcal=K*(Tcal**2.5)
    Lcal=L_cal_nucleo(i, P, T, r, L, X, Y, M, h, Pcal, Tcal, Mcal)
    if i==j:
     decision = "True"
    else:
        decision="False"
    if decision == "True":
        loop1 = False
    else:
        P.append(Pcal)
        T.append(Tcal)
        L.append(Lcal)
        M.append(Mcal)
        E.append(ritmo_generacion_energia(X, Y, Mcal, r[i+1], Tcal)[4])
        i += 1   # pasamos a la siguiente capa
        
 return (P, T, L, M, E,i)


#Calculamos los valores en la frontera con la envoltura radiativa
#para comparar con los obtenidos anteriormente integrando desde la
#superficie

def valores_en_frontera_centro( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N):
    
 i=capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[5]
 P=capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[0]
 T=capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[1]
 L=capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[2]
 M=capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[3]
 r=np.flip(radio_capa(R_ini,N))
 r_frontera=valores_en_frontera_superficie( X, Y, R_tot, L_tot, M_tot, R_ini, N)[0]
 
 P_1=[P[i-1],P[i]]
 T_1=[T[i-1],T[i]]
 L_1=[L[i-1],L[i]]
 M_1=[M[i-1],M[i]]
 r_1=[r[i-1],r[i]]
 

 #y=mx+c
 
 
 #presion
 [m,c]=np.polyfit(P_1,r_1,deg=1)
 P_frontera=(r_frontera-c)/m
 
 #temperatura
 [m,c]=np.polyfit(T_1,r_1,deg=1)
 T_frontera=(r_frontera-c)/m
 
 #luminosidad
 [m,c]=np.polyfit(L_1,r_1,deg=1)
 L_frontera=(r_frontera-c)/m
 
 #Masa
 [m,c]=np.polyfit(M_1,r_1,deg=1)
 M_frontera=(r_frontera-c)/m
 
  
 
 return (P_frontera,T_frontera,L_frontera,M_frontera)




#Comparación de las soluciones obtenidas en la frontera para 
#la integración desde la superficie y la integración desde el
#centro

def error_relativo_total( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N):
    
    #valores en la frontera para integración desde la superficie
    P_d=valores_en_frontera_superficie( X, Y, R_tot, L_tot, M_tot, R_ini, N)[1]
    T_d=valores_en_frontera_superficie( X, Y, R_tot, L_tot, M_tot, R_ini, N)[2]
    L_d=valores_en_frontera_superficie( X, Y, R_tot, L_tot, M_tot, R_ini, N)[3]
    M_d=valores_en_frontera_superficie( X, Y, R_tot, L_tot, M_tot, R_ini, N)[4]
    
    #valores en la frontera para integración desde el centro
    
    P_u=valores_en_frontera_centro( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[0]
    T_u=valores_en_frontera_centro( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[1]
    L_u=valores_en_frontera_centro( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[2]
    M_u=valores_en_frontera_centro( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[3]
    
    #calculamos el error relativo
    e_rel_P=abs(P_u-P_d)*100/P_d
    e_rel_T=abs(T_u-T_d)*100/T_d
    e_rel_L=abs(L_u-L_d)*100/L_d
    e_rel_M=abs(M_u-M_d)*100/M_d
    
    #error relativo total
    e_rel_total=np.sqrt((e_rel_P**2)+(e_rel_T**2)+(e_rel_L**2)+(e_rel_M**2))
    return(e_rel_P,e_rel_T,e_rel_L,e_rel_M,e_rel_total)



#Buscamos el valor de Tc que minimiza el error relativo total
# T_ini y T_final son los extremos del itervalo de temperaturas
#n es el paso entre las temperaturas del intervalo

def min_error_total( X, Y, R_tot, L_tot, M_tot, R_ini, N,n,T_ini,T_final):
    Tc=T_ini
    T=[]
    error=[]
    while Tc<=T_final:
        T.append(Tc)
        error.append(error_relativo_total( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[4])
        Tc=Tc+n
    if Tc>T_final:
        T.append(T_final)
        error.append(error_relativo_total( X, Y,T_final, R_tot, L_tot, M_tot, R_ini, N)[4])
 
    e_min=min(error)
    i=error.index(e_min)
    T_min=T[i]
    e_min=float(e_min)
    T_min=float(T_min)
    return(T_min,e_min)
    

        
        
#MODELO COMPLETO
#Modelo completo obtenido uniendo las soluciones de la integración
#desde la superficie y la integración desde el centro

def modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N):
    fase=[]
    R=R_ini
    r1=[]
    r2 = radio_capa(R_ini, N)
    while R<=R_tot:
        r1.append(R)
        R=R+R_ini/N
    r1=np.flip(r1)
    r1=r1[0:len(r1)-1]
    s=len(r1)
    
    #Tres primeras capas
    T1 = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,s,r1)[0]
    P1 = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,s,r1)[1]
    L1 = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,s,r1)[2]
    M1 = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,s,r1)[3]
    E1 = primeras_capas( X, Y, R_tot, L_tot, M_tot, R_ini, N,s,r1)[4]
    n1=np.zeros(s)
    E1=list(pd.core.common.flatten(E1))
    for i in range(s):
        fase.append("^^^^^^")
    
    #Envoltura radiativa
    P2=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[0]
    T2=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[1]
    L2=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[2]
    M2=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[3]
    n2=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[4]
    E2=envoltura_radiativa( X, Y, R_tot, L_tot, M_tot, R_ini, N)[7]
    E2=list(pd.core.common.flatten(E2))
    for i in range(3):
        fase.append("INICIO")
    for i in range(len(P2)-3):
        fase.append("RADIAT")
    
    #Núcleo convectivo
    P3=np.flip(capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[0])
    T3=np.flip(capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[1])
    L3=np.flip(capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[2])
    M3=np.flip(capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[3])
    n3=np.zeros(len(P3))
    E3=np.flip(capas_posteriores( X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[4])
    E3=list(pd.core.common.flatten(E3))
    for i in range(len(P3)-4):
        fase.append("CONVEC")
    for i in range(3):
        fase.append("CENTRO")
        
    #Unimos todas las capas
    E=np.concatenate((E1,E2))
    E=np.concatenate((E,E3[1:]))
    P=np.concatenate((P1,P2))
    P=np.concatenate((P,P3[1:]))
    T=np.concatenate((T1,T2))
    T=np.concatenate((T,T3[1:]))
    L=np.concatenate((L1,L2))
    L=np.concatenate((L,L3[1:]))
    M=np.concatenate((M1,M2))
    M=np.concatenate((M,M3[1:]))
    n=np.concatenate((n1,n2[0:len(P2)]))
    n=np.concatenate((n,n3[1:]))
    r=np.concatenate((r1,r2))
    fase=np.array(fase)
    
    
    T = ['%.7f' % elem for elem in T]
    P = ['%.7f' % elem for elem in P]
    

    #Representación en forma de tabla
    Datos = {'E':E,'fase':fase,'r': r, 'P': P, 'T': T, 'L': L, 'M': M,'n+1':n}
    datos = pd.DataFrame(data=Datos)
    i=np.arange(-s,N+1,1)
    datos=datos.set_index(i)
    
    tabla=open("Modelo_Interior.txt","w")
    datos=pd.DataFrame.to_string(datos)
    tabla.write(datos)
    
    T = [float(t) for t in T]
    P = [float(t) for t in P]
    L = [float(t) for t in L]
    M = [float(t) for t in M]
    r= [float(t) for t in r]

    
    return(T, P, L, M, E,r)

    

#Buscamos los valores para el radio y la luminosidad total
#que proporcionan el modelo con el menor error relativo total

def minimizar_error_relativo(X,Y,M_tot,R_tot,L_tot,delta_R,delta_L,T_ini,T_final,n,N):
    x=np.linspace(R_tot-2*delta_R,R_tot+2*delta_R,5)
    y=np.linspace(L_tot-2*delta_L,L_tot+2*delta_L,5)
    R_ini=[0.9*elem for elem in x]
  
    errores=np.zeros((5,5))
    Tc=np.zeros((5,5))
    #Valores obtenidos para el error relativo total
    
    for i in range(5):
        for j in range(5):
          R_tot=x[i]
          L_tot=y[j]
          r_ini=R_ini[i]
          errores[i,j]=min_error_total( X, Y, R_tot, L_tot, M_tot, r_ini, N,n,T_ini,T_final)[1]
          Tc[i,j]=min_error_total( X, Y, R_tot, L_tot, M_tot, r_ini, N,n,T_ini,T_final)[0]
    
    plt.subplots(figsize=(7,5))
    XX,YY=np.meshgrid(x,y)
    U=errores.transpose()
    plt.title("Minimización error relativo total ",fontsize=14)
    plt.xlabel("$R_{total}$",fontsize=20, labelpad=15)
    plt.ylabel("$L_{total}$",fontsize=20,labelpad=15)
    plt.xticks(fontsize=14)
    plt.xticks(rotation=45)
    plt.yticks(fontsize=14)
    plt.ticklabel_format(style='plain')
    plt.pcolormesh(XX, YY, U, cmap='viridis', shading='nearest') 
    plt.gca().invert_yaxis()          
    cbar = plt.colorbar()
    cbar.set_label('Error relativo total', rotation=270,fontsize=13,labelpad=15)
    plt.margins(0.5)
    cbar.ax.tick_params(labelsize=14)
    print(errores)
    print(Tc)
    return (errores,Tc)


#gráficas de los distintos parametros
def grafica_modelo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N):
    
    T = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[0]
    P = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[1]
    L = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[2]
    M = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[3]
    r = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[5]
    
    T_max=max(T)
    P_max=max(P)
    L_max=max(L)
    M_max=max(M)
    r_max=max(r)
    
    T=[i/T_max for i in T]
    P=[i/P_max for i in P]
    L=[i/L_max for i in L]
    M=[i/M_max for i in M]
    r=[i/r_max for i in r]
    
    plt.subplots(figsize=(7,5))
    plt.title("Parámetros físicos del modelo ",fontsize=13)
    plt.xlabel("$\mathrm{R/R_{tot}}$",fontsize=20)
    plt.ylabel(r"$\mathrm{M/M_{tot}}$,  $\mathrm{L/L_{tot}}$,  $\mathrm{P/P_{c}}$,  $\mathrm{T/T_{c}}$",fontsize=13)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.plot(r,T,'g',label='T(r)')
    plt.legend(loc='upper right')
    plt.plot(r,P,'r',label='P(r)')
    plt.legend(loc='upper right')
    plt.plot(r,L,'b',label='L(r)')
    plt.legend(loc='upper right')
    plt.plot(r,M,label='M(r)')
    plt.legend(loc='upper right')
    plt.plot(0.144*np.ones(20),np.linspace(0,1,20),'--')
    plt.show()
 
#gráfica del ritmo de generacion de energia y luminosidad
def grafica_generacion_energia(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N):
    r = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[5]
    T = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[0]
    M = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[3]
    L = modelo_completo(X, Y,Tc, R_tot, L_tot, M_tot, R_ini, N)[2]
    r_max=max(r)
    L_max=max(L)
    M_max=max(M)
    
    e=np.zeros(len(T))
    for i in range(len(T)-3):
        e[i]=ritmo_generacion_energia(X, Y, M[i], r[i], T[i])[5]
 
    
    e_max=max(e)
    
    R=[i/r_max for i in r]
    Mmax=[i/M_max for i in M]
    m_max=[i/M_max for i in M][0:len(T)-3]
    r=[i/r_max for i in r][0:len(T)-3]
    e=[i/e_max for i in e][0:len(T)-3]
    L=[i/L_max for i in L]
    
    plt.subplots(figsize=(7,5))
    plt.title("Ritmo generación de energía y luminosidad ",fontsize=14)
    plt.xlabel("$\mathrm{M/M_{tot}}$",fontsize=20)
    plt.ylabel(r"$\mathrm{\varepsilon/\varepsilon_{max}}$,  $\mathrm{L/L_{tot}}$",fontsize=13)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=12)
    plt.plot(m_max,e,label=r"$\varepsilon (M)$")
    plt.legend(fontsize=13)
    plt.plot(Mmax,L,label='L(M)')
    plt.legend(fontsize=13)
    plt.plot(0.104*np.ones(20),np.linspace(0,1,20),'--')
    plt.show()
    
    
    
    
modelo_completo(0.85, 0.1,1.7852, 11.5665, 28.63, 5.1, 10.40985, 100)