# -*- coding: utf-8 -*-
"""
Created on Fri May 13 08:39:29 2022

@author: Andrés Coto y James Ruiz
"""

# Se inicializan las bibliotecas necesarias
import numpy as np
import matplotlib.pyplot as plt

# Simulación del modelo de Ising

# Se inicializan las variables necesarias
J = 1 # Comportamiento ferromagnético
kBT = 0.5 # Constante de Boltzmann por la temperatura
nEspines = 100
nPasos = 1000
B = 0.0005 # Se asume que el campo magnético es muy pequeño
kB = 1.380649*10**(-23) # Constante de Boltzmann

# Se realiza una codificación de la configuración de los espines
frios_arriba = 1
frios_abajo = 2
calientes = 3

# Se define una rutina que simule el modelo de Ising en 1D
def ModeloIsing1D(j, kbt, nespines, tipo, npasos):
    """
    La función realiza la simulación del modelo de Ising para n cantidad de espines
    unidimensional, para una n cantidad de pasos. Esta simulación se realiza a partir
    del algoritmo Metrópolis.
    Entradas:
        j : Comportamiento ferromagnético del material
        kbt : constante de Boltzmann por la tempertura del material
        nespines: cantidad de espines a simular
        tipo : tipo de configuración de los espines(fríos-arriba = 1, fríos-abajo = 2, calientes = 3)
        npasos : cantidad de pasos del tiempo a simular
    Salidas:
        La función retorna un array de tamaño nespines x npasos con las configuraciones 
        de los nespines para npasos.
    """
    # Se inicializa un arreglo de nespines con condiciones de contorno periódicas
    arregloEspines = np.zeros(nespines)
    arregloEspines[0] = arregloEspines[nespines-1]
    arregloEspines[nespines-1] = arregloEspines[0]
    
    # Se genera un array 2D donde se guarde la configuración de espines para cada paso
    array_config_espines = np.zeros((nespines, npasos))   
    
    # Según el tipo de configuración de los espines se generan en el array
    if tipo == 1:
        arregloEspines[:] = 1/2
    elif tipo == 2:
        arregloEspines[:] = -1/2
    elif tipo == 3:
        # Se genera un array con números aleatorios enteros con posibles valores de 0 o 1
        n_aleatorios = np.random.randint(2, size=nespines)
        for i in range(nespines):
            if n_aleatorios[i] == 0:
                arregloEspines[i] = -1/2
            else:
                arregloEspines[i] = 1/2
    else:
        print("Ingrese un tipo de configuración válida")
    
    # Se guarda la configuración del estado en el array 2D
    array_config_espines[:,0] = arregloEspines
    
    for n in range(1, npasos):
        
        # Se guarda el estado actual del arreglo de los espines
        arregloEspinesPrueba = array_config_espines[:,n-1].copy()
        
        #Se determina la energía del sistema de espines
        suma_espines = 0
        for i in range(nespines-1):
            suma_espines += array_config_espines[i,n-1] * array_config_espines[i+1,n-1]
        energiaSistema = -j * suma_espines
        
        # Se elige un espín aleatoriamente mediante una distribución uniforme
        posicion_aleatoria = np.random.randint(nespines)
        arregloEspinesPrueba[posicion_aleatoria] = - arregloEspinesPrueba[posicion_aleatoria]
        
        # Se calcula la nueva energía del sistema
        suma_espines = 0
        for i in range(nespines-1):
            suma_espines += arregloEspinesPrueba[i] * arregloEspinesPrueba[i+1]
        nueva_energiaSistema = -j * suma_espines
        
        # De calcula el cambio de energía
        deltaE = nueva_energiaSistema - energiaSistema
        
        # Se calcula la probabilidad de aceptación de la transición
        if nueva_energiaSistema <= energiaSistema:
            array_config_espines[:,n] = arregloEspinesPrueba
        else:
            probabilidad = np.exp(-deltaE/kbt)
            prob_aleatorio = np.random.rand()
            
            # Se acepta o se rechaza la transición
            if probabilidad >= prob_aleatorio:
                array_config_espines[:,n] = arregloEspinesPrueba
            else:
                array_config_espines[:,n] = array_config_espines[:,n-1]      
        
    return array_config_espines
            
#config_espines1 = ModeloIsing1D(J, kBT, nEspines, frios_arriba, nPasos)
#config_espines2 = ModeloIsing1D(J, kBT, nEspines, frios_abajo, nPasos)
#config_espines3 = ModeloIsing1D(J, kBT, nEspines, calientes, nPasos)


# Se genera el gráfico del comportamiento de los espines mediante una rutina
def GraficoComportamientoEspines(configespines):
    """
    La función se encarga de realizar la gráfica con el comportamiento de los espines
    a través del tiempo (pasos)
    Entradas:
        configespines : array con las configuraciones de los espines con el tiempo (pasos)
    Salidas:
        La función grafica las configuraciones de los espines
    """
    fig, ax = plt.subplots(dpi=120)
    ax.patch.set_facecolor("gray")
    for (y, x), w in np.ndenumerate(configespines):
        color = "white" if w == 1/2 else "black"
        size = 1
        rect = plt.Rectangle([x-size/2, y-size/2], size, size, facecolor=color, edgecolor=color)
        
        ax.add_patch(rect)
    ax.autoscale_view() 
    ax.set_title('Simulación del Modelo de Ising')
    ax.set_xlabel('pasos')
    ax.set_ylabel('Espines')
    plt.show()
    
#GraficoComportamientoEspines(config_espines1)
#GraficoComportamientoEspines(config_espines2)
#GraficoComportamientoEspines(config_espines3)


# Propiedades termodinámicas del modelo de Ising

# Se crea una rutina que permita medir el comportamiento de la energía interna, la magnetización y el calor específico
# del sistema de espines en el tiempo

# Se define la cantidad de simulaciones
nSimulaciones = 20

# Se va a usar un rango para el kbt de entre 0 y 5 con 10 puntos
array_kBT = np.linspace(0.00001, 5, 50)

def ModeloIsingPropiedades(j, array_kbt, nespines, tipo, npasos, nsimulaciones):
    """
    La función realiza la simulación del modelo de Ising para n cantidad de espines
    unidimensional, para una n cantidad de pasos. Esta simulación se realiza a partir
    del algoritmo Metrópolis.
    Entradas:
        j : Comportamiento ferromagnético del material
        array_kbt : array con los valores de la constante de Boltzmann por la tempertura del material
        nespines: cantidad de espines a simular
        tipo : tipo de configuración de los espines(fríos-arriba = 1, fríos-abajo = 2, calientes = 3)
        npasos : cantidad de pasos del tiempo a simular
        nsimulaciones : cantidad de simulaciones a realizar para determinar U, C y la magnetización
    Salidas:
        Un array con los valores 
    """
    # Se inicializan arrays donde se guadarán los valores para U, C y la magnetización
    arrayU = np.zeros(np.size(array_kbt))
    arrayC = np.zeros(np.size(array_kbt))
    arrayMag = np.zeros(np.size(array_kbt))
    
    for o in range(np.size(array_kbt)):
               
        # Se inicializan arrays donde se guarden U, C y la magnetización para cada simulación
        # Luego se determina el promedio para cada parámetro
        
        arrayUsim = np.zeros(nsimulaciones)
        arrayMagsim = np.zeros(nsimulaciones)
        
        for p in range(nsimulaciones):
            # Se inicializa un arreglo de nespines con condiciones de contorno periódicas
            arregloEspines = np.zeros(nespines)
            arregloEspines[0] = arregloEspines[nespines-1]
            
            # Se inicializa una variable donde se sume las energías para determinar la energía interna y el calor específico
            sumaenergias = 0
            
            # Se genera un array 2D donde se guarde la configuración de espines para cada paso
            array_config_espines = np.zeros((nespines, npasos))   
            
            # Según el tipo de configuración de los espines se generan en el array
            if tipo == 1:
                arregloEspines[:] = 1/2
            elif tipo == 2:
                arregloEspines[:] = -1/2
            elif tipo == 3:
                # Se genera un array con números aleatorios enteros con posibles valores de 0 o 1
                n_aleatorios = np.random.randint(2, size=nespines)
                for i in range(nespines):
                    if n_aleatorios[i] == 0:
                        arregloEspines[i] = -1/2
                    else:
                        arregloEspines[i] = 1/2
            else:
                print("Ingrese un tipo de configuración válida")
            
            
            # Se guarda la configuración del estado en el array 2D
            array_config_espines[:,0] = arregloEspines
            
            # Se crea un array para determinar la magnetización
            magnetizacion = np.zeros(npasos)
            magnetizacion[0] = np.abs(np.sum(array_config_espines[:,0]))
            
            for n in range(1, npasos):
                
                # Se guarda el estado actual del arreglo de los espines
                arregloEspinesPrueba = array_config_espines[:,n-1].copy()

                #Se determina la energía del sistema de espines
                suma_espines = 0
                for i in range(nespines-1):
                    suma_espines += array_config_espines[i,n-1] * array_config_espines[i+1,n-1]
                energiaSistema = -j * suma_espines
                
                # Se elige un espín aleatoriamente mediante una distribución uniforme
                posicion_aleatoria = np.random.randint(nespines)
                arregloEspinesPrueba[posicion_aleatoria] = - arregloEspinesPrueba[posicion_aleatoria]
                
                # Se calcula la nueva energía del sistema
                suma_espines = 0
                for i in range(nespines-1):
                    suma_espines += arregloEspinesPrueba[i] * arregloEspinesPrueba[i+1]
                nueva_energiaSistema = -j * suma_espines
                
                # De calcula el cambio de energía
                deltaE = nueva_energiaSistema - energiaSistema
                
                # Se calcula la probabilidad de aceptación de la transición
                if nueva_energiaSistema <= energiaSistema:
                    array_config_espines[:,n] = arregloEspinesPrueba
                else:
                    probabilidad = np.exp(-deltaE/array_kbt[o])
                    prob_aleatorio = np.random.rand()
                    
                    # Se acepta o se rechaza la transición
                    if probabilidad >= prob_aleatorio:
                        array_config_espines[:,n] = arregloEspinesPrueba
                    else:
                        array_config_espines[:,n] = array_config_espines[:,n-1]
                
                # Se calcula la suma de la energía del sistema en el estado de "equilibrio"
                if n >= npasos/5:  # Se asume que el equilibrio se alcanza a partir de npasos/5
                    sumaenergias += energiaSistema
                    magnetizacion[int(n-npasos/5)] = np.abs(np.sum(array_config_espines[:,int(n-npasos/5)]))
            
            # Se determina el promedio de las energías para la simulación
            arrayUsim[p] = sumaenergias/npasos
            
            # Se determina la magnetización de la última configuración de espines
            arrayMagsim[p] = np.mean(magnetizacion)
            
        # Se determina U2 para el cálculo del calor específico
        u2 = 0
        u2 = (1/nsimulaciones)*np.sum(arrayUsim**2)
        
        # Se determina el promedio de la energía interna para las nsimulaciones
        arrayU[o] = np.mean(arrayUsim)
                    
        # Se determina el promedio de la magnetización para las nsimulaciones
        arrayMag[o] = np.mean(arrayMagsim)      
            
        # Se calcula el calor específico
        arrayC[o] = (1/nespines**2)*((u2-(arrayU[o])**2)/((array_kbt[o])**2))
         
        
    return arrayU, arrayMag, arrayC

arrayU1, arrayMag1, arrayC1 = ModeloIsingPropiedades(J, array_kBT, nEspines, frios_arriba, nPasos, nSimulaciones)
arrayU2, arrayMag2, arrayC2 = ModeloIsingPropiedades(J, array_kBT, nEspines, frios_abajo, nPasos, nSimulaciones)
arrayU3, arrayMag3, arrayC3 = ModeloIsingPropiedades(J, array_kBT, nEspines, calientes, nPasos, nSimulaciones)

# Las soluciones analíticas son las siguientes(con B=0.05):
uAnalitica = -J*nEspines*np.tanh(J/array_kBT)
cAnalitica = (J/array_kBT)**2/(np.cosh(J/array_kBT)**2)
mAnalitica = (nEspines*np.exp(J/array_kBT)*np.sinh(B/array_kBT))/np.sqrt(np.exp((2*J)/array_kBT)*(np.sinh(B/array_kBT))**2 + np.exp((-2*J)/array_kBT))


# Se realizan las gráficas de energía interna, magnetización y calor específico
# Para la configuración incial de espines fríos arriba

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayU1, label="Solución numérica")
ax.plot(array_kBT, uAnalitica, label="Solución analítica")
ax.set_title('Energía interna del modelo de Ising 1-D con espines fríos (arriba)')
ax.set_xlabel('kBT')
ax.set_ylabel('Energía interna')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayMag1, label="Solución numérica")
ax.plot(array_kBT, mAnalitica, label="Solución analítica")
ax.set_title('Magnetización del modelo de Ising 1-D con espines fríos (arriba)')
ax.set_xlabel('kBT')
ax.set_ylabel('Magnetización')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayC1, label="Solución numérica")
ax.plot(array_kBT, cAnalitica, label="Solución analítica")
ax.set_title('Calor específico del modelo de Ising 1-D con espines fríos (arriba)')
ax.set_xlabel('kBT')
ax.set_ylabel('Calor específico')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayC1, label="Solución numérica")
ax.set_title('Calor específico del modelo de Ising 1-D con espines fríos (arriba)')
ax.set_xlabel('kBT')
ax.set_ylabel('Calor específico')
ax.legend()
plt.show()


# Para la configuración incial de espines fríos abajo

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayU2, label="Solución numérica")
ax.plot(array_kBT, uAnalitica, label="Solución analítica")
ax.set_title('Energía interna del modelo de Ising 1-D con espines fríos (abajo)')
ax.set_xlabel('kBT')
ax.set_ylabel('Energía interna')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayMag2, label="Solución numérica")
ax.plot(array_kBT, mAnalitica, label="Solución analítica")
ax.set_title('Magnetización del modelo de Ising 1-D con espines fríos (abajo)')
ax.set_xlabel('kBT')
ax.set_ylabel('Magnetización')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayC2, label="Solución numérica")
ax.plot(array_kBT, cAnalitica, label="Solución analítica")
ax.set_title('Calor específico del modelo de Ising 1-D con espines fríos (abajo)')
ax.set_xlabel('kBT')
ax.set_ylabel('Calor específico')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayC2, label="Solución numérica")
ax.set_title('Calor específico del modelo de Ising 1-D con espines fríos (abajo)')
ax.set_xlabel('kBT')
ax.set_ylabel('Calor específico')
ax.legend()
plt.show()


# Para la configuración incial de espines calientes

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayU3, label="Solución numérica")
ax.plot(array_kBT, uAnalitica, label="Solución analítica")
ax.set_title('Energía interna del modelo de Ising 1-D con espines calientes')
ax.set_xlabel('kBT')
ax.set_ylabel('Energía interna')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayMag3, label="Solución numérica")
ax.plot(array_kBT, mAnalitica, label="Solución analítica")
ax.set_title('Magnetización del modelo de Ising 1-D con espines calientes')
ax.set_xlabel('kBT')
ax.set_ylabel('Magnetización')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayC3, label="Solución numérica")
ax.plot(array_kBT, cAnalitica, label="Solución analítica")
ax.set_title('Calor específico del modelo de Ising 1-D con espines calientes')
ax.set_xlabel('kBT')
ax.set_ylabel('Calor específico')
ax.legend()
plt.show()

fig, ax = plt.subplots(dpi=120)
ax.plot(array_kBT, arrayC3, label="Solución numérica")
ax.set_title('Calor específico del modelo de Ising 1-D con espines calientes')
ax.set_xlabel('kBT')
ax.set_ylabel('Calor específico')
ax.legend()
plt.show()
