#Compare integrations

import numpy as np

#%%
"""computamos el analitico usando SAM con dt=10^-18"""
#!!!!!!!!!!!!!!!!!!!!! revisar, no esta hecho todavia

from SAMcomentado import SAM

rAnaly = np.array([])
aux = SAM()
RT = float(aux[0])
"""colapsa a r en un array con las filas pegadas una tras otra"""
r = np.array(aux[1]).flatten()
print(r)

rAnaly = np.append(rAnaly, r)
print(rAnaly)

#%%
"""computamos RT y r para multiples dt de cada metodo, lo guardamos en un csv"""

from SAMcomentado import SAM
from EM import EM
from RK4M import RK4M



"""esto controla el orden de los dt y la cantidad de puntos computados"""
dt = np.array([10**-i for i in np.arange(2,3.1,0.3)])

#method 0 (SAM) method 1 (EM) method 2 (RK4M)
data = np.zeros((len(dt)*3,67))

rSAM = np.zeros((len(dt),64))
rtSAM = np.zeros(len(dt))
rEM = np.zeros((len(dt),64))
rtEM = np.zeros(len(dt))
rRK4M = np.zeros((len(dt),64))
rtRK4M = np.zeros(len(dt))

for i in range(len(dt)):
    #SAM computation
    aux = SAM(dt=dt[i])
    rt = float(aux[0])
    """colapsa a la matriz de r para un dado dt en un vector con una fila tras otra"""
    r = np.array(aux[1]).flatten()
    rtSAM[i] = rt
    rSAM[i] = r
    """agrego a la matriz de data la fila correspondiente a esta iteracion"""
    data[i,0] = 0
    data[i,1] = dt[i]
    data[i,2:] = rt
    
    
    #EM computation
    aux = EM(dt=dt[i])
    rt = float(aux[0])
    r = np.array(aux[1]).flatten()
    rtEM[i] = rt
    rEM[i] = r
    data[i+1,0] = 1
    data[i+1,1] = dt[i]
    data[i+1,2:] = rt
    
    #RK4M computation
    aux = RK4M(dt=dt[i])
    rt = float(aux[0])
    r = np.array(aux[1]).flatten()
    rtRK4M[i] = rt
    rRK4M[i] = r
    data[i+2,0] = 2
    data[i+2,1] = dt[i]
    data[i+2,2:] = rt
    
    #percentage completed
    print(str(round(i*100/len(dt),1))+'%')
    
    

"""agregue una matriz con primera columna metodo (0 para SAM 1 para EM y 2 para
 RK4M), segunda columna dt, tercera columna RT y siguientes 64 columnas el 
response field, revisar porque no se esta llenando bien"""
