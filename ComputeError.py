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
dt = np.array([10**-i for i in np.arange(2,7.1,0.3)])

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
    
    #EM computation
    aux = EM(dt=dt[i])
    rt = float(aux[0])
    r = np.array(aux[1]).flatten()
    rtEM[i] = rt
    rEM[i] = r
    
    #RK4M computation
    aux = RK4M(dt=dt[i])
    rt = float(aux[0])
    r = np.array(aux[1]).flatten()
    rtRK4M[i] = rt
    rRK4M[i] = r
    
    #percentage completed
    print(str(round(i*100/len(dt),1))+'%')
    
print(rSAM)

