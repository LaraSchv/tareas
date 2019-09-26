"""This function solves the IN model ecuations 

INPUT
op1: first operand (integer)
op2: second operand (integer)

OUTPUT
RT: Response Time
prod
Decade and unit parts of the aswers
"""
import numpy as np


def EM(op1=8,op2=7,alpha=0.75,C=0.5,dt=0.001):
    """Solves the IN model ecuations 

    INPUT
    op1: first operand, maximum (integer)
    op2: second operand, minimum (integer)
    
    OUTPUT
    RT: Response Time
    r: Response field values at RT
    """
    
    Dec = np.array([[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,1,1,0,0,0,0,0],
                    [1,1,2,2,0,0,0,0],[1,1,2,3,3,0,0,0],[1,2,2,3,4,4,0,0],
                    [1,2,3,4,4,5,6,0],[1,2,3,4,5,6,7,8]])
    Uni = np.array([[4,0,0,0,0,0,0,0],[6,9,0,0,0,0,0,0],[8,2,6,0,0,0,0,0],
                    [0,5,0,5,0,0,0,0],[2,8,4,0,6,0,0,0],[4,1,8,5,2,9,0,0],
                    [6,4,2,0,8,6,4,0],[8,7,6,5,4,3,2,1]])
    
    #input parameters
    B = 20
    M = 10
    S = 2
    tau = 19.5
    
    F = np.zeros((8,8)) #semantic field
    T = np.zeros(9) #tens field
    U = np.zeros(10) #unit field
    r = np.zeros((8,8)) #responce field
    maxi = np.zeros(8) #op1 field
    mini = np.zeros(8) #op2 field
    
    maxi[op1-1] = B
    mini[op2-1] = B

    test = 0
    t = 0
    while t<0.1:
        test += 1
        
        #response
        for aux1 in range(8):
            for aux2 in range(aux1):
                sum1 = 0
                for i in range(9):
                    if Dec[aux1,aux2]==i:
                        sum1 += S*T[i]
                for i in range(10):
                    if Uni[aux1,aux2]==i:
                        sum1 += S*U[i]
                r[aux1,aux2] -= r[aux1,aux2]*dt + dt*(B-r[aux1,aux2])*sum1
                if r[aux1,aux2]<0:
                    r[aux1,aux2] = 0
        
        #decades field
        sum2 = np.zeros(9)
        for aux in range(9):
            sum2[aux] = sum(T)- T[aux]
        for aux in range(9):
            sum1 = 0
            for aux1 in range(8):
                for aux2 in range(aux1):            
                    if Dec[aux1,aux2]==aux:
                        sum1 += F[aux1,aux2]
            T[aux] -= T[aux]*dt + dt*M*sum1*(B-T[aux]) - dt*T[aux]*sum2[aux]
            if T[aux]<0:
                T[aux] = 0
        
        #units field
        sum3 = np.zeros(10)
        for aux in range(10):
            sum3[aux] = sum(U)-U[aux]
        for aux in range(10):
            sum1 = 0
            for aux1 in range(8):
                for aux2 in range(aux1):            
                    if Uni[aux1,aux2]==aux:
                        sum1 += F[aux1,aux2]
            U[aux] -= U[aux]*dt + dt*M*sum1*(B-U[aux]) - dt*U[aux]*sum3[aux]
            if U[aux]<0:
                U[aux] = 0 
        
        #semantic field
        for aux1 in range(8):
            for aux2 in range (aux1): #aux1 y aux2 son los índices de cada campo semántico          
                aux3 = 0
                for i in range(8):
                    aux3 += maxi[i]*(np.exp(-alpha*np.absolute(aux1-i))-C)
                for i in range(8):
                    aux3 += mini[i]*(np.exp(-alpha*np.absolute(aux2-i))-C)
                F[aux1,aux2] -= F[aux1,aux2]*dt + dt*aux3*(B-F[aux1,aux2])
                if F[aux1,aux2]<0:
                    F[aux1,aux2] = 0
    
        # =====================================================================
        #     maxi *= np.exp(-dt)
        #     mini *= np.exp(-dt)
        #     analytic solution for the input field could be used
        # =====================================================================
    
        #input fields
        maxi -= maxi*dt
        mini -= mini*dt

        if np.amax(r) > tau:
                    break
                
        t = t + dt
        
    RT = t
    return RT, r

#!!!!!!!!!!!!!!!!!!!! revisar porque no esta computando nada
print(EM())
