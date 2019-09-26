import numpy as np

def SAM(op1=8,op2=7,alpha=0.75,C=0.5,dt=0.001):
    """Solves the IN model ecuations 

    INPUT
    op1: first operand, maximum (integer)
    op2: second operand, minimum (integer)
    
    OUTPUT
    RT: Response Time
    r: Response field values at RT
    """
    
    

    """corresponde a media tabla de multiplicar, desde el 2 hasta el 9"""
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
    r = np.zeros((8,8)) #response field
    maxi = np.zeros(8) #op1 field
    mini = np.zeros(8) #op2 field
    
    """al lugar del array que corresponde a op1 le coloca un B=20 (supone op1 
    como el operando maximo) ya que a t=0 el coeficiente b de este cambio es BI
    con I=1 solo si se presento el numero, y para los demas t I=0 para todo 
    nodo del campo""" 
    maxi[op1-2] = B
    mini[op2-2] = B
    
    test = 0
    t = 0
    while t < 0.1:
        test += 1

        #response
        """aux1 y aux2 se recorren de manera tal de atravesar la media 
        tabla de multiplicar a la que se refieren Dec y Uni, que a su vez son 
        matrices llenas de ceros en la mitad que no será recorrida"""
        """cada uno de los i recorren los z, es decir, los del campo semantico 
        (tanto las decenas como las unidades, evaluando si el correspondiente
        del campo de respuesta coincide con el numero al que debiera estar 
        asociado, si esto se cumple entonces se suma a la sumatoria sum1"""
        """multiplica por S los valores de la sumatoria ya que los pesos en el 
        response field son 0 o S"""
        for aux1 in range(8):
            for aux2 in range(aux1):
                sum1 = 0
                for i in range(9):
                    if Dec[aux1,aux2]==i:
                        sum1 += S*T[i]
                for i in range(10):
                    if Uni[aux1,aux2]==i:
                        sum1 += S*U[i]
                coef_a = -1 - sum1
                coef_b = B*sum1
                r[aux1,aux2] =  (coef_b/coef_a + r[aux1,aux2]) * np.exp(coef_a*dt) - coef_b/coef_a
                
                if r[aux1,aux2]<0: 
                    r[aux1,aux2] = 0
        
        #decades field
        """sum2 guarda en un array los distintos valores de la sumatoria sobre
        z (nodos del decades field) que no incluyen al i en cuestion, aux1 y 
        aux2 cumplen la misma funcion que en el response field, recorriendo
        ahora los posibles valores del semantic field, analogo al response"""
        sum2 = np.zeros(9)
        for aux in range(9):
            sum2[aux] = sum(T)- T[aux]
        for aux in range(9):
            sum1 = 0
            for aux1 in range(8):
                for aux2 in range(aux1):            
                    if Dec[aux1,aux2]==aux:
                        sum1 += F[aux1,aux2]
            coef_a = -1 - M*sum1 - sum2[aux]
            coef_b = M*B*sum1
            T[aux] =  (coef_b/coef_a + T[aux]) * np.exp(coef_a*dt) - coef_b/coef_a
            if T[aux]<0:
                T[aux]  = 0

        #units field
        sum2 = np.zeros(10)
        for aux in range(10):
            sum2[aux] = sum(U)- U[aux]
        for aux in range(10):
            sum1 = 0
            for aux1 in range(8):
                for aux2 in range(aux1):            
                    if Uni[aux1,aux2]==aux:
                        sum1 += F[aux1,aux2]
            coef_a = -1 - M*sum1 - sum2[aux]
            coef_b = M*B*sum1
            U[aux] =  (coef_b/coef_a + U[aux]) * np.exp(coef_a*dt) - coef_b/coef_a
            if U[aux]<0:
                U[aux] = 0
        
        """aca aux1 y aux2 recorren la estructura del semantic field como antes
        mientras que i recorre los posibles valores de inputs"""
        #semantic field
        for aux1 in range(8):
            for aux2 in range (aux1):         
                sum1 = 0
                for i in range(8):
                    sum1 += maxi[i]*(np.exp(-alpha*np.absolute(aux1-i))-C)
                for i in range(8):
                    sum1 += mini[i]*(np.exp(-alpha*np.absolute(aux2-i))-C)
                coef_a = -1 - sum1
                coef_b = B*sum1
                F[aux1,aux2] =  (coef_b/coef_a + F[aux1,aux2]) * np.exp(coef_a*dt) - coef_b/coef_a
                if F[aux1,aux2]<0:
                    F[aux1,aux2] = 0
    
        """revisar condicion inicial del input field y su evolucion en el 
        tiempo dada dicha condicion inicial (que incluye las lineas previas al
        while en donde se coloca una B en los nodos activados)"""    
        #input fields
        maxi *= np.exp(-dt)
        mini *= np.exp(-dt)
        
        """cuando el valor en el response field llega a tau se da una 
        respuesta, la correspondiente a ese tau"""
        if np.amax(r) > tau:
                    break
            
        t = t + dt

    """ encuentro el indice para el maximo valor de r y lo guardo en un array
    de numpy, calculo el producto prod
    
    index = np.array([int(np.where(r==np.amax(r))[0]),int(np.where(r==np.amax(r))[1])])
    prod = np.array([Dec[index[0],index[1]],Uni[index[0],index[1]]])"""
    RT = t
    return RT, r



print(SAM())





