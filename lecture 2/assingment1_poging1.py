import numpy as np
import matplotlib.pyplot as plt


def rad2cart(data, r, theta, x, y):
    hr = r[1] - r[0] #gelijke afstand aangenomen
    ht = theta[1] - theta[0]
    result = []
    for xcord in x:
        values = []
        for ycord in y:
            
            #omrekening naar poolcoördinaten
            rval = np.sqrt(xcord**2 + ycord**2)
            #ro en r1 zijn de punten waar (x,y ) in poolcoördinaten tussen ligt
            r0 = int((rval - r[0]) // hr)
            r1 = r0 + 1
            
            tval = np.arctan(ycord/xcord)
            t0 = int((tval - theta[0]) // ht)
            t1 = t0 + 1
            
            if r1 < len(r) and t1 < len(theta):
                s = 0            
                #feitelijk worden de ingesloten oppervlaktes hier berekend, als in 
                #het verslag uitgelegd is dat hier niet helemaal de juiste methode
                s += data[r1][t1]*((rval - r[r0])/hr)*((tval - theta[t0])/ht)
                s += data[r1][t0]*((rval - r[r0])/hr)*((theta[t1] - tval)/ht)
                s += data[r0][t1]*((r[r1] - rval)/hr)*((tval - theta[t0])/ht)
                s += data[r0][t0]*((r[r1] - rval)/hr)*((theta[t1] - tval)/ht)
                values.append(s)
            else: 
                #programma kan enkel interpoleren, niet extrapoleren. Om in dat geval
                #0 neer te zetten is een arbitraire keuze.
                values.append(0) 
                
        result.append(values)
    return result

#een simpele test voor de werking met specifieke waardes
def testrad2cart():
    data = [[0,1], [2,3]]
    r = [1,2]
    theta = [0,np.pi/2]
    x = [0.001,0.999]
    y = [1,1.5]
    print(rad2cart(data, r, theta, x, y))

