import random as rd
from re import I
import numpy as np
import matplotlib.pyplot as plt

def rad2cart(data, r, theta, x, y):
    nx = len(x)
    ny = len(y)
    output = np.zeros(nx,ny)
    #For each point x,y
    for i in range(nx):
        for j in range(ny):
            ix = x[i]
            iy = y[j]

            #Find the polar coordinates around x,y and translate to cartesian
            points = findSurPoints(data,r,theta,ix,iy)
            if (len(points) == 0):
                output[j][i] = 0
                continue

            #find a,b,c,d from axy + bx + cy + d = p(x,y)
            # (a,b,c,d) = x in Ax = b waarin b de waarden van de punten is
            A = np.zeros(4,4)
            b = [[points[0][2]],[points[1][2]],points[2][2],points[3][2]]
            for k in range(4):
                A[k][0] = points[k][0]*points[k][1]
                A[k][1] = points[k][0] 
                A[k][2] = points[k][1]
                A[k][3] = 1
            constants = np.linalg.solve(A,b)

            #value of x,y is now p(x,y)
            p = constants[0]*ix*iy + constants[1]*ix + constants[2]*iy + constants[3]
            output[j][i] = p
    return output

def rad2cartCoords(r,theta):
    return (r*np.cos(theta),r*np.sin(theta))

def cart2radCoords(x,y):
    if (x != 0):
        theta = np.arctan(y/x)
        if (x < 0):
            theta *= -1
    else:
        theta = np.pi
        if (y < 0):
            theta *= -1
    return (x**2 + y**2, theta)


def findSurPoints(data, r, theta, x, y):
    (ir,it) = cart2radCoords(x,y)
    
    #find the first instance where r[i] > ir
    i = 0
    while (i < len(r) and r[i] < ir):
        i += 1
    
    #find the first instance where theta[j] > it
    j = 0
    while (j < len(theta) and theta[j]< it):
        j += 1
    
    #Check of (ir,it) buiten onze grenzen ligt
    if (j == 0 or j == len(theta) or i == 0 or i == len(r)):
        return []

    
    points = [[i-1,j-1], [i - 1,j], [i, j - 1], [i,j]]#indexes van onze punten
    for point in points:
        (ix,iy) = rad2cartCoords(r[point[0]], theta[point[1]])
        ig = data[point[0],point[1]]
        point = [ix,iy,ig]#returnen van alle punten eromheen de x,y coordinaten samen met de waarde op dat punt
    return points

