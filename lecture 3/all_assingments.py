import numpy as np


def RadonSparseMatrix(n, sn, thetan):
    result = np.zeros(len(sn)*len(thetan), n*n)
    idxPtr = [0 for i in range(len(sn)*len(thetan) + 1)]
    indices = []
    data = []
    for idxT in range(len(thetan)):
        for idxS in range(len(sn)):
            idxPtr[idxS + len(sn)*idxT + 1] = idxPtr[idxS + len(sn)*idxT] #by default it's the previous value, increase by 1 for every value we add
            s = sn[idxS]
            theta = thetan[idxT]

            #calculate intersections with [-1,1]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

            #start at the smallest t
            #which one is smaller decides if we start at the top or bottom of the cube and which way we move
            sx = 0
            ix = 0
            dix = 0
            if (ax < bx):
                sx = ax
                ix = -1
                dix = 1
            else:
                sx = bx
                ix = n
                dix = -1
    
            sy = 0
            iy = 0
            diy = 0
            if (ay < by):
                sy = ay
                iy = -2
                diy = 1
            else:
                sy = by
                iy = n
                diy = -1
            end = max([ax,bx,ay,by]) #after the end we dont need to know if it is an x or y intersection

            #exact distance between each x and y intersection
            dx = abs((2/n)/np.sin(theta)) 
            dy = abs((2/n)/np.cos(theta))

            #Start outside of the cube
            t = min(sy,sx)
            while (t < end):
                #move out of current pixel
                dt = min(sx,sy) - t
                t = min(sx,sy)
                if (ix >= 0 and ix < n and iy >= 0 and iy < n):
                    idxPtr[idxS + len(sn)*idxT + 1] += 1

                    indices.append(ix + n*iy)
                    data.append(dt)
                    #at location (ix,iy) we want to put value dt
                    #dt is how long we were in the pixel
                #find if we moved in the x or y direction
                if (sy < sx):
                    sy += dy
                    iy += diy #y pixel index
                else:
                    sx += dx
                    ix += dix #x pixel index
    return result

def RadonMatrix(n, sn, thetan):
    result = np.zeros(len(sn)*len(thetan), n*n)
    for idxS in range(len(sn)):
        for idxT in range(len(thetan)):
            s = sn[idxS]
            theta = thetan[idxT]
            #calculate intersections with [-1,1]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

            #start at the smallest t
            #which one is smaller decides if we start at the top or bottom of the cube and which way we move
            sx = 0
            ix = 0
            dix = 0
            if (ax < bx):
                sx = ax
                ix = -1
                dix = 1
            else:
                sx = bx
                ix = n
                dix = -1
    
            sy = 0
            iy = 0
            diy = 0
            if (ay < by):
                sy = ay
                iy = -2
                diy = 1
            else:
                sy = by
                iy = n
                diy = -1
            end = max([ax,bx,ay,by]) #after the end we dont need to know if it is an x or y intersection

            #exact distance between each x and y intersection
            dx = abs((2/n)/np.sin(theta)) 
            dy = abs((2/n)/np.cos(theta))

            #Start outside of the cube
            t = min(sy,sx)
            while (t < end):
                #move out of current pixel
                dt = min(sx,sy) - t
                t = min(sx,sy)
                if (ix >= 0 and ix < n and iy >= 0 and iy < n):
                    result[ix + iy*n][idxS + len(sn)*idxT] = dt
                    #at location (ix,iy) we want to put value dt
                    #dt is how long we were in the pixel
                #find if we moved in the x or y direction
                if (sy < sx):
                    sy += dy
                    iy += diy #y pixel index
                else:
                    sx += dx
                    ix += dix #x pixel index
    return result

def Radon(u, n, sn, thetan):
    result = np.zeros(len(sn), len(thetan))
    for idxS in range(len(sn)):
        for idxT in range(len(thetan)):
            s = sn[idxS]
            theta = thetan[idxT]
            #calculate intersections with [-1,1]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

            #start at the smallest t
            #which one is smaller decides if we start at the top or bottom of the cube and which way we move
            sx = 0
            ix = 0
            dix = 0
            if (ax < bx):
                sx = ax
                ix = -1
                dix = 1
            else:
                sx = bx
                ix = n
                dix = -1
    
            sy = 0
            iy = 0
            diy = 0
            if (ay < by):
                sy = ay
                iy = -2
                diy = 1
            else:
                sy = by
                iy = n
                diy = -1
            end = max([ax,bx,ay,by]) #after the end we dont need to know if it is an x or y intersection

            #exact distance between each x and y intersection
            dx = abs((2/n)/np.sin(theta)) 
            dy = abs((2/n)/np.cos(theta))

            #Start outside of the cube
            t = min(sy,sx)
            r = 0 #value of this line
            while (t < end):
                #move out of current pixel
                dt = min(sx,sy) - t
                t = min(sx,sy)
                if (ix >= 0 and ix < n and iy >= 0 and iy < n):
                    r += u[ix][iy]*dt #add value of current pixel times length

                #find if we moved in the x or y direction
                if (sy < sx):
                    sy += dy
                    iy += diy #y pixel index
                else:
                    sx += dx
                    ix += dix #x pixel index
            result[idxT][idxS] = r
    return result


