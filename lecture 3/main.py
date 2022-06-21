import numpy as np


def RadonSparseMatrix(n, sn, thetan):
    result = np.zeros(len(sn)*len(thetan), n*n)
    idxPtr = [0 for i in range(len(sn)*len(thetan) + 1)]
    indices = []
    data = []
    for idxT in range(len(thetan)):
        for idxS in range(len(sn)):
            idxPtr[idxS + len(sn)*idxT + 1] = idxPtr[idxS + len(sn)*idxT] 
            s = sn[idxS]
            theta = thetan[idxT]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

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
            end = max([ax,bx,ay,by])

            dx = abs((2/n)/np.sin(theta))
            dy = abs((2/n)/np.cos(theta))

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
                    iy += diy
                else:
                    sx += dx
                    ix += dix
    return result

def RadonMatrix(n, sn, thetan):
    result = np.zeros(len(sn)*len(thetan), n*n)
    for idxS in range(len(sn)):
        for idxT in range(len(thetan)):
            s = sn[idxS]
            theta = thetan[idxT]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

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
            end = max([ax,bx,ay,by])

            dx = abs((2/n)/np.sin(theta))
            dy = abs((2/n)/np.cos(theta))

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
                    iy += diy
                else:
                    sx += dx
                    ix += dix
    return result

def Radon(u, n, sn, thetan):
    result = np.zeros(len(sn), len(thetan))
    for idxS in range(len(sn)):
        for idxT in range(len(thetan)):
            s = sn[idxS]
            theta = thetan[idxT]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

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
            end = max([ax,bx,ay,by])

            dx = abs((2/n)/np.sin(theta))
            dy = abs((2/n)/np.cos(theta))

            t = min(sy,sx)
            r = 0
            while (t < end):
                #move out of current pixel
                dt = min(sx,sy) - t
                t = min(sx,sy)
                if (ix >= 0 and ix < n and iy >= 0 and iy < n):
                    r += u[ix][iy]

                #find if we moved in the x or y direction
                if (sy < sx):
                    sy += dy
                    iy += diy
                else:
                    sx += dx
                    ix += dix
            result[idxT][idxS] = r
    return result


def RadonBasic(n,sn,thetan):
    for idxS in range(len(sn)):
        for idxT in range(len(thetan)):
            s = sn[idxS]
            theta = thetan[idxT]
            ax = (-1 - s*np.cos(theta))/np.sin(theta)
            bx = (1 - s*np.cos(theta))/np.sin(theta)

            ay = -(-1 - s*np.sin(theta))/np.cos(theta)
            by = -(1 - s*np.sin(theta))/np.cos(theta)

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
            end = max([ax,bx,ay,by])

            dx = abs((2/n)/np.sin(theta))
            dy = abs((2/n)/np.cos(theta))

            t = min(sy,sx)
            while (t < end):
                #move out of current pixel
                dt = min(sx,sy) - t
                t = min(sx,sy)
                if (ix >= 0 and ix < n and iy >= 0 and iy < n):
                    #at location (ix,iy) we want to put value dt
                    #dt is how long we were in the pixel
                    dt = 0
                #find if we moved in the x or y direction
                if (sy < sx):
                    sy += dy
                    iy += diy
                else:
                    sx += dx
                    ix += dix
    
        
        
