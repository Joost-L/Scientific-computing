import numpy as np

x = np.array([-1, -.5, 0, .5, 1])
y = np.array([-1, -.5, 0, .5, 1])
h = 0.5
s = np.array([-1.25, -.75, -.25, .25, .75, 1.25])

theta = 1
si = -.5

sx = min((x[0] - si*np.cos(theta))/np.sin(theta),(x[-1] - si*np.cos(theta))/np.sin(theta)) 
sy = min((y[0] - si*np.cos(theta))/np.sin(theta), (y[-1] - si*np.cos(theta)/np.sin(theta)))



tx = (x - si*np.cos(theta)) / np.sin(theta)
ty = -(y - si*np.sin(theta)) / np.cos(theta)
t = np.sort(np.concatenate((tx,ty)))

xi = si*np.cos(theta) + t*np.sin(theta)
yi = si*np.sin(theta) - t*np.cos(theta)

xc = h*((si*np.cos(theta) + 0.5*(t[0:-1]+t[1:])*np.sin(theta))//h) + h/2
yc = h*((si*np.sin(theta) - 0.5*(t[0:-1]+t[1:])*np.cos(theta))//h) + h/2

def RadonMatrix(n,sn,thetan):
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
    
        
        
