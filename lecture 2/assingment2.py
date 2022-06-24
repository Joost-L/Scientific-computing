import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import radon, iradon
from skimage.transform.radon_transform import _get_fourier_filter
from scipy.fft import fft, ifft, fftfreq, fftshift,ifftshift, fft2, ifft2
from skimage.data import shepp_logan_phantom

def rad2cart(data, r, theta, x, y):
    """
    Bi-linear interpolation from radial grid (r,theta) to cartesian grid (x,y).
    Polar coordinates are represented with r > 0 and theta \in [0,2*pi).
    All grids are assumed to be regularly spaced.
    
    Input: 
        data     - 2d array of size (nr,ntheta) containing values on the radial grid
        r, theta - 1d arrays of size nr and ntheta containing the polar gridpoints
        x,y      - 1d arrays of size nx and ny containing the cartesian gridpoints
        
    Output:
        result   - 2d array of size (nx,ny) containing the interpolated values
    """

    # get some grid related quantities (assumes regularly spaced grids)
    nr = len(r)
    dr = r[1] - r[0]
    nt = len(theta)
    dt = theta[1] - theta[0]    
    nx = len(x)
    dx = x[1] - x[0]
    ny = len(y)
    dy = y[1] - y[0]
    
    
    # initialise result
    result = np.zeros((nx,ny), dtype=type(data[0,0]))
    
    # loop over output grid
    for i in range(nx):
        for j in range(ny):
            # polar coordinates of current point (x[i], y[j])
            theta_ij = np.arctan2(x[i], y[j]) + np.pi # to map from [-pi,pi) to [0,2*pi)
            r_ij     = np.sqrt(x[i]**2 + y[j]**2)
            
            # look for lower-left gridpoint
            k = int((r_ij - r[0]) // dr)
            l = int((theta_ij - theta[0]) // dt)
            
            # check if i,j is an internal grid point and apply interpolation formula, 
            # see e.g., https://en.wikipedia.org/wiki/Bilinear_interpolation
            if (0 <= k < nr-1) and (0 <= l < nt-1):
                A = (r[k+1] - r_ij) * (theta[l+1] - theta_ij) / (dr * dt)
                B = (r_ij - r[k  ]) * (theta[l+1] - theta_ij) / (dr * dt)
                C = (r[k+1] - r_ij) * (theta_ij - theta[l  ]) / (dr * dt)
                D = (r_ij - r[k  ]) * (theta_ij - theta[l  ]) / (dr * dt)
                
                result[i,j] = A*data[k,l] + B*data[k+1,l] + C*data[k,l+1] + D*data[k+1,l+1]
    return result

#check
nr = 400
nt = 400

# radial grid
r     = np.linspace(0,1,nr)
theta = np.linspace(0,2*np.pi,nt)

#cartesian grid
x = np.linspace(-1,1,nr)
y = np.linspace(-1,1,nr)
        
# settings
nx = 400
na = 360
theta = np.linspace(0., 360, na)
sigmas = [0, 1, 5, 10, 50]

# phantom
u = shepp_logan_phantom()[::1, ::1]
u_rec = []

for i in range(len(sigmas)):  
    f = radon(u, theta = theta)
    f_noisy = f + sigmas[i] * np.random.randn(nx,na)

    # 1D Fourier and corresponding coordinates
    f_hat = fftshift(fft(ifftshift(f_noisy,axes=0), axis = 0), axes = 0)
    ks = fftshift(fftfreq(nx, 2/nx))

    # 2D frequency grid
    kx = fftshift(fftfreq(nx, 2/nx))
    ky = fftshift(fftfreq(nx, 2/nx))
    
    # interpolate; note that we use a sinogram for theta \in (0,2*pi) because we chose to implement
    # the rad2cart function to work only with a positive radius and angles in (0,2*pi). 
    # This means that currently only half of the sinogram is used; the other half is redundant due to symmetry.
    u_hat = rad2cart(f_hat, -ks, 2*np.pi*np.flip(theta) / 360, kx, ky)
    
    # inverse fft
    u_rec.append(np.real(fftshift(ifft2(ifftshift(u_hat)))))

def teken_verstoringen():
    
    fig, ax = plt.subplots(1,len(sigmas))
    for i in range(len(sigmas)):     
        ax[i].imshow(u_rec[i],vmax=1, vmin=0,extent=(x[0],x[-1],y[-1],y[0]))
        ax[i].set_aspect(1)
        ax[i].set_title(r'$\sigma$ = ' + str(sigmas[i]))

def bereken_verstoring(u1, u2): #kan gebruikt worden om sinogrammen te vergelijken of om plaatjes te vergelijken
     #Zoals in het verslag genoemd was er eerst een absolute in plaats van relatieve aanpak
    return (np.linalg.norm(u1 - u2))/(np.linalg.norm(u1))

def plot_verstoringen():    
    fouten = []
    for i in range(len(sigmas)):
        fouten.append(bereken_verstoring(u, u_rec[i]))
    
    plt.xlabel(r'$\sigma$', fontsize = 14)
    plt.ylabel('relatieve fout', fontsize=14)
    plt.plot(sigmas, fouten)
        

teken_verstoringen()