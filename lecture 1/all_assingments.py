import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import radon, iradon
from skimage.data import shepp_logan_phantom

nx = 400 #dimensie van de afbeelding
na = 360 #aantal hoeken waaronder wordt gemeten
theta = np.linspace(0., 180., na)
sigma = [0, 1, 5, 10, 50]
 
def slp():
    return shepp_logan_phantom()

def cirkel(dim):
    cirkel = np.zeros((dim, dim))
    for x in range(dim):
        for y in range(dim):
            #De telling begint bij 0. Dat zorgt dat bij bv. 400 pixels het midden tussen pixel 199 en 200 ligt
            if (x-(dim-1)/2)**2 + (y-(dim-1)/2)**2 < 10000: #de waarde 10.000 is arbitrair
                cirkel[x][y] = 1    
    return cirkel
    
def vierkant(dim):
    vierkant = np.zeros((dim, dim))
    for x in range(dim):
        for y in range(dim):
            if (dim-1)/2-(dim/4) < min(x,y) <= max(x,y) < (dim-1)/2+(dim/4):
                vierkant[x][y] = 1    
    return vierkant

    
def teken_sinogram(u): #tekent oorspronkelijke plaatje en sinogram  
    f = radon(u, theta=theta)
    fig,ax = plt.subplots(1,2)

    ax[0].imshow(u,extent=(-1,1,-1,1),vmin=0)
    ax[0].set_xlabel(r'$x$')
    ax[0].set_ylabel(r'$y$')
    #in dit geval hebben de waarden op de y-as een concrete betekenis en niet enkel relatief ten opzichte van elkaar, daarom is de extent anders
    ax[1].imshow(f,extent=(0,180, -1, 1), vmin = 0)
    ax[1].set_xlabel(r'$\theta$ in graden')
    ax[1].set_ylabel(r'positie ten opzichte van midden')
    ax[1].set_aspect(90)

    fig.tight_layout()

def func05(theta): #geeft de functie van de intensiteit op hoogte 0.5
    theta_rad = (theta*np.pi)/180
    #de motivatie achter deze functie staat in het verslag. Het nut van deze specifieke is beperkt.
    return ((1 - np.tan(theta_rad))/(2) + 1/(2*np.sin(theta_rad)*np.cos(theta_rad)) - 1/(2*np.sin(theta_rad)))/(np.sqrt(2)*np.cos(theta_rad))
    
def vergelijk(): #maakt een vergelijking tussen het exacte sinogram op hoogte 0.5 en die door de computer berekend worden
    u1 = vierkant(400)
    u2 = vierkant(32)

    f1 = radon(u1, theta=theta)
    f2 = radon(u2, theta=theta)
    
    f1max = np.matrix(f1).max()
    f2max = np.matrix(f2).max()
    
    #we duiden de eerste met boven aan, aangezien deze in de matrix van het sinogram hoger ligt en visueel in de getoonde sinogrammen dus ook
    f1_boven = [element/f1max for element in f1[299][0:46]]
    f1_onder = [element/f1max for element in f1[300][0:46]]
    
    f2_boven = [element/f2max for element in f2[23][0:46]]
    f2_onder = [element/f2max for element in f2[24][0:46]]
    
    exact = np.zeros(46)
    exact[0] = 1/np.sqrt(2) #de eerste waarde moet exact worden berekend aangezien de meetkunde uitleg er niet voor werkt, puur tussen de 0 en 45 graden
    for i in range(1,46):
        exact[i] = func05(i)
    
    waarden = list(range(46))
    
    plt.plot(waarden, exact, label='exact')
    plt.plot(waarden, f1_onder, label = '400 waarden, direct onder 0.5')
    plt.plot(waarden, f1_boven, label = '400 waarden, direct boven 0.5')
    plt.plot(waarden, f2_onder, label = '32 waarden, direct onder 0.5')
    plt.plot(waarden, f2_boven, label = '32 waarden, direct boven 0.5')    
    plt.xlabel('hoek in graden', fontsize = 14)
    plt.ylabel('geschaalde intensiteit', fontsize = 14)
    plt.legend(loc="upper right", fontsize = 14)
    plt.show()

def teken_verstoringen(sigma, u): #tekent oorspronkelijk plaatje en verstoorde versie
    f = radon(u, theta=theta)    

    aantal = len(sigma)
    fig,ax = plt.subplots(1, aantal)
    
    for i in range(aantal):
        verstoring = sigma[i]
        f_noisy = f + verstoring * np.random.randn(nx, na)
        u_noisy = iradon(f_noisy, theta=theta)
        ax[i].imshow(u_noisy, extent = (-1,1,-1,1), aspect = 1, vmin = 0)
        ax[i].set_xlabel(r'$x$')
        ax[i].set_ylabel(r'$y$')
        ax[i].set_title(r'$\sigma =$' + str(verstoring))

def bereken_verstoring(u1, u2): #kan gebruikt worden om sinogrammen te vergelijken of om plaatjes te vergelijken
     #Zoals in het verslag genoemd was er eerst een absolute in plaats van relatieve aanpak
    return (np.linalg.norm(u1 - u2))/(np.linalg.norm(u1))

def plot_verstoringen(sigma, u):
    f = radon(u, theta=theta)    
    aantal = len(sigma)
    fouten = np.zeros(aantal)
    
    for i in range(aantal):
        verstoring = sigma[i]
        f_noisy = f + verstoring * np.random.randn(nx, na)
        u_noisy = iradon(f_noisy, theta=theta)
        fouten[i] = bereken_verstoring(u, u_noisy)
    
    plt.xlabel(r'$\sigma$', fontsize = 14)
    plt.ylabel('relatieve fout', fontsize=14)
    plt.plot(sigma, fouten)
    
    
plot_verstoringen(sigma, slp())