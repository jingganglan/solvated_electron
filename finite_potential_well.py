%matplotlib inline
from IPython.display import HTML
import matplotlib as mpl # matplotlib library for plotting and visualization
import matplotlib.pylab as plt # matplotlib library for plotting and visualization
import numpy as np #numpy library for numerical manipulation, especially suited for data arrays
import warnings
warnings.filterwarnings('ignore')

# Reading the input variables from the user
def fintite_potential(vbe,gyr):
    Vo = vbe
    L =  gyr

    hbar = 1.05457180013e-34 # planck's constant divided by 2*pi, in Joules*second
    melec = 9.10938356e-31 # mass of an electron in kg
    eVtoJ = 1.60217662e-19 # conversion factor from eV to Joules
    AngstromtoMeter = 1e-10 # conversion factor from Angstroms to meters
    val = np.sqrt(2.0*melec*eVtoJ)*AngstromtoMeter/(2.0*hbar) # prefactor sqrt(2*melec)/(2*hbar), for when L is in angstroms and Vo is in eV

    # Generating the graph
    plt.rcParams.update({'font.size': 18, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
    #fig, axes = plt.subplots(1, 2, figsize=(13,4))

    E = np.linspace(0.0, Vo, 10000)
    num = int(round((L*np.sqrt(Vo)*val-np.pi/2.0)/np.pi))
    # Removing discontinuity points
    for n in range(10000):
        for m in range(num+2):
            if abs(E[n]-((2.0*float(m)+1.0)*np.pi/(2.0*L*val))**2)<0.01: E[n] = np.nan
            if abs(E[n]-(float(m)*np.pi/(L*val))**2)<0.01: E[n] = np.nan



    #print ("The allowed bounded energies are:")
    # We want to find the values of E in which f_even and f_odd are zero
    f_even = lambda E : np.sqrt(Vo-E)-np.sqrt(E)*np.tan(L*np.sqrt(E)*val)
    f_odd = lambda E : np.sqrt(Vo-E)+np.sqrt(E)/np.tan(L*np.sqrt(E)*val)
    E_old = 0.0
    f_even_old = f_even(0.0)
    f_odd_old = f_odd(0.0)
    n = 1
    E_vals = np.zeros(999)
    # Here we loop from E = 0 to E = Vo seeking roots
    for E in np.linspace(0.0, Vo, 2000):
        f_even_now = f_even(E)
        # If the difference is zero or if it changes sign then we might have passed through a root
        if f_even_now == 0.0 or f_even_now/f_even_old < 0.0:
            # If the old values of f are not close to zero, this means we didn't pass through a root but
            # through a discontinuity point
            if (abs(f_even_now)<1.0 and abs(f_even_old)<1.0):
                E_vals[n-1] = (E+E_old)/2.0
                #print "  State #%3d (Even wavefunction): %9.4f eV, %13.6g J" % (n,E_vals[n-1],E_vals[n-1]*1.60217662e-19)
                n += 1
        f_odd_now = f_odd(E)
        # If the difference is zero or if it changes sign then we might have passed through a root
        if f_odd_now == 0.0 or f_odd_now/f_odd_old < 0.0:
            # If the old values of f are not close to zero, this means we didn't pass through a root but
            # through a discontinuity point
            if (abs(f_odd_now)<1.0 and abs(f_odd_old)<1.0):
                E_vals[n-1] = (E+E_old)/2.0
                #print "  State #%3d  (Odd wavefunction): %9.4f eV, %13.6g J" % (n,E_vals[n-1],E_vals[n-1]*1.60217662e-19)
                n += 1
        E_old = E
        f_even_old = f_even_now
        f_odd_old = f_odd_now
    nstates = n-1
    #print ("\nTHERE ARE %3d POSSIBLE BOUNDED ENERGIES" % nstates)
    if (E_vals[1]-E_vals[0] > 0):
        ex=E_vals[1]-E_vals[0]
    else:
        ex=0
    # Generating the energy diagram
    return gyr,vbe,ex,E_vals
x=np.empty([20,20])
y=np.empty([20,20])
z=np.empty([20,20])
i=0
j=0
for gyr in np.linspace(2,3,20):
    j=0
    for vbe in np.linspace(3,4,20):
        #print(i,j)
        x[i,j],y[i,j],z[i,j]=fintite_potential(vbe,gyr*2)
        j=j+1
    i=i+1

        
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
#matplotlib notebook



plt.figure(figsize = (5,5),dpi=100)
plt.imshow(z.transpose(), cmap='jet', interpolation = 'bicubic',origin='lower',aspect = 1.2)
plt.colorbar()
plt.xlabel('Gyration Radius (Å)',fontsize=14)
plt.ylabel('VBE',fontsize=14)
plt.xticks(np.arange(0,20,5),np.arange(2,3,1/4).round(2),fontsize=14)
plt.yticks(np.arange(0,20,5),np.arange(3,4,1/4).round(2),fontsize=14)

