# -------------------------Assumptions/changes---------------------
# correction on R near poles for globe to flat map
# added correction on R for J2 effect (11 degree shift per orbit)
# albedo data file scaled to 0-0.35 matching found albedo
# 

from progress.bar import Bar
import Interpolation
import numpy as np
import matplotlib.pyplot as plt
from math import *
from matplotlib import gridspec
from PIL import Image
import os, sys, time
import image
from tqdm import tqdm_gui
import warnings
# ------------------------- Derivatives -------------------------------

# Forward, backward euler and central difference derivatives
def fd_deriv(x, f): return (f[1:] - f[:-1])/(x[1:] - x[:-1]), x[:-1]
def bd_deriv(x, f): return (f[1:] - f[:-1])/(x[1:] - x[:-1]), x[1:]
def cd_deriv(x, f): return (f[2:] - f[:-2])/(x[2:] - x[:-2]), x[1:-1]

def plot_derivs():
    f_fd, x_fd = fd_deriv(xx, spline)
    f_bd, x_bd = bd_deriv(xx, spline)
    f_cd, x_cd = cd_deriv(xx, spline)

    plt.suptitle('Element ' + str(element_number) + ' in ' + filename)
    plt.subplot(211); plt.title("Spline reconstruction"); plt.xlim(t_start, t_end); plt.xlabel("Time [s]"); plt.ylabel("Temperature [°C]")
    plt.plot(xx, spline, label='Spline reconstruction')
    plt.subplot(212); plt.title("Numerical derivative");  plt.xlim(t_start, t_end); plt.xlabel("Time [s]"); plt.ylabel(r"$ {\delta T}/{\delta t}\  [K / t] $")
    plt.plot(x_cd, f_cd, label='Central difference', linestyle='-'); plt.plot(x_bd, f_bd, label='Backward difference derivative', linestyle='-')
    plt.plot(x_fd, f_fd, label='Forward difference', linestyle='-'); plt.legend(loc='best')
    plt.axhline(y=0, color = 'k')
    plt.show()


# ------------------------- Incoming heat flux calculations -------------------------------

t_start, t_end, dt = 0, 15000, 0.708047411                       # Set boundary for graph and spline evaluation
t_steps = int(t_end/dt) + 1    
xx = np.linspace(t_start, t_end, t_steps)                     # Set datapoints to evaluate splne at


splines = np.zeros(shape = (8, t_steps))

for i in range(0, 8):

    filename, element_number = 'Raw_data/Mirror_segments.csv', i  #CHANGE THIS BASED ON THE LOCATION ON YOUR DEVICE!
    spline = Interpolation.spline(filename, element_number, xx, t_start, t_end, plotting = False)
    splines[i-2] = spline

T_front = (splines[0] + splines[2] + splines[4] + splines[6])/4 + 273.15
T_back = (splines[1] + splines[3] + splines[5] + splines[7])/4 + 273.15
T_avg = (splines[0] + splines[2] + splines[4] + splines[1] + splines[3] + splines[5] + splines[6] + splines[7])/8 + 273.15
dT_dt = cd_deriv(xx, T_avg)[0]

T_front, T_back, T_avg, xx = T_front[1:-1], T_back[1:-1], T_avg[1:-1], xx[1:-1]
sigma = 5.670374*10**(-8)
e_front, e_back, A, m, C = 0.035, 0.05, 1.0751, 6.5, 690

Q_out , Q_balance = e_front * sigma * (T_front*T_front*T_front*T_front) + e_back * sigma * (T_back*T_back*T_back*T_back), 1/A * m * C * dT_dt
Q_in = (Q_out + Q_balance)/e_front

avg = 717.5

Albedo = 0.35*((Q_in - avg)/356.5929343)
for q in range(len(Albedo)):
    if Albedo[q] < 0:
        Albedo[q] = 0
    
# -------------------------  World map -------------------------------
# initial values
dogleg = True
ktab = []
sumtab = []
plotsum = []
kplottab = []
plot2sum = []
kplot2tab = []

#load map and its data
img = Image.open('Earthalbedo.png').convert('LA')
img.save('Earthalbedogrey.png')
FILENAME='Earthalbedogrey.png' #image can be in gif jpeg or png format 
im=Image.open(FILENAME).convert('RGB')
pix=im.load()

k, k_final, dk = -180, 180, 1
longrun = False
zeroiteration = True
firstiteration = False
seconditeration = False

#loop to obtain all data
while dogleg == True:
    warnings.simplefilter("ignore")
    for i in tqdm_gui(range(388)):
        time.sleep(0.0001)
     

        #intial values
        xtab = []
        ytab = []
        x2tab = []
        y2tab = []
        albedotab = []
        ttab = []
        difftab = []

        #obtain orbit of best fit
        if k > k_final and zeroiteration == True:
            zeroiteration = False
            firstiteration = True
            for i in range(len(sumtab)-90):
                plotsum.append(sumtab[i])
                kplottab.append(ktab[i] + 90)
            for i in range(len(sumtab)- 270):
                i = i + 270
                plot2sum.append(sumtab[i])
                kplot2tab.append(ktab[i] - 270)
            k = ktab[sumtab.index(min(sumtab))] - 10
            k_final = ktab[sumtab.index(min(sumtab))] + 10
            longrun = True
            n = 10
        if k > k_final and firstiteration == True:
            firstiteration = False
            seconditeration = True
            k = ktab[sumtab.index(min(sumtab))] - 2
            k_final = ktab[sumtab.index(min(sumtab))] + 2
            Longrun = True
            n = 40
        if k > k_final and seconditeration == True:
            dogleg = False
            seconditeration == False
            k = ktab[sumtab.index(min(sumtab))] 
            longrun = True
            n = 119

        #initial valus
        jump = False
        orbit = True
        V_y = 0
        i = 0
        t = 0
        pxold = 1000
        pyold = 1000
        v_absground = 7.062206527 #km/s
        inclination_orbit = 98.31743089052 #degrees
        R_e = 6378.137
        longitude = k
        x = (longitude/180)*(20037.50834)
        x_old = x
        y = 18048.35
        y_old = y
        y_new = y

        #running loop for the orbit
        while y_new <= y:
            if abs(x - x_old) <= 20038:
                fact = (abs(1-(abs(x-x_old)/10019)))*(-inclination_orbit+180)
            else:
                fact = (abs(3-(abs(x-x_old)/10019)))*(-inclination_orbit+180)
            if V_y <=0:
                latitude = (y_new/111.318845)-90
                R = 6378.137*cos(radians(latitude)) #km
                v_earth = (2*pi*R)/(24*60*60)
                v_abs = sqrt(v_absground**2-v_earth**2)
                v_angle = degrees(atan((-v_earth+v_absground*sin(radians(-inclination_orbit+90-fact)))/(v_absground*cos(radians(-inclination_orbit+90-fact)))))
                V_y = v_abs*cos(radians(-180+v_angle))
                V_x = v_abs*sin(radians(-180+v_angle))*(6378.137/R)
                x_new = x_old + (5/abs(v_abs))*V_x
                y_new = y_old + (5/abs(v_abs))*V_y
                if y_new <= 20038 - y:
                    V_y = 0.01
            if V_y > 0:
                latitude = (y_new/111.318845)-90
                R = 6378.137*cos(radians(latitude)) #km
                v_earth = (2*pi*R)/(24*60*60)
                v_abs = sqrt(v_absground**2+v_earth**2)
                v_angle = degrees(atan((v_earth+v_absground*sin(radians(inclination_orbit-90+fact)))/(v_absground*cos(radians(inclination_orbit-90+fact)))))
                V_y = v_abs*cos(radians(v_angle))
                V_x = v_abs*sin(radians(v_angle))*(6378.137/R) 
                x_new = x_old + (5/abs(v_abs))*V_x
                y_new = y_old + (5/abs(v_abs))*V_y
            #add data of position and albedo
            if x_new >= 20038:
                dub = x_new
                jump = True
                x_new = dub - 40076
            x_old = x_new
            y_old = y_new
            if jump == False and dogleg == False:
                ytab.append(y_new)
                xtab.append(x_new)
            if jump == True and dogleg == False:
                y2tab.append(y_new)
                x2tab.append(x_new)
            pixel_x = int(2159*((x_new + 20038)/40076))
            pixel_y = int(1079-1079*(y_new/20038))
            if longrun == True:
                if pixel_x == pxold and pixel_y == pyold:
                    albedotab.append(albedotab[-1])
                else:
                    m =  n + 1
                    dattab = []
                    angtab = []
                    for f in range(-n, m):
                        for g in range( -abs(n - abs(f)), abs(m - abs(f) + 1)):
                            ang = cos(atan((sqrt(((18.554*g)**2)+((18.554*f)**2)))/505.1964))
                            angtab.append(ang)
                            if pixel_x + f <= 2159:
                                Pf = pixel_x + f
                            else:
                                pf = pixel_x + f - 2160

                            if pixel_y + g <= 1079:
                                Py = pixel_y + g
                            else:
                                Py = pixel_y + g - 1080
                            dattab.append(ang*0.35*((pix[Pf, Py][0])/255))
                    albedotab.append(sum(dattab)/sum(angtab))
            else:
                albedotab.append(0.35*((pix[pixel_x, pixel_y][0])/255))

            #update time and y for next loop            
            ttab.append(t)
            t = t + (5/abs(v_abs))
            pxold = pixel_x
            pyold = pixel_y
        
        # plot graph 3 and 4
        for w in range(len(albedotab)):
            difftab.append(abs(Albedo[w+2119]-albedotab[w]))
        sumtab.append((sum(difftab)/len(difftab)))
        
        ktab.append(k)
        print(k)
        k = k + dk
 
#delete final k as it is a run for the plot
del plotsum[-1]
del kplottab[-1]
del plotsum[-1]
del kplottab[-1]

print("Orbital period = ", ttab[-1], "[s]")
print("Average error = ", sumtab[-1], "[-]")
print("Right ascending node = ", ktab[-1] + 90 , "[degrees]")

#----------------------Plotting all graphs----------------------------------------------
im = plt.imread("EarthAlbedo.png")
fig = plt.figure(figsize=(15, 5)) 
gs = gridspec.GridSpec(1, 5, width_ratios=[10, 1, 1, 1, 3]) 
ax0 = plt.subplot(gs[0])
plt.title("ground track")
plt.axis('off')
ax0.plot(xtab, ytab, linewidth = 2, color="red")
ax0.plot(x2tab, y2tab, linewidth = 2, color="red")
plt.xlim(-20038, 20038)
plt.ylim(0,20038)
plt.imshow(im, aspect="auto", extent=(-20038, 20038, 0, 20038))
ax1 = plt.subplot(gs[1])
plt.title("albedo dataset")
plt.xlabel("Albedo [-]")
plt.ylabel("Time [s]")
ax1.plot(Albedo, xx - 1500)
plt.ylim(0, 5675)
ax2 = plt.subplot(gs[2])
plt.title("albedo map")
plt.xlabel("Albedo [-]")
plt.ylabel("Time [s]")
ax2.plot(albedotab, ttab)
plt.ylim(0, 5675)
ax3 = plt.subplot(gs[3])
plt.title("difference")
plt.xlabel("Albedo [-]")
plt.ylabel("Time [s]")
ax3.plot(difftab, ttab)
plt.ylim(0, 5675)
ax4 = plt.subplot(gs[4])
plt.title("difference per orbit")
plt.xlabel("right ascending node [degrees]")
plt.ylabel("average error [-]")
ax4.plot(kplottab, plotsum, color = "blue")
ax4.plot(kplot2tab, plot2sum, color = "blue")
plt.show()



