# Runtime: about 3 minutes
# -------------------------Assumptions/changes---------------------
# correction on R near poles for globe to flat map
# no orbit shift (neglect j2) as data is 10 times the same

# Import modules
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
from Flux_intensity import xx, Albedo_true

# Data tables
dogleg = True
ktab = []
sumtab = []
plotsum = []
kplottab = []
plot2sum = []
kplot2tab = []

#load map and its data
img = Image.open('Albedomap_nov.png').convert('LA')
img.save('Albedomap_grey_nov.png')
FILENAME='Albedomap_grey_nov.png' #image can be in gif jpeg or png format 
im=Image.open(FILENAME).convert('RGB')
pix=im.load()

# Define longitude of ascending node and loop statements
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
        lattab = []
        difftab = []

        #Several runs over range of possible k at different precisions
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
            n = 88 #maximum visible range

        #initial values
        jump = False
        orbit = True
        V_y = 0
        i = 0
        t = 0
        pxold = 1000
        pyold = 1000
        v_absground = 7.062206527 #km/s
        inclination_orbit = 97.42 #degrees
        R_e = 6378.137
        longitude = k
        x = (longitude/180)*(20037.50834)
        x_old = x
        y = 18183.91
        y_old = y
        y_new = y

        #running loop for the orbit shape
        while y_new <= y:
            # Angle shift over earth map
            if abs(x - x_old) <= 20038:
                fact = (abs(1-(abs(x-x_old)/10019)))*(-inclination_orbit+180)
            else:
                fact = (abs(3-(abs(x-x_old)/10019)))*(-inclination_orbit+180)
            # Physical model for descending flight
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
            #Physical model for ascending flight
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
            pixel_x = int(1439*((x_new + 20038)/40076))
            pixel_y = int(719-719*(y_new/20038))

            # Read data from greyscale picture 
            if longrun == True:
                if pixel_x == pxold and pixel_y == pyold:
                    albedotab.append(albedotab[-1])
                else:
                    m =  n + 1
                    dattab = []
                    angtab = []
                    for f in range(-n, m):
                        for g in range( -abs(n - abs(f)), abs(m - abs(f) + 1)):
                            ang = cos(atan((sqrt(((27.831*g)**2)+((27.831*f)**2)))/505.1964)) + radians((sqrt(((27.831*g)**2)+((27.831*f)**2))/(2*pi*R_e))*360) 
                            angtab.append(ang)
                            if pixel_x + f <= 1439:
                                Pf = pixel_x + f
                            else:
                                pf = pixel_x + f - 1440

                            if pixel_y + g <= 719:
                                Py = pixel_y + g
                            else:
                                Py = pixel_y + g - 720
                            dattab.append(ang*1371*((pix[Pf, Py][0])/255))
                    albedotab.append(sum(dattab)/sum(angtab))
            else:
                albedotab.append(1372*((pix[pixel_x, pixel_y][0])/255))
            if t >= 3120 and t <= 5195:
                albedotab[-1] = 0

            #update time and y for next loop            
            ttab.append(t)
            lattab.append(latitude)
            t = t + (5/abs(v_abs))
            pxold = pixel_x
            pyold = pixel_y
        
        # Difference between data and this model
        for w in range(len(albedotab)):
            difftab.append(abs(Albedo_true[w+2119]-albedotab[w]))
        sumtab.append((sum(difftab)/len(difftab)))
        
        ktab.append(k)
        print(k)
        k = k + dk


print("---------------------------------------------------------")
print("Orbital period = ", ttab[-1], "[s]")
print("Average error = ", sumtab[-1], "[W/m^2]")
print("Right ascending node = ", ktab[-1] , "[degrees]")

#----------------------Plotting all graphs----------------------------------------------
im = plt.imread("Albedomap_nov.png")
fig = plt.figure(figsize=(10, 5)) 
gs = gridspec.GridSpec(3, 1, height_ratios=[10, 2, 4]) 
ax0 = plt.subplot(gs[0])
plt.title("Ground Track of the DST Satellite")
plt.axis('off')
ax0.plot(xtab, ytab, linewidth = 2, color="red")
ax0.plot(x2tab, y2tab, linewidth = 2, color="red")
plt.xlim(-20038, 20038)
plt.ylim(0,20038)
plt.imshow(im, aspect="auto", extent=(-20038, 20038, 0, 20038))
ax1 = plt.subplot(gs[2])
plt.title("Albedo Data Over One Orbit")
plt.ylabel(r"$Albedo\ [W/m^2]$")
plt.xlabel(r"$Time\ [s]$")
ax1.plot(xx - 1500, Albedo_true, color = "blue", label='Albedo from dataset', linestyle='dotted')
ax1.plot(ttab, albedotab, color = "red", label='Albedo from albedomap', linestyle='dashed')
plt.legend(loc='best')
plt.minorticks_on()
plt.grid(which='both', axis = 'both', color='#999999', linestyle='-', alpha=5)
plt.xlim(0, 5675)
"""
ax2 = plt.subplot(gs[4])
plt.title("Albedo retrieved from the orbit over albedo map")
plt.ylabel("Albedo [W/m^2]")
plt.xlabel("Time [s]")
ax2.plot(ttab, albedotab, color = "blue")
plt.xlim(0, 5675)
ax3 = plt.subplot(gs[6])
plt.title("Difference between the albedo of the dataset and map")
plt.ylabel("Albedo [W/m^2]")
plt.xlabel("Time [s]")
ax3.plot(ttab, difftab, color = "blue")
plt.xlim(0, 5675)
ax4 = plt.subplot(gs[8])
plt.title("difference per orbit")
plt.xlabel("right ascending node [degrees]")
plt.ylabel("average error [W/m^2]")
ax4.plot(kplottab, plotsum, color = "blue")
ax4.plot(kplot2tab, plot2sum, color = "blue")
"""
plt.show()



