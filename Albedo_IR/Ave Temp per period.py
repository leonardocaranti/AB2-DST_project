from TAS_Data_reading import choose_dataset
import matplotlib.pyplot as plt
import numpy as np


# ------------ Data handling ------------

data = choose_dataset('AB2_data/AB2/Mirror_segments.csv')

period = np.array([0,0])
for i in range(1,len(data[1])):
    if data[1][i][0] < 5674.57: # first period
        newrow = np.array([data[1][i][0],data[1][i][1]])
        period = np.vstack([period, newrow])
    else:
        points_T = i #data points per period
        break

period1 = data[1][1:points_T]
period2 = data[1][points_T:(2*points_T)]
period3 = data[1][(2*points_T):(3*points_T)]
period4 = data[1][(3*points_T):(4*points_T)]
period5 = data[1][(4*points_T):(5*points_T)]
period6 = data[1][(5*points_T):(6*points_T)]
period7 = data[1][(6*points_T):(7*points_T)]
period8 = data[1][(7*points_T):(8*points_T)]
period9 = data[1][(8*points_T):(9*points_T)]
period10 = data[1][(9*points_T):]

def find_Tave(set): #to find ave T
    sum = 0
    for i in range(0,len(set)):
        sum = sum + set[i][1]
        points = len(set)-1
    return (sum/points)

T_ave_points = np.array([find_Tave(period1),find_Tave(period2),find_Tave(period3),find_Tave(period4),find_Tave(period5),find_Tave(period6),find_Tave(period7), find_Tave(period8), find_Tave(period9), find_Tave(period10)])


N=0
f = []
grid = []
for i in range(0,len(T_ave_points)):
    T = T_ave_points[i]
    time = 2837 + (5674.57*N)
    #plt.plot(time,T,'or')
    f.append(T)
    grid.append(time)
    N = N+1
#plt.show()

f = np.array(f) + 273.15
x = np.array(grid)


# ------------ Regression: decaying exponential (y = C*exp(-lambda*x))------------

A = np.ones((len(x), 3))
A[:,1] = -x
A[:,2] = 1/2*x**2
a = np.dot(np.linalg.inv(np.dot(A.T, A)), np.dot(A.T, f))

_lambda = -a[2]/a[1]
C = -a[1]/_lambda
B = a[0] - C

def regress_func(y): return B + C*np.exp(_lambda*y)

rel_error, abs_error = np.sum(np.absolute(regress_func(x)-f) / f ), np.sum(np.absolute(regress_func(x)-f))

# ------------ Plotting ------------

xx = np.linspace(0, x[-1]*1.01, 1000)
plt.plot(x, f, marker='s', color = 'red', label='Average temperature per period')
plt.plot(xx, regress_func(xx), color = 'black', label=r'$ a_0 + a_1 \cdot e^{- \lambda t}$')
plt.legend(loc='best')
plt.title("Decaying exponential regression for average temperature per period")
plt.xlabel("Time [s]"); plt.ylabel("Temperature [K]")
plt.text(0, 272.4, r'$a_1 = $' + str(np.round(B,2)) + '\n' + r'$a_2 = $' + str(np.round(C,2)) + '\n' + \
    r'$\lambda$ = ' + str(np.round(-_lambda,7)) + '\n' + 'Relative error: ' + str(np.round(rel_error*100,2)) \
         + '%' + '\n' + 'Absolute error: ' + str(np.round(abs_error,2)))
plt.show()
