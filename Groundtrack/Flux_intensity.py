import Interpolation
import numpy as np
import matplotlib.pyplot as plt

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

t_start, t_end, dt = 0, 25000, 0.70804741                           # Set boundary for graph and spline evaluation
t_steps = int(t_end/dt) + 1    
xx = np.linspace(t_start, t_end, t_steps)                     # Set datapoints to evaluate splne at


splines = np.zeros(shape = (8, t_steps))

for i in range(0, 8):
    filename, element_number = 'Raw_data/Mirror_segments.csv', i  #CHANGE THIS BASED ON THE LOCATION ON YOUR DEVICE!
    spline = Interpolation.spline(filename, element_number, xx, t_start, t_end, plotting = False)
    splines[i-2] = spline

T_front = (splines[0])/1 + 273.15
T_back = (splines[1])/1 + 273.15
T_avg = (splines[0] + splines[1])/2 + 273.15
T_E = np.mean(T_avg)
dT_dt = cd_deriv(xx, T_avg)[0]

T_front, T_back, T_avg, xx = T_front[1:-1], T_back[1:-1], T_avg[1:-1], xx[1:-1]
sigma = 5.670374*10**(-8)
e_front, e_back, A, m, C = 0.035, 0.05, 1.076/4, 6. , 690

Q_out , Q_balance = e_front * sigma * (T_front*T_front*T_front*T_front) + e_back * sigma * (T_back*T_back*T_back*T_back), 1/A * m * C * dT_dt
Q_in = (Q_out + Q_balance)/e_front

t_eclipse_1, t_eclipse_2 = 5000, 6500
avg = np.mean(Q_in[int(t_eclipse_1/dt):int(t_eclipse_2/dt)])
Albedo_true = (Q_in - avg)*(1360/1883.861)



def make_f_positive(f): #input: array of numbers
    for i in range(0,len(f)):
        if f[i] < 0: #changes negative values to positive
            f[i] = 0
    return f  #returns array with only postive or 0 values


Albedo_true = make_f_positive(Albedo_true)





