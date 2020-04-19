import Interpolation
import numpy as np
import matplotlib.pyplot as plt

# ------------------------- Period calculations -------------------------------

for i in range(8):

    filename, element_number = 'AB2_data/AB2/Baffle_surfaces.csv', i

    t_start, t_end, dt = 0, 20000, 0.5                            # Set boundary for graph and spline evaluation
    t_steps = int(t_end/dt) + 1    
    xx = np.linspace(t_start, t_end, t_steps)                     # Set datapoints to evaluate splne at
    spline = Interpolation.spline(filename, element_number, xx, t_start, t_end, plotting = False)


    spline_copy = list(spline)
    plt.title("Spline reconstruction"); plt.xlim(t_start, t_end); plt.xlabel("Time [s]"); plt.ylabel("Temperature [°C]")
    plt.plot(xx, spline, label='Spline reconstruction'); plt.show()

    print('Inputs for period (maximum 1 boundaries)')
    s_low  = int(input("Lower boundary [s]:"))
    s_high = int(input("Upper boundary [s]:"))
    i_low, i_high = int(s_low/dt), int(s_high/dt)
    i_T_max_1 = spline_copy.index(max(spline[i_low:i_high]))
    print()

    print('Inputs for eclipse period (minimum boundaries)')
    s_low  = int(input("Lower boundary [s]:"))
    s_high = int(input("Upper boundary [s]:"))
    i_low, i_high = int(s_low/dt), int(s_high/dt)
    i_T_min = spline_copy.index(min(spline[i_low:i_high]))
    print()

    print('Inputs for period (maximum 2 boundaries)')
    s_low  = int(input("Lower boundary [s]:"))
    s_high = int(input("Upper boundary [s]:"))
    i_low, i_high = int(s_low/dt), int(s_high/dt)
    i_T_max_2 = spline_copy.index(max(spline[i_low:i_high]))
    print()

    period = (i_T_max_2 - i_T_max_1)*dt
    eclipse_per = (i_T_min - i_T_max_1)*dt
    sunlight_per = period - eclipse_per
    print('Period [s]:', period)
    print('Eclipse period [s]:', eclipse_per)
    print('Sunlight period [s]:', sunlight_per)
    print()
    print()
