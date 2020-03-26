from TAS_Data_reading import choose_dataset
import matplotlib.pyplot as plt
import numpy as np


# -------------------------- Spline definition --------------------------
def S(x, xj, a, b, c, d):
    return a + b * (x-xj) + c * (x-xj)**2 + d * (x-xj)**3

def system_matrix(grid):
    N = len(grid)-1
    h = grid[1:] - grid[:-1]   # width of interval i - vector size N
    A = np.zeros((N+1,N+1))
    ### Equations describing conditions at interior nodes
    for i in range(1, N):
        A[i,i]   = 2*(h[i-1]+h[i])
        A[i,i-1] = h[i-1]
        A[i,i+1] = h[i]
    return A

def rhs(grid, fi):           ### RHS - interior points
    N = len(grid)-1
    h = grid[1:] - grid[:-1]   # width of interval i - vector size N
    rhs = np.zeros(N+1)
    for i in range(1, N):
        rhs[i] = 3.*(fi[i+1] - fi[i])/h[i] - 3.*(fi[i] - fi[i-1])/h[i-1]
    return rhs

def spline_natural(xi, fi, xx):
    """
    One-shot function for spline interpolation (with natural BCs).
    
    Args:
      xi (array, n+1): Sample locations
      fi (array, n+1): Sample values
      xx (array, M):   Reconstuction locations
    Return:
      ff (array, M): Reconstructed values at xx
    """
    h = xi[1:] - xi[:-1]       # Interval width
    N = len(h)
                               ### Setup system
    A = system_matrix(xi)      # Left-hand side 
    frhs = rhs(xi, fi)         # Right-hand side
    A[0,0] = A[N,N] = 1        # BC for LHS (natural)
    frhs[0] = 0                # BC for RHS (natural)
    frhs[-1] = 0 
                               ### Solve system for coefficients
    c = np.linalg.solve(A, frhs)
    a = fi[:]                                                  # N+1
    b = (a[1:] - a[:-1]) / h[:] - h[:]/3. * (2*c[:-1] + c[1:]) # N
    d = (c[1:] - c[:-1]) / (3. * h[:])                         # N
                               ### Reconstuct spline at locations xx  
    ii = np.digitize(xx, xi)   # Find to which interval each xx belongs
    ii = np.fmin(np.fmax(ii-1,0),N-1) 
    ff = np.zeros(xx.shape)
    for j, i in enumerate(ii): # Compute spline for each x
        ff[j] = S(xx[j], xi[i], a[i], b[i], c[i], d[i])
    return ff


# -------------------------- Spline implementation and plotting --------------------------
def spline():
    # Possible filenames: 'AB2_data/AB2/Baffle_surfaces.csv', 'AB2_data/AB2/Mirror_segments.csv', 'AB2_data/AB2/Secondary_mirror.csv'
    element_number = 0
    data = choose_dataset('AB2_data/AB2/Secondary_mirror.csv')[0]       # Might have to change directory based on where you saved file on own your computer
    xi, fi = np.array([data[i][element_number] for i in range(1, len(data))]), np.array([data[i][1] for i in range(1, len(data))])

    t_start, t_end, t_steps = 0, 15000, 10001         # Set boundary for graph and spline evaluation
    xx = np.linspace(t_start, t_end, t_steps)         # Set datapoints to evaluate splne at

    plt.title('Element ' + data[0]); plt.xlabel("Time [s]"); plt.ylabel("Temperature [°C]")
    plt.plot(xi, fi, label='Experimental data', marker='o'); plt.plot(xx, spline_natural(xi, fi, xx), label='Spline reconstruction')
    plt.legend(loc='best'); plt.xlim(t_start, t_end)
    plt.show()

# The following command just makes sure that whetever is inside of this is run only when this file is run. 
# For example, if I was to import this file (Interpolation.py) into another file, whetever is inside of this loop would not be run. 
if __name__ == '__main__': 
    spline()