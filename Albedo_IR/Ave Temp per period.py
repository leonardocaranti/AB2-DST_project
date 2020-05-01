from TAS_Data_reading import choose_dataset
import matplotlib.pyplot as plt
import numpy as np

data = choose_dataset('Raw_data/Mirror_segments.csv')
#print(data[1])
#I want element [1][i][1]
period = np.array([0,0])
for i in range(1,len(data[1])):
    if data[1][i][0] < 5674.57: # first period
        newrow = np.array([data[1][i][0],data[1][i][1]])
        period = np.vstack([period, newrow])
    else:
        points_T = i #data points per period
        break
#print(period)
#print(i)
#print(data[1][points_T:(2*points_T)])
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
#print(T_ave_points)
#print(type(T_ave_points))

for i in range(1,len(data[1])):
    plt.plot(data[1][i][0],data[1][i][1],'.k')

N=0
f = []
grid = []
for i in range(0,len(T_ave_points)):
    T = T_ave_points[i]
    time = 2837 + (5674.57*N)
    plt.plot(time,T,'or')
    f.append(T)
    grid.append(time)
    N = N+1
plt.show()

f = np.array(f)
grid = np.array(grid)
print()
print(f)
print(grid)

a = -75000
b = 135000
n = len(grid) -1
xx = np.linspace(a,b,101)

#now our data grid is the time/x-axis, f is data points

#POLY INTERPOLATION
'''
def basis_lagrange(x, grid):
    phi = np.ones(len(grid))
    for i, xi in enumerate(grid):
        for j,xj in enumerate(grid):
            if i != j:
                phi[i] *= (x - xj) / (xi - xj)
    return phi

def plot_basis(basisfn):
    xx = np.linspace(a, b, 101)
    phi = np.zeros((101, n+1))
    for i,x in enumerate(xx):
        phi[i] = basisfn(x,grid)
    for j in range(n+1):
        plt.plot(xx, phi[:,j])
    #plt.plot(grid, np.zeros(grid.shape), 'ok', label='samples')
    plt.xlabel(r'$x$'); plt.ylabel(r'$l_i$')
    plt.title('Lagrange basis')
    plt.legend()
    plt.show()
plot_basis(basis_lagrange)
def basis_newton(x, grid): #newton phi = PI (x - xj)
    phi = np.ones(len(grid))
    for i, xi in enumerate(grid):  #sets i to index and xi to item at that index in grid
        for j in range(i): #for j in range 0 to (index - 1)
            phi[i] *= (x - grid[j])
    return phi
def basis_monomial(x, grid): #monomial phi = x^i
    phi = np.zeros(len(grid))
    for i in range(len(grid)): phi[i] = x**i
    return phi


#left hand side of A*a = f
def interpolation_matrix(grid, basisfn):
    n = len(grid)-1
    A = np.zeros((n+1,n+1))
    for i in range(n+1):
        A[i,:] = basisfn(grid[i],grid)
    return A

A = interpolation_matrix(grid, basis_lagrange)
print(A)
plt.imshow(A, interpolation='none')
plt.xlabel(r'$j$'); plt.ylabel(r'$i$')
plt.title ('Visualization A')
plt.show()

#solve matrix
aa = np.linalg.solve(A, f)
print(aa)

def reconstruct(x, grid):
    return np.sum(aa * basis_lagrange(x, grid))

xx2 = np.linspace(a,b,101)
reconstruction = np.zeros(len(xx2))
for i,x in enumerate(xx):
    reconstruction[i] = reconstruct(x, grid)

#plt.plot(xx, f(xx), label='original')
plt.plot(grid, f, 'ok', label='samples')
plt.plot(xx, reconstruction, label='interpolation')
plt.xlabel(r'$x$')
plt.legend(loc='lower right')
plt.show()
#with poly interpolation i get a increasing function (eg cubic), same with lagrange, monomial and newton
'''

#RADIAL INTERPOLATION
def phi_invquad(r,l): return 1./(1 + (l*r)**2)
def phi_invmultiquad(r,l): return 1./np.sqrt(1 + (l*r)**2)
def phi_multiquad(r,l): return np.sqrt(1 + (l*r)**2)
def phi_linear(r,l): return r
def phi_gaussian(r,l): return np.exp(-(l*r*l*r))

N_1d = 10                         ### Number of samples
a,b = -50000,100000                        ### Interval of interest
#xi = np.linspace(a,b,N_1d+1)      ### Uniform sample locs
xi = grid
phi,l = phi_gaussian,0.00001  #choose basis funct and l
xx = np.linspace(a,b,501)         ### Sampling of x, for plotting

#euclidean distance (in 1D), sample/distance between two points
def dist_1d(x1, x2): return np.abs(x1-x2)  # this is r

for x in xi:
    plt.plot(xx, phi(dist_1d(x,xx),l))
plt.plot(xi,0.*xi,'ok')
plt.title('Basis functions based on data points')
plt.show()

#interpolation matrix A
def interpolation_matrix_1d(xi,phi,l):
    N = len(xi)-1
    A = np.zeros((N+1,N+1))
    for i in range(N+1):
        A[i,:] = phi(dist_1d(xi[i],xi),l)
    return A

A = interpolation_matrix_1d(grid,phi,l)
plt.imshow(A, interpolation='none')
plt.show()

#solve
coeffs = np.linalg.solve(A,f)
#print(coeffs)

def reconstruct(xx):
    out = np.zeros(xx.shape)
    for i,x in enumerate(xi):
         out += coeffs[i] * phi(dist_1d(x,xx),l)
    return out

plt.plot(xx, reconstruct(xx), label='RBF')
#plt.plot(xx, f_1d(xx), label='f exact')
plt.plot(xi, f, 'ok', label='samples')
plt.legend()
plt.title('Interpolation')
plt.show()


