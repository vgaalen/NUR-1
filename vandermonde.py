import numpy as np
import matplotlib.pyplot as plt
from timeit import timeit

def LU_decomp(mat):
    """
    Perform a LU decomposition of a matrix.
    Both the L and U matrices are stored in the same matrix, with all values below the diagonal representing the U matrix
    and the diagonal and above representing the L matrix (with the U matrix's diagonal consisting of 1's and the other missing values being 0's).
    """
    mat[1:,0] = mat[1:,0] / mat[0,0]
    for j in range(1,mat.shape[1]):
        for i in range(1,mat.shape[0]):
            if i <= j:
                mat[i,j] = mat[i,j] - np.sum(np.multiply(mat[i,:i],mat[:i,j]))
            else:
                mat[i,j] = (mat[i,j] - np.sum(np.multiply(mat[i,:j],mat[:j,j]))) / mat[j,j]
    return mat

def LU_solve(LU, y):
    """
    Solve a system of linear equations using the LU decomposition of a matrix.
    """
    #y[0] = y[0]
    for i in range(1,len(y)):
        y[i] = (y[i] - np.sum(np.multiply(LU[i,:i],y[:i])))
    
    y[-1] = y[-1]/LU[-1,-1]
    for i in range(len(y)-2, -1, -1):
        y[i] = 1/LU[i,i]*(y[i] - np.sum(np.multiply(LU[i,i+1:],y[i+1:])))
    
    return y

def TwoA(data, silent=True):
    ##############
    ##### 2a #####
    ##############
    
    # Convert the dataset into a Vandermonde matrix
    data.astype(np.float128)
    matrix = np.repeat([data[:,0]], data.shape[0], axis=0)
    power = np.arange(data.shape[0], dtype=np.float128)
    powers = np.repeat([power], data.shape[0], axis=0)
    matrix = np.power(matrix, powers.T).astype(np.float128)

    # Perform a LU decomposition and use it to solve the system of linear equations
    LU = LU_decomp(matrix.copy().T)
    c = LU_solve(LU, data[:,1].astype(np.float128).copy())

    # Print (if desired) and write to file
    if not silent:
        print("Solution found (using LU Decomposition):")
        print(c)
    np.savetxt('vandermonde_coefficients.txt', c, delimiter=' ')

    # Determine the polynomial found at 1000 points between the minimum and maximum x values
    x = np.linspace(np.min(data[:,0]), np.max(data[:,0]), 1000)
    powers = np.repeat([power], 1000, axis=0)
    xx = np.repeat([x], data.shape[0], axis=0)
    cc = np.repeat([c], 1000, axis=0)
    y = np.sum( np.multiply(cc.T, np.power(xx, powers.T)), axis=0 )

    # Plot the data and the polynomial
    fig, ax = plt.subplots(2, 2, figsize=(10,5), dpi=192, gridspec_kw={'height_ratios': [6, 3]})
    ax[0,0].plot(data[:,0], data[:,1], 'ro', label="datapoints")
    ax[0,1].plot(data[:,0], data[:,1], 'ro', label="datapoints")
    ax[0,0].plot(x, y, 'b', label='LU decomposition')
    ax[0,1].plot(x, y, 'b', label='LU decomposition')

    # Plot the deviation of the polynomial from the datapoints
    dy = np.zeros((data.shape[0]))
    for i in range(data.shape[0]):
        dy[i] = np.abs(data[i,1]-np.sum( np.multiply(c.T, np.power(data[i,0], power))))
    ax[1,0].plot(data[:,0], dy, label='LU')
    ax[1,1].plot(data[:,0], dy, label='LU')

    ax[1,0].set_xlabel('x')
    ax[1,1].set_xlabel('x')
    ax[0,0].set_ylabel('y')
    ax[1,0].set_ylabel('|y(x)-y_0|')
    fig.suptitle('Solution to Vandermonde Matrix')
    ax[0,1].set_ylim([np.min(data[:,1])-50, np.max(data[:,1])+50])
    ax[1,1].set_yscale('log')
    ax[0,0].legend()
    ax[0,1].legend()
    #plt.yscale('log')
    fig.tight_layout()
    fig.savefig('plots/vandermonde.png')
    return fig, ax



def bisection(x, x0, num=2):
    """
    Bisection algorithm for finding the index of the reference data point that x is between

    Parameters
    ----------
    x : float
        The x value to interpolate at
    x0 : float
        The x values of the reference data points
    """
    if x < x0[0]:
        return [0]
    elif x > x0[-1]:
        return [len(x0)-1]
    indices = range(len(x0))
    while len(indices) > 2:
        #print(len(indices))
        if x > x0[indices[len(indices)//2]]:
            indices = indices[len(indices)//2:]
        else:
            indices = indices[:len(indices)//2+1]
    if num==2:
        return indices
    elif num==1:
        closest = np.argmin(np.abs(x-x0[indices]))
        return indices[closest]

def neville(x, x0, y0, ord=None, fill_value=None):
    """
    Neville's algorithm for polynomial interpolation

    Parameters
    ----------
    x : float
        The x value to interpolate at
    x0 : float
        The x values of the reference data points
    y0 : float
        The y values of the reference data points
    ord : int, optional
        The order of the polynomial to interpolate with. If None, the order is set to len(x0)-1
    fill_value : float, optional
        The value to return if x is outside the range of x0. If None, an error is raised.
    """

    if ord==None:
        ord = len(x0)-1
    # if x < np.min(x0) or x > np.max(x0):
    #     if fill_value is not None:
    #         return fill_value
    #     else:
    #         raise ValueError('x is out of range')
    
    p = np.zeros((ord,ord))
    x_i = np.zeros((ord))

    for i in range(ord):
        closest = bisection(x, x0, num=1)
        x_i[i] = x0[closest]
        x0 = np.delete(x0,closest,0)
        p[i,i] = y0[closest]
        y0 = np.delete(y0,closest,0)
    
    for k in range(1,ord):
        for i in range(ord-k):
            j = i+k
            p[i,j] = ( (x-x_i[i])*p[i+1,j] - (x-x_i[j])*p[i,j-1] ) / (x_i[j]-x_i[i])
    return p[0,ord-1]



def TwoB(data, fig=None, ax=None):
    ##############
    ##### 2b #####
    ##############

    # Use nevilles algorithm to interpolate the data
    x = np.linspace(np.min(data[:,0]), np.max(data[:,0]), 1000)
    y_neville = np.zeros((len(x)))
    for i in range(len(x)):
        y_neville[i] = neville(x[i], data[:,0], data[:,1])

    if fig is not None and ax is not None:
        ax[0,0].plot(x, y_neville, 'g', label="Neville's Algorithm")
        ax[0,1].plot(x, y_neville, 'g', label="Neville's Algorithm")

        # Plot the deviation from the datapoints
        dy = np.zeros((data.shape[0]))
        for i in range(data.shape[0]):
            dy[i] = np.abs(data[i,1]-neville(data[i,0], data[:,0], data[:,1]))
        ax[1,0].plot(data[:,0], dy, label='Neville')
        ax[1,1].plot(data[:,0], dy, label='Neville')

        fig.suptitle("Neville's Algorithm")
        ax[0,0].legend()
        ax[0,1].legend()
        ax[1,0].legend()
        ax[1,1].legend()
        fig.tight_layout()
        plt.savefig('plots/neville_compare.png')
        plt.close()

def TwoC(data):
    ##############
    ##### 2c #####
    ##############

    # Iterative Improvements on LU decomposition
    # first itteration gives c' in V@c=y with c'=c+dc, where dc is the error
    # second itteration we solve for V@dc=y-V@c', subtracting dc from c' then gives c''
    # and so on

    matrix = np.repeat([data[:,0]], data.shape[0], axis=0)
    power = np.arange(data.shape[0], dtype=np.float128)
    powers = np.repeat([power], data.shape[0], axis=0)
    matrix = np.power(matrix, powers.T).astype(np.float128)

    LU = LU_decomp(matrix.copy().T)
    c = LU_solve(LU, data[:,1].astype(np.float128).copy())

    LU1 = c - LU_solve(LU, matrix.T@c-data[:,1].copy())

    LU10 = LU1
    for i in range(9):
        LU10 = LU10 - LU_solve(LU, matrix.T@LU10-data[:,1].copy())

    fig, ax = plt.subplots(2, 2, figsize=(10,5), dpi=192)
    ax[0,0].plot(data[:,0], data[:,1], 'ro', label="datapoints")
    ax[0,1].plot(data[:,0], data[:,1], 'ro', label="datapoints")
    
    x = np.linspace(np.min(data[:,0]), np.max(data[:,0]), 1000)
    x_matrix = np.repeat([x], data.shape[0], axis=0)

    power = np.arange(data.shape[0], dtype=np.float128)
    powers = np.repeat([power], x.shape[0], axis=0)
    x_matrix = np.power(x_matrix, powers.T).astype(np.float128)

    y = x_matrix.T@c
    ax[0,0].plot(x,y, label='LU0')
    ax[0,1].plot(x,y, label='LU0')

    y = x_matrix.T@LU1
    ax[0,0].plot(x, y, label='LU1')
    ax[0,1].plot(x, y, label='LU1')

    y = x_matrix.T@LU10
    ax[0,0].plot(x, y, label='LU10')
    ax[0,1].plot(x, y, label='LU10')

    val_LU0 = matrix.T@c
    val_LU1 = matrix.T@LU1
    val_LU10 = matrix.T@LU10
    text = []

    dy_LU0 = np.zeros((data.shape[0]))
    dy_LU1 = np.zeros((data.shape[0]))
    dy_LU10 = np.zeros((data.shape[0]))
    for i in range(data.shape[0]):
        dy_LU0[i] = np.abs(data[i,1]-val_LU0[i])
        dy_LU1[i] = np.abs(data[i,1]-val_LU1[i])
        dy_LU10[i] = np.abs(data[i,1]-val_LU10[i])
    ax[1,0].plot(data[:,0], dy_LU0, label='LU0')
    ax[1,1].plot(data[:,0], dy_LU0, label='LU0')
    ax[1,0].plot(data[:,0], dy_LU1, label='LU1')
    ax[1,1].plot(data[:,0], dy_LU1, label='LU1')
    ax[1,0].plot(data[:,0], dy_LU10, label='LU10')
    ax[1,1].plot(data[:,0], dy_LU10, label='LU10')

    ax[1,0].set_xlabel('x')
    ax[1,1].set_xlabel('x')
    ax[0,0].set_ylabel('y')
    ax[1,0].set_ylabel('y')
    fig.suptitle('Itterating on Vandermonde Matrix Solution')
    ax[0,1].set_ylim([np.min(data[:,1])-50, np.max(data[:,1])+50])
    ax[0,0].legend()
    ax[0,1].legend()
    ax[1,0].legend()
    ax[1,1].legend()
    ax[1,0].set_yscale('log')
    ax[1,1].set_yscale('log')
    fig.tight_layout()
    fig.savefig('plots/vandermonde_itt.png')

    plt.close()


if __name__=='__main__':
    data = np.loadtxt("Vandermonde.txt", dtype=np.float64)
    fig, ax = TwoA(data, silent=False)
    TwoB(data, fig=fig, ax=ax)
    TwoC(data)

    ##############
    ##### 2d #####
    ##############
    t_2a = timeit('TwoA(data)', globals=globals(), number=10)
    t_2b = timeit('TwoB(data)', globals=globals(), number=10)
    t_2c = timeit('TwoC(data)', globals=globals(), number=10)

    with open('vandermonde_timing.txt', 'w') as f:
        f.write(f"2a: {t_2a:.2e} s\n2b: {t_2b:.2e} s\n2c: {t_2c:.2e} s\n")