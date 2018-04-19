import copy, math, cmath
import matplotlib.pyplot as plt
import numpy as np

def TiagonalMatrixSolve(A, B):
    n = len(B) 
    L = [0 for _ in range(n)]
    F = [0 for _ in range(n)]
    L[0] = A[0][0]
    F[0] = B[0]
    for i in range(1, n):
        L[i] = A[i][i] - A[i][i-1]*A[i-1][i]/L[i-1]
        F[i] = B[i] - A[i][i-1]*F[i-1]/L[i-1]
    x = [0 for _ in range(n)]
    x[n-1] = F[n-1] / L[n-1]
    
    for i in range(n-2, -1, -1):
        x[i] = (F[i]-A[i][i+1]*x[i+1])/L[i]

    return x
    
def spline(x,y,x0):
    n=len(x)
    a = y[:]
    b = [0.0]*(n-1)
    d = [0.0]*(n-1)
    h = [x[i+1]-x[i] for i in range(n-1)]
    C = [[2*(h[0]+h[1]), h[1], 0],
         [h[1], 2*(h[1]+h[2]), h[3]],
         [0, h[2], 2*(h[2]+h[3])]]
    F = [3*((a[2]-a[1])/h[1] - (a[1]-a[0])/h[0]), 3*((a[3]-a[2])/h[2] - (a[2]-a[1])/h[1]), 3*((a[4]-a[3])/h[3] - (a[3]-a[2])/h[2])]
    c = [0] + TiagonalMatrixSolve(C, F)
    for j in range(n-2):
        b[j] = (a[j+1]-a[j])/h[j] - (h[j]*(c[j+1]+2*c[j]))/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
    b[n-2] = (a[n-1] - a[n-2])/h[n-2] - 2*h[n-2]*c[n-2]/3
    d[n-2] = -c[n-2]/(3*h[n-2])
    for i in range(n-1):
        if x0 >= x[i] and x0 <= x[i+1]:
            return a[i]+b[i]*(x0-x[i]) + c[i]*((x0-x[i])**2) + d[i]*((x0-x[i])**3)

def main():

    x = [0.0, 0.9, 1.8, 2.7, 3.6]
    y = [0.0, 0.36892, 0.85408, 1.7856, 6.3138]
    
    xnew = np.linspace(min(x),max(x),100)
    ynew = [spline(x,y,i) for i in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    print("f(x*):")
    print(spline(x,y,1.5))
    


if __name__ == '__main__':
		main()
