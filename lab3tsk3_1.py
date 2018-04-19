import copy, math, cmath
import matplotlib.pyplot as plt
import numpy as np

def matrixmult(A, B):
    C = [[0.0 for col in range(len(B[0]))] for row in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                C[i][j] += A[i][k]*B[k][j]
    return C

def pivot_matrix(M):
    m = len(M)
    MCopy = copy.deepcopy(M)
                                                                                                                                                                                               
    id_mat = [[float(i==j) for i in range(m)] for j in range(m)]
    row_exchanges = 0
                                                                                                                                                                                              
    for i in range(m):
        maxElem = abs(MCopy[i][i])
        maxRow = i
        for k in range(i+1, m):
            if(abs(MCopy[k][i]) > maxElem):
                maxElem = abs(MCopy[k][i]) 
                maxRow = k
        if i != maxRow:                                                                                                                                                                                                                          
            id_mat[i], id_mat[maxRow] = id_mat[maxRow], id_mat[i]
            MCopy[i], MCopy[maxRow] = MCopy[maxRow], MCopy[i]
            row_exchanges += 1

    return id_mat, row_exchanges

def lup_decomposition(A):
    
    n = len(A)
                                                                                                                                                                                                               
    L = [[0.0] * n for i in range(n)]
    U = [[0.0] * n for i in range(n)]

                                                                                                                                                                                           
    P, rowExc = pivot_matrix(A)
    PA = matrixmult(P, A)
                                                                                                                                                                                                                     
    for j in range(n):
                                                                                                                                                                                                
        L[j][j] = 1.0

                                                                                                                                                                                      
        for i in range(j+1):
            s = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = PA[i][j] - s

        for i in range(j, n):
            s = sum(L[i][k] * U[k][j] for k in range(j))
            L[i][j] = (PA[i][j] - s) / U[j][j]

    return (P, L, U, rowExc)

def lup_solve(P,L,U,B):
    n = len(P)

    Bt = matrixmult(P, [[i] for i in B])

    Y = [0.0 for i in range(n)]
    for i in range(n):
        Y[i] = Bt[i][0]/L[i][i]
        for k in range(i):
            Y[i] -= Y[k]*L[i][k]

    X = [0.0 for i in range(n)]
    for i in range(n-1,-1,-1):
        s = sum(X[k]*U[i][k] for k in range(i+1,n))
        X[i] = (Y[i] - s)/U[i][i]

    return X
    
def lsm1(x,y):
    n = len(x)
    A = [[n, sum(x)],
         [sum(x), sum(el**2 for el in x)]]
    P, L, U, _ = lup_decomposition(A)
    
    ynew = [sum(y),sum(x[i]*y[i] for i in range(n))]
    a = lup_solve(P,L,U,ynew)
    
    return lambda x: a[0] + a[1]*x
    
def main():

    x = [-0.9, 0.0, 0.9, 1.8, 2.7, 3.6]
    y = [-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.3138]
    
    func = lsm1(x,y)
    xnew = np.linspace(min(x),max(x),100)
    ynew = [func(el) for el in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    print("squared errors:")
    print(sum((func(x[i])-y[i])**2 for i in range(len(x))))
    


if __name__ == '__main__':
		main()
