import copy

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

def lup_invert(P,L,U):
    n = len(P)

    IA = [[float(i==j) for i in range(n)] for j in range(n)]
    for i in range(n) :
        b = [IA[i][k] for k in range(n)]
        IA[i] = lup_solve(P,L,U,b)
    for i in range(n):
        for j in range(i):
            IA[i][j],IA[j][i] = IA[j][i],IA[i][j]
	
    return IA

def lup_determinant(U, rowExc):
    n = len(U[0])

    det = U[0][0]
    for j in range(1,n):
        det *= U[j][j]

    det *= (-1)**rowExc

    return det
	
def main():
    print("Matrix A:")
    A = []
    A.append([int(j) for j in raw_input().strip().split(" ")])
    for i in range(1,len(A[0])) :
        A.append([int(j) for j in raw_input().strip().split(" ")])

    print("Vector b:")
    b = [int(j) for j in raw_input().strip().split(" ")]



    P, L, U, rowExc = lup_decomposition(A)
    print("Matrix P:")
    for elem in P:
        print(elem)
    print("Matrix L:")
    for elem in L:
        print(elem)
    print("Matrix U:")
    for elem in U:
        print(elem)
        
    x = lup_solve(P,L,U,b)
    print("X:")
    print(x)

    IA = lup_invert(P,L,U)
    print("Inverted A:")
    for elem in IA:
        print(elem)

    detA = lup_determinant(U, rowExc)
    print("Determinant of A: %d" % detA)
	
if __name__ == '__main__':
		main()


    
