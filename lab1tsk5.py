import copy, math, cmath

def mult_matrix(A, B):
	C = [[0.0 for col in range(len(B[0]))] for row in range(len(A))]
	for i in range(len(A)):
		for j in range(len(B[0])):
			for k in range(len(B)):
				C[i][j] += A[i][k]*B[k][j]
	return C
	
def trans_matrix(M):
	n = len(M)
	return [[ M[i][j] for i in range(n)] for j in range(n)]

def norm(x):
    return math.sqrt(sum([x_i**2 for x_i in x]))

def Q_i(Q_min, i, j, k):
    if i < k or j < k:
        return float(i == j)
    else:
        return Q_min[i-k][j-k]

def QR_decomp(A):
    n = len(A)

    R = copy.deepcopy(A)
    Q = [[0.0] * n for i in range(n)]

    for k in range(n-1):                                                                    
        I = [[float(i == j) for i in range(n)] for j in range(n)]

		
        x = [row[k] for row in R[k:]]
        e = [row[k] for row in I[k:]]
        alpha = -cmp(x[0],0) * norm(x)

        u = map(lambda p,q: p + alpha * q, x, e)
        norm_u = norm(u)
        v = map(lambda p: p/norm_u, u)

        Q_min = [ [float(i==j) - 2.0 * v[i] * v[j] for i in range(n-k)] for j in range(n-k) ]

        Q_t = [[ Q_i(Q_min,i,j,k) for i in range(n)] for j in range(n)]

        if k == 0:
            Q = Q_t
            R = mult_matrix(Q_t,A)
        else:
            Q = mult_matrix(Q_t,Q)
            R = mult_matrix(Q_t,R)

    return trans_matrix(Q), R

def QR_solve(A, eps):
	n = len(A)
	M = copy.deepcopy(A)
	it = 0
	while True:
		it += 1
		Q, R = QR_decomp(M)
		M = mult_matrix(R, Q)
		La = [M[l][l] for l in range(n)]
		t = True
		i = -1
		while i < n-2:
			i += 1
			if math.sqrt(sum(M[i][j]**2 for j in range(i+1, n))) > eps:
				if math.sqrt(sum(M[i][j]**2 for j in range(i+2, n))) <= eps:
					a = 1
					b = -(M[i][i]+M[i+1][i+1])
					c = M[i][i]*M[i+1][i+1] - M[i][i+1]*M[i+1][i]
					disc = b**2 - 4 * a * c
					l1 = (-b + cmath.sqrt(disc))/ (2 * a)
					l2 = (-b - cmath.sqrt(disc))/ (2 * a)
					La[i] = l1
					La[i+1] = l2
					i += 1
				else:
					t = False
					break
		if t:
			return La, it
	
def main():
	print("Matrix A:")
	A = []
	A.append([int(j) for j in raw_input().strip().split(" ")])
	for i in range(1,len(A[0])) :
		A.append([int(j) for j in raw_input().strip().split(" ")])
		
	print("Epsilon:")
	eps = float(input())
	
	Lambdas, it = QR_solve(A,eps)
	print("Lambdas:")
	print(Lambdas) 
	print("Number of iterations:")
	print(it)

if __name__ == '__main__':
		main()
