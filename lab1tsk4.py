import copy, math

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
	
def check_symmetry(M):
	return M == trans_matrix(M)
	
def givens_rot(M, eps):
	n = len(M)
	A = copy.deepcopy(M)
	it = 0
	Ux = [[float(i == j) for j in range(n)] for i in range(n)]
	while True:
		it += 1
		a, i, j = 0.0, 0, 0
		for l in range(n):
			for m in range(l):
				if abs(A[l][m]) > a:
					a, i, j = abs(A[l][m]), l, m
		phi = 0.5 * (math.atan2(2*A[i][j], (A[i][i]-A[j][j]))) if A[i][i]-A[j][j] != 0 else math.pi / 4
		U = [[float(l == m) for m in range(n)] for l in range(n)]
		U[i][i] = math.cos(phi)
		U[j][j] = U[i][i]
		U[j][i] = math.sin(phi)
		U[i][j] = - U[j][i]
		A = mult_matrix(trans_matrix(U), mult_matrix(A, U))
		Ux = mult_matrix(Ux, U)
		t = math.sqrt(sum(sum(A[l][m]**2 for m in range(l)) for l in range(n)))
		if t <= eps :
			return [A[l][l] for l in range(n)], trans_matrix(Ux), it

def main():
	print("Matrix A:")
	A = []
	A.append([int(j) for j in raw_input().strip().split(" ")])
	for i in range(1,len(A[0])) :
		A.append([int(j) for j in raw_input().strip().split(" ")])
	
	if not(check_symmetry(A)):
		print("not solvable")
		return
		
	print("Epsilon:")
	eps = float(input())
	
	Lambdas, Ux, it = givens_rot(A,eps)
	print("Lambdas:")
	print(Lambdas)
	print("All X:")
	for elem in Ux:
		print(elem) 
	print("Number of iterations:")
	print(it)

if __name__ == '__main__':
		main()
