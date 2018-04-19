import copy

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
	
def main():
	print("Dimension of A:")
	n = int(input())
	print("Diagonal elements of A:")
	A = [[0 for _ in range(n)] for _ in range(n)]
	tmplist = [float(j) for j in raw_input().strip().split(" ")]
	A[0][0] = tmplist[0]
	A[0][1] = tmplist[1]
	for i in range(1,n-1) :
		tmplist = [float(j) for j in raw_input().strip().split(" ")]
		A[i][i-1] = tmplist[0]
		A[i][i] = tmplist[1]
		A[i][i+1] = tmplist[2]
	tmplist = [float(j) for j in raw_input().strip().split(" ")]
	A[n-1][n-2] = tmplist[0]
	A[n-1][n-1] = tmplist[1]

	print("Vector b:")
	b = [float(j) for j in raw_input().strip().split(" ")]

	x = TiagonalMatrixSolve(A,b)
	print("X:")
	print(x)
	
if __name__ == '__main__':
		main()
