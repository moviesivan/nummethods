import copy, math, cmath

def Newton(f1, f2, df1dx1, df2dx1, df1dx2, df2dx2, x1, x2, eps):
	it = 0
	while True:
		prev_x1 = x1
		prev_x2 = x2
		
		detJ = df1dx1(prev_x1, prev_x2) * df2dx2(prev_x1, prev_x2) - df1dx2(prev_x1, prev_x2) * df2dx1(prev_x1, prev_x2)
		detA1 = f1(prev_x1, prev_x2) * df2dx2(prev_x1, prev_x2) - df1dx2(prev_x1, prev_x2) * f2(prev_x1, prev_x2)
		detA2 = df1dx1(prev_x1, prev_x2) * f2(prev_x1, prev_x2) - f1(prev_x1, prev_x2) * df2dx1(prev_x1, prev_x2)
		
		x1 = prev_x1 - float(detA1)/detJ
		x2 = prev_x2 - float(detA2)/detJ

		it += 1
		
		if max(abs(x1 - prev_x1), abs(x2 - prev_x2))  < eps:
			return x1, x2, it
	

def f1(x1, x2):
	return (x1**2 + 4**2) * x2 - 4**3
	
def f2(x1, x2):
	return (x1 - 2)**2 + (x2 - 2)**2 - 4**2

def df1dx1(x1, x2):
	return 2*x1*x2

def df2dx1(x1, x2):
	return 2*(x1-2)

def df1dx2(x1, x2):
	return x1**2 + 16

def df2dx2(x1, x2):
	return 2*(x2-2)

	
def main():
		
	print("Epsilon:")
	eps = float(input())
	
	x1, x2, it = Newton(f1, f2, df1dx1, df2dx1, df1dx2, df2dx2, 6.0, 2.0, eps)
	print("X1: ")
	print(x1) 
	print("X2: ")
	print(x2) 
	print("Number of iterations:")
	print(it)

if __name__ == '__main__':
		main()
