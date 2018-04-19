import copy, math, cmath

def Newton(f, dfdx, x, eps):
	f_value = f(x)
	it = 0
	while True:
		prev_x = x
		x = x - float(f_value)/dfdx(x)

		f_value = f(x)
		it += 1
		
		if abs(x - prev_x) < eps:
			return x, it
	

def f(x):
	return x**3 + x**2 - x - 0.5

def dfdx(x):
	return 3*x**2 + 2*x - 1

	
def main():
		
	print("Epsilon:")
	eps = float(input())
	
	x, it = Newton(f, dfdx, 0.5, eps)
	print("X:")
	print(x) 
	print("Number of iterations:")
	print(it)

if __name__ == '__main__':
		main()
