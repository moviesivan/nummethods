import copy, math, cmath

def cubicroot(x):
	if x > 0:
		return math.pow(x, float(1)/3)
	elif x < 0:
		return -math.pow(abs(x), float(1)/3)
	else:
		return 0

def SimpleIter(f, x, q, eps):
	it = 0
	while True:
		prev_x = x
		x = f(x)

		it += 1

		if float(q) * abs(x - prev_x) / (1 - q) <= eps:
			return x, it
	

def phi(x):
	return cubicroot(-x**2 + x + 0.5)

def dphidx(x):
	return (-2.0/3.0 * x + 1.0 / 3) / cubicroot((-x**2 + x + 0.5)**2)

	
def main():
		
	print("Epsilon:")
	eps = float(input())
	q = abs(dphidx(1))
	
	x, it = SimpleIter(phi, 0.5, q, eps)
	print("X:")
	print(x) 
	print("Number of iterations:")
	print(it)

if __name__ == '__main__':
		main()
