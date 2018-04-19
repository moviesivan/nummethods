import copy, math, cmath

def SimpleIter(phi1, phi2, x1, x2, q, eps):
	it = 0
	while True:
		prev_x1 = x1
		prev_x2 = x2
		x1 = phi1(x1, x2)
		x2 = phi2(x1, x2)

		it += 1

		if max(float(q) * abs(x1 - prev_x1) / (1 - q), float(q) * abs(x2 - prev_x2) / (1 - q)) <= eps:
			return x1, x2, it
	

def phi1(x1, x2):
	return math.sqrt(16 - (x2 - 2)**2) + 2

def phi2(x1, x2):
	return 64.0 / (x1**2 + 16)
	
def main():
		
	print("Epsilon:")
	eps = float(input())
	
	x1, x2, it = SimpleIter(phi1, phi2, 6.0, 2.0, 0.381, eps)
	print("X1: ")
	print(x1) 
	print("X2: ")
	print(x2) 
	print("Number of iterations:")
	print(it)

if __name__ == '__main__':
		main()
