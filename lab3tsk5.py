import math

def rectangle(f, x0, x1, h):
    x = x0
    res = 0
    while x < x1:
        x += h
        res += h * f((2*x - h)/2)
    return res

def  trapeze(f, x0, x1, h):
    x = x0
    res = 0
    while x < x1:
        x += h
        res += h/2 * (f(x) + f(x-h))
    return res
    
def simson(f, x0, x1, h):
    xprev = x0
    xcur = x0 + h
    res = h/6*(f(xprev) + 4 *f((xprev + xcur)/2) + f(xcur))
    while xcur < x1:
        xprev = xcur
        xcur += h
        res += h/6*(f(xprev) + 4 *f((xprev + xcur)/2) + f(xcur))
    return res
    
def main():

    x0 = -1.0
    x1 = 1.0
    h1 = 0.5
    h2 = 0.25
    f = lambda x: x/((3*x+4)**3)
    
    rectres1 = rectangle(f,x0,x1,h1)
    recterr1 = (x1-x0)/24*h1**2*126 #126 - max|f''(x)| for x from -1 to 1
    print("rectangle method(h=0.5):")
    print(rectres1)
    print("error:")
    print(recterr1)
    
    rectres2 = rectangle(f,x0,x1,h2)
    recterr2 = (x1-x0)/24*h2**2*126 #126 - max|f''(x)| for x from -1 to 1
    print("rectangle method(h=0.25):")
    print(rectres2)
    print("error:")
    print(recterr2)
    
    trapres1 = trapeze(f,x0,x1,h1)
    traperr1 = (x1-x0)/12*h1**2*126 #126 - max|f''(x)| for x from -1 to 1
    print("trapeze method(h=0.5):")
    print(trapres1)
    print("error:")
    print(traperr1)
    
    trapres2 = trapeze(f,x0,x1,h2)
    traperr2 = (x1-x0)/12*h2**2*126 #126 - max|f''(x)| for x from -1 to 1
    print("trapeze method(h=0.5):")
    print(trapres2)
    print("error:")
    print(traperr2)
    
    simres1 = simson(f,x0,x1,h1)
    simerr1 = (x1-x0)/180*h1**4*35640 #35640 - max|f''''(x)| for x from -1 to 1
    print("simson method(h=0.5):")
    print(simres1)
    print("error:")
    print(simerr1)
    
    simres2 = simson(f,x0,x1,h2)
    simerr2 = (x1-x0)/180*h2**4*35640 #35640 - max|f''''(x)| for x from -1 to 1
    print("simson method(h=0.5):")
    print(simres2)
    print("error:")
    print(simerr2)
    
    


if __name__ == '__main__':
		main()
