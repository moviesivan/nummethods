import copy, math, cmath
import matplotlib.pyplot as plt
import numpy as np

def n(j,xc,x):
    n = 1
    for i in range(j):
        n *= (xc-x[i])

    return n
 
def a(i,j,x,y):
    if j==0:
        return y[0]
    elif j-i==1 :
        return (y[i]-y[j])/(x[i]-x[j])
    else:
        return (a(i,j-1,x,y)-a(i+1,j,x,y))/(x[i]-x[j])

def newton(x,y,x0):
    z = 0
    for j in range(len(x)):
        z += a(0,j,x,y)*n(j,x0,x)
    
    return z
   
def main():
		
    print("task a)")
    x = [0,  math.pi/6, math.pi*2/6, math.pi*3/6]
    y = [math.cos(k) for k in x]
    
    xnew = np.linspace(min(x),max(x),100)
    ynew = [newton(x,y,i) for i in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    print("error:")
    print(abs(newton(x,y,math.pi*3/16) - math.cos(math.pi*3/16)))
    
    print("task b)")
    x = [0,  math.pi/8, math.pi/3, math.pi*3/8]
    y = [math.cos(k) for k in x]
    
    xnew = np.linspace(min(x),max(x),100)
    ynew = [newton(x,y,i) for i in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    print("error:")
    print(abs(newton(x,y,math.pi*3/16) - math.cos(math.pi*3/16)))


if __name__ == '__main__':
		main()
