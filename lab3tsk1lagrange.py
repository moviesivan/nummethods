import copy, math, cmath
import matplotlib.pyplot as plt
import numpy as np

def lagrange(x,y,x0):
    z=0
    for j in range(len(y)):
        w1=1.0; w2=1.0
        for i in range(len(x)):  
            if i != j: 
                w1*=(x0-x[i])
                w2*=(x[j]-x[i])
        z+=y[j]*w1/w2
    return z
	
def main():
		
    print("task a)")
    x = [0,  math.pi/6, math.pi*2/6, math.pi*3/6]
    y = [math.cos(k) for k in x]
    
    xnew = np.linspace(min(x),max(x),100)
    ynew = [lagrange(x,y,i) for i in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    print("error:")
    print(abs(lagrange(x,y,math.pi*3/16) - math.cos(math.pi*3/16)))
    
    print("task b)")
    x = [0,  math.pi/8, math.pi/3, math.pi*3/8]
    y = [math.cos(k) for k in x]
    
    xnew = np.linspace(min(x),max(x),100)
    ynew = [lagrange(x,y,i) for i in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    print("error:")
    print(abs(lagrange(x,y,math.pi*3/16) - math.cos(math.pi*3/16)))


if __name__ == '__main__':
		main()
