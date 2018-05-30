import copy, math, cmath
import matplotlib.pyplot as plt
import numpy as np

def calcExpSpline(x0,x,A,B,C,D,p):
    for i in range(len(x)-1):
        if x0 >= x[i] and x0 <= x[i+1]:
            if p[i] == 0:
                return A[i]+B[i]*(x0-x[i]) + C[i]*((x0-x[i])**2) + D[i]*((x0-x[i])**3)
            else:
                return A[i]+B[i]*(x0-x[i]) + C[i]*math.exp(p[i]*(x0-x[i])) + D[i]*math.exp(-p[i]*(x0-x[i]))
            

def expSplineCoefs(x,y,dy0,dyn, p):
    n = len(x)
    A = [0.0]*(n-1)
    B = [0.0]*(n-1)
    C = [0.0]*(n-1)
    D = [0.0]*(n-1)
    dy = [0.0]*n; dy[0] = dy0; dy[-1] = dyn
    dk = [0.0]*n
    f = [0.0]*(n-1)
   
    if (n - 2) > 0:
        tmp1 = (y[1]-y[0])/(x[1]-x[0])
        for i in range(1,(n-1)):
            tmp2 = (y[i+1]-y[i])/(x[i+1]-x[i])
            dk[i] = tmp2 - tmp1
            tmp1 = tmp2
        dk[0] = (y[1]-y[0])/(x[1]-x[0]) - dy0
        dk[n-1] = dyn - (y[n-1]-y[n-2])/(x[n-1]-x[n-2])
        for i in range(0,n):
            if dk[i] == 0:
                return A,B,C,D,f
        
        L = [0.0]*(n-1)
        R = [0.0]*(n-1)
        l1 = 0; m1 = 0; r1 = 0
        l2 = 0; m2 = 0; r2 = 0
        for i in range(0,n-1):
            
            if p[i] == 0:
                l2 = 2.0/(x[i+1]-x[i])
                m2 = 2*l2
                r2 = (m2 + l2) * (y[i+1]-y[i])/(x[i+1]-x[i])
            else :
                zk = math.sinh(p[i]*(x[i+1]-x[i]))/(x[i+1]-x[i])
                vk = (p[i]*math.cosh(p[i]*(x[i+1]-x[i])) - zk)/(zk-p[i])
                tk = p[i]**2 * zk * (x[i+1]-x[i])/((vk**2-1)*(zk-p[i]))
                l2 = tk; m2 = tk*vk
                r2 = (vk + 1)*tk*(y[i+1]-y[i])/(x[i+1]-x[i])
            if i == 0:
                l1 = l2; m1 = m2; r1 = r2
                continue
            summ = m1 + m2; sumr = r1 + r2
            if i == 1:
                sumr -= l1 * dy0
            elif i == n-2:
                sumr -= l2 * dyn
            frac = 1.0/(summ - l1 * L[i-1])
            L[i] = - frac*l2
            R[i] = frac * (sumr - l1 * R[i-1])
            l1 = l2; m1 = m2; r1 = r2
            
        dy[n-2] = R[n-2]
        for i in range(n-3, 0, -1):
            dy[i] = dy[i+1] * L[i] + R[i]
        
    for i in range(0,n-1):
        alpha = y[i+1] / (x[i+1]-x[i])
        beta = y[i] / (x[i+1]-x[i])
        if p[i] == 0:
            gamma = (dy[i] + 2*dy[i+1] - 3.0*(y[i+1]-y[i])/(x[i+1]-x[i]))/3
            delta = -(2*dy[i] + dy[i+1] - 3.0*(y[i+1]-y[i])/(x[i+1]-x[i]))/3
            A[i] = (x[i+1]-x[i]) * beta
            B[i] = alpha - gamma - beta - 2*delta
            C[i] = 3*delta/(x[i+1]-x[i])
            D[i] = (gamma - delta)/((x[i+1]-x[i])**2)
        else:
            zk = math.sinh(p[i]*(x[i+1]-x[i]))/(x[i+1]-x[i])
            wk = -zk/(zk-p[i])
            vk = p[i]*math.cosh(p[i]*(x[i+1]-x[i]))/(zk-p[i]) + wk 
            gamma = (dy[i] + vk * dy[i+1] - (vk+1)*(y[i+1]-y[i])/(x[i+1]-x[i]))/(vk**2-1)
            delta = (-vk*dy[i] - dy[i+1] + (vk+1)*(y[i+1]-y[i])/(x[i+1]-x[i]))/(vk**2-1)
            A[i] = (x[i+1]-x[i]) * (beta + wk*delta)
            B[i] = alpha - beta + wk*(gamma - delta)
            C[i] = (gamma - delta*math.exp(-p[i]*(x[i+1]-x[i])))/(2*(zk-p[i]))
            D[i] = (delta*math.exp(p[i]*(x[i+1]-x[i])) - gamma )/(2*(zk-p[i]))
        
    for i in range(0,n-1):
        if dk[i]*dk[i+1] > 0:
            if D[i] == 0:
                f[i] = 0
                continue
            extr = 0
            if p[i] == 0:
                extr = -C[i]/(3*D[i]) + x[i]
            else:
                if (-C[i]/D[i]) <= 0:
                    f[i] = 0
                    continue
                extr = -math.log(-C[i]/D[i])/(2*p[i])+x[i]
                   
            if extr <= x[i+1] and extr >= x[i]:
                f[i] = 1
            else:
                f[i] = 0
        else:
            f[i] = 0
                
                
    return A,B,C,D,f
    
    
def main():
    
    
    x = [0.0, 0.9, 1.8, 2.7, 3.6]
    y = [0.0, 0.36892, 0.85408, 1.7856, 6.3138]
    dy0 = 0
    dyn = 10
    p = [5,6,3,0]
    '''
    x = [0.0, 1.0, 1.5, 2.5, 4.0, 4.5, 5.5, 6, 8, 10]
    y = [10, 8, 5, 4, 3.5, 3.4, 6, 7.1, 8, 8.5]
    dy0 = -1
    dyn = 0.5
    p = [4,4,6,12,12,0,1,4,8]
    '''
    A,B,C,D,f = expSplineCoefs(x,y,dy0,dyn, p)
    print(A)
    print(B)
    print(C)
    print(D)
    print(f)
    
    
    xnew = np.linspace(min(x),max(x),10000)
    ynew = [calcExpSpline(i,x,A,B,C,D,p) for i in xnew]
    plt.plot(x,y,'o',xnew,ynew)
    plt.grid(True)
    plt.show()
    #print("f(x*):")
    #print(calcExpSpline(1.5,x,A,B,C,D,p))
    


if __name__ == '__main__':
		main()
