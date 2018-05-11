import math

def shooting(dy, dz, z0, fy1, x0, x1, h, eps):
    nu1,nu2 = 1.5,1.7
    xlist, ylist, zlist = rugne_kutta(dy, dz, nu1, z0, x0, x1, h)
    y1 = fy1(zlist[-1])
    phi1 = abs(ylist[-1] - y1)
    if phi1 <= eps:
        return xlist, ylist
    xlist, ylist, zlist = rugne_kutta(dy, dz, nu2, z0, x0, x1, h)
    y1 = fy1(zlist[-1])
    phi2 = abs(ylist[-1] - y1)
    while phi2 > eps:
        nu1, nu2 = nu2, nu2 - phi2*(nu2-nu1)/(phi2-phi1)
        xlist, ylist, zlist = rugne_kutta(dy, dz, nu2, z0, x0, x1, h)
        y1 = fy1(zlist[-1])
        phi1, phi2 = phi2, abs(ylist[-1] - y1)
    return xlist, ylist
    
    

def  rugne_kutta(dy, dz, y0, z0, x0, x1, h):
    x = xprev = x0
    xlist = [x]
    yprev = y = y0
    ylist = [y]
    z = zprev = z0
    zlist = [z]
    while x < x1:
        xprev = x
        x = round(x+h,6)
        xlist.append(x)
        yprev = y
        zprev = z
        
        k1 = h*dy(zprev)
        l1 = h*dz(xprev,yprev)
        k2 = h*dy(zprev+1.0/2*l1)
        l2 = h*dz(xprev+1.0/2*h,yprev+1.0/2*k1)
        k3 = h*dy(zprev+1.0/2*l2)
        l3 = h*dz(xprev+1.0/2*h,yprev+1.0/2*k2)
        k4 = h*dy(zprev+l3)
        l4 = h*dz(xprev+h,yprev+k3)
        
        y = yprev + 1.0/6*(k1+2*k2+2*k3+k4)
        ylist.append(y)
        z = zprev + 1.0/6*(l1+2*l2+2*l3+l4)
        zlist.append(z)
    return xlist, ylist, zlist
    
    
    
def runge_romberg(h1data, h2data, p):
    h1x, h1y = h1data
    h2x, h2y = h2data
    res = 0
    i2 = 0
    for i in range(len(h1x)):
        tmp = abs((h1y[i]-h2y[i2])/(2**p-1))
        i2 += 2
        if tmp > res:
            res = tmp
    return res
    
    
def abserr(data, y):
    xlist, ylist = data
    res = 0
    for i in range(len(xlist)):
        tmp = abs(ylist[i]-y(xlist[i]))
        if tmp > res:
            res = tmp
    return res
    
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

    dy = (lambda z: z)
    dz = (lambda x,y: 2.0*y/(x**2*(x+1)))
    fy1 = (lambda z: 2+2*z)
    z0 = -1.0
    x0 = 1.0
    x1 = 2.0
    h1 = 0.1
    h2 = h1/2
    eps = 1e-5
    fy = (lambda x: 1/x+1)
    
    shootingdata = shooting(dy, dz, z0, fy1, x0, x1, h1, eps)
    print("shooting method:")
    for i in range(len(shootingdata[0])):
        print("x:   ", shootingdata[0][i], " y:  ", shootingdata[1][i])
    shootingdata2 = shooting(dy, dz, z0, fy1, x0, x1, h2, eps)
    shootingrungeerr = runge_romberg(shootingdata, shootingdata2, 1)
    print("calculated error: ",shootingrungeerr)
    shootingabserr = abserr(shootingdata, fy)
    print("absolute error: ",shootingabserr)
    
    xlist = shootingdata[0]
    qx = (lambda x: -2/(x**2*(x+1)))
    n = len(xlist)
    h = h1
    A = [[0 for _ in range(n)] for _ in range(n)]
    b = [0 for _ in range(n)]
    A[0][0] = -2 + h**2*qx(xlist[1]) + 1
    A[0][1] = 1
    blist = [-h]
    for i in range(1,n-1) :
        A[i][i-1] = 1
        A[i][i] = -2 + h**2*qx(xlist[i+1])
        A[i][i+1] = 1
        blist.append(0)
    A[n-1][n-2] = (2-h)
    A[n-1][n-1] = (2-h)*(-2+h**2*qx(xlist[-1])) + 2
    blist.append(2*h)
    ylist = TiagonalMatrixSolve(A, blist)
    print("mesh method:")
    for i in range(n):
        print("x:   ", xlist[i], " y:  ", ylist[i])
        
    xlist2 = shootingdata2[0]
    qx = (lambda x: -2/(x**2*(x+1)))
    n = len(xlist2)
    h = h2
    A = [[0 for _ in range(n)] for _ in range(n)]
    b = [0 for _ in range(n)]
    A[0][0] = -2 + h**2*qx(xlist2[1]) + 1
    A[0][1] = 1
    blist = [-h]
    for i in range(1,n-1) :
        A[i][i-1] = 1
        A[i][i] = -2 + h**2*qx(xlist2[i+1])
        A[i][i+1] = 1
        blist.append(0)
    A[n-1][n-2] = (2-h)
    A[n-1][n-1] = (2-h)*(-2+h**2*qx(xlist2[-1])) + 2
    blist.append(2*h)
    ylist2 = TiagonalMatrixSolve(A, blist)
    rungeerr = runge_romberg((xlist, ylist), (xlist2, ylist2), 1)
    print("calculated error: ",rungeerr)
    absolerr = abserr((xlist, ylist), fy)
    print("absolute error: ",absolerr)
        
        
    
    
    


if __name__ == '__main__':
		main()
