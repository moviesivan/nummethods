import math

def euler(dy, dz, y0, z0, x0, x1, h):
    x = xprev = x0
    xlist = [x]
    yprev = y = y0
    ylist = [y]
    z = zprev = z0
    while x < x1:
        xprev = x
        x = round(x+h,6)
        xlist.append(x)
        yprev = y
        zprev = z
        y = yprev + h*dy(zprev)
        ylist.append(y)
        z = zprev + h*dz(xprev,yprev)
    return xlist, ylist

def  rugne_kutta(dy, dz, y0, z0, x0, x1, h):
    x = xprev = x0
    xlist = [x]
    yprev = y = y0
    ylist = [y]
    z = zprev = z0
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
    return xlist, ylist
    
def adams(dy, dz, y0, z0, x0, x1, h):
    x = xprev = x0
    xlist = [x]
    yprev = y = y0
    ylist = [y]
    z = zprev = z0
    zlist = [z]
    while x < x1 and x < (x0+h*3):
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
        
    while x < x1:
        xprev = x
        x = round(x+h,6)
        xlist.append(x)
        yprev = y
        zprev = z
        
        y = yprev + h*(55*ylist[-1] - 59*ylist[-2] + 37*ylist[-3] - 9*ylist[-4])/24
        ylist.append(y)
        z = zprev + h*(55*zlist[-1] - 59*zlist[-2] + 37*zlist[-3] - 9*zlist[-4])/24
        zlist.append(z)
    return xlist, ylist
    
    
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
    
    
    # y' = z
    # z' = 2y+4x^2e^x^2
def main():

    dy = (lambda z: z)
    dz = (lambda x,y: 2*y+4*(x**2)*(math.exp(x**2)))
    y0 = 3.0
    z0 = 0.0
    x0 = 0.0
    x1 = 1.0
    h1 = 0.1
    h2 = h1/2
    fy = (lambda x: math.exp(x**2)+math.exp(x*math.sqrt(2))+math.exp(-x*math.sqrt(2)))
    
    eulerdata = euler(dy, dz, y0, z0, x0, x1, h1)
    print("euler method:")
    for i in range(len(eulerdata[0])):
        print("x:   ", eulerdata[0][i], " y:  ", eulerdata[1][i])
    eulrungeerr = runge_romberg(eulerdata, euler(dy, dz, y0, z0, x0, x1, h2), 1)
    print("calculated error: ",eulrungeerr)
    eulabserr = abserr(eulerdata, fy)
    print("absolute error: ",eulabserr)
    
    rungekuttadata = rugne_kutta(dy, dz, y0, z0, x0, x1, h1)
    print("rugne-kutta method:")
    for i in range(len(rungekuttadata[0])):
        print("x:   ", rungekuttadata[0][i], " y:  ", rungekuttadata[1][i])
    rungekuttarungeerr = runge_romberg(rungekuttadata, euler(dy, dz, y0, z0, x0, x1, h2), 2)
    print("calculated error: ",rungekuttarungeerr)
    rungekuttaabserr = abserr(rungekuttadata, fy)
    print("absolute error: ",rungekuttaabserr)
        
    adamsdata = adams(dy, dz, y0, z0, x0, x1, h1)
    print("adams method:")
    for i in range(len(adamsdata[0])):
        print("x:   ", adamsdata[0][i], " y:  ", adamsdata[1][i])
    adamsrungeerr = runge_romberg(adamsdata, adams(dy, dz, y0, z0, x0, x1, h2), 1)
    print("calculated error: ",adamsrungeerr)
    adamsabserr = abserr(adamsdata, fy)
    print("absolute error: ",adamsabserr)
    


if __name__ == '__main__':
		main()
