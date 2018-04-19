def main():

    x = [1.0, 1.5, 2.0, 2.5, 3.0]
    y = [0.0, 0.40547, 0.69315, 0.91629, 1.0986]
    xp = 2.0
    
    yd = (y[2] - y[1])/(x[2]-x[1])
    print("left-hand derivative 1:")
    print(yd)
    
    yd = (y[3] - y[2])/(x[3]-x[2])
    print("right-hand derivative 1:")
    print(yd)
    
    yd = (y[2] - y[1])/(x[2]-x[1]) + (((y[3]-y[2])/(x[3]-x[2])-(y[2]-y[1])/(x[2]-x[1]))/(x[3]-x[1]))*(2*xp-x[1]-x[2])
    print("derivative 2:")
    print(yd)
    
    yd = 2*(((y[3]-y[2])/(x[3]-x[2])-(y[2]-y[1])/(x[2]-x[1]))/(x[3]-x[1]))
    print("second derivative:")
    print(yd)
    


if __name__ == '__main__':
		main()
