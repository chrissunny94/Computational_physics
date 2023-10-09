import numpy as np
from numpy.random import rand
from LUdecomp import *
#from LUpivot import *
#swap function , arguments A-matrix ,i-ith row , k-kth row , n -size of matrix
import math

def inversePower(a,s,tol=1.0e-6):
    n = len(a)
    aStar = a - np.identity(n)*s # Form [a*] = [a] - s[I]
    aStar = LUdecomp(aStar) # Decompose [a*]
    x = np.zeros(n)
    for i in range(n):
        # Seed [x] with random numbers
        x[i] = np.random()
        xMag = math.sqrt(np.dot(x,x)) # Normalize [x]
        x =x/xMag
    for i in range(50): # Begin iterations
        xOld = x.copy() # Save current [x]
    x = LUsolve(aStar,x)
    # Solve [a*][x] = [xOld]
    xMag = math.sqrt(np.dot(x,x)) # Normalize [x]
    x = x/xMag
    if np.dot(xOld,x) < 0.0:
        # Detect change in sign of [x]
        sign = -1.0
        x = -x
    else: sign = 1.0
    if math.sqrt(np.dot(xOld - x,xOld - x)) < tol:
        return s + sign/xMag,x
    print('Inverse power method did not converge')
        
def inversePower5(Bv,d,e,f,tol=1.0e-6):
    n = len(d)
    d,e,f = LUdecomp5(d,e,f)
    x = rand(n) # Seed x with random numbers
    xMag = math.sqrt(np.dot(x,x)) # Normalize {x}
    x = x/xMag
    for i in range(30):
        # Begin iterations
        xOld = x.copy() # Save current {x}
        x = Bv(xOld) # Compute [B]{x}
        x = LUsolve5(d,e,f,x) # Solve [A]{z} = [B]{x}
        xMag = math.sqrt(np.dot(x,x))
        # Normalize {z}
        x = x/xMag
        if np.dot(xOld,x) < 0.0:
            sign = -1.0
            x = -x
        else: sign = 1.0
    if math.sqrt(np.dot(xOld - x,xOld - x)) < tol:
        return sign/xMag,x
    print('Inverse power method did not converge')  

def Bv(v):
    # Compute {z} = [B]{v}
    n = len(v)
    z = np.zeros(n)
    z[0] = 2.0*v[0] - v[1]
    for i in range(1,n-1):
        z[i] = -v[i-1] + 2.0*v[i] - v[i+1]
    z[n-1] = -v[n-2] + 2.0*v[n-1]
    return z

def pprint(A):
    n = len(A)
    for i in range(0, n):
        line = ""
        for j in range(0, n+1):
            line += str(A[i][j]) + "\t"
            if j == n-1:
                line += "| "
        print(line)
    print("")

if __name__ == "__main__":
    A = np.array([\
    [1/3,-1/3,0],\
    [-1/3,4/3,-1],\
    [0,-1,2]])
    
    
    n = 3 # Number of interior nodes
    d = np.ones(n)*6.0 # Specify diagonals of [A] = [f\e\d\e\f]
    d[0] = 5.0
    d[n-1] = 7.0
    e = np.ones(n-1)*(-4.0)
    f = np.ones(n-2)*1.0
    print(inversePower(A,1,))    
    lam,x = inversePower5(Bv,d,e,f)
    print("PL",lam*(n+1)**2)
    input("\nPress return to exit")
    
    print ("\n\nPROGRAM ENDED")    
