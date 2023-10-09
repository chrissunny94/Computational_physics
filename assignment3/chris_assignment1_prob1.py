import numpy as np
import swap
import error
#from LUpivot import *
#swap function , arguments A-matrix ,i-ith row , k-kth row , n -size of matrix
def matInv(a):
	n = len(a[0])
	aInv = np.identity(n)
	a,seq = LUdecomp(a)
	for i in range(n):
		aInv[:,i] = LUsolve(a,aInv[:,i],seq)
	return aInv

def LUsolve(a,b):
    n = len(a)
    for k in range(1,n):
        b[k] = b[k] - np.dot(a[k,0:k],b[0:k])
        b[n-1] = b[n-1]/a[n-1,n-1]
        for k in range(n-2,-1,-1):
            b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    return b

def LUsolve(a,b,seq):
	n = len(a)
	# Rearrange constant vector; store it in [x]
	x = b.copy()
	for i in range(n):
		x[i] = b[seq[i]]
	# Solution
	for k in range(1,n):
		x[k] = x[k] - np.dot(a[k,0:k],x[0:k])
	x[n-1] = x[n-1]/a[n-1,n-1]
	for k in range(n-2,-1,-1):
		x[k] = (x[k] - np.dot(a[k,k+1:n],x[k+1:n]))/a[k,k]
	return x

    
def LUdecomp(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a [i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
    return a

def LUdecomp(a,tol=1.0e-100):
	n = len(a)
	seq = np.array(range(n))
	# Set up scale factors
	s = np.zeros((n))
	for i in range(n):
		s[i] = max(abs(a[i,:]))
	for k in range(0,n-1):
		# Row interchange, if needed
		p = np.argmax(np.abs(a[k:n,k])/s[k:n]) + k
		if abs(a[p,k]) < tol: error.err('Matrix is singular')
		if p != k:
			swap.swapRows(s,k,p)
			swap.swapRows(a,k,p)
			swap.swapRows(seq,k,p)
	# Elimination
		for i in range(k+1,n):
			if a[i,k] != 0.0:
				lam = a[i,k]/a[k,k]
				a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
				a[i,k] = lam
	return a,seq




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
    a = np.array([\
    [6,-4,1],\
    [-3,2,5],\
    [6,-1,5]])
    A=[[6,-4,1],[-3,2,5],[6,-1,5]]
    A=np.array(A)
    print("A=")
    print(A)
    
    B=[[1,2,3],[4,5,6],[3,6,9]]
    B=np.array(B)
    #print("B=")
    #print(B)
    
    C=[[1,4],[9,1],[6,2]]
    C=np.array(C)
    print("C=")
    print(C)
    
    D,seq=LUdecomp(A,)
    print("Decomposition result ,\n")
    print(D)
    print("\nSeq")
    print(seq)
    
    print("\nSolution ,\n")
    print(LUsolve(D,C,seq))
    
    aInv = matInv(A)
    print("\naInv =\n")
    print(aInv)
    print("\nCheck: a*aInv =\n") 
    print(np.dot(A,aInv))
    
    print ("\n\nPROGRAM ENDED")    
