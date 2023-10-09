## module LUdecomp3
''' c,d,e = LUdecomp3(c,d,e).
LU decomposition of tridiagonal matrix [c\d\e]. On output
{c},{d} and {e} are the diagonals of the decomposed matrix.
x = LUsolve(c,d,e,b).
Solves [c\d\e]{x} = {b}, where {c}, {d} and {e} are the
vectors returned from LUdecomp3.
'''



import swap
import error
import numpy as np

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

 
def LUdecomp3(c,d,e):
    n = len(d)
    for k in range(1,n):
        lam = c[k-1]/d[k-1]
        d[k] = d[k] - lam*e[k-1]
        c[k-1] = lam
    return c,d,e
def LUsolve3(c,d,e,b):
    n = len(d)
    for k in range(1,n):
        b[k] = b[k] - c[k-1]*b[k-1]
    b[n-1] = b[n-1]/d[n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - e[k]*b[k+1])/d[k]
    return b
 
 ## module LUdecomp5
''' d,e,f = LUdecomp5(d,e,f).
LU decomposition of symmetric pentadiagonal matrix [a], where
{f}, {e} and {d} are the diagonals of [a]. On output
{d},{e} and {f} are the diagonals of the decomposed matrix.
x = LUsolve5(d,e,f,b).
Solves [a]{x} = {b}, where {d}, {e} and {f} are the vectors
returned from LUdecomp5.
'''
def LUdecomp5(d,e,f):
    n = len(d)
    for k in range(n-2):
        lam = e[k]/d[k]
        d[k+1] = d[k+1] - lam*e[k]
        e[k+1] = e[k+1] - lam*f[k]
        e[k] = lam
        lam = f[k]/d[k]
        d[k+2] = d[k+2] - lam*f[k]
        f[k] = lam
    lam = e[n-2]/d[n-2]
    d[n-1] = d[n-1] - lam*e[n-2]
    e[n-2] = lam
    return d,e,f
def LUsolve5(d,e,f,b):
    n = len(d)
    b[1] = b[1] - e[0]*b[0]
    for k in range(2,n):
        b[k] = b[k] - e[k-1]*b[k-1] - f[k-2]*b[k-2]
    b[n-1] = b[n-1]/d[n-1]
    b[n-2] = b[n-2]/d[n-2] - e[n-2]*b[n-1]
    for k in range(n-3,-1,-1):
        b[k] = b[k]/d[k] - e[k]*b[k+1] - f[k]*b[k+2]
    return b


