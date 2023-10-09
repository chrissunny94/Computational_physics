#swap function , arguments A-matrix ,i-ith row , k-kth row , n -size of matrix
def swap(A,i,k):
        n = len(A)
        temp =0
        for j in range(n):
                temp=A[i][j]
                A[i][j]=A[k][j]
                A[k][j]=temp
        
        

def gaussElimin(a,b):
        n = len(b)
        for k in range(0,n-1):
                for i in range(k+1,n):
                        if a[i,k] != 0.0:
                                lam = a [i,k]/a[k,k]
                                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                                b[i] = b[i] - lam*b[k]
        for k in range(n-1,-1,-1):
                b[k] = (b[k] - dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
        return b



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


def gauss(A):
    n = len(A)

    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k][i]/A[i][i]
            for j in range(i, n+1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x


if __name__ == "__main__":
    #from fractions import Fraction
    #n = input()
    #A = [[0 for j in range(n+1)] for i in range(n)]
    # Read input data
    #for i in range(0, n):
    #    line = map(Fraction, raw_input().split(" "))
    #    for j, el in enumerate(line):
    #        A[i][j] = el
    #raw_input()
    #line = raw_input().split(" ")
    #lastLine = map(Fraction, line)
    #for i in range(0, n):
    #    A[i][n] = lastLine[i]
    # Print input
    n = 3
    count = 100
    eps = -100.000000000000
    A=[[eps,-1,1,0],[-1,2,-1,0],[2,-2,0,0]]
    pprint(A)
    for i in range(count):	
            
            A=[[eps,-1,1,0],[-1,2,-1,0],[2,-2,0,0]]
            swap(A,0,2)
            
            # Calculate solution
            x = gauss(A)
            #pprint()    

            # Print result
            line = "Result:\t"
            for i in range(0, n):
                line += str(x[i]) + "\t"
            print(line)
            print ("eps=",eps)    
            eps = (eps - (eps/10))

    print ("PROGRAM ENDED")    
