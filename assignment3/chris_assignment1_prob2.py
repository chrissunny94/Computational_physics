#swap function , arguments A-matrix ,i-ith row , k-kth row , n -size of matrix
def swap(A,i,k):
        n = len(A)
        temp =0
        for j in range(n):
                temp=A[i][j]
                A[i][j]=A[k][j]
                A[k][j]=temp

def mult_matrix(M, N):
    """Multiply square matrices of same dimension M and N"""

    # Converts N into a list of tuples of columns                                                                                                                                                                                                      
    tuple_N = zip(*N)

    # Nested list comprehension to calculate matrix multiplication                                                                                                                                                                                     
    return [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in tuple_N] for row_m in M]

def pivot_matrix(M):
    """Returns the pivoting matrix for M, used in Doolittle's method."""
    m = len(M)

    # Create an identity matrix, with floating point values                                                                                                                                                                                            
    id_mat = [[float(i ==j) for i in xrange(m)] for j in xrange(m)]

    # Rearrange the identity matrix such that the largest element of                                                                                                                                                                                   
    # each column of M is placed on the diagonal of of M                                                                                                                                                                                               
    for j in xrange(m):
        row = max(xrange(j, m), key=lambda i: abs(M[i][j]))
        if j != row:
            # Swap the rows                                                                                                                                                                                                                            
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat

def lu_decomposition(A):
    """Performs an LU Decomposition of A (which must be square)                                                                                                                                                                                        
    into PA = LU. The function returns P, L and U."""
    n = len(A)

    # Create zero matrices for L and U                                                                                                                                                                                                                 
    L = [[0.0] * n for i in xrange(n)]
    U = [[0.0] * n for i in xrange(n)]

    # Create the pivot matrix P and the multipled matrix PA                                                                                                                                                                                            
    P = pivot_matrix(A)
    PA = mult_matrix(P, A)

    # Perform the LU Decomposition                                                                                                                                                                                                                     
    for j in xrange(n):
        # All diagonal entries of L are set to unity                                                                                                                                                                                                   
        L[j][j] = 1.0

        # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}                                                                                                                                                                                      
        for i in xrange(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
            U[i][j] = PA[i][j] - s1

        # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik} )                                                                                                                                                                  
        for i in xrange(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in xrange(j))
            L[i][j] = (PA[i][j] - s2) / U[j][j]

    return (P, L, U)

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
    from fractions import Fraction
    n = input()
    A = [[0 for j in range(n+1)] for i in range(n)]
    # Read input data
    for i in range(0, n):
        line = map(Fraction, raw_input().split(" "))
        for j, el in enumerate(line):
            A[i][j] = el
    raw_input()
    line = raw_input().split(" ")
    lastLine = map(Fraction, line)
    for i in range(0, n):
        A[i][n] = lastLine[i]
    # Print input
    pprint(A)
    # Calculate solution
    
    P,L,U = lu_decomposition(A)
    print "P:"
    pprint(P)
    print "L:"
    pprint(L)
    print "U"
    pprint(U)
    # Print result
    line = "Result:\t"
    for i in range(0, n):
        line += str(x[i]) + "\t"
    print(line)
