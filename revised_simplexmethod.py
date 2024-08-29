import sys
import numpy as np
from numpy.linalg import inv

# Function to take the system as input
def take_input():

    # taking the number of variables and equations
    n = int(input("Enter the number of variables : "))
    m = int(input("Enter the number of equations : "))

    # the user enters the coefficient matrix A as an mxn matrix
    print("Enter the coefficient matrix A :")

    A = np.zeros((m, n))
    # takes each equation individually
    for i in range(m):
        A[i] = np.array([float(j) for j in input().split()])

    # the matrix with slack variables included
    AI = np.append(A, np.identity(m), axis=1)

    # taking the RHS values 
    print("Enter the RHS matrix b :")
    b = np.array([float(i) for i in input().split()])

    # taking the coefficients of actual variables
    print("Enter the coefficients of the objective function :")
    c = np.array([float(i) for i in input().split()])

    # C contains coefficients of actual and slack variables
    C = np.append(c, np.zeros(m))

    # bv stores the indices of basic variables
    # nbv stores the indices of non basic variables
    bv = np.array([int(i) for i in range(n, n+m)])
    nbv = np.array([int(i) for i in range(n)])

    # the constant value in Z
    const = int(input("Enter the constant value in the objective function : "))

    return A, AI, b, c, C, bv, nbv, const

# function to find the entering variable index and also check the 
# optimality of the solution at the present iteration 
def entering(A, c, C, bv, nbv, B):
    
    # C_b is the coefficient matrix of basic variables
    C_b = np.array([C[index] for index in bv])

    # mat1 is C_bxB (slack variable coefficients after an iteration) 
    # , mat2 is C_bxBxA - c (actual variable coefficients after the iteration) 
    mat_1 = np.matmul(C_b, B)
    mat_2 = np.subtract(np.matmul(mat_1, A), c)
    
    # we take the coefficients of non basic variables and take 
    # the entering variable as the one with the least coefficient
    nbv_coef = np.append(mat_2, mat_1)
    nbv_coef = np.array([nbv_coef[index] for index in nbv])
    
    entering_index = np.argmin(nbv_coef)
    
    # if the least coefficient is non negative, there can be no 
    # further optimisation and so is_optimal = 1
    is_optimal = 0
    if nbv_coef[entering_index] >= 0:
        is_optimal = 1
    else:
        is_optimal = 0

    return is_optimal, entering_index
    
# function to find the leaving index based on the fastest to reach zero 
def leaving(AI, bv, nbv, b, B, entering_index):

    # q stores the increase in the entering variable for the basic variable 
    # to reach zero, and we take the minimum of that as leaving
    q = []

    # finding the coefficients of the non basic variables in the 
    # updated set of equations and storing them in nbv_coefs
    A_nbv = np.transpose(np.array([AI[:, index] for index in nbv]))
    
    nbv_coefs = np.matmul(B, A_nbv)

    for i in range(len(bv)):
        if nbv_coefs[i][entering_index] <= 0:
            q.append(sys.maxsize)
        else:
            q.append(b[i]/nbv_coefs[i][entering_index])
            
    leaving_index = q.index(min(q))

    return leaving_index

# update the inverse of the basis using product form of inverse
def update_basis_inv(Basis_inv, AI, nbv, entering, leaving):

    m = AI.shape[0]

    r = leaving

    # Get the coefficient matrix on the last iteration, divider is
    # coefficient of entering variable at the place it is entering
    A_new = np.matmul(Basis_inv, AI)
    divider = A_new[r][nbv[entering]]
    
    # E is an identity matrix with the column at the leaving position changed
    E = np.identity(m)

    # E_ir : -a_ie/div for i !+ r
    # E_rr : 1/div
    for i in range(m):
        E[i][r] = -((A_new[i][nbv[entering]])/divider)

    E[r][r] = (1/divider)

    # inverse of basis is product of E and previous inverse
    Basis_inv = np.matmul(E, Basis_inv)

    return Basis_inv

# function to print the problem statement 
def print_problem(A, b, c, const):
    
    print("\n\nProblem statement :\n")

    # printing the objective function to optimise
    Z = "Z = "
    for i in range(len(c)):
        z = str(c[i]) + " x_" + str(i+1) + " + "
        Z = Z + z
    
    Z = Z + str(const)

    print("Maximise " + Z)
    print("Subject to constraints")

    # printing the constraint equations one by one
    for i in range(len(b)):
        eq = ""
        for j in range(len(c)-1):
            z = str(A[i][j]) + " x_" + str(j+1) + " + "
            eq = eq + z
        
        j = j + 1
        z = str(A[i][j]) + " x_" + str(j+1) + " "

        eq = eq + z + "<= "
        eq = eq + str(b[i])

        print(eq)
    
    print("\n")

# function to print the basis and the basic variables in each iteration
def print_iter(Basis_inv, X_b, bv, i):

    # Basis is inverse of Basis_inv
    Basis = inv(Basis_inv)

    print("Iteration", i, ":\n")

    # printing the basis matrix
    print("The basis matrix is :")
    for items in Basis:
        print("| ", end="")
        for item in items:
            print(item, end=" ")
        print(" |")
    print("\n")

    # printing the basic variables
    print("The basic variables are :")
    for i in range(len(bv)):
        print("X_" + str(bv[i]+1) + " = " + str(X_b[i]))
    
    print("\n")

# function to print the optimum solution 
def print_solution(c, X_b, bv, nbv, const):

    print("Optimal solution :\n")

    # X is the matrix [x1, x2, ....., x(n+m)] 
    X = np.zeros(len(bv)+len(nbv))
    
    # since X_b is the basic variables and everything 
    # else is zero, we set the bvs in X
    for i in range(len(bv)):
        X[bv[i]] = X_b[i]
    
    # change X to [x1, x2, ...., x(n)] (actual variables)
    X = np.array(X[0:len(nbv)])
    
    # objective function is multiplication of c, X and adding constant
    Z = np.matmul(c, X)
    
    Z = Z + const

    # print the optimum value of objective function 
    print("The maximum value of Z is " + str(Z))
    print("The values of x at which maximum occurs :")
    for i in range(len(X)):
        print("X_" + str(i+1) + " = " + str(X[i]))

def Lpp_solve(A, AI, b, c, C, bv, nbv, const):

    # finding the entering and leaving index after the basic solution is found
    is_optimal, entering_index = entering(A, c, C, bv, nbv, np.identity(len(bv)))
    leaving_index = leaving(AI, bv, nbv, b, np.identity(len(bv)), entering_index)

    X_b = b
    Basis_inv = np.identity(len(bv))

    i = 1
    while i:


        # we try to find the inverse of the basis at each step
        # if the basis does not have an inverse, the problem is unbounded
        try:

            # updating the Basis inverse using the product form of inverse
            Basis_inv = update_basis_inv(Basis_inv, AI, nbv, entering_index, leaving_index)
            
            # basic solution is matrix multiplication of B inverse and b at any step
            X_b = np.matmul(Basis_inv, b)

            # at each iteration swap the entering and leaving 
            # variables between bv and nbv
            temp = bv[leaving_index]
            bv[leaving_index] = nbv[entering_index]
            nbv[entering_index] = temp

            is_optimal, entering_index = entering(A, c, C, bv, nbv, Basis_inv)
            leaving_index = leaving(AI, bv, nbv, X_b, Basis_inv, entering_index)
        
            # if the solution is optimal, we print the solution and break the loop
            if is_optimal:
                print_iter(Basis_inv, X_b, bv, i)
                print_solution(c, X_b, bv, nbv, const)
                break
        
            # else continue the loop
            else:
                print_iter(Basis_inv, X_b, bv, i)
                i = i+1
                continue
        # in case of the inverse throwing an exception (meaning unbounded solution)
        except:
            print("This Lpp is unbounded")
            break
        

if __name__ == "__main__":

    # set the variables by taking input
    A, AI, b, c, C, bv, nbv, const = take_input()

    # print the problem
    print_problem(A, b, c, const)

    # pass the problem into the Lpp solver function
    Lpp_solve(A, AI, b, c, C, bv, nbv, const)

