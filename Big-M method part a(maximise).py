# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import sys
import numpy as np

M = 1000
alternate = 0
unbounded = 0

def TakeInput():
    # taking the number of variables and equations
    n = int(input("Enter the No. variables : "))
    m = int(input("Enter the No. constraints : "))
    l = int(input("Enter No. '<=' constraints: "))
    g = int(input("Enter No. '>=' constraints: "))

    # the user enters the coefficient matrix A as an mxn matrix
    print("Enter the coefficient matrix A :")

    A = np.zeros((m, n))
    # takes each equation individually
    for i in range(m):
        A[i] = np.array([float(j) for j in input().split()])

    # initialise I with coefficients of slack, surplus and additional variables for each constraint
    I = np.zeros((m, g))
    for i in range(m):
        for j in range(g):
            if i >= m - g and i == m-g+j:
                I[i][j] = -1

    # append I to A
    AI = np.append(A, I, axis=1)
    AI = np.append(AI, np.identity(m), axis=1)

    # taking constraint types for each constraint as ('>=': 1, '=': 0, '<=': -1)
    print("enter constraint types (for '>=': 1, '=': 0, '<=': -1): ")
    CT = np.array([int(i) for i in input().split()])

    # taking the RHS values
    print("Enter the RHS matrix b :")
    b = np.array([float(i) for i in input().split()])

    # taking the coefficients of actual variables
    print("Enter the coefficients of the objective function :")
    c = np.array([float(i) for i in input().split()])
    c = np.array([i for i in c])

    # the constant value in Z
    const = int(input("Enter the constant value in the objective function : "))
    const = const

    # bv stores the indices of basic variables
    bv = np.array([int(i) for i in range(n+g, n+m+g)])

    # C contains coefficients of actual and slack variables
    C = np.append(c, np.zeros(l+g))
    temp = np.zeros(m-l)
    for i in range(m-l):
        temp[i] = float(-M)
    C = np.append(C, temp)

    return A, AI, CT, b, c, C, bv, const


def print_problem(A, CT, b, c, const):
    print("\n\nProblem statement")

    # printing the objective function to optimise
    Z = "Z = "
    for i in range(len(c)):
        z = str(c[i]) + " x_" + str(i + 1) + " + "
        Z = Z + z

    Z = Z + str(const)

    print("Maximise " + Z)
    print("Subject to constraints")

    # printing the constraint equations one by one
    for i in range(len(b)):
        eq = ""
        for j in range(len(c) - 1):
            z = str(A[i][j]) + " x_" + str(j + 1) + " + "
            eq = eq + z

        j = j + 1
        z = str(A[i][j]) + " x_" + str(j + 1) + " "

        eq = eq + z
        switcher = {
            -1: "<=",
            0: "=",
            1: ">="
        }
        eq = eq + switcher.get(CT[i])
        eq = eq + str(b[i])

        print(eq)

    print("\n")

def opt_check(AI, C, bv, b):
    global alternate, unbounded
    # C_b is the coefficient matrix of basic variables
    C_b = np.array([C[index] for index in bv])

    mat = np.subtract(np.matmul(C_b,AI), C)

    entering_index = np.argmin(mat)

    # if the least coefficient is non negative, there can be no
    # further optimisation and so is_optimal = 1
    is_optimal = 0
    if mat[entering_index] >= 0:
        is_optimal = 1
    else:
        is_optimal = 0

    q = []
    min = sys.maxsize
    leaving_index = -1
    # finding the leaving variable coefficient by checking
    # the least b[i]/a[i][e]
    for i in range(len(bv)):
        if b[i] > 0 and AI[i][entering_index] > 0:
            val = b[i]/AI[i][entering_index]
            if val < min:
                min = val
                leaving_index = i

    if leaving_index == -1:
        unbounded = 1
        print("LPP is UNBOUNDED")

    return is_optimal, entering_index, leaving_index


def print_solution(c, C, X_b, bv, const):
    # X is the matrix [x1, x2, ....., x(n+m)]
    X = np.zeros(len(C))

    # since X_b is the basic variables and everything
    # else is zero, we set the bvs in X
    for i in range(len(bv)):
        X[bv[i]] = X_b[i]

    # change X to [x1, x2, ...., x(n)] (actual variables)
    X = np.array(X[0:len(c)])

    # objective function is multiplication of c, X and adding constant
    Z = np.matmul(c, X)

    Z = Z + const

    # print the optimum value of objective function
    print("The maximumum value of Z is " + str(Z))
    print("The values of x at which minimum occurs :")
    print(X)

def Big_M(AI, b, c, C, bv, const):
    # finding the entering and leaving index after the basic solution is found
    is_optimal, entering_index, leaving_index = opt_check(AI, C, bv, b)

    while is_optimal == 0 and unbounded == 0:
        k_e = AI[leaving_index][entering_index]
        for i in range(len(C)):
            AI[leaving_index][i] /= k_e
        b[leaving_index] /= k_e

        for i in range(len(bv)):
            if i != leaving_index:
                b[i] -= AI[i][entering_index] * b[leaving_index]
                AI[i, :] -= (AI[i, entering_index] * AI[leaving_index, :])

        bv[leaving_index] = entering_index
        is_optimal, entering_index, leaving_index = opt_check(AI, C, bv, b)

    if is_optimal:
        print_solution(c, C, b, bv, const)


if __name__ == "__main__":

    # set the variables by taking input
    A, AI, CT, b, c, C, bv, const = TakeInput()

    # print the problem
    print_problem(A, CT, b, c, const)

    # pass the problem into the Lpp solver function
    Big_M(AI, b, c, C, bv, const)

