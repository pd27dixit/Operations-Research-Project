import numpy as np
import math

INFINITY = 10e9

M = 10000

def getInput():
    print("Enter the optimization type(min/max):")
    optimization_type = str(input())
    print("Enter the number of variables in the objective function:")
    n = int(input())
    print("Enter the number of constraints:")
    m = int(input())
    print("Enter the coefficients of the Objective Function:")
    c = [float(i) for i in input().split(" ")]
    if optimization_type == 'min':
        c = [-z for z in c]
    print("Enter the value of the constant in the objective function:")
    d = float(input())
    print("Enter the matrix A which is the coefficient of the constraints row by row:")
    a = []
    for i in range(m):
        a.append([float(j) for j in input().split(" ")])
    print("Enter the constants/RHS bi for the constraint equations:")
    b = [float(i) for i in input().split(" ")]
    print("Enter the equation or in-equation type fo the variable(=/<=/>=):")
    s = input().split(" ")
    for i in range(len(b)):
        if b[i] >= 0:
            continue
        b[i] *= -1
        if s[i] == '>=':
            s[i] = '<='
        elif s[i] == '<=':
            s[i] = '>='
        for j in range(n):
            a[i][j] *= -1
    return n, m, a, np.array(c), d, b, s, optimization_type


def print_standard_form(n, m, a, c, d, b, s, M, optimization_type):
    No_artf_var = 0
    No_surplus_var = 0
    No_slack_var = 0
    k = n + 1
    x = []
    add_coloum = []
    total_var = []
    nbv = []
    bv = []
    sv = []
    av = []
    suv = []
    a = np.array(a)
    b = np.array(b)
    objective_fun = ""
    for i, j in enumerate(c):
        objective_fun += str(j) + "x" + str(i + 1)
        objective_fun += " + "
    size = len(c)
    objective_fun += str(d)
    print("The objective function to maximize is: \n" + objective_fun)
    print("\nThe standard form of the constraints using slack/surplus variables is:")
    k = n + 1
    x = []
    add_coloum = []
    for i in range(n):
        nbv.append("x" + str(i + 1))
        x.append(0)
    for j, row in enumerate(a):
        constraint = ""
        for i, q in enumerate(row):
            constraint += str(q) + "x" + str(i + 1) + " + " * (i != len(row) - 1)
        if s[j] == "=":
            constraint += " + x" + str(k)
            bv.append("x" + str(k))
            x.append(b[j])
            k += 1
            No_artf_var += 1
        elif s[j] == "<=":
            constraint += " + x" + str(k)
            bv.append("x" + str(k))
            x.append(b[j])
            k += 1
            No_slack_var += 1
        else:
            constraint += " - x" + str(k) + " + x" + str(k + 1)
            bv.append("x" + str(k + 1))
            x.append(0)
            nbv.append("x" + str(k))
            x.append(b[j])
            No_surplus_var += 1
            temp_col = np.zeros((m,))
            temp_col[j] = -1
            add_coloum.append(temp_col)
            No_artf_var += 1
            k += 2
        constraint += " = " + str(b[j])
        print(constraint)

    if len(nbv) - size > 0:
        additional = np.array(add_coloum).transpose()
        a = np.concatenate([a, additional], axis=1)

    total_var = nbv + bv
    cn = np.zeros((len(nbv)))
    cb = np.zeros((len(bv)))
    for i in range(len(c)):
        cn[i] = c[i]
    for i, var in enumerate(bv):
        if s[i] == '<=':
            cb[i] = 0
        else:
            cb[i] = -M
    zj_cj = np.empty((len(nbv)))
    print(a)
    print(cb)
    print(cn)
    print(nbv)
    print(bv)
    for i in range(len(nbv)):
        zj_cj[i] = 0
        zj_cj[i] += (np.dot(np.transpose(cb), a[:, i]) - cn[i])
    print("The objective function to maximize is")
    objective_fun = ""
    for i, var in enumerate(nbv):
        objective_fun += "-" + str(-cn[i]) + str(var) + " "
    print(objective_fun + " + ", d)
    return a, b, cb, cn, zj_cj, np.array(
        x), nbv, bv, total_var, av, sv, suv


def valueObjectiveFunction(cb, b, d, optimization_type):
    if optimization_type == 'max':
        return np.dot(np.transpose(cb), b) + d
    else:
        return -np.dot(np.transpose(cb), b) + d


def simplex(a, c, cn, cb, d, b, nbv, bv, total_var,
            x, optimization_type, No_surplus_var, No_slack_var, No_artf_var):
    count = 0
    N = a.shape[1]
    M = a.shape[0]
    while True:
        count += 1
        print("-" * 60)
        print("Iteration ", count)
        print("The value of Objective function for this iteration is ",
              round(valueObjectiveFunction(cb, b, d, optimization_type), 5))
        print("The values of x are:")
        l_ = []
        for i in range(len(x)):
            l_.append("x" + str(i + 1) + '=' + str(x[i]))
        print(l_)
        print("list of basic variables : ", bv)
        print("list of non-basic variables : ", nbv)
        print("The matrix A is")
        print(a)
        print("The value of Zj-Cj are")
        print(c)
        print("The values of basic solutions X_b are")
        print(b)
        v = np.argmin(c)
        cv = c[v]
        if cv >= 0:
            if count > No_artf_var:
                print("-" * 60)
                print("The iterations have ended")
                print("list of basic variables : ", nbv)
                print("list of non-basic variables : ", bv)
                print("The values of x are:")
                li = []
                for i in range(len(x)):
                    li.append("x" + str(i + 1) + '=' + str(round(x[i], 5)))
                print(li)
                print("Final value of objective function is: ",
                      round(valueObjectiveFunction(cb, b, d, optimization_type), 5))
                return a, b, cb, cn, c, x, nbv, bv, total_var, sv
            else:
                print("-" * 60)
                print("Infeasible solution because all artificial variables are not zero")
                print("list of basic variables : ", nbv)
                print("list of non-basic variables : ", bv)
                print("The values of x are:")
                li = []
                for i in range(len(x)):
                    li.append("x" + str(i + 1) + '=' + str(round(x[i], 5)))
                print(li)
                print("value of objective function is :", round(valueObjectiveFunction(cb, b, d, optimization_type), 5))
                return a, b, cb, cn, c, x, nbv, bv, total_var, sv
        else:
            print("The value of most negative c is", cv, " of respective column", v + 1)
            rs = np.zeros((M,))
            u = 0
            pivot = -1
            min_ratio = INFINITY
            for i in range(M):
                if a[i][v] == 0:
                    continue
                rs[i] = b[i] / a[i][v]
                if min_ratio > rs[i] >= 0:
                    u = i
                    min_ratio = rs[i]
                    pivot = a[u][v]
            print("ratios are for corresponding column", rs)
            print("minimum ratio is:", min_ratio)
            print("pivot element is ", pivot, " and respective coordinates(1 based indexing) is", u + 1, " ",
                  v + 1)
            if min_ratio == INFINITY:
                print("-" * 60)
                print("problem is unbounded")
                print("value of objective function is:", round(valueObjectiveFunction(cb, b, d, optimization_type), 5))
                print("values of x are:", x)
                return a, b, cb, cn, c, x, nbv, bv, total_var, sv
            A = np.empty((M, N))
            for i in range(M):
                for j in range(N):
                    if i == u and j == v:
                        A[i][j] = 1 / a[i][j]
                    elif i == u:
                        A[i][j] = a[i][j] / pivot
                    elif j == v:
                        A[i][j] = -a[i][j] / pivot
                    else:
                        A[i][j] = (pivot * a[i][j] - a[i][v] * a[u][j]) / pivot
            C = np.copy(c)
            for j in range(N):
                if j == v:
                    C[j] = round(-c[j] / pivot, 6)
                else:
                    C[j] = round((pivot * c[j] - c[v] * a[u][j]) / pivot, 6)
            B = np.copy(b)
            for j in range(M):
                if j == u:
                    B[j] = b[j] / pivot
                else:
                    B[j] = (pivot * b[j] - a[j][v] * b[u]) / pivot
            a = np.copy(A)
            c = np.copy(C)
            b = np.copy(B)
            temp2 = cb[u]
            cb[u] = cn[v]
            cn[v] = temp2
            temp1 = nbv[v]
            nbv[v] = bv[u]
            bv[u] = temp1
            for i in range(len(bv)):
                s = bv[i]
                x[int(s[1:]) - 1] = b[i]
            for i in range(len(nbv)):
                s = nbv[i]
                x[int(s[1:]) - 1] = 0


def getFractionalPart(x):
    if x == 0:
        return 0
    elif x > 0:
        return x - int(x)
    else:
        return 1 + x - int(x)


def cuttingPlane(a, b, cb, cn, c, x, nbv, bv, total_var,
                  sv):
    n = len(nbv)
    m = len(bv)
    sv.append('x'+str(n+m+1))
    total_var.append('x'+str(n+m+1))
    bv.append('x'+str(n+m+1))
    i = -1
    beta_i = -INFINITY
    for j in range(len(b)):
        f = b[j] - int(b[j])
        if beta_i < f:
            beta_i = f
            i = j
    add_row = np.empty((1, n))
    for item in range(n):
        f_i = getFractionalPart(a[i][item])
        add_row[0][item] = -f_i
    B = np.empty((len(bv)))
    for item in range(len(b)):
        B[item] = b[item]
    B[len(b)] = -getFractionalPart(b[i])
    X = list(x)
    X.append(-getFractionalPart(b[i]))
    x = np.array(X)
    b = np.copy(B)
    cb_list = cb.tolist()
    cb_list.append(0)
    cb = np.array(cb_list)
    a = np.concatenate([a, add_row], axis=0)
    return a, b, cb, cn, c, x, nbv, bv, total_var, sv


def print_table(a, c, d, b, cb, x, nbv, bv, optimization_type):
    l_ = []
    for basic in nbv:
        l_.append(basic)
    print("Non Basic Variables = ", l_)
    l_ = []
    for non_basic in bv:
        l_.append(non_basic)
    print("Basic Variables = ", l_)
    print("matrix A: ")
    print(a)
    print("Xb: ")
    print(b)
    print("'zj-cj' : ")
    print(c)
    print("Cb (coefficients of basic variable): ")
    print(cb)
    print("Cn (coefficients of non-basic variable): ")
    print(cn)
    print("So the optimal value of objective function is:", round(valueObjectiveFunction(cb, b, d, optimization_type), 5))


def print_pivot(cv, v, u, rs, min_ratio, pivot):
    print("The value of c is", -cv, " corresponding to column", v + 1)
    print("The ratios for corresponding rows", rs)
    print("minimum ratio is:", min_ratio)
    print("pivot element is ", pivot, " and corresponding coordinates is", u + 1, " ", v + 1)


def dualSimplex(a, c, cn, cb, d, b, nbv, bv, total_var,
                 x, optimization_type, No_surplus_var, No_slack_var, No_artf_var):
    N = a.shape[1]
    M = a.shape[0]
    count = 0
    while True:
        count += 1
        if count > 10:
            print(
                "Table has repeated. Due to this there is infinite iterations. Hence cannot solve by dual simplex.")
            exit(0)
        print("-" * 60)
        print("Iteration ", count)
        print_table(a, c, d, b, cb, x, nbv, bv,
                    optimization_type)
        print("The values of x are:")
        l_ = []
        for i in range(len(x)):
            l_.append("x" + str(i + 1) + '=' + str(x[i]))
        print(l_)
        u = np.argmin(b)
        cv = b[u]
        if cv >= 0:
            if count > No_artf_var:
                print("*" * 60)
                print("The iterations has ended because no negative values of Xb are present")
                print("list of non-basic variables : ", nbv)
                print("list of basic variables : ", bv)
                print("The values of x are:")
                l_ = []
                for i in range(len(x)):
                    l_.append("x" + str(i + 1) + '=' + str(x[i]))
                print(l_)
                print("So the Final value of objective function is:",
                      round(valueObjectiveFunction(cb, b, d, optimization_type), 5))
                return a, b, cb, cn, c, x, nbv, bv, total_var, sv
            else:
                print("*" * 60)
                print("Infeasible solution because all artificial variables are not zero")
                print("list of non-basic variables : ", nbv)
                print("list of basic variables : ", bv)
                print("The values for x are:")
                li = []
                for i in range(len(x)):
                    li.append("x" + str(i + 1))
                print(li)
                print("As you can see the the value of objective function is :",
                      round(valueObjectiveFunction(cb, b, d, optimization_type), 5))
                return a, b, cb, cn, c, x, nbv, bv, total_var, sv
        else:
            print("The value of most negative Xb is", cv, " Corresponding to row", u + 1)
            rs = np.empty((N,))
            v = 0
            pivot = -1
            min_ratio = INFINITY
            for i in range(N):
                if a[u][i] == 0:
                    continue
                rs[i] = abs(c[i] / a[u][i])
                if min_ratio > rs[i] >= 0:
                    v = i
                    min_ratio = rs[i]
                    pivot = a[u][v]
            print("The rs are for corresponding row", rs)
            print("The minimum ratio is:", min_ratio)
            print("The pivot element is ", pivot, " and corresponding coordinates(1 based indexing) is", u + 1, " ",
                  v + 1)
            if min_ratio == INFINITY:
                print("-" * 60)
                print("The problem is unbounded")
                print("The value of objective function is:", valueObjectiveFunction(cb, b, d, optimization_type))
                print("The values for x are:", x)
                return a, b, cb, cn, c, x, nbv, bv, total_var, sv
            A = np.empty((M, N))
            for i in range(M):
                for j in range(N):
                    if i == u and j == v:
                        A[i][j] = 1
                    elif i == u:
                        A[i][j] = round(a[i][j] / pivot, 6)
                    elif j == v:
                        A[i][j] = 0
                    else:
                        A[i][j] = round((pivot * a[i][j] - a[i][v] * a[u][j]) / pivot, 6)
            C = np.copy(c)
            for j in range(N):
                if j == v:
                    C[j] = 0
                else:
                    C[j] = round((pivot * c[j] - c[v] * a[u][j]) / pivot, 6)
            B = np.copy(b)
            for j in range(M):
                if j == u:
                    B[j] = round(b[j] / pivot)
                else:
                    B[j] = round((pivot * b[j] - a[j][v] * b[u]) / pivot, 6)
            a = np.copy(A)
            c = np.copy(C)
            b = np.copy(B)
            temp2 = cb[u]
            cb[u] = cn[v]
            cn[v] = temp2
            temp1 = nbv[v]
            nbv[v] = bv[u]
            bv[u] = temp1
            for i in range(len(bv)):
                s = bv[i]
                x[int(s[1:]) - 1] = b[i]
            for i in range(len(nbv)):
                s = nbv[i]
                x[int(s[1:]) - 1] = 0


def fraction(x, n):
    for i in range(n):
        if abs(math.ceil(x[i]) - x[i]) > 10e-4 and abs(math.floor(x[i]) - x[i]) > 10e-4:
            return False
    return True


if __name__ == '__main__':

    n, m, a, c, d, b, s, optimization_type = getInput()
    # print standard form
    a, b, cb, cn, c, x, nbv, bv, total_var, av, sv, suv = print_standard_form(
        n, m, a, c, d, b, s, M, optimization_type)
    # simplex method

    a, b, cb, cn, c, x, nbv, bv, total_var, sv = simplex(
        a, c, cn, cb, d, b, nbv, bv, total_var, x, optimization_type,
        len(suv), len(sv), len(av))
    print_table(a, c, d, b, cb, x, nbv, bv,
                optimization_type)
    print("The values for x are:")
    l_ = []
    for i in range(len(x)):
        l_.append("x" + str(i + 1) + '=' + str(x[i]))
    print(l_)
    count = 0
    while not fraction(x,n) and count < 4:
        count+=1
        print("^*"*45)
        print("Dual Simplex Iteration: ", count)
        # Introduce fractional cutting plane
        a, b, cb, cn, c, x, nbv, bv, total_var, sv = cuttingPlane(a, b, cb, cn, c, x, nbv, bv, total_var, sv)
        # Dual Simplex Method
        a, b, cb, cn, c, x, nbv, bv, total_var, sv = dualSimplex(a, c, cn, cb, d, b, nbv, bv, total_var, x, optimization_type, len(suv), len(sv), len(list_of_artificial_variables))
    # End