from pysat.card import *
from pysat.solvers import Glucose3
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
import pysat.pb as pb
import numpy as np
from pysat.examples.lsu import LSU

n = 0
d = np.zeros(21)
r = np.zeros(21)
maxxD = 0


def sgn(x):
    if x < 0:
        return -1
    return 1


def yy(x, y, ml):
    return int(n * ml + (x - 1) * ml + y)


def xx(x, y, ml):
    return int((x - 1) * ml + y)


def tt(x, y, ml):
    return int((x - 1) * ml + y + yy(n + 1, ml, ml))


def create_formula(maxL, weighted):
    print(maxL)
    # if it is not possible, return an unsatisfiable formula
    if maxL < maxxD:
        return "1 0\n-1 0", 2, 1000
    new_form = ""
    new_clauses = 0
    for i in range(1, n + 1):
        lst = []

        for j in range(1, int(maxL) + 2 - int(d[i])):
            lst.append(int(xx(i, j, maxL)))
        formula = CardEnc.equals(lst, 1, encoding=EncType.pairwise)
        for c in formula:
            for x in c:
                new_form += str(x) + " "
            new_form += "0\n"
            new_clauses += 1
    maxL = int(maxL)

    # add a constraint for project i to take d_{i} units of time
    for i in range(1, n + 1):
        for j in range(1, int(maxL + 2 - d[i])):
            for k in range(j, min(maxL + 1, int(j + d[i]))):
                new_form += str(-xx(i, j, maxL)) + " " + str(yy(i, k, maxL)) + " 0\n"
                new_clauses += 1

    # dependencies encoding
    # next n * maxL variables are t_{x}{y}
    for i in range(1, n + 1):
        for j in range(1, maxL + 1):
            for k in range(int(j + d[i]), int(maxL + 1)):
                new_form += str(-xx(i, j, maxL)) + " " + str(tt(i, k, maxL)) + " 0\n"
                new_clauses += 1
            # make sure t is correctly formed
            for k in range(1, min(maxL + 1, int(j + d[i]))):
                new_form += str(-xx(i, j, maxL)) + " " + str(-tt(i, k, maxL)) + " 0\n"
                new_clauses += 1
    # encoding the actual dependencies
    for i in range(1, n + 1):
        for j in range(1, maxL + 1):
            for x in s[i]:
                new_form += str(-xx(i, j, maxL)) + " " + str(tt(x, j, maxL)) + " 0\n"
                new_clauses += 1

    # resource overloading constraints
    lst = []
    w = []
    maxx = tt(n, maxL, maxL)
    for i in range(1, n + 1):
        lst.append(i)
        w.append(int(r[i]))
    maxV = 0
    formula = pb.PBEnc.atmost(lits=lst, weights=w, bound=capacity, encoding=pb.EncType.best)
    for i in formula:
        for x in i:
            maxV = max(maxV, abs(x))
    for t in range(1, maxL + 1):
        for c in formula:
            for x in c:
                val = abs(x)
                sign = sgn(x)
                if val > n:
                    val = val - n + tt(n, maxL, maxL) + t * maxV
                else:
                    val = yy(val, t, maxL)
                new_form += str(sign * val) + " "
                maxx = max(maxx, val)
            new_form += "0\n"
            new_clauses += 1

    return new_form, new_clauses, maxx


if __name__ == '__main__':
    # Start parsing
    file = open(r"C:\Users\Traian\PycharmProjects\RCPSPBinarySearch - Copy\rcpsp-instances\Q3_1.dzn").read()
    ind = 0
    while not file[ind].isdigit():
        ind += 1
    capacity = 0
    while file[ind].isdigit():
        capacity = capacity * 10 + int(file[ind])
        ind += 1
    ind += 1
    while not file[ind].isdigit():
        ind += 1
    n = 0
    while file[ind].isdigit():
        n = n * 10 + int(file[ind])
        ind += 1
    ind += 2
    i = 1
    D = 0
    while file[ind] != ';':
        while not file[ind].isdigit():
            ind += 1
        while file[ind].isdigit():
            d[i] = d[i] * 10 + int(file[ind])
            ind += 1
        d[i] = int(d[i])
        maxxD = max(int(maxxD), int(d[i]))
        i += 1
        D += d[i - 1]
        ind += 2
    ind += 1
    i = 1
    while file[ind] != ';':
        while not file[ind].isdigit():
            ind += 1
        while file[ind].isdigit():
            r[i] = r[i] * 10 + int(file[ind])
            ind += 1
        i += 1
        ind += 2
    i = 1

    s = [[] for _ in range(n + 1)]
    while ind < len(file):
        while ind < len(file) and file[ind] != '{':
            ind += 1
        ind += 2
        while ind < len(file) and file[ind] != ',' and file[ind] != ' ':
            nr = 0
            while ind < len(file) and file[ind].isdigit():
                nr = nr * 10 + int(file[ind])
                ind += 1
            s[i].append(nr)
            ind += 2
        i += 1
    # End parsing
    # First n * maxL variables are x_{i}{t}
    # Next n * maxL are y_{i}{t}
    # Then come t_{i}{t}
    # binary search for lowest possible finish time

    msk = 1 << (int(np.log2(D)) + 1)
    cost = 0

    for _ in range(int(np.log2(D)) + 2):
        if msk == 0:
            break
        result = create_formula(cost + msk, False)
        form = result[0]
        clauses = result[1]
        variables = result[2]
        header = "p cnf {0} {1} \n".format(result, clauses)
        form = header + form
        with open(r"C:\Users\Traian\PycharmProjects\RCPSPBinarySearch - Copy\output.cnf", 'w') as file:
            file.write(form)
        cnf = CNF()
        cnf.from_file(r"C:\Users\Traian\PycharmProjects\RCPSPBinarySearch - Copy\output.cnf")
        solver = Glucose3()
        var = solver.append_formula(cnf.clauses, no_return=False)
        if not solver.solve():
            cost += msk
        msk /= 2
        if msk == 4:
            check = form
    cost += 1
    # end binary search, found the lowest finish time possible
    # create the formula for the found cost
    result = create_formula(cost, False)
    form = result[0]
    clauses = result[1]
    variables = result[2]
    header = "p cnf {0} {1} \n".format(variables, clauses)
    form = header + form
    with open(r"C:\Users\Traian\PycharmProjects\pythonProject10\output.cnf", 'w') as file:
        file.write(form)
    cnf = CNF()
    cnf.from_file(r"C:\Users\Traian\PycharmProjects\pythonProject10\output.cnf")
    solver = Glucose3()
    solver.append_formula(cnf.clauses, no_return=False)
    solver.solve()
    print(cost)
    model = solver.get_model()
    for i in range(1, n + 1):
        start = 0
        for j in range(1, int(cost) + 1):  # here we will go up to cost
            if model[xx(i, j, cost) - 1] > 0:
                start = j
                break
        print("project {0} starts at time {1}".format(i, start))
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
