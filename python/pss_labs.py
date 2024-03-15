from gurobipy import *
from itertools import *
from math import fabs
import sys

n = int(sys.argv[1])
M = 1000

def f(x):
  N = len(x)
  R = N
  result = 0
  for i in range(1,N-R+1 + 1):
#    for d in range(1,R-1 + 1):
    for d in range(1,2):
      factor = sum( (2*x[j-1]-1)*(2*x[j+d-1] - 1) for j in range(i,i+R-1-d + 1) )
      result += factor * factor
  return result

zeroVector = tuple([0] * n)
unitVectors = []
for j in range(n):
  vector = [0] * n
  vector[j] = 1
  unitVectors.append(tuple(vector))

def fNonlinear(x):
  global zeroVector
  global unitVectors

  return f(x) - f(zeroVector) - sum( (f(unitVectors[j]) - f(zeroVector)) * x[j] for j in range(n) )

# Compute a list of all binary vectors of length n.
binaryVectors = []
for j in range(2**n):
  x = [ 0 ] * n
  for i in range(n):
    if (j >> (n-i-1)) & 1:
      x[i] = 1
  binaryVectors.append(tuple(x))
binaryVectors = tuple(binaryVectors)


Y = set()
for vector in binaryVectors:
  y = fNonlinear(vector)
  Y.add(y)
  print(f'{vector} -> {f(vector)} and {y}')

print(f'|Y| = {len(Y)}; Y = {Y}')

# w_i in \R
# x_i = 0 => w_i = 0
# z^y indicator vector for y
# sum_i w_i z^y_i 

sys.exit(1)

# Consider all subsets.
if False:
  subsets = list(map(tuple, chain.from_iterable(combinations( binaryVectors , r) for r in range(len(binaryVectors)+1)) ))

# Consider only the subsets stemming from products of (complemented) variables.
if True:
  subsets = []
  for vector in product([0, 1, -1], repeat=n):
    X = []
    for x in binaryVectors:
      add = True
      for i in range(n):
        if vector[i] == 1 and x[i] == 0:
          add = False
          break
        elif vector[i] == -1 and x[i] == 1:
          add = False
          break
      if add:
        X.append(tuple(x))
    subsets.append(tuple(X))
  print(f'Generated {len(subsets)} subsets corresponding to products of (complemented) variables.')

# Consider only monomials without negations.
if False:
  subsets = []
  for vector in product([0, 1], repeat=n):
    X = []
    for x in binaryVectors:
      add = True
      for i in range(n):
        if vector[i] == 1 and x[i] == 0:
          add = False
          break
      if add:
        X.append(tuple(x))
    subsets.append(tuple(X))
  print(f'Generated {len(subsets)} subsets corresponding to products of variables.')

#if not binaryVectors in subsets:
#  subsets.append(binaryVectors)

#print(subsets)

model = Model("blc")
model.modelSense = GRB.MINIMIZE

var_a = {}
for i in range(n):
  var_a[i] = model.addVar(name=f'a#{i}', lb=-GRB.INFINITY)

var_beta = model.addVar(name='beta', lb=-GRB.INFINITY)

var_b = {}
for X in subsets:
  var_b[X] = model.addVar(lb=-M, ub=M)

var_x = {}
for X in subsets:
  var_x[X] = model.addVar(vtype=GRB.BINARY, obj=1.0)

model.update()

for X in subsets:
  X = tuple(X)
  model.addConstr( var_b[X] <= M*var_x[X] )
  model.addConstr( var_b[X] >= -M*var_x[X] )

for x in binaryVectors:
  lhs = quicksum( x[i] * var_a[i] for i in range(n) ) + var_beta
  for X in subsets:
    if x in X:
       lhs += var_b[X]
  model.addConstr( lhs == f(x) )

model.write('debug.lp')

model.optimize()

print(f'\n\nPrinting solution of size {len([0 for X in subsets if var_x[X].x > 0.5])}.')

for i in range(n):
  if fabs(var_a[i].x) > 1.0e-5:
    print(f'Nonzero coefficient {var_a[i].x} of a_{i}.')
print(f'Value of beta: {var_beta.x}\n')
for X in subsets:
  if var_x[X].x > 0.5:
    print(f'Nonzero coefficient {var_b[X].x} of set {X}.')

