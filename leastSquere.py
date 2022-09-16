import numpy as np
from sympy import symbols
### ----- linear problem -------------
# equation : y = a.x + b  true values for simulation -> a=2, b =1
# Do the simulation to get the observations
observation_x = [0, 1, 2]
observation_y = list(map(lambda x: 2*x+1, observation_x))
print(f" linear simulation x:{observation_x}, y:{observation_y}")
# Use least square for linear problems
# X = inv(A.T * A) * (A.T * L)
x, y, a, b = symbols('x y a b')
y = a * x + b
print("Polynomial function, f(x):\n", y)
diffs_wrt_unknowns = [y.diff(a), y.diff(b)]
print(f"derivative w.r.t unknowns-> {diffs_wrt_unknowns}")
design_matrix = []
L = []
for item1, item2 in zip(observation_x, observation_y):
    design_matrix.append([each.subs(x, item1).evalf() for each in diffs_wrt_unknowns])
    L.append(item2)

design_matrix = np.array(design_matrix, dtype=float)
print(f"A(design_matrix) = {design_matrix}")
print(f"L(observations) = {L}")
# Solve it
X_cap = np.linalg.inv(design_matrix.T @ design_matrix) @ (design_matrix.T @ L)
print(f"calculated X : {X_cap}")

# Gradient descent:
alpha = 0.5  #'learning rate'
a_init = 5
b_init = 3
for jj in range(20):
    residual = list()
    for i in range(3):
        r = float(y.subs(a, a_init).subs(b, b_init).subs(x, observation_x[i]).evalf() - observation_y[i])
        residual.append(r)

    a_init = a_init - alpha * (1/3)*sum([item1*item2 for item1, item2 in zip(residual, observation_x)])
    b_init = b_init - alpha * (1/3)*sum(residual)
    print(f'a = {a_init}, b= {b_init}')





print("_"*20)
### ----- Non-linear problem -------------
# equation : y = a.(b^2).x + b  true values for simulation -> a=3, b=2
# Do the simulation to get the observations
observation_x = [0, 1, 2]
observation_y = list(map(lambda x: 3*4*x + 2, observation_x))
print(f" non_linear simulation x:{observation_x}, y:{observation_y}")

# Use least square for non_linear problems
# delta_L = A.X_prev - L
# delta_X = inv(A.T * A) * (A.T * delta_L)
x, y, a, b = symbols('x y a b')
y = a * (b**2) * x + b
print("Polynomial function, f(x):\n", y)

diffs_wrt_unknowns = [y.diff(a), y.diff(b)]
print(f"derivative w.r.t unknowns-> {diffs_wrt_unknowns}")

# we need initial values for unknowns (non zero values)

a_init = 1
b_init = 1
for jj in range(5):
    design_matrix = []
    L = []
    for item1, item2 in zip(observation_x, observation_y):
        design_matrix.append([each.subs(x, item1).subs(a, a_init).subs(b, b_init).evalf() for each in diffs_wrt_unknowns])
        calculate_L_with_initial = y.subs(x, item1).subs(a, a_init).subs(b, b_init).evalf()
        L.append(calculate_L_with_initial - item2)


    design_matrix = np.array(design_matrix, dtype=float)
    initial_x = np.array([a_init, b_init]).T
    print(f"A(design_matrix) = {design_matrix}")
    print(f"L(observations) = {L}")
    print(f"init initial_x {initial_x}")
    delta_L = L
    print(f"delta_L {delta_L}")
    # Solve it
    # delta_X = inv(A.T * A) * (A.T * delta_L)
    delta_X_cap = np.linalg.inv(design_matrix.T @ design_matrix) @ (design_matrix.T @ delta_L)
    print(f"calculated delta X : {delta_X_cap}")
    a_init -= delta_X_cap[0]
    b_init -= delta_X_cap[1]
    print(a_init)
    print(b_init)
    print('*'*10)


# Gradient descent:
alpha = 0.0006  #'learning rate'
a_init = 5
b_init = 3
for jj in range(20):
    residual = list()
    for i in range(3):
        r = float(y.subs(a, a_init).subs(b, b_init).subs(x, observation_x[i]).evalf() - observation_y[i])
        residual.append(r)

    a_init = a_init - alpha * (1/3)*sum([item1*item2*(b_init**2) for item1, item2 in zip(residual, observation_x)])
    b_init = b_init - alpha * (1/3)*sum([item1*item2*2*a_init*b_init+1 for item1, item2 in zip(residual, observation_x)])
    print(f'a = {a_init}, b= {b_init}')