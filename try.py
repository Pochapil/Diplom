import numpy as np

cos = np.dot([1, 0, 1], [0, 1, 1])
print(cos)

from scipy.optimize import fsolve

R_e = 2

x_origin = 0
y_origin = 1
z_origin = 0

x_direction = 1
y_direction = 1
z_direction = 1


def equations(variables):
    theta, phi, t = variables
    eqn = np.empty(3)

    x = x_origin + x_direction * t
    y = y_origin + y_direction * t
    z = z_origin + z_direction * t

    eqn[0] = R_e * np.sin(theta) ** 3 * np.cos(phi) - x
    eqn[1] = R_e * np.sin(theta) ** 3 * np.sin(phi) - y
    eqn[2] = R_e * np.sin(theta) ** 2 * np.cos(theta) - z

    return eqn


# result_theta, result_phi, result_t = fsolve(equations, (0, 0, 0))
result = fsolve(equations, (0, 0, 0))
for i in result:
    print(i)
# print(result_theta, result_phi, result_t)

print(equations(result))

print(x_origin + x_direction * result[2])
print(y_origin + y_direction * result[2])
print(z_origin + z_direction * result[2])
