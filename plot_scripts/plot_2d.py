# from tkinter import *
import tkinter as Tk
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np

import config

fi_0 = 0
betta_mu = 70 * config.grad_to_rad

step_phi_accretion = 2 * np.pi / (config.N_phi_accretion - 1)
phi_range = np.array([fi_0 * config.grad_to_rad + step_phi_accretion * i for i in range(config.N_phi_accretion)])

theta_range = np.zeros_like(phi_range)

for i in range(len(phi_range)):
    theta_range[i] = np.pi / 2 - betta_mu * np.cos(phi_range[i])

R_alfven = (config.mu ** 2 / (2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
R_alfven = R_alfven / config.R_ns

R_e = R_alfven * config.ksi_param

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)

# внутренний радиус колонки
r = R_e * np.sin(theta_range) ** 2
r1 = r * np.sin(theta_range)
x = r1 * np.cos(phi_range)
y = r1 * np.sin(phi_range)
ax.plot(x, y, 'blue', label='inner R')

# внешний радиус колонки
r = (1 + config.dRe_div_Re) * R_e * np.sin(theta_range) ** 2
r1 = r * np.sin(theta_range)
x = r1 * np.cos(phi_range)
y = r1 * np.sin(phi_range)
ax.plot(x, y, 'green', label='outer R')

# круг R_e
x = R_e * np.cos(phi_range)
y = R_e * np.sin(phi_range)
ax.plot(x, y, 'red', label='R_e circle')

# попытка менять R_e - 69 стр - там ksi
'''  The process also requires selecting a
value for the normalization constant ξ appearing in Equation (22), such that the minimum
value of the Alfv´en radius equals the maximum value of R2, disk, corresponding to α = 0'''

c = 1
r = ((3 * np.cos(theta_range) ** 2 + 1) / c) ** (2 / 7)
ksi = max(R_e * np.sin((theta_range)) ** 2) / min(r)
r *= ksi

# ksi = 2 * max(np.sin((theta_range)) ** 2)
# print(ksi)
# r = ksi * R_e * 1 / 2 * (1 + 3 * np.cos(theta_range) ** 2) ** (1 / 2) / (R_e * np.sin(theta_range) ** 3)
rotate_angle = 0
x = r * np.cos(phi_range + rotate_angle)
y = r * np.sin(phi_range + rotate_angle)
ax.plot(x, y, 'black', label='R_e eq22')

# вектор магнитного поля
origin = np.array([0, 0])
vector = np.array([R_e * 1.25, 0])
# plt.quiver(origin, vector, angles='xy', scale_units='xy',color='black', scale=10)
ax.arrow(0, 0, *vector, head_width=0.05, head_length=1, color='purple')

ax.legend()
ax.set_xlabel('r/R_ns', fontsize=24)
ax.set_ylabel('r/R_ns', fontsize=24)

plt.title(r'$\beta_{\mu}$' + f'={betta_mu / config.grad_to_rad}')
plt.axis('equal')
plt.show()

crit_betta_mu = np.arctan((1 / 12) ** (1 / 2)) / config.grad_to_rad
betta_mu_arr = np.array([crit_betta_mu + 5 * i for i in range(18)])

phi_arr = np.arccos(1 / (12 * np.tan(betta_mu_arr * config.grad_to_rad) ** 2)) / config.grad_to_rad
a_arr = 2 * phi_arr / 360

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.plot(betta_mu_arr, a_arr)

ax.set_xlabel(r'$\beta_{\mu}$', fontsize=24)
ax.set_ylabel('a', fontsize=24)

plt.show()
