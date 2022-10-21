import numpy as np
import matplotlib.pyplot as plt

import geometricTask.matrix as matrix
import config
import accretionColumnService
from accretionColumn import AccretionColumn
import vectors

R_alfven = (config.mu ** 2 / (2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
print('R_e = %f' % (R_e / config.R_ns))
R_e_outer_surface, R_e_inner_surface = R_e, R_e
# вектор на наблюдателя в системе координат двойной системы
e_obs = np.array([0, np.sin(config.i_angle), np.cos(config.i_angle)])
file_name_variables = "betta_omega=%d betta_mu=%d a_portion=%f M_rate_c2_Led=%d" \
                      % (config.betta_rotate, config.betta_mu, config.a_portion, config.M_rate_c2_Led)
approx_method = accretionColumnService.approx_method

file_folder = 'figs/'

# от поверхности NS - угол при котором радиус = радиусу НЗ
# ----------------- начало инициализации верхней колонки ------------------------
theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(R_e_outer_surface)
theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(R_e_inner_surface)

top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                             theta_accretion_begin_inner_surface, True)
# ----------------- конец инициализации верхней колонки ------------------------

# ----------------- начало инициализации нижней колонки ------------------------
theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(top_column.outer_surface.theta_range, top_column.outer_surface.T_eff)
plt.show()

bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                             theta_accretion_begin_inner_surface, False)
# ----------------- конец инициализации нижней колонки ------------------------

surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface, 2: bot_column.outer_surface,
            3: bot_column.inner_surface}

print('phi_theta_range saved')
file_name = "save_phi_range.txt"
full_file_name = file_folder + file_name
np.savetxt(full_file_name, top_column.outer_surface.phi_range)
file_name = "save_theta_range.txt"
full_file_name = file_folder + file_name
np.savetxt(full_file_name, top_column.outer_surface.theta_range)

# print('T_eff:')
# print(top_column.outer_surface.T_eff)

# ----------------- углы для нахождения пересечений -------------------------
theta_accretion_begin = top_column.outer_surface.theta_range[0]
theta_accretion_end = top_column.outer_surface.theta_range[-1]
# ---------------------------------------------------------------------------

# ------------------ начало заполнения матриц косинусов ---------------------------
for key, surface in surfaces.items():
    surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end, top_column.outer_surface.phi_range,
                               bot_column.outer_surface.phi_range, e_obs)
# ------------------ конец заполнения матриц косинусов ---------------------------

arr_simps_integrate = [0] * 4
sum_simps_integrate = 0
for key, surface in surfaces.items():
    arr_simps_integrate[key] = surface.calculate_integral_distribution()
    # file_name = "%s %s %d.txt" % (file_name_variables, approx_method, key)
    # np.savetxt(file_name, arr_simps_integrate[key])
    sum_simps_integrate += np.array(arr_simps_integrate[key])
print('ksi_shock = %f' % bot_column.outer_surface.ksi_shock)

fig = plt.figure(figsize=(8, 8))
phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

append_index = config.t_max_for_plot - config.t_max
for i in range(4):
    arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

ax = fig.add_subplot(111)
ax.plot(phi_for_plot, arr_simps_integrate[0],
        label='top outer')
ax.plot(phi_for_plot, arr_simps_integrate[1],
        label='top inner')
ax.plot(phi_for_plot, arr_simps_integrate[2],
        label='bot outer', marker='*')
ax.plot(phi_for_plot, arr_simps_integrate[3],
        label='bot inner')
ax.plot(phi_for_plot, sum_simps_integrate,
        label='sum')
ax.legend()
fig.suptitle('total luminosity of surfaces', fontsize=14)
# plt.yscale('log')
file_name = 'total_luminosity_of_surfaces.png'
full_file_name = file_folder + file_name
fig.savefig(full_file_name + 'total_luminosity_of_surfaces.png', dpi=fig.dpi)
plt.show()

observer_theta = [0] * config.t_max_for_plot
observer_phi = [0] * config.t_max_for_plot

for t in range(config.t_max_for_plot):
    phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
    # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
    A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
    observer_phi[t], observer_theta[t] = vectors.get_angles_from_vector(e_obs_mu)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(phi_for_plot, observer_theta, label=r'$\theta_{observer}$')
ax.plot(phi_for_plot, observer_phi, label=r'$\phi_{observer}$')
ax.legend()
fig.suptitle('Observer angles', fontsize=14)
file_name = 'Observer_angles.png'
full_file_name = file_folder + file_name
fig.savefig(full_file_name + 'Observer_angles.png', dpi=fig.dpi)
plt.show()

while True:
    energy_bot = float(input('введите нижний предел в КэВ: '))
    energy_top = float(input('введите верхний предел в КэВ: '))

    arr_simps_integrate = [0] * 4
    sum_simps_integrate = 0
    for key, surface in surfaces.items():
        arr_simps_integrate[key] = surface.calculate_integral_distribution_in_range(energy_bot, energy_top)
        sum_simps_integrate += np.array(arr_simps_integrate[key])

    PF = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

    fig = plt.figure(figsize=(8, 8))
    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    append_index = config.t_max_for_plot - config.t_max
    for i in range(4):
        arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
    sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

    ax = fig.add_subplot(111)
    # ax.plot(phi_for_plot, arr_simps_integrate[0],
    #         label='top outer')
    # ax.plot(phi_for_plot, arr_simps_integrate[1],
    #         label='top inner')
    # ax.plot(phi_for_plot, arr_simps_integrate[2],
    #         label='bot outer', marker='*')
    # ax.plot(phi_for_plot, arr_simps_integrate[3],
    #         label='bot inner')
    ax.plot(phi_for_plot, sum_simps_integrate,
            label='sum')
    ax.legend()
    file_name = "sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (energy_bot, energy_top)
    full_file_name = file_folder + file_name
    np.savetxt(full_file_name, sum_simps_integrate)

    fig_title = 'luminosity in range %0.2f - %0.2f KeV of surfaces, PF = %0.3f' % (energy_bot, energy_top, PF)
    file_name = 'luminosity_in_range%0.2f_-_%0.2f_KeV_of_surfaces.png' % (energy_bot, energy_top)
    full_file_name = file_folder + file_name
    fig.suptitle(fig_title, fontsize=14)
    fig.savefig(full_file_name, dpi=fig.dpi)
    # plt.yscale('log')
    plt.show()
