import numpy as np
import matplotlib.pyplot as plt

import accretionColumnService
import config
import main_service

plt.style.use(['science', 'notebook', 'grid'])
# plt.style.use(['default'])

folder = 'L_nu/'
working_folder = config.full_file_folder + folder

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

y_axis_label = r'$Spectrum \: L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
# -------------------------------------------------------------------------------------------
x_axis_label = r'$h \nu$' + ' [KeV]'

file_name = "PF.txt"
PF = main_service.load_arr_from_txt(working_folder, file_name)

file_name = "energy.txt"
energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)
N_energy = config.N_energy

fig = main_service.create_figure(energy_arr, PF, x_axis_label=x_axis_label, y_axis_label='PF', is_y_2d=False)

file_name = "PF.png"
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
x_axis_label = 'Phase'

file_name = "L_nu.txt"
data_array = main_service.load_arr_from_txt(working_folder, file_name)

arr_to_plt = [0] * len(data_array)
for i in range(len(data_array)):
    arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

# -------------------------------------------------------------------------------------------
for energy_i in range(N_energy):
    fig_title = 'Spectrum of energy %0.2f KeV of surfaces, PF = %0.3f' % (energy_arr[energy_i], PF[energy_i])
    fig = main_service.create_figure(phi_for_plot, arr_to_plt[energy_i], labels_arr=r'$L_{\nu}(phase)$',
                                     x_axis_label=x_axis_label, y_axis_label=y_axis_label, figure_title=fig_title,
                                     is_y_2d=False)

    file_name = 'L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy_arr[energy_i]
    main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
# по идее переписать!!
labels_arr = [''] * N_energy
for i in range(N_energy):
    labels_arr[i] = '%0.2f KeV' % energy_arr[i]
fig_title = r'$L_{\nu}$'
fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr=labels_arr, x_axis_label=x_axis_label,
                                 y_axis_label=y_axis_label, figure_title=fig_title)

file_name = 'L_nu' + '.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
x_axis_label = r'$h \nu$' + ' [KeV]'
phase_index = 0  # индекс фазы для L_nu(nu)

L_nu = [0] * N_energy
for i in range(N_energy):
    L_nu[i] = arr_to_plt[i][phase_index]

fig_title = r'$L_{\nu}$'
fig = main_service.create_figure(energy_arr, L_nu, figure_title=fig_title, x_axis_label=x_axis_label,
                                 y_axis_label=y_axis_label, is_y_2d=False)

file_name = 'L_nu(nu)' + '.png'
main_service.save_figure(fig, working_folder, file_name)

# ----------------------------------- L_nu(nu)_avg ----------------------------------------
L_nu_avg_on_phase = [0] * N_energy
for i in range(N_energy):
    L_nu_avg_on_phase[i] = np.mean(arr_to_plt[i])

# fig_title = r'$L_{\nu}$'
# fig = main_service.create_figure(energy_arr, L_nu_avg_on_phase, x_axis_label=x_axis_label,
#                                  y_axis_label=y_axis_label, figure_title=fig_title, is_y_2d=False)
#
# file_name = 'L_nu(nu)_avg' + '.png'
# main_service.save_figure(fig, working_folder, file_name)

fig_title = r'$L_{\nu}$'
fig = main_service.create_figure(energy_arr, L_nu_avg_on_phase, figure_title=fig_title, x_axis_label=x_axis_label,
                                 y_axis_label=y_axis_label, is_y_2d=False, is_x_log_scale=True, is_y_log_scale=True)

file_name = 'L_nu(nu)_avg_log_log' + '.png'
main_service.save_figure(fig, working_folder, file_name)

# ------------------------ phases and avg ------------------------------
# plt.style.use(['science', 'notebook', 'grid'])
N_phases = 2
# phase_indexes = [0 + i * 7 for i in range(N_phases)]  # индекс фазы для L_nu(nu)
phase_indexes = [4, 15, 33]
L_nu_phases = [0] * N_phases
L_nu = [0] * N_energy

for j in range(N_phases):
    for i in range(N_energy):
        L_nu[i] = arr_to_plt[i][phase_indexes[j]]
    L_nu_phases[j] = L_nu.copy()

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)
ax.plot(energy_arr, L_nu_avg_on_phase, label=r'$L_{\nu} \, avg$', color='red')
for i in range(N_phases):
    ax.plot(energy_arr, L_nu_phases[i], label=r'$L_{\nu}$' + ' on phase %0.2f' % phi_for_plot[phase_indexes[i]])
plt.xscale('log')
plt.yscale('log')
ax.legend()
ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)

file_name = 'L_nu(nu)_avg_and_phases' + '.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
N_column_plot = config.N_column_plot
energy_indexes = config.energy_indexes
fig, axes = plt.subplots(N_column_plot, 1, figsize=(12, 3 * N_column_plot), sharex=True)
for i in range(N_column_plot):
    ax = axes[i]
    label = "%0.1f KeV\n PF=%0.3f" % (energy_arr[energy_indexes[i]], PF[energy_indexes[i]])
    ax.tick_params(axis='both', labelsize=12)
    ax.plot(phi_for_plot, arr_to_plt[energy_indexes[i]], color='black', lw=0.8)
    # ax.plot(phi_for_plot, arr_to_plt[i], color='black', lw=0.8, label=label)
    ax.text(0.98, 0.87, label, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black'), ha='right',
            va='top')
    # ax.legend(loc='upper right')

# fig.add_subplot(111, frameon=False)
# plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# plt.xlabel('phase')
# plt.ylabel('luminosity [erg/s]')
plt.rc('font', size=24)
fig.text(0.5, 0.07, 'Phase', ha='center')
fig.text(0.06, 0.5, r'$L_{\nu} \, [erg \cdot s^{-1} \cdot hz^{-1}]$', va='center', rotation='vertical')
file_name = 'pretty_fig.png'
main_service.save_figure(fig, working_folder, file_name)
plt.rc('font', size=10)

# -------------------------------------------------------------------------------------------

freq_arr = accretionColumnService.get_frequency_from_energy(energy_arr)
black_body = accretionColumnService.plank_energy_on_frequency(freq_arr, 3 * 10 ** 7)
# тут засунуть функцию get_black_body_approximation(self, energy, T_eff)

file_name = "Black_body.txt"
black_body_arr = main_service.load_arr_from_txt(working_folder, file_name)

black_body_avg = [0] * N_energy
for i in range(N_energy):
    black_body_avg[i] = np.mean(black_body_arr[i])

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)

ax.plot(energy_arr, black_body_avg, label='black body', color='black')
ax.plot(energy_arr, L_nu_avg_on_phase, label=r'$L_{\nu} \, avg$', color='red')

plt.xscale('log')
plt.yscale('log')
ax.legend()

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)

file_name = 'L_nu(nu)_avg_and_black_body' + '.png'
main_service.save_figure(fig, working_folder, file_name)
