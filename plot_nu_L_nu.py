import numpy as np
import matplotlib.pyplot as plt

import config
import main_service

plt.style.use(['science', 'notebook', 'grid'])
# plt.style.use(['default'])

folder = 'nu_L_nu/'
working_folder = config.full_file_folder + folder

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

x_axis_label = 'Phase'
y_axis_label = r'$Spectral \: energy \: \nu \cdot L_{\nu} \: [erg/s]$'
# -------------------------------------------------------------------------------------------
file_name = "PF.txt"
PF = main_service.load_arr_from_txt(working_folder, file_name)

file_name = "energy.txt"
energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)
N_energy = config.N_energy

fig = main_service.create_figure(energy_arr, PF, x_axis_label='Phase', y_axis_label='PF', is_y_2d=False)

file_name = "PF.png"
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
file_name = "nu_L_nu.txt"
data_array = main_service.load_arr_from_txt(working_folder, file_name)

arr_to_plt = [0] * len(data_array)
for i in range(len(data_array)):
    arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

# -------------------------------------------------------------------------------------------
for energy_i in range(N_energy):
    fig_title = 'Spectral Energy Distribution of energy %0.2f KeV of surfaces, PF = %0.3f' % (
        energy_arr[energy_i], PF[energy_i])
    fig = main_service.create_figure(phi_for_plot, arr_to_plt[energy_i], labels_arr=r'$\nu \cdot L_{\nu}(\nu)$',
                                     x_axis_label=x_axis_label, y_axis_label=y_axis_label, figure_title=fig_title,
                                     is_y_2d=False)

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy_arr[energy_i]
    main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
# по идее переписать!!
labels_arr = [''] * N_energy
for i in range(N_energy):
    labels_arr[i] = '%0.2f KeV' % energy_arr[i]

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr=labels_arr, x_axis_label=x_axis_label,
                                 y_axis_label=y_axis_label, figure_title=fig_title)

# -------------------------------------------------------------------------------------------
phase_index = 0  # индекс фазы для nu_L_nu(nu)

nu_L_nu = [0] * N_energy
for i in range(N_energy):
    nu_L_nu[i] = arr_to_plt[i][phase_index]

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(energy_arr, nu_L_nu, x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                                 figure_title=fig_title, is_y_2d=False)

file_name = 'nu_L_nu(nu)' + '.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
nu_L_nu_avg_on_phase = [0] * N_energy
for i in range(N_energy):
    nu_L_nu_avg_on_phase[i] = np.mean(arr_to_plt[i])

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(energy_arr, nu_L_nu_avg_on_phase, x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                                 figure_title=fig_title, is_y_2d=False)

file_name = 'nu_L_nu(nu)_avg' + '.png'
main_service.save_figure(fig, working_folder, file_name)

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(energy_arr, nu_L_nu_avg_on_phase, x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                                 figure_title=fig_title, is_y_2d=False, is_x_log_scale=True, is_y_log_scale=True)

file_name = 'nu_L_nu(nu)_avg_log_log' + '.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
plt.style.use(['science', 'notebook', 'grid'])

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
fig.text(0.5, 0.08, 'Phase', ha='center')
fig.text(0.04, 0.5, r'$\nu \cdot L_{\nu} \: [erg/s]$', va='center', rotation='vertical')
file_name = 'pretty_fig.png'
main_service.save_figure(fig, working_folder, file_name)
plt.rc('font', size=10)
