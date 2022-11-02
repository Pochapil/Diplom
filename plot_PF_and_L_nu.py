import numpy as np
import matplotlib.pyplot as plt

import config
import main_service

file_folder = 'figs/'
args_folder = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)
full_file_folder = file_folder + args_folder
full_file_folder = config.full_file_folder


N_energy = 10
phase_index = 0  # индекс фазы для nu_L_nu(nu), L_nu(nu)
# --------------------------- PF --------------------------
folder = 'luminosity_in_range/'
# folder = 'L_nu/'
# folder = 'nu_L_nu/'

file_name = "PF.txt"
PF = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

energies = [0] * N_energy
for i in range(N_energy):
    if i == 0:
        energies[i] = 1
    else:
        energies[i] = i * 4

fig = main_service.create_figure(energies, PF, is_y_2d=False)

file_name = "PF.png"
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- PF --------------------------

# --------------------------- L_nu(phase) --------------------------
folder = 'L_nu/'
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    arr[i] = main_service.load_arr_from_txt(full_file_folder + folder + 'txt/', file_name)

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
labels_arr = [''] * N_energy
for i in range(N_energy):
    labels_arr[i] = '%0.2f KeV' % energies[i]
fig_title = r'$L_{\nu}$'
fig = main_service.create_figure(phi_for_plot, arr, labels_arr=labels_arr, figure_title=fig_title)

file_name = 'L_nu' + '.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- L_nu(phase) --------------------------

# --------------------------- nu_L_nu(phase) --------------------------
folder = 'nu_L_nu/'
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    arr[i] = main_service.load_arr_from_txt(full_file_folder + folder + 'txt/', file_name)

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
labels_arr = [''] * N_energy
for i in range(N_energy):
    labels_arr[i] = '%0.2f KeV' % energies[i]

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(phi_for_plot, arr, labels_arr=labels_arr, figure_title=fig_title)

file_name = 'nu_L_nu' + '.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- nu_L_nu(phase) --------------------------

# --------------------------- L_nu(nu) --------------------------
folder = 'L_nu/'
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    arr[i] = main_service.load_arr_from_txt(full_file_folder + folder + 'txt/', file_name)

L_nu = []
for i in range(N_energy):
    L_nu.append(arr[i][phase_index])

energy = [1] * N_energy
for i in range(1, N_energy):
    energy[i] = i * 4

fig_title = r'$L_{\nu}$'
fig = main_service.create_figure(energy, L_nu, figure_title=fig_title, is_y_2d=False)

file_name = 'L_nu(nu)' + '.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- L_nu(nu) --------------------------

# --------------------------- L_nu(nu)_avg --------------------------
L_nu_avg_on_phase = [0] * N_energy
for i in range(N_energy):
    sum = 0
    for j in range(config.t_max):
        sum += arr[i][j]
    avg = sum / config.t_max
    L_nu_avg_on_phase[i] = avg

fig_title = r'$L_{\nu}$'
fig = main_service.create_figure(energy, L_nu_avg_on_phase, figure_title=fig_title, is_y_2d=False)

file_name = 'L_nu(nu)_avg' + '.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- L_nu(nu)_avg --------------------------

# --------------------------- nu_L_nu(nu) --------------------------
folder = 'nu_L_nu/'
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    arr[i] = main_service.load_arr_from_txt(full_file_folder + folder + 'txt/', file_name)

nu_L_nu = []
for i in range(N_energy):
    nu_L_nu.append(arr[i][phase_index])

energy = [1] * N_energy
for i in range(1, N_energy):
    energy[i] = i * 4

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(energy, nu_L_nu, figure_title=fig_title, is_y_2d=False)

file_name = 'nu_L_nu(nu)' + '.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- nu_L_nu(nu) --------------------------

# --------------------------- nu_L_nu(nu)_avg --------------------------
nu_L_nu_avg_on_phase = [0] * N_energy
for i in range(N_energy):
    sum = 0
    for j in range(config.t_max):
        sum += arr[i][j]
    avg = sum / config.t_max
    nu_L_nu_avg_on_phase[i] = avg

fig_title = r'$\nu L_{\nu}$'
fig = main_service.create_figure(energy, nu_L_nu_avg_on_phase, figure_title=fig_title, is_y_2d=False)

file_name = 'nu_L_nu(nu)_avg' + '.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
# --------------------------- nu_L_nu(nu)_avg --------------------------
