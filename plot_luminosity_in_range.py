import config
import numpy as np
import matplotlib.pyplot as plt

import main_service

full_file_folder = config.full_file_folder
folder = 'luminosity_in_range/'

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
N_energy = 10

# -------------------------------------------------------------------------------------------
file_name = "PF.txt"
PF = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

energy_arr = [1] * N_energy
for i in range(1, N_energy):
    energy_arr[i] = i * 4

fig = main_service.create_figure(energy_arr, PF, is_y_2d=False)

file_name = "PF.png"
main_service.save_figure(fig, full_file_folder + folder, file_name)

# -------------------------------------------------------------------------------------------
file_name = "luminosity_in_range.txt"
data_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

arr_to_plt = [0] * len(data_array)
for i in range(len(data_array)):
    arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

# -------------------------------------------------------------------------------------------
energy_min = [0] * N_energy
energy_max = [0] * N_energy
energy_min[0] = 1
energy_max[0] = 4
for i in range(1, N_energy):
    energy_min[i] = 4 * i
    energy_max[i] = 4 * (i + 1)

for energy_i in range(N_energy):
    fig_title = 'luminosity in range %0.2f - %0.2f KeV of surfaces, PF = %0.3f' % (
        energy_min[energy_i], energy_max[energy_i], PF[energy_i])
    fig = main_service.create_figure(phi_for_plot, arr_to_plt[energy_i], labels_arr='sum', figure_title=fig_title,
                                     is_y_2d=False)

    file_name = 'luminosity_in_range%0.2f_-_%0.2f_KeV_of_surfaces.png' % (energy_min[energy_i], energy_max[energy_i])
    main_service.save_figure(fig, full_file_folder + folder, file_name)

# -------------------------------------------------------------------------------------------
N_energy = 6
labels_arr = [''] * N_energy

for i in range(N_energy):
    labels_arr[i] = "%0.2f - %0.2f KeV" % (energy_min[i], energy_max[i])

fig = main_service.create_figure(phi_for_plot, arr_to_plt[:N_energy], labels_arr=labels_arr)

file_name = 'sum_of_luminosity_in_range.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)

# -------------------------------------------------------------------------------------------
# plt.style.use(['science', 'notebook', 'grid'])
