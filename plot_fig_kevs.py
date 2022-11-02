import numpy as np
import matplotlib.pyplot as plt

import config
import main_service

N_energy = 6

energy_min = [0] * N_energy
energy_max = [0] * N_energy

energy_min[0] = 1
energy_max[0] = 4

for i in range(1, N_energy):
    energy_min[i] = 4 * i
    energy_max[i] = 4 * (i + 1)

file_folder = 'figs/'
args_folder = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)
full_file_folder = file_folder + args_folder
full_file_folder = config.full_file_folder

folder = 'luminosity_in_range/'
arr_to_plt = [0] * N_energy
for i in range(N_energy):
    file_name = "sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (energy_min[i], energy_max[i])
    arr_to_plt[i] = main_service.load_arr_from_txt(full_file_folder + folder + 'txt/', file_name)

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

labels_arr = [''] * N_energy

for i in range(N_energy):
    labels_arr[i] = "%0.2f - %0.2f KeV" % (energy_min[i], energy_max[i])

fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr=labels_arr)
plt.show()
file_name = 'sum_of_luminosity_in_range.png'
main_service.save_figure(fig, full_file_folder + folder, file_name)
