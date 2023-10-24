import numpy as np
import matplotlib.pyplot as plt

import main_service
import config

parent_folder = config.PROJECT_DIR + 'figs/difference_N/'
curr_folder = 'luminosity_in_range/'
file_folder_accretion_args = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)

first_folder = 'Th 100 Phi 100/'
full_file_folder = parent_folder + first_folder + file_folder_accretion_args + curr_folder

file_name = "luminosity_in_range.txt"
first_data_array = main_service.load_arr_from_txt(full_file_folder, file_name)

second_folder = 'Th 300 Phi 300 Om 3/'
full_file_folder = parent_folder + second_folder + file_folder_accretion_args + curr_folder

second_data_array = main_service.load_arr_from_txt(full_file_folder, file_name)

first_phase_index = 3  # 8 * (2+1)
second_phase_index = 8  # 3 * (7+1)
arr_difference = [0] * len(first_data_array)
for i in range(len(first_data_array)):
    arr_difference[i] = np.abs(first_data_array[i][first_phase_index] - second_data_array[i][second_phase_index]) / \
                        first_data_array[i][first_phase_index]

print(max(arr_difference))

first_phi_for_plot = list(8 * config.grad_to_rad * i / (2 * np.pi) for i in range(360//8))
second_phi_for_plot = list(3 * config.grad_to_rad * i / (2 * np.pi) for i in range(360//3))

plt.style.use(['science', 'notebook', 'grid'])
N_column_plot = 2
fig, axes = plt.subplots(N_column_plot, 1, figsize=(12, 6 * N_column_plot))

energy_index = 17
ax = axes[0]
label = first_folder[:-1]
ax.plot(first_phi_for_plot, first_data_array[energy_index], color='black', lw=0.8)
ax.text(0.98, 0.87, label, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black'), ha='right',
        va='top')
ax = axes[1]
label = second_folder[:-1]
ax.plot(second_phi_for_plot, second_data_array[energy_index], color='black', lw=0.8)
ax.text(0.98, 0.87, label, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black'), ha='right',
        va='top')
plt.show()

# phase_index = 0
# arr_difference = [0] * len(data_array_100[phase_index])
# for i in range(len(data_array_100[phase_index])):
#     arr_difference[i] = np.abs(data_array_600[phase_index][i] - data_array_100[phase_index][i]) / \
#                         data_array_600[phase_index][i]
#
# print(max(arr_difference))
