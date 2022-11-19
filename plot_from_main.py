import config
import numpy as np
import matplotlib.pyplot as plt

import main_service

plt.style.use(['science', 'notebook', 'grid'])
working_folder = config.full_file_folder

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

# -------------------------------------------------------------------------------------------
file_name = 'total_luminosity_of_surfaces.txt'
data_array = main_service.load_arr_from_txt(working_folder, file_name)

arr_to_plt = [0] * len(data_array)
# нужно расширить массивы, чтобы покрыть фазу [0,2]
for i in range(len(data_array)):
    arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

labels_arr = ['top outer', 'top inner', 'bot outer', 'bot inner', 'sum']
fig_title = 'total luminosity of surfaces'
# fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr, x_axis_label='Phase',
#                                  y_axis_label=r'$Luminosity \: [erg/s]$', figure_title=fig_title)

fig = plt.figure(figsize=(21, 10))
ax = fig.add_subplot(111)
line = [0] * 5
i = 0
ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='red', marker='.', linestyle=':')
i += 1
ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='red', marker='*')
i += 1
ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='green', marker='+', linestyle=':')
i += 1
ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='green', marker='^')
i += 1
ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='black')

# for i in range(len(data_array)):
#     line[i], = ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i])
# ax.set_xlabel('Phase', fontsize=24)
# ax.set_ylabel('L [erg/s]', fontsize=24)

ax.legend()

file_name = 'total_luminosity_of_surfaces.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
file_name = 'observer_angles.txt'
data_array = main_service.load_arr_from_txt(working_folder, file_name)
observer_phi = data_array[0]
observer_theta = data_array[1]

observer_phi = main_service.extend_arr_for_phase(observer_phi)
observer_theta = main_service.extend_arr_for_phase(observer_theta)

labels_arr = [r'$\theta_{observer}$', r'$\phi_{observer}$']
fig_title = 'Observer angles'
combined_arrays_for_plot = np.append([observer_theta], [observer_phi], 0)
fig = main_service.create_figure(phi_for_plot, combined_arrays_for_plot, labels_arr, x_axis_label='Phase',
                                 y_axis_label=r'$Angle \: [rad]$', figure_title=fig_title)

file_name = 'Observer_angles.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
file_name = 'surfaces_T_eff.txt'
data_array = main_service.load_arr_from_txt(working_folder, file_name)

file_name = "save_theta_range.txt"
theta_range = main_service.load_arr_from_txt(working_folder, file_name)

fig = main_service.create_figure(theta_range, data_array[0], x_axis_label='Phase', y_axis_label=r'$T_{eff} \: [K]$',
                                 figure_title=r'$T_{eff}$', is_y_2d=False)
file_name = 'T_eff.png'
main_service.save_figure(fig, working_folder, file_name)
