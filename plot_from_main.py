import config
import numpy as np

import main_service

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
fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr, x_axis_label='phase',
                                 y_axis_label='luminosity, erg/s', figure_title=fig_title)
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
fig = main_service.create_figure(phi_for_plot, combined_arrays_for_plot, labels_arr, x_axis_label='phase',
                                 figure_title=fig_title)

file_name = 'Observer_angles.png'
main_service.save_figure(fig, working_folder, file_name)

# -------------------------------------------------------------------------------------------
file_name = 'surfaces_T_eff.txt'
data_array = main_service.load_arr_from_txt(working_folder, file_name)

file_name = "save_theta_range.txt"
theta_range = main_service.load_arr_from_txt(working_folder, file_name)

fig = main_service.create_figure(theta_range, data_array[0], x_axis_label='phase', y_axis_label=r'$T_{eff}, K$',
                                 figure_title=r'$T_{eff}$', is_y_2d=False)
file_name = 'T_eff.png'
main_service.save_figure(fig, working_folder, file_name)

