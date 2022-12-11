import config
import numpy as np
import matplotlib.pyplot as plt

import main_service


def plot_save_sky_map(obs_i_angle_arr):
    file_folder = 'figs/sky_map/'
    config.full_file_folder = file_folder + config.file_folder_args

    file_name = "L_x.txt"
    L_x = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    data_array = [0] * len(obs_i_angle_arr)
    for i in range(len(obs_i_angle_arr)):
        file_name = "i=%d.txt" % obs_i_angle_arr[i]
        data_array[i] = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    phase = np.linspace(0, 1, config.t_max)

    # fig, ax = plt.subplots()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # нормировка на L_nu_avg
    data_to_plot = []
    for arr in data_array:
        data_to_plot.append(arr / L_x)

    im = ax.pcolormesh(phase, obs_i_angle_arr, data_to_plot)

    x_axis_label = 'Phase'
    y_axis_label = r'$i_{obs}$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    fig_title = 'L_iso/L_x'
    fig.suptitle(fig_title, fontsize=14)

    plt.colorbar(im)

    file_name = 'map' + '.png'
    main_service.save_figure(fig, config.full_file_folder, file_name)


if __name__ == '__main__':
    obs_i_angle_arr = np.linspace(0, 180, 19)
    plot_save_sky_map(obs_i_angle_arr)
