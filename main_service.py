import numpy as np
import matplotlib.pyplot as plt
import pathlib
# import matplotlib
# matplotlib.use('Agg')

import config


def create_file_path(file_path):
    pathlib.Path(file_path).mkdir(parents=True, exist_ok=True)


def fill_arr_with_func(func, surface, energy):
    return func(surface, energy)


def create_figure(x, y_arr, labels_arr='', x_axis_label='', y_axis_label='', figure_title='', is_y_2d=True,
                  is_x_log_scale=False, is_y_log_scale=False):
    #plt.style.use(['science', 'notebook', 'grid'])
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    if is_y_2d:
        if labels_arr == '':
            for i in range(len(y_arr)):
                ax.plot(x, y_arr[i])
        else:
            for i in range(len(y_arr)):
                ax.plot(x, y_arr[i], label=labels_arr[i])
    else:
        ax.plot(x, y_arr, label=labels_arr)

    ax.set_xlabel(x_axis_label, fontsize=24)
    if is_x_log_scale:
        plt.xscale('log')
    ax.set_ylabel(y_axis_label, fontsize=24)
    if is_y_log_scale:
        plt.yscale('log')
    if not (labels_arr == ''):
        ax.legend()  # frameon=False - без рамки
    fig.suptitle(figure_title, fontsize=14)
    return fig


def save_figure(fig, file_path, file_name):
    create_file_path(file_path)
    full_file_name = file_path + file_name
    fig.savefig(full_file_name, dpi=fig.dpi)
    plt.close(fig)


def extend_arr_for_phase(arr):
    # нужно расширить массивы, чтобы покрыть фазу [0,2]
    append_index = config.t_max_for_plot - config.t_max
    array_to_plot = np.append(arr, arr[0:append_index])
    return array_to_plot


def prepare_phase_and_extend_arr_for_phase(combined_arrays):
    phi_for_plot = list(
        config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    append_index = config.t_max_for_plot - config.t_max
    arr_dimensions = combined_arrays.ndim
    if arr_dimensions > 1:
        arrays_to_plot = [0] * len(combined_arrays)
        for i in range(len(combined_arrays)):
            arrays_to_plot[i] = np.append(combined_arrays[i], combined_arrays[i][0:append_index])
            return phi_for_plot, arrays_to_plot
    else:
        array_to_plot = np.append(combined_arrays, combined_arrays[0:append_index])
        return phi_for_plot, array_to_plot


def save_arr_as_txt(arr, file_folder, file_name):
    create_file_path(file_folder)
    full_file_name = file_folder + file_name
    np.savetxt(full_file_name, arr)


def load_arr_from_txt(file_folder, file_name):
    return np.loadtxt(file_folder + file_name)
