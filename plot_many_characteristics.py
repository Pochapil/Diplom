import numpy as np
import multiprocessing as mp
from itertools import repeat
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import scienceplots

import geometricTask.matrix as matrix
import config
import accretionColumnService
from accretionColumn import AccretionColumn
import vectors
import main_service

import plot_from_main
import plot_luminosity_in_range
import plot_L_nu
import plot_nu_L_nu

plt.style.use(['science', 'notebook', 'grid'])  # для красивых графиков

# jet =  mpl.colormaps['jet']
# energy_index = 8 - 4.7 keV


marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}


def get_folder_with_args(obs_i_angle_deg, betta_mu_deg, mc2, a_portion, fi_0):
    file_folder = 'figs/loop/'
    file_folder_angle = 'i=%d betta_mu=%d/' % (obs_i_angle_deg, betta_mu_deg)
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (mc2, a_portion, fi_0)
    full_file_folder = file_folder + file_folder_angle + file_folder_args

    return full_file_folder


def get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr, energy_index=8):
    folder = 'nu_L_nu/'
    file_name = 'PF.txt'
    # такой порядок: a, mc2, i, beta, fi - для того, чтобы делать colormesh x = beta, y=i; по fi делаем усреднение
    final_final_array = np.zeros((len(a_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr), len(fi_0_arr)))

    for a_index in range(len(a_arr)):
        for mc2_index in range(len(mc2_arr)):
            for i_angle_index in range(len(i_angle_arr)):
                for betta_mu_index in range(len(betta_mu_arr)):
                    final_array = []
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = get_folder_with_args(i_angle_arr[i_angle_index],
                                                                betta_mu_arr[betta_mu_index],
                                                                mc2_arr[mc2_index], a_arr[a_index],
                                                                fi_0_arr[fi_0_index])

                        PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
                        final_array.append(PF_array[energy_index])

                    # label = r'$\theta{obs}$' + ('=%d' % obs_i_angle_deg) + (r'$\beta_{\mu}$') + ('=%d' % betta_mu_deg) + \
                    #         ('a=%0.2f' % a_portion) + (r'$\dot{m}$') + ('=%d' % M_rate_c2_Led)

                    # label = 'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led)
                    final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index] = final_array

    return final_final_array


def get_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr, energy_index=8,
                      final_final_array=np.array([])):
    ''' Дисперсия среднне и max - min colormesh '''

    # if not isinstance(final_final_array, np.ndarray):
    if not final_final_array.size:
        final_final_array = get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr)

    dispersion_arr = np.zeros((len(a_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr)))
    mean_arr = np.zeros((len(a_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr)))
    max_min_arr = np.zeros((len(a_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr)))

    for a_index in range(len(a_arr)):
        for mc2_index in range(len(mc2_arr)):
            for i_angle_index in range(len(i_angle_arr)):
                for betta_mu_index in range(len(betta_mu_arr)):
                    dispersion_arr[a_index][mc2_index][i_angle_index][betta_mu_index] = np.var(
                        final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index]) ** (1 / 2)

                    mean_arr[a_index][mc2_index][i_angle_index][betta_mu_index] = np.mean(
                        final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index])

                    max_min_arr[a_index][mc2_index][i_angle_index][betta_mu_index] = np.max(
                        final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index]) - np.min(
                        final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index])

    return dispersion_arr, mean_arr, max_min_arr


def plot_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr, energy_index=8,
                       final_final_array=np.array([])):
    # функция для сохранения colormesh
    def plot_and_save_col_bar(n_rows, n_cols, save_file_name, limits_arr=[]):
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 6))
        for i in range(len(axes)):
            if limits_arr:
                max_value, min_value = limits_arr[i]
                col_bar = axes[i].pcolormesh(betta_mu_arr, i_angle_arr, data_dict[i][a_index][mc2_index],
                                             vmax=max_value, vmin=min_value)
            else:
                col_bar = axes[i].pcolormesh(betta_mu_arr, i_angle_arr, data_dict[i][a_index][mc2_index])

            axes[i].set_xlabel(x_axis_label, fontsize=24)
            axes[i].set_ylabel(y_axis_label, fontsize=24)
            axes[i].set_title(title_dict[i], fontsize=24)

            plt.colorbar(col_bar, ax=axes[i])

        fig.tight_layout()
        # plt.show()
        main_service.save_figure(fig, save_folder, save_file_name)

    dispersion_arr, mean_arr, max_min_arr = get_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr,
                                                              final_final_array=final_final_array)

    # словари для циклов
    title_dict = {0: 'mean(PF)', 1: r'$\sigma = \sqrt{D[PF]}$', 2: 'max(PF) - min(PF)'}
    data_dict = {0: mean_arr, 1: dispersion_arr, 2: max_min_arr}

    # находим min and max для pcolormesh с vmax и vmin - в относительной мере
    limits_arr = []

    max_value_mean = np.max(mean_arr)
    min_value_mean = np.min(mean_arr)
    limits_arr.append((max_value_mean, min_value_mean))

    max_value_disp = np.max(dispersion_arr)
    min_value_disp = np.min(dispersion_arr)
    limits_arr.append((max_value_disp, min_value_disp))

    max_value_max_min = np.max(max_min_arr)
    min_value_max_min = np.min(max_min_arr)
    limits_arr.append((max_value_max_min, min_value_max_min))

    for a_index in range(len(a_arr)):
        for mc2_index in range(len(mc2_arr)):
            x_axis_label = r'$\beta_{\mu}$'
            y_axis_label = r'$\theta_{obs}$'

            save_folder = 'figs/colormesh_mean_disp/mc2=%d/a=%0.2f/' % (mc2_arr[mc2_index], a_arr[a_index])

            save_file_name = 'colormesh_mean_disp_fig.png'
            plot_and_save_col_bar(1, 3, save_file_name)

            save_file_name = 'colormesh_mean_disp_without_fig.png'
            plot_and_save_col_bar(1, 2, save_file_name)

            # тут с лимитами по max min
            save_file_name = 'colormesh_mean_disp_fig_fixed.png'
            plot_and_save_col_bar(1, 3, save_file_name, limits_arr)


def print_L_nu(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr, energy_index=8, final_final_array=np.array([])):
    #  PF(L_nu) - pf берутся среднее по фазе, из-за fi_0 - множество точек

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    # если не передали - то создадим
    # if not isinstance(final_final_array, np.ndarray):
    if not final_final_array.size:
        final_final_array = get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_arr, fi_0_arr)

    # x = np.zeros((len(a_portion_arr), len(mc2)))
    # y = np.zeros((len(a_portion_arr), len(mc2)))

    markers_dict = {0: 'x', 1: '+', 2: 'o', 3: '.'}

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):

            fig = plt.figure(figsize=(12, 6))
            ax = fig.add_subplot(111)

            fig1 = plt.figure(figsize=(12, 6))
            ax1 = fig1.add_subplot(111)

            for a_index in range(len(a_arr)):
                for mc2_index in range(len(mc2_arr)):
                    L_nu_data_dict = {}
                    color_dict = {}
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = get_folder_with_args(i_angle_arr[i_angle_index],
                                                                betta_mu_arr[betta_mu_index], mc2_arr[mc2_index],
                                                                a_arr[a_index], fi_0_arr[fi_0_index])

                        L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                        L_nu_data_dict[np.mean(L_nu_array)] = \
                            final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index][fi_0_index]

                        color_dict[np.mean(L_nu_array)] = fi_0_arr[fi_0_index]

                    # print(full_file_folder)
                    lists = sorted(L_nu_data_dict.items())  # sorted by key, return a list of tuples
                    x, y = zip(*lists)  # unpack a list of pairs into two tuples

                    lists = sorted(color_dict.items())
                    buffer, colors = zip(*lists)
                    colors = np.array(colors) / 340

                    # fig = main_service.create_figure(x, y, x_axis_label=r'$L_{\nu}$', y_axis_label='PF',
                    #                                  figure_title=title, is_y_2d=False)

                    label = 'a=%0.2f m=%d' % (a_arr[a_index], mc2_arr[mc2_index])
                    ax.scatter(x, y, marker=marker_dict[a_index * len(mc2_arr) + mc2_index], color=cm.jet(colors),
                               label=label)

                    ax1.scatter(x, y, marker=marker_dict[a_index * len(mc2_arr) + mc2_index], label=label)

                save_folder = 'figs/PF_to_L_nu/mc_a/colors/'
                save_file_name = 'i=%d betta_mu=%d PF_to_L_nu_color.png' % (
                    i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                x_axis_label = r'$L_{\nu}$'
                y_axis_label = r'$PF_{' + '%.2f' % energy_arr[energy_index] + r'}$'
                figure_title = r'$\theta{obs}$' + ('=%d' % i_angle_arr[i_angle_index]) + (r'$^\circ$') + ' ' + \
                               (r'$\beta_{\mu}$') + ('=%d' % betta_mu_arr[betta_mu_index]) + (r'$^\circ$')

                ax.set_xlabel(x_axis_label, fontsize=24)
                ax.set_ylabel(y_axis_label, fontsize=24)
                fig.suptitle(figure_title, fontsize=24)

                ax.set_xscale('log')
                ax.legend()

                main_service.save_figure(fig, save_folder, save_file_name)

                save_folder = 'figs/PF_to_L_nu/mc_a/noncolors/'
                save_file_name = 'i=%d betta_mu=%d PF_to_L_nu.png' % (
                    i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                # figure_title = 'i=%d betta_mu=%d' % (obs_i_angle_deg, betta_mu_deg)

                ax1.set_xlabel(x_axis_label, fontsize=24)
                ax1.set_ylabel(y_axis_label, fontsize=24)
                fig1.suptitle(figure_title, fontsize=24)

                ax1.set_xscale('log')
                ax1.legend()

                main_service.save_figure(fig1, save_folder, save_file_name)


def plot_L_nu_iso(i_angle_arr, betta_mu_arr, a_arr, mc2_arr, fi_0_arr, energy_index=8):
    """2 мерный график nu L nu от phi_0, Phi"""

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    phase = np.linspace(0, 1, config.t_max)

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):
            for a_index in range(len(a_arr)):
                for mc2_index in range(len(mc2_arr)):

                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_subplot(111)

                    L_nu_data = []
                    color_dict = {}
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = get_folder_with_args(i_angle_arr[i_angle_index],
                                                                betta_mu_arr[betta_mu_index],
                                                                mc2_arr[mc2_index], a_arr[a_index],
                                                                fi_0_arr[fi_0_index])

                        L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                        L_nu_data.append(L_nu_array)

                    im = ax.contourf(phase, fi_0_arr, L_nu_data, levels=30)
                    clb = plt.colorbar(im)
                    clb.set_label(r'$\nu L_{\nu}$' + ' [erg/s]', fontsize=26)

                    figure_title = r'$L_{\nu}$ ' + 'i=%d betta_mu=%d a=%0.2f m=%d' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index], a_arr[a_index], mc2_arr[mc2_index])

                    save_folder = 'figs/L_nu_to_phi_to_phase/mc2=%d a=%0.2f/' % (mc2_arr[mc2_index], a_arr[a_index])
                    save_file_name = 'i=%d betta_mu=%d L_nu_to_phi_to_phase.png' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                    x_axis_label = r'$\Phi$'
                    y_axis_label = r'$\phi_0$'

                    ax.set_xlabel(x_axis_label, fontsize=26)
                    ax.set_ylabel(y_axis_label, fontsize=26)
                    # fig.suptitle(figure_title, fontsize=22)

                    main_service.save_figure(fig, save_folder, save_file_name)


def plot_L_nu_avg_on_fi(i_angle_arr, betta_mu_arr, a_arr, mc2_arr, fi_0_arr, energy_index=8):
    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    phase = np.linspace(0, 1, config.t_max)

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):
            for a_index in range(len(a_arr)):
                for mc2_index in range(len(mc2_arr)):

                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_subplot(111)

                    L_nu_data = []

                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = get_folder_with_args(i_angle_arr[i_angle_index],
                                                                betta_mu_arr[betta_mu_index],
                                                                mc2_arr[mc2_index], a_arr[a_index],
                                                                fi_0_arr[fi_0_index])

                        L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                        L_nu_data.append(np.mean(L_nu_array))

                    ax.scatter(fi_0_arr, L_nu_data)

                    figure_title = r'$\theta{obs}$' + ('=%d' % i_angle_arr[i_angle_index]) + (r'$^\circ$') + ' ' + \
                                   (r'$\beta_{\mu}$') + ('=%d' % betta_mu_arr[betta_mu_index]) + (r'$^\circ$') + \
                                   'a=%0.2f m=%d' % (a_arr[a_index], mc2_arr[mc2_index])

                    save_folder = 'figs/L_nu_avg_on_phase_to_phi/mc2=%d a=%0.2f/' % (mc2_arr[mc2_index], a_arr[a_index])
                    save_file_name = 'i=%d betta_mu=%d L_nu_avg_on_phase_to_phi.png' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                    x_axis_label = r'$\phi_0$'
                    y_axis_label = r'$L_{\nu}^{avg}$ '

                    ax.set_xlabel(x_axis_label, fontsize=22)
                    ax.set_ylabel(y_axis_label, fontsize=22)
                    fig.suptitle(figure_title, fontsize=16)

                    main_service.save_figure(fig, save_folder, save_file_name)





i_angle_arr = [10 * i for i in range(1, 10)]
betta_mu_arr = [10 * i for i in range(1, 10)]
a_arr = [0.25, 0.65]
mc2_arr = [30, 100]
fi_0_arr = [20 * i for i in range(18)]

energy_step = (config.energy_max / config.energy_min) ** (1 / (config.N_energy - 1))
energy_arr = list(config.energy_min * energy_step ** i for i in range(config.N_energy - 1))
energy_arr.append(config.energy_max)

final_final_array = get_PF_data(i_angle_arr, betta_mu_arr, a_arr, mc2_arr, fi_0_arr)
plot_disp_max_mean(i_angle_arr, betta_mu_arr, a_arr, mc2_arr, fi_0_arr, final_final_array=final_final_array)
