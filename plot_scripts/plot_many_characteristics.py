import numpy as np
import multiprocessing as mp
from itertools import repeat
import time
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import scienceplots
from matplotlib import ticker

import geometricTask.matrix as matrix
import config
import accretionColumnService
from accretionColumn import AccretionColumn
import vectors
import main_service

plt.style.use(['science', 'notebook', 'grid'])  # для красивых графиков

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')


# jet =  mpl.colormaps['jet']
# energy_index = 8 - 4.7 keV


def make_new_phi(a_portion, fi_0_arr):
    # new_phi_0 - перевод из старых в новые (чтобы 0 соответствовал старым)
    for k in range(len(fi_0_arr)):
        fi_0_arr[k] = (config.fi_0_dict[a_portion] + 20 * (k)) % 360
        # fi_0_arr[k] = (config.fi_0_dict[a_portion] + fi_0_arr[k]) % 360


def make_new_phi_0(a_portion, fi_0):
    return (config.fi_0_dict[a_portion] + fi_0) % 360


marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}


def get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8):
    '''достаю массивы PF для указанных входных аргументов'''

    folder = 'L_nu/'
    if config.new_magnet_lines_flag:
        folder = 'scattered_on_magnet_lines/' + folder

    file_name = 'PF.txt'
    # такой порядок: a, mc2, i, beta, fi - для того, чтобы делать colormesh x = beta, y=i; по fi делаем усреднение
    final_final_array = np.zeros((len(a_portion_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr), len(fi_0_arr)))

    for a_index in range(len(a_portion_arr)):
        for mc2_index in range(len(mc2_arr)):
            for i_angle_index in range(len(i_angle_arr)):
                for betta_mu_index in range(len(betta_mu_arr)):
                    final_array = []
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = config.get_folder_with_args(i_angle_arr[i_angle_index],
                                                                       betta_mu_arr[betta_mu_index],
                                                                       mc2_arr[mc2_index], a_portion_arr[a_index],
                                                                       fi_0_arr[fi_0_index])

                        PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
                        final_array.append(PF_array[energy_index])

                    # label = r'$\theta{obs}$' + ('=%d' % obs_i_angle_deg) + (r'$\beta_{\mu}$') + ('=%d' % betta_mu_deg) + \
                    #         ('a=%0.2f' % a_portion) + (r'$\dot{m}$') + ('=%d' % M_rate_c2_Led)

                    # label = 'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led)
                    final_final_array[a_index][mc2_index][i_angle_index][betta_mu_index] = final_array

    return final_final_array


def get_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8,
                      final_final_array=np.array([])):
    ''' Дисперсия среднне и max - min colormesh - получаю распределения - усреднены по fi_0'''

    # if not isinstance(final_final_array, np.ndarray):
    if not final_final_array.size:
        final_final_array = get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr)

    dispersion_arr = np.zeros((len(a_portion_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr)))
    mean_arr = np.zeros((len(a_portion_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr)))
    max_min_arr = np.zeros((len(a_portion_arr), len(mc2_arr), len(i_angle_arr), len(betta_mu_arr)))

    for a_index in range(len(a_portion_arr)):
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


def plot_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8,
                       final_final_array=np.array([])):
    # функция для сохранения colormesh
    '''рисую mean, disp, max(PF)-min(PF)
    2D рисунки - по ox - beta_mu, oy - i_obs, цвет - значение
    1 график - их абсолютные
    2 график - с границами max min для каждого рисунка
    3 график - без картинки с max-min'''

    def plot_and_save_col_bar(n_rows, n_cols, save_file_name, limits_arr=[]):
        # вспомогательная функция для отрисовки и сохранения графиков
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

            plt.colorbar(col_bar, ax=axes[i], pad=0.01)

        fig.tight_layout()
        # plt.show()
        main_service.save_figure(fig, save_folder, save_file_name)

    # получаем массивы для отрисовки
    dispersion_arr, mean_arr, max_min_arr = get_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr,
                                                              fi_0_arr, final_final_array=final_final_array)
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

    for a_index in range(len(a_portion_arr)):
        for mc2_index in range(len(mc2_arr)):
            x_axis_label = r'$\beta_{\mu}$'
            y_axis_label = r'$\theta_{obs}$'

            save_folder = config.PROJECT_DIR + 'figs/colormesh_mean_disp/mc2=%d/a=%0.2f/' % (
                mc2_arr[mc2_index], a_portion_arr[a_index])

            save_file_name = 'colormesh_mean_disp_fig.png'
            plot_and_save_col_bar(1, 3, save_file_name)

            save_file_name = 'colormesh_mean_disp_without_fig.png'
            plot_and_save_col_bar(1, 2, save_file_name)

            # тут с лимитами по max min
            save_file_name = 'colormesh_mean_disp_fig_fixed.png'
            plot_and_save_col_bar(1, 3, save_file_name, limits_arr)


def plot_L_nu(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8,
              final_final_array=np.array([])):
    #  PF(L_nu) - pf берутся среднее по фазе, из-за fi_0 - множество точек
    ''' добавляю в мапу: ключ - mean(L_nu_array) - усреднение по фазе, значения - PF
     1 график - точки раскрашены по fi_0
     2 график - комбинация a и mc2'''

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    # если не передали - то создадим
    # if not isinstance(final_final_array, np.ndarray):
    if not final_final_array.size:
        final_final_array = get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr)

    # x = np.zeros((len(a_portion_arr), len(mc2)))
    # y = np.zeros((len(a_portion_arr), len(mc2)))

    markers_dict = {0: 'x', 1: '+', 2: 'o', 3: '.'}

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):

            fig = plt.figure(figsize=(12, 6))
            ax = fig.add_subplot(111)

            fig1 = plt.figure(figsize=(12, 6))
            ax1 = fig1.add_subplot(111)

            for a_index in range(len(a_portion_arr)):
                for mc2_index in range(len(mc2_arr)):
                    L_nu_data_dict = {}
                    color_dict = {}
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = config.get_folder_with_args(i_angle_arr[i_angle_index],
                                                                       betta_mu_arr[betta_mu_index], mc2_arr[mc2_index],
                                                                       a_portion_arr[a_index], fi_0_arr[fi_0_index])

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

                    label = 'a=%0.2f m=%d' % (a_portion_arr[a_index], mc2_arr[mc2_index])
                    ax.scatter(x, y, marker=marker_dict[a_index * len(mc2_arr) + mc2_index], color=cm.jet(colors),
                               label=label)

                    ax1.scatter(x, y, marker=marker_dict[a_index * len(mc2_arr) + mc2_index], label=label)

                x_axis_label = r'$L_{\nu}$'
                y_axis_label = r'$PF_{' + '%.2f' % energy_arr[energy_index] + r'}$'
                figure_title = r'$\theta{obs}$' + ('=%d' % i_angle_arr[i_angle_index]) + (r'$^\circ$') + ' ' + \
                               (r'$\beta_{\mu}$') + ('=%d' % betta_mu_arr[betta_mu_index]) + (r'$^\circ$')

                ax.set_xlabel(x_axis_label, fontsize=24)
                ax.set_ylabel(y_axis_label, fontsize=24)
                fig.suptitle(figure_title, fontsize=24)

                ax.set_xscale('log')
                ax.legend()

                save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a/colors/'
                save_file_name = 'i=%d betta_mu=%d PF_to_L_nu_color.png' % (
                    i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])
                main_service.save_figure(fig, save_folder, save_file_name)

                # figure_title = 'i=%d betta_mu=%d' % (obs_i_angle_deg, betta_mu_deg)

                ax1.set_xlabel(x_axis_label, fontsize=24)
                ax1.set_ylabel(y_axis_label, fontsize=24)
                fig1.suptitle(figure_title, fontsize=24)

                ax1.set_xscale('log')
                ax1.legend()

                save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a/noncolors/'
                save_file_name = 'i=%d betta_mu=%d PF_to_L_nu.png' % (
                    i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])
                main_service.save_figure(fig1, save_folder, save_file_name)


def plot_L_nu_iso(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8):
    """2 мерный график nu_L_nu от fi_0, Phi (фазы)"""

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    phase = np.linspace(0, 1, config.t_max)

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):
            for a_index in range(len(a_portion_arr)):
                for mc2_index in range(len(mc2_arr)):

                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_subplot(111)

                    L_nu_data = []
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = config.get_folder_with_args(i_angle_arr[i_angle_index],
                                                                       betta_mu_arr[betta_mu_index],
                                                                       mc2_arr[mc2_index], a_portion_arr[a_index],
                                                                       fi_0_arr[fi_0_index])

                        L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                        L_nu_data.append(L_nu_array)

                    im = ax.contourf(phase, fi_0_arr, L_nu_data, levels=30)
                    clb = plt.colorbar(im, pad=0.01)
                    clb.set_label(r'$\nu L_{\nu}$' + ' [erg/s]', fontsize=26)

                    figure_title = r'$L_{\nu}$ ' + 'i=%d betta_mu=%d a=%0.2f m=%d' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index], a_portion_arr[a_index],
                        mc2_arr[mc2_index])

                    save_folder = config.PROJECT_DIR + 'figs/L_nu_to_phi_to_phase/mc2=%d a=%0.2f/' % (
                        mc2_arr[mc2_index], a_portion_arr[a_index])
                    save_file_name = 'i=%d betta_mu=%d L_nu_to_phi_to_phase.png' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                    x_axis_label = r'$\Phi$'
                    y_axis_label = r'$\phi_0$'

                    ax.set_xlabel(x_axis_label, fontsize=26)
                    ax.set_ylabel(y_axis_label, fontsize=26)
                    # fig.suptitle(figure_title, fontsize=22)

                    main_service.save_figure(fig, save_folder, save_file_name)


def plot_L_nu_avg_on_fi(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8):
    # L_nu_avg_on_phase(fi_0)
    '''Иду в цикле, добавляю в массив среднее по фазе nu_L_nu
    рисую - по ох fi_0; по оу - nu_L_nu - этот массив'''

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    phase = np.linspace(0, 1, config.t_max)

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):
            for a_index in range(len(a_portion_arr)):
                for mc2_index in range(len(mc2_arr)):

                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_subplot(111)

                    L_nu_data = []

                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = config.get_folder_with_args(i_angle_arr[i_angle_index],
                                                                       betta_mu_arr[betta_mu_index],
                                                                       mc2_arr[mc2_index], a_portion_arr[a_index],
                                                                       fi_0_arr[fi_0_index])

                        L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                        L_nu_data.append(np.mean(L_nu_array))

                    ax.scatter(fi_0_arr, L_nu_data)

                    figure_title = r'$\theta{obs}$' + ('=%d' % i_angle_arr[i_angle_index]) + (r'$^\circ$') + ' ' + \
                                   (r'$\beta_{\mu}$') + ('=%d' % betta_mu_arr[betta_mu_index]) + (r'$^\circ$') + \
                                   'a=%0.2f m=%d' % (a_portion_arr[a_index], mc2_arr[mc2_index])

                    save_folder = config.PROJECT_DIR + 'figs/L_nu_avg_on_phase_to_phi/mc2=%d a=%0.2f/' % (
                        mc2_arr[mc2_index], a_portion_arr[a_index])
                    save_file_name = 'i=%d betta_mu=%d L_nu_avg_on_phase_to_phi.png' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                    x_axis_label = r'$\phi_0$'
                    y_axis_label = r'$L_{\nu}^{avg}$ '

                    ax.set_xlabel(x_axis_label, fontsize=22)
                    ax.set_ylabel(y_axis_label, fontsize=22)
                    fig.suptitle(figure_title, fontsize=16)

                    main_service.save_figure(fig, save_folder, save_file_name)


def plot_masses_PF_L_nu(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8):
    # PF(L_nu) много точек, берутся линии по mc, a
    ''' рисует графики PF от nu_L_nu (nu_L_nu усреднили по фазе)

    сначала в цикле читаю PF, заношу в 2D массив
    потом в цикле по fi_0 в мапу - ключ среднее по фазе nu_Lnu значение - значение PF массив
    одновременно для L_nu запоминаю fi_0, чтобы окрасить в цвет (1 график)

    рисую 3 графика ox - nu_L_nu, oy - PF:
    1 - множество точек для разных mc2, a - папки по a - точки окрашены по fi_0
    2 - множество точек для разных mc2, a - папки по a - точки окрашены для комбинаций mc2, a
    3 - множество точек для разных mc2, a - все точки - точки окрашены для комбинаций mc2, a'''

    final_final_array = np.zeros((len(a_portion_arr), len(mc2_arr), len(fi_0_arr)))

    folder = 'L_nu/'
    file_name = 'PF.txt'
    bounds = fi_0_arr.copy()
    for a_index in range(len(a_portion_arr)):
        # новые fi_0
        make_new_phi(a_portion_arr[a_index], fi_0_arr)
        # new_phi_0
        for mc2_index in range(len(mc2_arr)):
            final_array = []
            for fi_0_index in range(len(fi_0_arr)):
                full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[mc2_index],
                                                               a_portion_arr[a_index],
                                                               fi_0_arr[fi_0_index])

                PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
                final_array.append(PF_array[energy_index])

            final_final_array[a_index][mc2_index] = final_array

    line_color_dict = {0: 'blue', 1: 'green', 2: 'orange', 3: 'red', 4: 'purple', 5: 'black', 6: 'yellow'}
    marker_index = 0

    line_style = ['-', '--']

    folder = 'nu_L_nu/txt/'
    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]

    fig2 = plt.figure(figsize=(12, 6))
    ax2 = fig2.add_subplot(111)

    # layout='compressed'
    fig3 = plt.figure(figsize=(12, 6))
    ax3 = fig3.add_subplot(111)
    for a_index in range(len(a_portion_arr)):
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

        fig1 = plt.figure(figsize=(12, 6))
        ax1 = fig1.add_subplot(111)

        line_color_index = 0
        for mc2_index in range(len(mc2_arr)):
            L_nu_data_dict = {}
            color_dict = {}

            full_dict = {}

            # new_phi_0
            make_new_phi(a_portion_arr[a_index], fi_0_arr)

            for fi_0_index in range(len(fi_0_arr)):
                full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[mc2_index],
                                                               a_portion_arr[a_index],
                                                               fi_0_arr[fi_0_index])

                L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                buf = main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/nu_L_nu/',
                                                     'bot_column_scatter_nu_L_nu.txt')[energy_index]
                if not np.isnan(buf).any():
                    L_nu_array += buf

                buf = main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/nu_L_nu/',
                                                     'top_column_scatter_nu_L_nu.txt')[energy_index]
                if not np.isnan(buf).any():
                    L_nu_array += buf

                L_nu_data_dict[np.mean(L_nu_array)] = final_final_array[a_index][mc2_index][fi_0_index]

                color_dict[np.mean(L_nu_array)] = fi_0_arr[fi_0_index]

                full_dict[np.mean(L_nu_array)] = (final_final_array[a_index][mc2_index][fi_0_index], bounds[fi_0_index])

                # L_nu_data_dict[L_nu_array[L_nu_array_index]] = final_final_array[i_angle_index][betta_mu_index][
                #     i]
            lists = sorted(full_dict.items(), key=lambda item: item[1])
            x_sort, y_sort = zip(*lists)
            y_sort, colors_sort = zip(*y_sort)

            lists = sorted(full_dict.items(), key=lambda item: item[1][1])
            x_sort_fi_0, y_sort_fi_0 = zip(*lists)
            y_sort_fi_0, _ = zip(*y_sort_fi_0)

            lists = sorted(L_nu_data_dict.items())  # sorted by key, return a list of tuples
            x, y = zip(*lists)  # unpack a list of pairs into two tuples

            lists = sorted(color_dict.items())
            buffer, colors = zip(*lists)  # L_nu, fi_0 =
            colors_deg = colors
            x = x_sort
            y = y_sort
            # colors = (np.array(colors) % 180) / 180
            colors = (np.array(colors_sort)) / 180

            # colors = (np.array(colors) % 360) / 360
            # fig = main_service.create_figure(x, y, x_axis_label=r'$L_{\nu}$', y_axis_label='PF',
            #                                  figure_title=title, is_y_2d=False)
            label = 'm=%d' % (mc2_arr[mc2_index])

            # ax.scatter(x, y, marker=marker_dict[marker_index % 4], color=cm.jet(colors), label=label)
            # ax1.scatter(x, y, marker=marker_dict[marker_index % 4], label=label)

            fillstyle = 'full'
            if marker_index == 0:
                fillstyle = 'none'

            ax.scatter(x, y, marker=marker_dict[marker_index % 4], color=cm.jet(colors), label=label)
            ax1.scatter(x, y, marker=marker_dict[marker_index % 4], color=line_color_dict[line_color_index % 6])

            if marker_index == 0:
                ax2.scatter(x, y, s=30, facecolors='none', edgecolors=line_color_dict[line_color_index % 6])
                ax3.plot(x_sort_fi_0, y_sort_fi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])
                ax3.scatter(x, y, s=30, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax2.scatter(x, y, marker=marker_dict[marker_index % 4], color=line_color_dict[line_color_index % 6])
                ax3.plot(x_sort_fi_0, y_sort_fi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])
                ax3.scatter(x, y, marker=marker_dict[marker_index % 4], color=cm.jet(colors))
            line_color_index += 1

        marker_index = 3

        x_axis_label = r'$\nu L_{\nu}$' + ' [erg/s]'
        y_axis_label = r'$PF_{' + '%.2f' % energy_arr[energy_index] + r'}$'
        # figure_title = 'i=%d betta_mu=%d a=%0.2f' % (i_angle, betta_mu, a_arr[a_index])
        ax.set_xlabel(x_axis_label, fontsize=24)
        ax.set_ylabel(y_axis_label, fontsize=24)

        ax.set_xscale('log')
        ax.legend()

        # fig.suptitle(figure_title, fontsize=14)

        save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/colors/' + 'a=%0.2f/' % a_portion_arr[a_index]
        save_file_name = 'i=%d betta_mu=%d PF_to_L_nu_color.png' % (i_angle, betta_mu)

        # cmap = mpl.cm.jet
        # norm = mpl.colors.Normalize(vmin=min(colors_deg), vmax=max(colors_deg))
        # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax)
        # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal', label='Some Units')

        cmap = mpl.cm.jet
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                     ax=ax, orientation='vertical', label=r'$\phi_0$', pad=0.01)

        # fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=180), cmap='jet'),
        #              ax=ax, orientation='vertical', label=r'$\phi_0$')

        main_service.save_figure(fig, save_folder, save_file_name)

        ax1.set_xlabel(x_axis_label, fontsize=24)
        ax1.set_ylabel(y_axis_label, fontsize=24)
        # fig1.suptitle(figure_title, fontsize=14)

        ax1.set_xscale('log')
        ax1.legend()

        save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/noncolors/' + 'a=%0.2f/' % a_portion_arr[a_index]
        save_file_name = 'i=%d betta_mu=%d PF_to_L_nu.png' % (i_angle, betta_mu)

        # fig1.colorbar(
        #     mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=180), cmap='jet'),
        #     ax=ax, orientation='vertical', label=r'$\phi_0$')

        main_service.save_figure(fig1, save_folder, save_file_name)

    ax2.set_xlabel(x_axis_label, fontsize=24)
    ax2.set_ylabel(y_axis_label, fontsize=24)
    # fig1.suptitle(figure_title, fontsize=14)

    ax2.set_xscale('log')

    save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/noncolors/'
    save_file_name = 'i=%d betta_mu=%d All_PF_to_L_nu.png' % (i_angle, betta_mu)
    main_service.save_figure(fig2, save_folder, save_file_name)

    ax3.set_xlabel(x_axis_label, fontsize=24)
    ax3.set_ylabel(y_axis_label, fontsize=24)
    ax3.set_xscale('log')

    # fig3.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=180), cmap='jet'), ax=ax3,
    #               orientation='vertical', label=r'$\phi_0$')

    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    fig3.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=ax3, orientation='vertical', label=r'$\phi_0$', pad=0.01)

    save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/colors/'
    save_file_name = 'i=%d betta_mu=%d All_PF_to_L_nu.png' % (i_angle, betta_mu)
    main_service.save_figure(fig3, save_folder, save_file_name)


def plot_masses_PF_L_nu_single_fi_0(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0, energy_index=8):
    # PF(L_nu) много точек, берутся линии по mc, a
    ''' рисует графики PF от nu_L_nu (nu_L_nu усреднили по фазе)

    сначала в цикле читаю PF, заношу в 2D массив
    потом в цикле по fi_0 в мапу - ключ среднее по фазе nu_Lnu значение - значение PF массив
    одновременно для L_nu запоминаю fi_0, чтобы окрасить в цвет (1 график)

    рисую 3 графика ox - nu_L_nu, oy - PF:
    1 - множество точек для разных mc2, a - папки по a - точки окрашены по fi_0
    2 - множество точек для разных mc2, a - папки по a - точки окрашены для комбинаций mc2, a
    3 - множество точек для разных mc2, a - все точки - точки окрашены для комбинаций mc2, a'''

    final_final_array = np.zeros((len(a_portion_arr), len(mc2_arr), fi_0))

    folder = 'L_nu/'
    file_name = 'PF.txt'

    for a_index in range(len(a_portion_arr)):
        # новые fi_0
        fi_0_new = make_new_phi_0(a_portion_arr[a_index], fi_0)
        # new_phi_0
        for mc2_index in range(len(mc2_arr)):
            final_array = []

            full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[mc2_index],
                                                           a_portion_arr[a_index],
                                                           fi_0_new)

            PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
            final_array.append(PF_array[energy_index])

            final_final_array[a_index][mc2_index] = final_array

    line_color_dict = {0: 'blue', 1: 'green', 2: 'orange', 3: 'red', 4: 'purple', 5: 'black', 6: 'yellow'}
    marker_index = 0

    folder = 'nu_L_nu/txt/'
    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]

    fig2 = plt.figure(figsize=(12, 6))
    ax2 = fig2.add_subplot(111)

    fig3 = plt.figure(figsize=(12, 6))
    ax3 = fig3.add_subplot(111)
    for a_index in range(len(a_portion_arr)):
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

        fig1 = plt.figure(figsize=(12, 6))
        ax1 = fig1.add_subplot(111)

        line_color_index = 0
        for mc2_index in range(len(mc2_arr)):
            L_nu_data_dict = {}
            color_dict = {}
            # new_phi_0
            make_new_phi(a_portion_arr[a_index], fi_0_arr)

            for fi_0_index in range(len(fi_0_arr)):
                full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[mc2_index],
                                                               a_portion_arr[a_index],
                                                               fi_0_arr[fi_0_index])

                L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                L_nu_array += main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/nu_L_nu/',
                                                             'bot_column_scatter_nu_L_nu.txt')[energy_index]
                L_nu_array += main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/nu_L_nu/',
                                                             'top_column_scatter_nu_L_nu.txt')[energy_index]

                L_nu_data_dict[np.mean(L_nu_array)] = final_final_array[a_index][mc2_index][fi_0_index]

                color_dict[np.mean(L_nu_array)] = fi_0_arr[fi_0_index]

                # L_nu_data_dict[L_nu_array[L_nu_array_index]] = final_final_array[i_angle_index][betta_mu_index][
                #     i]

            lists = sorted(L_nu_data_dict.items())  # sorted by key, return a list of tuples
            x, y = zip(*lists)  # unpack a list of pairs into two tuples

            lists = sorted(color_dict.items())
            buffer, colors = zip(*lists)  # L_nu, fi_0 =
            colors_deg = colors
            colors = (np.array(colors) % 180) / 180
            # colors = (np.array(colors) % 360) / 360
            # fig = main_service.create_figure(x, y, x_axis_label=r'$L_{\nu}$', y_axis_label='PF',
            #                                  figure_title=title, is_y_2d=False)

            label = 'm=%d' % (mc2_arr[mc2_index])

            # ax.scatter(x, y, marker=marker_dict[marker_index % 4], color=cm.jet(colors), label=label)
            # ax1.scatter(x, y, marker=marker_dict[marker_index % 4], label=label)

            fillstyle = 'full'
            if marker_index == 0:
                fillstyle = 'none'

            ax.scatter(x, y, marker=marker_dict[marker_index % 4], color=cm.jet(colors), label=label)
            ax1.scatter(x, y, marker=marker_dict[marker_index % 4], color=line_color_dict[line_color_index % 6])

            if marker_index == 0:
                ax2.scatter(x, y, s=30, facecolors='none', edgecolors=line_color_dict[line_color_index % 6])
                ax3.scatter(x, y, s=30, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax2.scatter(x, y, marker=marker_dict[marker_index % 4], color=line_color_dict[line_color_index % 6])
                ax3.scatter(x, y, marker=marker_dict[marker_index % 4], color=cm.jet(colors))
            line_color_index += 1

        marker_index = 3

        x_axis_label = r'$\nu L_{\nu}$' + ' [erg/s]'
        y_axis_label = r'$PF_{' + '%.2f' % energy_arr[energy_index] + r'}$'
        # figure_title = 'i=%d betta_mu=%d a=%0.2f' % (i_angle, betta_mu, a_arr[a_index])
        ax.set_xlabel(x_axis_label, fontsize=24)
        ax.set_ylabel(y_axis_label, fontsize=24)

        ax.set_xscale('log')
        ax.legend()

        # fig.suptitle(figure_title, fontsize=14)

        save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/colors/' + 'a=%0.2f/' % a_portion_arr[a_index]
        save_file_name = 'i=%d betta_mu=%d PF_to_L_nu_color.png' % (i_angle, betta_mu)

        # cmap = mpl.cm.jet
        # norm = mpl.colors.Normalize(vmin=min(colors_deg), vmax=max(colors_deg))
        # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax)
        # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal', label='Some Units')

        cmap = mpl.cm.jet
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                     ax=ax, orientation='vertical', label=r'$\phi_0$')

        # fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=180), cmap='jet'),
        #              ax=ax, orientation='vertical', label=r'$\phi_0$')

        main_service.save_figure(fig, save_folder, save_file_name)

        ax1.set_xlabel(x_axis_label, fontsize=24)
        ax1.set_ylabel(y_axis_label, fontsize=24)
        # fig1.suptitle(figure_title, fontsize=14)

        ax1.set_xscale('log')
        ax1.legend()

        save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/noncolors/' + 'a=%0.2f/' % a_portion_arr[a_index]
        save_file_name = 'i=%d betta_mu=%d PF_to_L_nu.png' % (i_angle, betta_mu)

        # fig1.colorbar(
        #     mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=180), cmap='jet'),
        #     ax=ax, orientation='vertical', label=r'$\phi_0$')

        main_service.save_figure(fig1, save_folder, save_file_name)

    ax2.set_xlabel(x_axis_label, fontsize=24)
    ax2.set_ylabel(y_axis_label, fontsize=24)
    # fig1.suptitle(figure_title, fontsize=14)

    ax2.set_xscale('log')

    save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/noncolors/'
    save_file_name = 'i=%d betta_mu=%d All_PF_to_L_nu.png' % (i_angle, betta_mu)
    main_service.save_figure(fig2, save_folder, save_file_name)

    ax3.set_xlabel(x_axis_label, fontsize=24)
    ax3.set_ylabel(y_axis_label, fontsize=24)
    ax3.set_xscale('log')

    # fig3.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=180), cmap='jet'), ax=ax3,
    #               orientation='vertical', label=r'$\phi_0$')

    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    fig3.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=ax3, orientation='vertical', label=r'$\phi_0$')

    save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/colors/'
    save_file_name = 'i=%d betta_mu=%d All_PF_to_L_nu.png' % (i_angle, betta_mu)
    main_service.save_figure(fig3, save_folder, save_file_name)


def plot_L_nu_flag_particular_fi_0(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8):
    '''считываем нужные PF -3х мерный массив [a][mc2][fi_0]
    fig1 - выводим все линии и mc2 и a - all fi_0
    fig2 - отдельно для каждого значения fi_0 мнжество точек при разных mc2, a
    '''

    final_final_array = np.zeros((len(a_portion_arr), len(mc2_arr), len(fi_0_arr)))

    folder = 'nu_L_nu/'
    file_name = 'PF.txt'

    for a_index in range(len(a_portion_arr)):
        for mc2_index in range(len(mc2_arr)):
            final_array = []
            for fi_0_index in range(len(fi_0_arr)):
                full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[mc2_index],
                                                               a_portion_arr[a_index],
                                                               fi_0_arr[fi_0_index])

                PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
                final_array.append(PF_array[energy_index])

            final_final_array[a_index][mc2_index] = final_array

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    # x = np.zeros((len(a_portion_arr), len(mc2)))
    # y = np.zeros((len(a_portion_arr), len(mc2)))

    # marker_dict = {0: 'x', 1: '+', 2: 'o', 3: '.'}
    # marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}
    line_color_dict = {0: 'blue', 1: 'green', 2: 'orange', 3: 'red', 4: 'purple', 5: 'black'}
    colors = np.array(fi_0_arr) / 340

    fig1 = plt.figure(figsize=(12, 6))
    ax1 = fig1.add_subplot(111)

    for fi_0_index in range(len(fi_0_arr)):

        fig2 = plt.figure(figsize=(12, 6))
        ax2 = fig2.add_subplot(111)

        marker_index = 0
        for a_index in range(len(a_portion_arr)):
            line_color_index = 0
            L_nu_data_dict_for_a = {}
            for mc2_index in range(len(mc2_arr)):
                L_nu_data_dict = {}

                full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[mc2_index],
                                                               a_portion_arr[a_index],
                                                               fi_0_arr[fi_0_index])

                L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                L_nu_data_dict[np.mean(L_nu_array)] = final_final_array[a_index][mc2_index][fi_0_index]

                # вроде как добавляю ключ+значение в dict  L_nu_data_dict_for_a - чтобы сохранить.
                # так как L_nu_data_dict рисуем в цикле по mc2
                L_nu_data_dict_for_a.update(L_nu_data_dict)

                # L_nu_data_dict[L_nu_array[L_nu_array_index]] = final_final_array[i_angle_index][betta_mu_index][
                #     i]

                print(full_file_folder)
                # print(L_nu_data_dict)
                # print(L_nu_data_dict)
                lists = sorted(L_nu_data_dict.items())  # sorted by key, return a list of tuples
                x, y = zip(*lists)  # unpack a list of pairs into two tuples

                fillstyle = 'full'
                if marker_index == 0:
                    fillstyle = 'none'

                if marker_index == 0:
                    ax2.scatter(x, y, s=30, facecolors='none', edgecolors=line_color_dict[line_color_index])
                else:
                    ax2.scatter(x, y, marker=marker_dict[marker_index % 4], color=line_color_dict[line_color_index])

                line_color_index += 1

            marker_index = 3

            lists = sorted(L_nu_data_dict_for_a.items())  # sorted by key, return a list of tuples
            x, y = zip(*lists)  # unpack a list of pairs into two tuples

            if a_index == 0:
                ax1.plot(x, y, label='fi_0 = %d' % fi_0_arr[fi_0_index], color=cm.jet(colors)[fi_0_index])
            else:
                ax1.plot(x, y, color=cm.jet(colors)[fi_0_index])

        x_axis_label = r'$\nu L_{\nu}$' + ' [erg/s]'
        y_axis_label = r'$PF_{' + '%.2f' % energy_arr[energy_index] + r'}$'

        figure_title = 'i=%d betta_mu=%d a=%0.2f' % (i_angle, betta_mu, a_portion_arr[a_index])

        ax2.set_xlabel(x_axis_label, fontsize=22)
        ax2.set_ylabel(y_axis_label, fontsize=22)

        ax2.set_xscale('log')

        save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/fi_0_particular/' \
                      + 'i=%d betta_mu=%d/' % (i_angle, betta_mu)
        save_file_name = 'fi_0=%d' % fi_0_arr[fi_0_index]
        main_service.save_figure(fig2, save_folder, save_file_name)

        x_axis_label = r'$\nu L_{\nu}$' + ' [erg/s]'
        y_axis_label = r'$PF_{' + '%.2f' % energy_arr[energy_index] + r'}$'

        ax1.set_xlabel(x_axis_label, fontsize=22)
        ax1.set_ylabel(y_axis_label, fontsize=22)

        ax1.set_xscale('log')
        ax1.legend(ncol=3, fontsize=10, framealpha=0.2)

        save_folder = config.PROJECT_DIR + 'figs/PF_to_L_nu/mc_a_many/fi_0_particular/' + \
                      'i=%d betta_mu=%d/' % (i_angle, betta_mu)
        save_file_name = 'all fi_0'
        main_service.save_figure(fig1, save_folder, save_file_name)


def plot_L_to_a_portion(i_angle, betta_mu, mc2, a_portion_arr, fi_0):
    save_folder = config.PROJECT_DIR + 'figs/L_to_a/' + f'i={i_angle} betta_mu={betta_mu}/' + f'mc2={mc2}/' f'fi_0={fi_0}/'

    data_array = [0] * len(a_portion_arr)
    buf_fi_0 = fi_0
    for i in range(len(a_portion_arr)):
        if buf_fi_0 == 'new':
            fi_0 = config.fi_0_dict[a_portion_arr[i]]
        else:
            fi_0 = (config.fi_0_dict[a_portion_arr[i]] + 90) % 360
        file_path = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion_arr[i], fi_0)

        file_name = "total_luminosity_of_surfaces.txt"
        data_array[i] = main_service.load_arr_from_txt(file_path, file_name)[4]

        file_path = file_path + 'scattered_on_magnet_lines/'

        file_name = "scattered_energy_bot.txt"
        buf = main_service.load_arr_from_txt(file_path, file_name)
        if not np.isnan(buf).any():
            data_array[i] += buf

        file_name = "scattered_energy_top.txt"
        buf = main_service.load_arr_from_txt(file_path, file_name)
        if not np.isnan(buf).any():
            data_array[i] += buf

    phase = np.linspace(0, 2, 2 * config.t_max)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # нормировка на L_nu_avg
    data_to_plot = [0] * len(a_portion_arr)

    for i, arr in enumerate(data_array):
        data_to_plot[i] = (arr / max(arr))
        data_to_plot[i] = main_service.extend_arr_for_phase(data_to_plot[i])

    im = ax.pcolormesh(phase, np.linspace(0, 1, len(a_portion_arr)), data_to_plot)
    ax.set_yticks(np.linspace(0, 1, len(a_portion_arr)), np.round(a_portion_arr, 2))

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$a$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{iso} \cdot max(L_{iso})^{-1}$', fontsize=24)

    file_name = 'map_contour' + '.png'
    main_service.save_figure(fig, save_folder, file_name)


def plot_L_to_accr_rate(i_angle, betta_mu, mc2_arr, a_portion, fi_0):
    save_folder = config.PROJECT_DIR + 'figs/L_to_mass/' + f'i={i_angle} betta_mu={betta_mu}/' + f'a={a_portion} fi_0={fi_0}/'

    working_folder = config.full_file_folder

    L_x_arr = [0] * len(mc2_arr)
    for i, mc2 in enumerate(mc2_arr):
        with open(working_folder + 'save_values.txt') as f:
            lines = f.readlines()
            L_x_arr[i] = float(lines[3][12:20]) * 10 ** float(lines[3][27:29])
            # total L_x = 4.396383 * 10**38 - 12 это индекс начала числа, 27-29 это степень 10

    data_array = [0] * len(mc2_arr)
    for i in range(len(mc2_arr)):
        file_path = config.get_folder_with_args(i_angle, betta_mu, mc2_arr[i], a_portion, fi_0)

        file_name = "total_luminosity_of_surfaces.txt"
        data_array[i] = main_service.load_arr_from_txt(file_path, file_name)[4]

        file_path = file_path + 'scattered_on_magnet_lines/'

        file_name = "scattered_energy_bot.txt"
        buf = main_service.load_arr_from_txt(file_path, file_name)
        if not np.isnan(buf).any():
            data_array[i] += buf

        file_name = "scattered_energy_top.txt"
        buf = main_service.load_arr_from_txt(file_path, file_name)
        if not np.isnan(buf).any():
            data_array[i] += buf

    phase = np.linspace(0, 2, 2 * config.t_max)

    # fig, ax = plt.subplots()

    # нормировка на L_nu_avg
    data_to_plot = [0] * len(mc2_arr)
    L_x = max(L_x_arr)
    for i, arr in enumerate(data_array):
        data_to_plot[i] = (arr / max(arr))
        data_to_plot[i] = main_service.extend_arr_for_phase(data_to_plot[i])

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)

    # data2 = [main_service.extend_arr_for_phase(item) for item in data_array]
    # im = ax.pcolormesh(phase, [np.mean(item) for item in data_array], data2, alpha=0)
    dummy = np.zeros(90)
    dummy1 = [np.mean(item) for item in data_array]
    for i,item in enumerate(dummy1):
        dummy[i] = item

    ax.plot(phase, dummy, alpha=0)

    # im = ax.pcolormesh(phase, np.linspace(0, 1, len(mc2_arr)), data_to_plot)

    def fmt_two_digits(x, pos):
        return f'{x:.2e}'

    ax.set_yscale('log')
    # ax.set_yticks(np.linspace(0, 1, len(mc2_arr)), labels=[fmt_two_digits(np.mean(item), 0) for item in data_array])

    secax1 = ax.twinx()
    im = secax1.pcolormesh(phase, np.linspace(0, 1, len(mc2_arr)), data_to_plot)
    secax1.set_yticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)

    # secax = ax.secondary_yaxis('right')
    # im = secax.pcolormesh(phase, np.linspace(0, 1, len(mc2_arr)), data_to_plot)
    # secax.set_yticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)

    # ax.plot_surface(phase, mc2_arr, data_to_plot)
    # triang = tri.Triangulation(phase, mc2_arr)
    # im = ax.tricontour(triang, data_to_plot, levels=14, linewidths=0.5, colors='k')

    # im = ax.contourf(phase, mc2_arr, data_to_plot, levels=30)
    # cmap = 'cividis', cmap='Spectral', cmap='Spectral_r'

    # ax.set_yticks([80, 100, 120], labels=['really, really, really', 'long', 'labels'])

    # plt.yticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)

    # ratio = 1.0
    # x_left, x_right = ax.get_xlim()
    # y_low, y_high = ax.get_ylim()
    # ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\dot{m}$'
    y_axis_label = r'$mean L_{iso} [erg/s]$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    y_axis_label = r'$\dot{m}$'
    secax1.set_ylabel(y_axis_label, fontsize=24)
    # secax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

    # secax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_two_digits))
    # secax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

    # secax.ticklabel_format()
    # secax.set_ylabel('radians')

    # ax2 = ax.twinx()
    # ax2.set_yticks(np.linspace(0, 1, len(mc2_arr)), labels=[np.mean(item) for item in data_to_plot])

    # fig_title = 'L_iso/L_x'
    # fig.suptitle(fig_title, fontsize=14)

    # fig_title = r'$L_{iso}/L_{x}$ ' + title
    # fig.suptitle(fig_title, fontsize=14)

    clb = plt.colorbar(im, pad=0.15, format="{x:.2}")
    clb.set_label(r'$L_{iso} \cdot max(L_{iso})^{-1}$', fontsize=24)

    # clb.ax.set_title(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)

    file_name = 'map_contour' + '.png'
    main_service.save_figure(fig, save_folder, file_name)

    L_x = max(L_x_arr)
    for i, arr in enumerate(data_array):
        data_to_plot[i] = (arr / L_x)
        data_to_plot[i] = main_service.extend_arr_for_phase(data_to_plot[i])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    im = ax.pcolormesh(phase, np.linspace(0, 1, len(mc2_arr)), data_to_plot)
    ax.set_yticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\dot{m}$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.02, format="{x:.2}")
    clb.set_label(r'$L_{iso} \cdot L_x^{-1}$', fontsize=24)
    # clb.ax.set_title(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)

    file_name = 'map_contour_L_x' + '.png'
    main_service.save_figure(fig, save_folder, file_name)


def plot_geometric_contribution(i_angle, betta_mu, mc2, a_portion, fi_0, mu_power=31, opacity_above_shock=0):
    var_path = f'i={i_angle} betta_mu={betta_mu}/' + f'mc2={mc2}/' + f'a={a_portion} fi_0={fi_0}/'
    save_folder = config.PROJECT_DIR + 'figs/geom_contr/' + var_path

    file_path = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion, fi_0)

    file_name = "total_luminosity_of_surfaces.txt"
    data_array_tau = main_service.load_arr_from_txt(file_path, file_name)[4]

    folder = 'scattered_on_magnet_lines/'

    file_name = 'scattered_energy_bot.txt'
    bot_scatter = main_service.load_arr_from_txt(file_path + folder, file_name)
    bot_scatter = main_service.extend_arr_for_phase(bot_scatter)

    file_name = 'scattered_energy_top.txt'
    top_scatter = main_service.load_arr_from_txt(file_path + folder, file_name)
    top_scatter = main_service.extend_arr_for_phase(top_scatter)

    file_path = config.PROJECT_DIR_ORIGIN + 'new_data/' + f'mu=0.1e{mu_power}/opacity={opacity_above_shock:.2f}/' + 'figs/loop/' + var_path
    file_name = "total_luminosity_of_surfaces.txt"
    data_array = main_service.load_arr_from_txt(file_path, file_name)[4]

    fig = plt.figure(figsize=(21, 10))
    ax = fig.add_subplot(111)

    phase = np.linspace(0, 2, 2 * config.t_max)

    arr_to_plt = [0] * 3
    arr_to_plt[0] = main_service.extend_arr_for_phase(data_array)
    arr_to_plt[1] = main_service.extend_arr_for_phase(data_array_tau)
    arr_to_plt[2] = main_service.extend_arr_for_phase(data_array_tau) + bot_scatter + top_scatter

    labels_arr = [r'$opacity=0$', r'$with \tau$', r'$with \tau \, and \, scatter$']

    for i in range(3):
        ax.plot(phase, arr_to_plt[i], label=labels_arr[i])

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$L_{\rm{iso}}$' + ' [erg/s]'
    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    ax.legend()

    file_name = 'geometric_contribution.png'
    main_service.save_figure(fig, save_folder, file_name)


def plot_L_nu_iso_for_new_phi(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8):
    """2 мерный график nu_L_nu от fi_0, Phi (фазы)"""

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_arr[energy_index]
    folder = 'nu_L_nu/txt/'

    phase = np.linspace(0, 1, config.t_max)

    for i_angle_index in range(len(i_angle_arr)):
        for betta_mu_index in range(len(betta_mu_arr)):
            for a_index in range(len(a_portion_arr)):

                make_new_phi(a_portion_arr[a_index], fi_0_arr)
                for mc2_index in range(len(mc2_arr)):
                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_subplot(111)

                    L_nu_data = [0] * 18
                    for fi_0_index in range(len(fi_0_arr)):
                        full_file_folder = config.get_folder_with_args(i_angle_arr[i_angle_index],
                                                                       betta_mu_arr[betta_mu_index],
                                                                       mc2_arr[mc2_index], a_portion_arr[a_index],
                                                                       fi_0_arr[fi_0_index])

                        L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                        L_nu_array += \
                            main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/nu_L_nu/',
                                                           'bot_column_scatter_nu_L_nu.txt')[energy_index]
                        L_nu_array += \
                            main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/nu_L_nu/',
                                                           'top_column_scatter_nu_L_nu.txt')[energy_index]

                        L_nu_data[fi_0_index] = L_nu_array

                    for i in range(8):
                        L_nu_data[9 + i + 1] = L_nu_data[9 - i - 1][::-1]

                    fi_0_arr_for_plot = [20 * i for i in range(18)]

                    for i in range(18):
                        L_nu_data[i] = main_service.extend_arr_for_phase(L_nu_data[i])

                    phase = np.linspace(0, 2, 2 * config.t_max)
                    im = ax.contourf(phase, fi_0_arr_for_plot, L_nu_data, levels=30)
                    clb = plt.colorbar(im, pad=0.01)
                    clb.set_label(r'$\nu L_{\nu}$' + ' [erg/s]', fontsize=26)

                    figure_title = r'$L_{\nu}$ ' + 'i=%d betta_mu=%d a=%0.2f m=%d' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index], a_portion_arr[a_index],
                        mc2_arr[mc2_index])

                    save_folder = config.PROJECT_DIR + 'figs/L_nu_to_phi_to_phase/mc2=%d a=%0.2f/' % (
                        mc2_arr[mc2_index], a_portion_arr[a_index])
                    save_file_name = 'i=%d betta_mu=%d L_nu_to_phi_to_phase.png' % (
                        i_angle_arr[i_angle_index], betta_mu_arr[betta_mu_index])

                    x_axis_label = r'$\Phi$'
                    y_axis_label = r'$\phi_0$'

                    ax.set_xlabel(x_axis_label, fontsize=26)
                    ax.set_ylabel(y_axis_label, fontsize=26)
                    # fig.suptitle(figure_title, fontsize=22)

                    main_service.save_figure(fig, save_folder, save_file_name)


def plot_L_to_new_fi_0(i_angle, betta_mu, mc2, a_portion, fi_0_arr):
    save_folder = config.PROJECT_DIR + 'figs/L_to_fi_0/' + f'i={i_angle} betta_mu={betta_mu}/' + f'mc2={mc2}/' f'a={a_portion}/'

    file_name = 'total_luminosity_of_surfaces.txt'
    # изменить данные - чтобы читать
    make_new_phi(a_portion, fi_0_arr)

    L_data = [0] * 18
    for fi_0_index in range(len(fi_0_arr)):
        full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion, fi_0_arr[fi_0_index])

        L_array = main_service.load_arr_from_txt(full_file_folder, file_name)[4]
        buf = main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/',
                                             'scattered_energy_top.txt')
        if not np.isnan(buf).any():
            L_array += buf
        buf = main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/',
                                             'scattered_energy_bot.txt')
        if not np.isnan(buf).any():
            L_array += buf

        L_data[fi_0_index] = L_array

    for i in range(8):
        L_data[9 + i + 1] = L_data[9 - i - 1][::-1]
        # симметрично для 180 - fi_0

    for i in range(18):
        L_data[i] = main_service.extend_arr_for_phase(L_data[i])

    fi_0_arr_for_plot = [20 * i for i in range(18)]

    phase = np.linspace(0, 2, 2 * config.t_max)

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    im = ax.contourf(phase, fi_0_arr_for_plot, L_data, levels=30)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{iso}$' + ' [erg/s]', fontsize=26)

    save_file_name = 'map_contour_L_iso' + '.png'

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\phi_0$'

    ax.set_xlabel(x_axis_label, fontsize=26)
    ax.set_ylabel(y_axis_label, fontsize=26)
    # fig.suptitle(figure_title, fontsize=22)

    main_service.save_figure(fig, save_folder, save_file_name)


def plot_PF_contour(mc2, a_portion, fi_0):
    # PF(L_nu) много точек, берутся линии по mc, a

    i_angle_arr = np.linspace(10, 90, 9)
    betta_mu_arr = np.linspace(10, 90, 9)

    final_final_array = np.zeros((len(i_angle_arr), len(betta_mu_arr)))

    folder = 'scattered_on_magnet_lines/L_nu/'
    file_name = 'PF.txt'

    save_file_name = f'a={a_portion}' + ' ' + f'fi_0={fi_0}.png'
    # new_phi_0
    fi_0 = make_new_phi_0(a_portion, fi_0)

    folder = 'scattered_on_magnet_lines/'
    file_name = 'total_luminosity_of_surfaces.txt'

    for i in range(len(i_angle_arr)):
        for j in range(len(betta_mu_arr)):
            full_file_folder = config.get_folder_with_args(i_angle_arr[i], betta_mu_arr[j], mc2, a_portion, fi_0)
            L = main_service.load_arr_from_txt(full_file_folder, file_name)[4]

            buf = main_service.load_arr_from_txt(full_file_folder + folder, 'scattered_energy_bot.txt')
            if not np.isnan(buf).any():
                L += buf

            buf = main_service.load_arr_from_txt(full_file_folder + folder, 'scattered_energy_top.txt')
            if not np.isnan(buf).any():
                L += buf

            PF = (max(L) - min(L)) / (max(L) + min(L))
            final_final_array[i][j] = PF

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)

    x, y = np.meshgrid(betta_mu_arr, i_angle_arr)
    z = final_final_array
    # z = np.transpose(final_final_array)

    # x, y = np.meshgrid(betta_mu_arr, i_angle_arr)
    # z = final_final_array

    cs = ax.contourf(x, y, z)

    ax.set_xlabel(r'$\chi$')
    ax.set_ylabel(r'$\theta_{obs}$')

    cbar = fig.colorbar(cs, pad=0.01)
    cbar.ax.set_ylabel('PF')

    # plt.show()

    save_folder = config.PROJECT_DIR + 'figs/PF_contour/' + f'mc2={mc2}/'
    main_service.save_figure(fig, save_folder, save_file_name)


energy_step = (config.energy_max / config.energy_min) ** (1 / (config.N_energy - 1))
energy_arr = list(config.energy_min * energy_step ** i for i in range(config.N_energy - 1))
energy_arr.append(config.energy_max)

all_flag = False
if all_flag:
    i_angle_arr = [10 * i for i in range(1, 10)]
    betta_mu_arr = [10 * i for i in range(1, 10)]
    mc2_arr = [30, 100]
    a_portion_arr = [0.25, 0.65]
    fi_0_arr = [20 * i for i in range(18)]

    energy_step = (config.energy_max / config.energy_min) ** (1 / (config.N_energy - 1))
    energy_arr = list(config.energy_min * energy_step ** i for i in range(config.N_energy - 1))
    energy_arr.append(config.energy_max)

    final_final_array = get_PF_data(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr)
    plot_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, final_final_array=final_final_array)

i_angle = 80
betta_mu = 80
mc2_arr = [10, 20, 30, 50, 100, 200]
a_portion_arr = [0.25, 0.65]
fi_0_arr = [20 * i for i in range(18)]

# i_angle = 60
# betta_mu = 40
# mc2_arr = [140]
# a_portion_arr = [0.25]
# fi_0_arr = [0]

# ---------------------------------------------------------
# i_angle_arr = [10 * i for i in range(1, 10)]
# betta_mu_arr = [10 * i for i in range(1, 10)]
# mc2_arr = [30, 100]
# a_portion_arr = [0.25, 0.65]
# fi_0_arr = [20 * i for i in range(18)]
#
# plot_disp_max_mean(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr)

# ---------------------------------------------------------
a_portion_arr = [0.22]
i_angle_arr = [20]
betta_mu_arr = [60]
# mc2_arr = np.linspace(10, 130, 13)
# mc2_arr = list(map(int, mc2_arr))
for i_angle in i_angle_arr:
    for betta_mu in betta_mu_arr:
        for a_portion in a_portion_arr:
            # fi_0 = (config.fi_0_dict[a_portion]) % 360
            # # fi_0 = config.fi_0_dict[a_portion]
            # plot_L_to_accr_rate(i_angle, betta_mu, mc2_arr, a_portion, fi_0)
            if a_portion == 0.11:
                mc2_arr = np.linspace(10, 60, 6)
                mc2_arr = list(map(int, mc2_arr))
            elif a_portion == 0.22:
                mc2_arr = np.linspace(10, 110, 11)
                mc2_arr = list(map(int, mc2_arr))
            else:
                mc2_arr = np.linspace(10, 130, 13)
                mc2_arr = list(map(int, mc2_arr))
            fi_0 = (config.fi_0_dict[a_portion]) % 360
            # fi_0 = config.fi_0_dict[a_portion]
            plot_L_to_accr_rate(i_angle, betta_mu, mc2_arr, a_portion, fi_0)

# ---------------------------------------------------------------
# # a_portion_arr = [0.22, 0.33, 0.44, 0.55, 0.66, 0.77]
# a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825]
# i_angle_arr = [20, 40, 60]
# betta_mu_arr = [20, 40, 60]
# mc2_arr = [30, 100]
# #
# # mc2_arr = np.linspace(10, 130, 13)
# # mc2_arr = list(map(int, mc2_arr))
#
# fi_0 = 'new'
#
# for i_angle in i_angle_arr:
#     for betta_mu in betta_mu_arr:
#         for mc2 in mc2_arr:
#             plot_L_to_a_portion(i_angle, betta_mu, mc2, a_portion_arr, fi_0)

# ---------------------------------------------------------
# i_angle_arr = [20, 40, 60]
# betta_mu_arr = [40, 60]
# mc2_arr = [30, 100]
# a_portion_arr = [0.25, 0.65]
# fi_0 = 0
#
# for i_angle in i_angle_arr:
#     for betta_mu in betta_mu_arr:
#         for mc2 in mc2_arr:
#             for a_portion in a_portion_arr:
#                 plot_geometric_contribution(i_angle, betta_mu, mc2, a_portion, fi_0)

# plot_geometric_contribution(i_angle, betta_mu, mc2, a_portion, fi_0)

# -------------------------------------------------------

# i_angle = 60
# betta_mu = 40
#
# mc2_arr = [30, 100]
# # mc2_arr = np.linspace(10, 130, 13)
#
# a_portion_arr = [0.22, 0.66]
#
# fi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
# # fi_0_arr = [0]
#
# plot_masses_PF_L_nu(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8)

# ---------------------------------------------------------
# i_angle_arr = [60]
# betta_mu_arr = [40]
# mc2_arr = [30, 100]
# a_portion_arr = [0.22, 0.66]
#
# fi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
# plot_L_nu_iso_for_new_phi(i_angle_arr, betta_mu_arr, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8)

# -------------------------------------------------------
# i_angle = 40
# betta_mu = 60
# mc2 = 100
# a_portion = 0.66
#
# fi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]


# plot_L_to_new_fi_0(i_angle, betta_mu, mc2, a_portion, fi_0_arr)

# ---------------------------------------------------------

# i_angle = 40
# betta_mu = 60
# mc2 = 100
# a_portion = 0.66
# fi_0 = 0
#
# a_portion_arr = [0.22, 0.44, 0.66]
# mc2_arr = [30, 100]
#
# for mc2 in mc2_arr:
#     for a_portion in a_portion_arr:
#         plot_PF_contour(mc2, a_portion, fi_0)

# plot_masses_PF_L_nu(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0_arr, energy_index=8)
# plot_masses_PF_L_nu(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0_arr)
# plot_L_nu_flag_particular_fi_0(i_angle, betta_mu, mc2_arr, a_portion_arr, fi_0_arr)
