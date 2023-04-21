import numpy as np
import multiprocessing as mp
from itertools import repeat
import time
import matplotlib.pyplot as plt
import pandas as pd

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

plt.style.use(['science', 'notebook', 'grid'])

mc2 = [10, 30, 100]
a_portion_arr = [0.1, 0.25, 0.65]
# a_portion_arr = [0.25]
fi_0 = [20 * i for i in range(18)]
# fi_0 = [0, 20, 40, 60, 100, 160, 200, 220, 260, 300]
i_angle = [10, 20, 30, 40, 50, 60, 70, 80, 90]
betta_mu = [10, 20, 30, 40, 50, 60, 70, 80, 90]
# betta_mu = [30]

obs_i_angle_deg = i_angle[0]
betta_mu_deg = betta_mu[0]
M_rate_c2_Led = mc2[2]
phi_accretion_begin_deg = fi_0[0]
a_portion = a_portion_arr[1]


def get_folder():
    file_folder = 'figs/loop/'
    file_folder_angle = 'i=%d betta_mu=%d/' % (obs_i_angle_deg, betta_mu_deg)
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg)
    full_file_folder = file_folder + file_folder_angle + file_folder_args

    return full_file_folder


full_file_folder = get_folder()

folder = 'nu_L_nu/'
working_folder = full_file_folder + folder

file_name = 'PF.txt'

PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

file_name = 'energy.txt'
energy_array = main_service.load_arr_from_txt(full_file_folder, file_name)

# for i, E in enumerate(energy_array):
#     print(i, E)

energy_index = 8  # 4.7 keV

marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}
folder = 'nu_L_nu/'
file_name = 'PF.txt'

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)

# меняем углы
a_portion = a_portion_arr[2]
M_rate_c2_Led = mc2[1]
buff_marker = 0

final_final_array = np.zeros((len(i_angle), len(betta_mu), len(fi_0)))

bins_size = 6
bins_arr = np.zeros((len(i_angle), len(betta_mu), bins_size))
position_in_bins_arr = np.zeros((len(i_angle), len(betta_mu), len(fi_0)))

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        final_array = []
        for i in range(len(fi_0)):
            obs_i_angle_deg = i_angle[i_angle_index]
            betta_mu_deg = betta_mu[betta_mu_index]

            phi_accretion_begin_deg = fi_0[i]
            if phi_accretion_begin_deg == 360:
                phi_accretion_begin_deg = 0

            full_file_folder = get_folder()
            working_folder = full_file_folder + folder

            PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
            final_array.append(PF_array[energy_index])

            label = 'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led)

        bins_arr[i_angle_index][betta_mu_index] = np.linspace(min(final_array), max(final_array), bins_size)
        position_in_bins_arr[i_angle_index][betta_mu_index] = np.digitize(final_array,
                                                                          bins_arr[i_angle_index][betta_mu_index])

        final_final_array[i_angle_index][betta_mu_index] = final_array
        marker = marker_dict[buff_marker % 4]
        buff_marker += 1
        ax.plot(fi_0, final_array, label=label, marker=marker)

# plt.legend(ncol=3, fontsize=10, framealpha=0.2)
plt.show()

'''Bars'''
bar_flag = False
if bar_flag:
    for i_angle_index in range(len(i_angle)):
        for betta_mu_index in range(len(betta_mu)):

            obs_i_angle_deg = i_angle[i_angle_index]
            betta_mu_deg = betta_mu[betta_mu_index]

            bins_dict = dict((i, 0) for i in range(1, bins_size + 1))

            unique, counts = np.unique(position_in_bins_arr[i_angle_index][betta_mu_index], return_counts=True)
            count_bins = dict(zip(unique, counts))

            for key, value in count_bins.items():
                bins_dict[key] = count_bins[key]

            count_bins = bins_dict

            count_bins[len(count_bins) - 1] = count_bins[len(count_bins) - 1] + count_bins[len(count_bins)]

            labels_arr = []
            for i in range(1, len(bins_arr[i_angle_index][betta_mu_index])):
                labels_arr.append('%0.3f - %0.3f' % (
                    bins_arr[i_angle_index][betta_mu_index][i - 1], bins_arr[i_angle_index][betta_mu_index][i]))

            # fig = plt.figure(figsize=(8, 6))

            plt.bar(range(len(count_bins) - 1), list(count_bins.values())[:-1], align='center', edgecolor='black')
            plt.xticks(range(len(count_bins) - 1), labels_arr, fontsize=12)

            plt.xlabel('PF')
            plt.ylabel('N')

            plt.title('i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led))

            save_folder = 'figs/bars/mc2=%d/a=%0.2f/' % (M_rate_c2_Led, a_portion)
            save_file_name = 'i=%d betta_mu=%d bar.png' % (obs_i_angle_deg, betta_mu_deg)

            main_service.create_file_path(save_folder)
            plt.savefig(save_folder + save_file_name)

            plt.cla()

dispersion_arr = np.zeros((len(i_angle), len(betta_mu)))
mean_arr = np.zeros((len(i_angle), len(betta_mu)))
max_min = np.zeros((len(i_angle), len(betta_mu)))

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        dispersion_arr[i_angle_index][betta_mu_index] = np.var(final_final_array[i_angle_index][betta_mu_index]) ** (
                1 / 2)
        mean_arr[i_angle_index][betta_mu_index] = np.mean(final_final_array[i_angle_index][betta_mu_index])
        max_min[i_angle_index][betta_mu_index] = np.max(final_final_array[i_angle_index][betta_mu_index]) - np.min(
            final_final_array[i_angle_index][betta_mu_index])

'''считываю L_iso'''

L_iso_data_dict = {}

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        obs_i_angle_deg = i_angle[i_angle_index]
        betta_mu_deg = betta_mu[betta_mu_index]

        full_file_folder = get_folder()

        with open(full_file_folder + 'save_values.txt') as f:
            lines = f.readlines()
            # print(lines[3][12:20])
            # print(lines[3][27:29])
            # print(lines)

            L_sio = float(lines[3][12:20]) * 10 ** float(lines[3][27:29])
            # print(L_sio)

        L_iso_data_dict[L_sio] = mean_arr[i_angle_index][betta_mu_index]

# print(L_iso_data_dict)
lists = sorted(L_iso_data_dict.items())  # sorted by key, return a list of tuples
x, y = zip(*lists)  # unpack a list of pairs into two tuples
plt.plot(x, y)
plt.show()

'''L_nu'''

L_nu_flag = True
L_nu_array_index = 0

if L_nu_flag:

    file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt' % energy_array[energy_index]
    folder = 'nu_L_nu/txt/'

    for i_angle_index in range(len(i_angle)):
        for betta_mu_index in range(len(betta_mu)):
            L_nu_data_dict = {}
            for i in range(len(fi_0)):
                obs_i_angle_deg = i_angle[i_angle_index]
                betta_mu_deg = betta_mu[betta_mu_index]

                phi_accretion_begin_deg = fi_0[i]

                full_file_folder = get_folder()

                L_nu_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)

                L_nu_data_dict[np.mean(L_nu_array)] = final_final_array[i_angle_index][betta_mu_index][i]

                # L_nu_data_dict[L_nu_array[L_nu_array_index]] = final_final_array[i_angle_index][betta_mu_index][
                #     i]

            print(full_file_folder)
            # print(L_nu_data_dict)
            # print(L_nu_data_dict)
            lists = sorted(L_nu_data_dict.items())  # sorted by key, return a list of tuples
            x, y = zip(*lists)  # unpack a list of pairs into two tuples

            title = 'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led)

            # fig = main_service.create_figure(x, y, x_axis_label=r'$L_{\nu}$', y_axis_label='PF',
            #                                  figure_title=title, is_y_2d=False)

            save_folder = 'figs/PF_to_L_nu/mc2=%d/a=%0.2f/' % (M_rate_c2_Led, a_portion)
            save_file_name = 'i=%d betta_mu=%d PF_to_L_nu.png' % (obs_i_angle_deg, betta_mu_deg)

            fig = plt.figure(figsize=(12, 6))
            ax = fig.add_subplot(111)

            ax.plot(x, y, color='black')

            x_axis_label = r'$L_{\nu}$'
            y_axis_label = 'PF'
            figure_title = title

            ax.set_xlabel(x_axis_label, fontsize=24)
            ax.set_ylabel(y_axis_label, fontsize=24)
            fig.suptitle(figure_title, fontsize=14)
            #
            # plt.show()

            main_service.save_figure(fig, save_folder, save_file_name)

            # plt.plot(x, y, color='black')
            #
            # plt.xlabel(r'$L_{\nu}$')
            # plt.ylabel('PF')

            # plt.title(
            #     'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led))
            #
            # main_service.create_file_path(save_folder)
            # plt.savefig(save_folder + save_file_name)
            #
            # plt.cla()

            # plt.show()

# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111)
#
# im = ax.pcolormesh(betta_mu, i_angle, mean_arr)
#
# x_axis_label = 'betta_mu'
# y_axis_label = 'i_angle'
#
# ax.set_xlabel(x_axis_label, fontsize=24)
# ax.set_ylabel(y_axis_label, fontsize=24)
#
# fig_title = 'mean'
# fig.suptitle(fig_title, fontsize=14)
#
# plt.colorbar(im)
# plt.show()
#
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111)
#
# im = ax.pcolormesh(betta_mu, i_angle, dispersion_arr)
#
# x_axis_label = 'betta_mu'
# y_axis_label = 'i_angle'
#
# ax.set_xlabel(x_axis_label, fontsize=24)
# ax.set_ylabel(y_axis_label, fontsize=24)
#
# fig_title = 'dispersion'
# fig.suptitle(fig_title, fontsize=14)
#
# plt.colorbar(im)
# plt.show()
#

''' Дисперсия среднне и max - min colormesh '''

title_dict = {0: 'mean', 1: 'dispersion', 2: 'max - min'}
data_dict = {0: mean_arr, 1: dispersion_arr, 2: max_min}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for i in range(len(axes)):
    col_bar = axes[i].pcolormesh(betta_mu, i_angle, data_dict[i])

    axes[i].set_xlabel('betta_mu')
    axes[i].set_ylabel('i_angle')
    axes[i].set_title(title_dict[i])

    plt.colorbar(col_bar, ax=axes[i])

fig_title = 'a = %0.2f, mc = %d' % (a_portion, M_rate_c2_Led)

# plt.subplots_adjust(hspace=0.5)

fig.suptitle(fig_title, fontsize=14)
fig.tight_layout()
# plt.show()

save_folder = 'figs/colormesh/mc2=%d/a=%0.2f/' % (M_rate_c2_Led, a_portion)
save_file_name = 'colormesh_fig.png'
main_service.save_figure(fig, save_folder, save_file_name)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for i in range(2, len(axes)):
    col_bar = axes[i].pcolormesh(betta_mu, i_angle, data_dict[i])

    axes[i].set_xlabel('betta_mu')
    axes[i].set_ylabel('i_angle')
    axes[i].set_title(title_dict[i])

    plt.colorbar(col_bar, ax=axes[i])

fig_title = 'a = %0.2f, mc = %d' % (a_portion, M_rate_c2_Led)
fig.suptitle(fig_title, fontsize=14)
fig.tight_layout()

max_value_mean = 0.7
min_value_mean = 0.1

col_bar = axes[0].pcolormesh(betta_mu, i_angle, data_dict[0], vmax=max_value_mean, vmin=min_value_mean)
axes[0].set_xlabel('betta_mu')
axes[0].set_ylabel('i_angle')
axes[0].set_title(title_dict[0])
plt.colorbar(col_bar, ax=axes[0])

max_value_disp = 0.2
min_value_disp = 0.02

col_bar = axes[1].pcolormesh(betta_mu, i_angle, data_dict[1], vmax=max_value_disp, vmin=min_value_disp)
axes[1].set_xlabel('betta_mu')
axes[1].set_ylabel('i_angle')
axes[1].set_title(title_dict[1])
plt.colorbar(col_bar, ax=axes[1])

# plt.show()

save_file_name = 'colormesh_fig_fixed.png'
main_service.save_figure(fig, save_folder, save_file_name)
# ---------------------------------------------------------------------------------------------------

flag_next = False
if flag_next:
    # меняем долю а
    M_rate_c2_Led = mc2[2]
    obs_i_angle_deg = 30
    betta_mu_deg = 60

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    buff_marker = 0
    for a_index in range(len(a_portion_arr)):
        final_array = []
        for i in range(len(fi_0)):
            a_portion = a_portion_arr[a_index]
            phi_accretion_begin_deg = fi_0[i]

            full_file_folder = get_folder()
            working_folder = full_file_folder + folder

            PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
            final_array.append(PF_array[energy_index])

            label = 'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led)
            buff_marker += 1

        ax.plot(fi_0, final_array, label=label)

    plt.legend(ncol=3, fontsize=10, framealpha=0.2)
    plt.show()

    # меняем mc
    a_portion = a_portion_arr[2]
    obs_i_angle_deg = 30
    betta_mu_deg = 60

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    buff_marker = 0
    for mc_index in range(len(mc2)):
        final_array = []
        for i in range(len(fi_0)):
            M_rate_c2_Led = mc2[mc_index]
            phi_accretion_begin_deg = fi_0[i]

            full_file_folder = get_folder()
            working_folder = full_file_folder + folder

            PF_array = main_service.load_arr_from_txt(full_file_folder + folder, file_name)
            final_array.append(PF_array[energy_index])

            label = 'i=%d betta_mu=%d a=%0.2f m=%d' % (obs_i_angle_deg, betta_mu_deg, a_portion, M_rate_c2_Led)
            buff_marker += 1

        ax.plot(fi_0, final_array, label=label)

    plt.legend(ncol=3, fontsize=10, framealpha=0.2)
    plt.show()
