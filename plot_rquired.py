import numpy as np
import multiprocessing as mp
from itertools import repeat
import time
import matplotlib.pyplot as plt

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
fi_0 = [20 * i for i in range(19)]
# fi_0 = [0, 20, 40, 60, 100, 160, 200, 220, 260, 300]
i_angle = [10, 20, 30, 40, 60, 90]
betta_mu = [30, 60, 90]
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
a_portion = a_portion_arr[1]
M_rate_c2_Led = mc2[1]
buff_marker = 0

final_final_array = np.zeros((len(i_angle), len(betta_mu), len(fi_0)))
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

        final_final_array[i_angle_index][betta_mu_index] = final_array
        marker = marker_dict[buff_marker % 4]
        buff_marker += 1
        ax.plot(fi_0, final_array, label=label, marker=marker)

plt.legend(ncol=3, fontsize=10, framealpha=0.2)
plt.show()

# print(final_array)


dispersion_arr = np.zeros((len(i_angle), len(betta_mu)))
mean_arr = np.zeros((len(i_angle), len(betta_mu)))

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        dispersion_arr[i_angle_index][betta_mu_index] = np.var(final_final_array[i_angle_index][betta_mu_index]) ** (
                    1 / 2)
        mean_arr[i_angle_index][betta_mu_index] = np.mean(final_final_array[i_angle_index][betta_mu_index])

_x, _y = np.meshgrid(i_angle, betta_mu)
x, y = _x.ravel(), _y.ravel()

dispersion_arr_rev = dispersion_arr.ravel()
mean_arr_rev = mean_arr.ravel()

print(mean_arr.shape)
print(mean_arr)
# print(dispersion_arr[i_angle_index][betta_mu_index])


# fig = plt.figure(figsize=(12, 6))
# ax1 = fig.add_subplot(121, projection='3d')
# ax2 = fig.add_subplot(122, projection='3d')
#
# ax1.bar3d(x, y, np.zeros_like(mean_arr_rev), 5, 5, mean_arr_rev, shade=True)
# ax1.set_title('mean')
# ax1.set_xlabel('i_angle')
# ax1.set_ylabel('betta_mu')
#
# ax2.bar3d(x, y, np.zeros_like(dispersion_arr_rev), 5, 5, dispersion_arr_rev, shade=True)
# ax2.set_title('dispersion')
# ax2.set_xlabel('i_angle')
# ax2.set_ylabel('betta_mu')
#
# plt.show()

# fig, axes = plt.subplots(2, 1, figsize=(12, 6))
#
# ax1 = axes[0]
# ax2 = axes[1]
#
# ax1.title.set_text('First Plot')



fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)

im = ax.pcolormesh(betta_mu, i_angle, mean_arr)

x_axis_label = 'betta_mu'
y_axis_label = 'i_angle'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)

fig_title = 'mean'
fig.suptitle(fig_title, fontsize=14)

plt.colorbar(im)
plt.show()




fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)

im = ax.pcolormesh(betta_mu, i_angle, dispersion_arr)

x_axis_label = 'betta_mu'
y_axis_label = 'i_angle'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)

fig_title = 'dispersion'
fig.suptitle(fig_title, fontsize=14)

plt.colorbar(im)
plt.show()


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
