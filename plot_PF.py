import numpy as np
import matplotlib.pyplot as plt

import config
import main_service

plt.style.use(['science', 'notebook', 'grid'])
# plt.style.use(['default'])

folder = 'nu_L_nu/'
working_folder = config.full_file_folder + folder

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

file_name = "energy.txt"
energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

N_deg = 3
phi_accretion_begin_deg = [0, 70, 140]

M_rate_c2_Led = 10
a_portion = 0.25

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)

for i in range(N_deg):
    file_folder = 'figs/new_energy/'
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg[i])
    full_file_folder = file_folder + file_folder_args

    working_folder = full_file_folder + folder

    file_name = "PF.txt"
    PF = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = "energy.txt"
    energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    ax.plot(energy_arr, PF, label='angle = %d' % phi_accretion_begin_deg[i])

x_axis_label = r'$h \nu$' + ' [KeV]'
y_axis_label = 'PF'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)
ax.legend()
# plt.show()

save_folder = 'figs/new_energy/PF/'
file_name = 'mc2=%d a=%0.2f.png' % (M_rate_c2_Led, a_portion)
main_service.save_figure(fig, save_folder, file_name)

M_rate_c2_Led = [10, 40]
a_portion = [0.25, 0.65]

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)

for j in range(2):
    for k in range(2):
        for i in range(N_deg):
            file_folder = 'figs/new_energy/'
            file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led[j], a_portion[k], phi_accretion_begin_deg[i])
            full_file_folder = file_folder + file_folder_args

            working_folder = full_file_folder + folder

            file_name = "PF.txt"
            PF = main_service.load_arr_from_txt(working_folder, file_name)

            file_name = "energy.txt"
            energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

            ax.plot(energy_arr, PF,
                    label='M=%d, a=%0.2f, angle=%d' % (M_rate_c2_Led[j], a_portion[k], phi_accretion_begin_deg[i]))

x_axis_label = r'$h \nu$' + ' [KeV]'
y_axis_label = 'PF'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)
ax.legend()
# plt.show()

file_name = 'all.png'
main_service.save_figure(fig, save_folder, file_name)