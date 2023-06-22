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

phi_accretion_begin_deg = [0, 20, 40, 70, 110, 140]
N_deg = len(phi_accretion_begin_deg)

M_rate_c2_Led = 10
a_portion = 0.65

fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot(111)

for i in range(N_deg):
    file_folder = config.file_folder
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg[i])
    full_file_folder = file_folder + file_folder_args

    working_folder = full_file_folder + folder

    file_name = "PF.txt"
    PF = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = "energy.txt"
    energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    label = r'$\varphi_0$' + ('=%d' % phi_accretion_begin_deg[i]) + r'$^{\circ}$'
    ax.plot(energy_arr, PF, label=label)

x_axis_label = r'$h \nu$' + ' [keV]'
y_axis_label = 'PF'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)
#ax.legend(loc='best', bbox_to_anchor=(0.27, 0.35), fontsize=14, ncol=3)
ax.legend(fontsize=14, ncol=3)
# plt.show()

save_folder = config.file_folder + 'PF/'
file_name = 'mc2=%d a=%0.2f.png' % (M_rate_c2_Led, a_portion)
main_service.save_figure(fig, save_folder, file_name)

# --------------------------------------------------------------------------------
M_rate_c2_Led = [10]
a_portion = [0.65]

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)

for j in range(len(M_rate_c2_Led)):
    for k in range(len(a_portion)):
        for i in range(N_deg):
            file_folder = config.file_folder
            file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led[j], a_portion[k], phi_accretion_begin_deg[i])
            full_file_folder = file_folder + file_folder_args

            working_folder = full_file_folder + folder

            file_name = "PF.txt"
            PF = main_service.load_arr_from_txt(working_folder, file_name)

            file_name = "energy.txt"
            energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

            ax.plot(energy_arr, PF,
                    label='M=%d, a=%0.2f, angle=%d' % (M_rate_c2_Led[j], a_portion[k], phi_accretion_begin_deg[i]))

x_axis_label = r'$h \nu$' + ' [keV]'
y_axis_label = 'PF'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)
ax.legend()
# plt.show()

file_name = 'all.png'
main_service.save_figure(fig, save_folder, file_name)

# ---------------------------------------

M_rate_c2_Led = [10, 20, 30, 40, 50, 70]
a_portion = 0.65
phi_accretion_begin_deg = 0


fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot(111)

for i in range(len(M_rate_c2_Led)):
    file_folder = config.file_folder
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led[i], a_portion, phi_accretion_begin_deg)
    full_file_folder = file_folder + file_folder_args

    working_folder = full_file_folder + folder

    file_name = "PF.txt"
    PF = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = "energy.txt"
    energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    label = r'$\dot{m}$' + ('=%d' % M_rate_c2_Led[i])
    ax.plot(energy_arr, PF, label=label)

x_axis_label = r'$h \nu$' + ' [KeV]'
y_axis_label = 'PF'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)
# ax.legend(loc='best', bbox_to_anchor=(0.3, 0.35), fontsize=14, ncol=3)
ax.legend(fontsize=14, ncol=3)
# plt.show()

save_folder = config.file_folder + 'PF/'
file_name = 'PF_on_energy'
main_service.save_figure(fig, save_folder, file_name)
