import numpy as np
import config
import matplotlib.pyplot as plt

arr_len = 6
arr = [0] * arr_len
energy_bot = [0] * arr_len
energy_top = [0] * arr_len
append_index = config.t_max_for_plot - config.t_max

energy_bot[0] = 1
energy_top[0] = 4
energy_bot[1] = 4
energy_top[1] = 8

for i in range(2, arr_len):
    energy_bot[i] = energy_bot[i - 1] + 4
    energy_top[i] = energy_top[i - 1] + 4

# energy_bot[0], energy_bot[1], energy_bot[2], energy_bot[3] = 1, 4, 8, 12
# energy_top[0], energy_top[1], energy_top[2], energy_top[3] = 4, 8, 12, 16
file_folder = 'figs/'
args_folder = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)
file_folder = file_folder + args_folder

folder = 'luminosity_in_range/'
for i in range(arr_len):
    file_name = file_folder + folder + 'txt/' + "sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (
        energy_bot[i], energy_top[i])
    arr[i] = np.loadtxt(file_name)

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
for i in range(arr_len):
    ax.plot(phi_for_plot, arr[i],
            label="%0.2f - %0.2f KeV" % (energy_bot[i], energy_top[i]))

# plt.yscale('log')
ax.legend()
plt.show()
file_name = file_folder + folder + 'sum_of_luminosity_in_range.png'
fig.savefig(file_name, dpi=fig.dpi)
