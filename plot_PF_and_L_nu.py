import config
import numpy as np
import matplotlib.pyplot as plt

file_folder = 'figs/'
args_folder = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)
file_folder = file_folder + args_folder

folder = 'luminosity_in_range/'
# folder = 'L_nu/'
# folder = 'nu_L_nu/'

file_name = "PF.txt"
full_file_name = file_folder + folder + file_name
PF = np.loadtxt(full_file_name)

N_energy = 10
energies = [0] * N_energy

for i in range(N_energy):
    if i == 0:
        energies[i] = 1
    else:
        energies[i] = i * 4

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(energies, PF)

file_name = "PF.png"
full_file_name = file_folder + folder + file_name
fig.savefig(full_file_name, dpi=fig.dpi)
plt.close()

folder = 'L_nu/'
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "txt/L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    full_file_name = file_folder + folder + file_name
    arr[i] = np.loadtxt(full_file_name)

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
for i in range(N_energy):
    ax.plot(phi_for_plot, arr[i], label='%0.2f KeV' % (energies[i]))
ax.legend()
fig_title = r'$L_{\nu}$'
fig.suptitle(fig_title, fontsize=14)

file_name = folder[:-1] + '.png'
full_file_name = file_folder + folder + file_name
fig.savefig(full_file_name, dpi=fig.dpi)
plt.close()

folder = 'nu_L_nu/'
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "txt/nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    full_file_name = file_folder + folder + file_name
    arr[i] = np.loadtxt(full_file_name)

phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
for i in range(N_energy):
    ax.plot(phi_for_plot, arr[i], label='%0.2f KeV' % (energies[i]))
ax.legend()
fig_title = r'$\nu L_{\nu}$'
fig.suptitle(fig_title, fontsize=14)

file_name = folder[:-1] + '.png'
full_file_name = file_folder + folder + file_name
fig.savefig(full_file_name, dpi=fig.dpi)
plt.close()

folder = 'L_nu/'
phase_index = 0
arr = [0] * N_energy
for i in range(N_energy):
    file_name = "txt/L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energies[i]
    full_file_name = file_folder + folder + file_name
    arr[i] = np.loadtxt(full_file_name)

L_nu = []
for i in range(N_energy):
    L_nu.append(arr[i][phase_index])
energy_i = 0
energy = [1] * N_energy
for i in range(1, N_energy):
    energy[i] = i * 4

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(energy, L_nu)
fig_title = r'$L_{\nu}$'
fig.suptitle(fig_title, fontsize=14)

file_name = 'L_nu(nu)' + '.png'
full_file_name = file_folder + folder + file_name
fig.savefig(full_file_name, dpi=fig.dpi)
plt.close()

L_nu_avg_on_phase = [0] * N_energy
for i in range(N_energy):
    sum = 0
    for j in range(config.t_max):
        sum += arr[i][j]
    avg = sum / config.t_max
    L_nu_avg_on_phase[i] = avg

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(energy, L_nu_avg_on_phase)
fig_title = r'$L_{\nu}$'
fig.suptitle(fig_title, fontsize=14)

file_name = 'L_nu(nu)_avg' + '.png'
full_file_name = file_folder + folder + file_name
fig.savefig(full_file_name, dpi=fig.dpi)
plt.close()