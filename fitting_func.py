import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import config
import main_service
import accretionColumnService

plt.style.use(['science', 'notebook', 'grid'])

i_angle = 20
betta_mu = 20
mc2 = 30
a_portion = 0.66

fi_0 = 0

config.set_e_obs(i_angle, 0)
config.set_betta_mu(betta_mu)

config.M_rate_c2_Led = mc2
config.a_portion = a_portion
config.phi_accretion_begin_deg = config.fi_0_dict[a_portion] + fi_0

config.update()

folder = 'L_nu/'
working_folder = config.full_file_folder + folder

file_name = "L_nu.txt"
data_array = main_service.load_arr_from_txt(working_folder, file_name)

N_energy = config.N_energy

L_nu_avg_on_phase = [0] * N_energy
for i in range(N_energy):
    L_nu_avg_on_phase[i] = np.mean(data_array[i])

file_name = 'energy.txt'
energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

fig_title = r'$L_{\nu}$'

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)

# ax.plot(energy_arr, L_nu_avg_on_phase, label='origin')
# plt.xscale('log')
# plt.yscale('log')

xdata = energy_arr
normalization_coeff = max(L_nu_avg_on_phase)

ydata = L_nu_avg_on_phase / energy_arr

y_to_plot = list(ydata).copy()

normalization_coeff = max(ydata)
ydata = ydata / normalization_coeff

nu_arr = accretionColumnService.get_frequency_from_energy(energy_arr)
nu_L_nu_avg_on_phase = L_nu_avg_on_phase * nu_arr

ydata = nu_L_nu_avg_on_phase / energy_arr

y_to_plot = list(ydata).copy()

normalization_coeff = max(ydata)
ydata = ydata / normalization_coeff

# normalization_coeff = max(nu_L_nu_avg_on_phase)
# ydata = nu_L_nu_avg_on_phase / normalization_coeff

ax.plot(energy_arr, y_to_plot, label='origin')
plt.xscale('log')
plt.yscale('log')


# cutoffpl
def cutoff(x, K, alpha, beta):
    y = K * x ** (-alpha) * np.exp(-x / beta)
    return y


def newfit(x, k, gamma, Ec, Ef):
    y = k * x ** (-gamma) * np.exp(-(x - Ec) / Ef)
    return y


parameters, covariance = curve_fit(cutoff, xdata, ydata)
fit_y = cutoff(xdata, *parameters)

print('cutoff')
names = ['K', 'alpha', 'beta']
dict_vals = dict(zip(names, parameters))
for key, val in dict_vals.items():
    print(f'{key} : {val:.2f}')

plt.plot(xdata, fit_y * normalization_coeff, label='cutoff')


parameters, covariance = curve_fit(newfit, xdata, ydata)
fit_y = newfit(xdata, *parameters)

print('new')
names = ['k', 'gamma', 'Ec', 'Ef']
dict_vals = dict(zip(names, parameters))
for key, val in dict_vals.items():
    print(f'{key} : {val:.2f}')

plt.plot(xdata, fit_y * normalization_coeff, label='new_func')

# fit_K = parameters[0]
# fit_alpha = parameters[1]
# fit_beta = parameters[2]

plt.xlabel(r'$h \nu$ [keV]', fontsize=24)
plt.ylabel(r'$\frac{L_{\nu}}{h \nu} $ [$erg/keV$]', fontsize=24)

plt.ylim(min(y_to_plot))

ax.legend()
plt.show()
