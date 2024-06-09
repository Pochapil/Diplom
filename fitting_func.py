import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import config
import main_service

i_angle = 20
betta_mu = 40
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

ax.plot(energy_arr, L_nu_avg_on_phase, label='origin')
plt.xscale('log')
plt.yscale('log')

xdata = energy_arr
normalization_coeff = max(L_nu_avg_on_phase)
ydata = L_nu_avg_on_phase / normalization_coeff

# cutoffpl
def cutoff(x, K, alpha, beta):
    y = K * x ** (-alpha) * np.exp(-x / beta)
    return y


parameters, covariance = curve_fit(cutoff, xdata, ydata)

fit_K = parameters[0]
fit_alpha = parameters[1]
fit_beta = parameters[2]

print(fit_K, fit_alpha, fit_beta)

fit_y = cutoff(xdata, *parameters)
plt.plot(xdata, fit_y * normalization_coeff, label='fit')

plt.ylim(min(L_nu_avg_on_phase))

ax.legend()
plt.show()


def func(x):
    x = 2
    y = 10


x = 10
func(x)
print(x)
