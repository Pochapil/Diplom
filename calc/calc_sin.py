import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.fft import fft

import config
import main_service


def make_new_phi(a_portion, fi_0_arr):
    # new_phi_0
    for k in range(len(fi_0_arr)):
        fi_0_arr[k] = (config.fi_0_dict[a_portion] + 20 * (k)) % 360


i_angle = 40
betta_mu = 60
mc2 = 100
a_portion = 0.22
fi_0 = 40

# fi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

file_name = 'total_luminosity_of_surfaces.txt'

# make_new_phi(a_portion, fi_0_arr)
# for fi_0_index in range(len(fi_0_arr)):
full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion, fi_0)

L_array = main_service.load_arr_from_txt(full_file_folder, file_name)[4]
L_array += \
    main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/', 'scattered_energy_top.txt')
L_array += \
    main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/', 'scattered_energy_bot.txt')

y_data = L_array
y_max = np.max(y_data)
y_data /= y_max
x_data = np.linspace(0, 2 * np.pi, 45)


def fit_func(x, a, b):
    return a * np.sin(x - b)


popt, pcov = curve_fit(fit_func, x_data, y_data)
fit_y_1 = fit_func(x_data, *popt)
fit_params_1 = popt

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.plot(x_data, y_data)
ax.plot(x_data, fit_y_1)
plt.show()


def fit_func2(x, a, b):
    return a * np.sin(2 * x - b)


popt, pcov = curve_fit(fit_func2, x_data, y_data - fit_y_1)
fit_y_2 = fit_func2(x_data, *popt)
fit_params_2 = popt

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.plot(x_data, y_data)
ax.plot(x_data, fit_y_1 + fit_y_2)
plt.show()


def fit_func3(x, a, b):
    return a * np.sin(3 * x - b)


popt, pcov = curve_fit(fit_func3, x_data, y_data - fit_y_1 - fit_y_2)
fit_y_3 = fit_func3(x_data, *popt)
fit_params_3 = popt

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.plot(x_data, y_data)
ax.plot(x_data, fit_y_1 + fit_y_2 + fit_y_3)
plt.show()

def fit_func4(x, c):
    return np.sin(0 * x) + c

popt, pcov = curve_fit(fit_func4, x_data, y_data)
fit_y_4 = fit_func4(x_data, *popt)
fit_params_4 = popt

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.plot(x_data, y_data)
ax.plot(x_data, fit_y_1 + fit_y_2 + fit_y_3 + fit_y_4)
plt.show()

print(fit_params_1)
print(fit_params_2)
print(fit_params_3)
print(fit_params_4)


n = len(x_data)
dx = np.diff(x_data)[0]
f_hat = np.fft.fft(y_data, n)
# print(f_hat)

PSD = f_hat * np.conj(f_hat)/n
freq = (1/(dx * n)) * np.arange(n)
L = np.arange(1, np.floor(n/2), dtype='int')

ind = np.argpartition(np.abs(f_hat.real), -4)[-4:]

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.plot(freq[L], PSD[L])
plt.show()

indecies = PSD > 0.1
indecies[-1] = False
indecies[-2] = False
indecies[0] = True
indecies[3:6] = True

PSDclean = PSD * indecies
f_hat = indecies * f_hat
print(indecies)
ffilt = np.fft.ifft(f_hat)

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
ax.plot(x_data, y_data)
ax.plot(x_data, ffilt, label='furie')
plt.legend()
plt.show()

print(PSD[indecies])

print(f_hat[ind])
print(ind)

