import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time

import config
import main_service

plt.style.use(['science', 'notebook', 'grid'])


# def make_new_phi(a_portion, fi_0_arr):
#     # new_phi_0
#     for k in range(len(fi_0_arr)):
#         fi_0_arr[k] = (config.fi_0_dict[a_portion] + 20 * (k)) % 360

def make_new_phi(a_portion, fi_0):
    # new_phi_0
    fi_0 = (config.fi_0_dict[a_portion] + fi_0) % 360
    return fi_0


check_flag = False

i_angle_arr = [10 * i for i in range(1, 10)]
betta_mu_arr = [10 * i for i in range(1, 10)]

# i_angle_arr = [20, 40, 60]
# betta_mu_arr = [20, 40, 60]
# a_portion_arr = [0.11]
a_portion_arr = [0.22, 0.44, 0.55, 0.66]
mc2_arr = [30, 100]
fi_0_arr = [0]

for i_angle in i_angle_arr:
    for betta_mu in betta_mu_arr:
        for mc2 in mc2_arr:
            for a_portion in a_portion_arr:
                for fi_0 in fi_0_arr:
                    fi_0 = make_new_phi(a_portion, fi_0)

                    print(
                        f'calculate for i_angle={i_angle} betta_mu={betta_mu} mc2={mc2} a_portion={a_portion} fi_0={fi_0}')

                    full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion, fi_0)

                    bot_scatter = main_service.load_arr_from_txt(
                        full_file_folder + 'scattered_on_magnet_lines/' + 'nu_L_nu/',
                        'bot_column_scatter_nu_L_nu.txt')
                    top_scatter = main_service.load_arr_from_txt(
                        full_file_folder + 'scattered_on_magnet_lines/' + 'nu_L_nu/',
                        'top_column_scatter_nu_L_nu.txt')

                    A1 = [0] * len(config.energy_arr)
                    A2 = [0] * len(config.energy_arr)
                    A3 = [0] * len(config.energy_arr)

                    x_data = np.linspace(0, 2 * np.pi, 45)

                    dict_vals = {'a1': [], 'b1': [], 'a2': [], 'b2': [], 'a3': [], 'b3': [], 'c': []}

                    for energy_index, current_energy in enumerate(config.energy_arr):
                        file_name = "nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % current_energy

                        nu_L_nu_array = main_service.load_arr_from_txt(full_file_folder + 'nu_L_nu/' + 'txt/',
                                                                       file_name)

                        if not (np.isnan(bot_scatter[energy_index]).any() or np.isinf(bot_scatter[energy_index]).any()):
                            nu_L_nu_array += bot_scatter[energy_index]
                        if not (np.isnan(top_scatter[energy_index]).any() or np.isinf(top_scatter[energy_index]).any()):
                            nu_L_nu_array += top_scatter[energy_index]

                        y_data = nu_L_nu_array
                        mean = np.mean(y_data)
                        std = np.std(y_data)

                        y_data = (y_data - mean) / std


                        def fit_func_all(x, a_1, b_1, a_2, b_2, a_3, b_3, c):
                            return a_1 * np.sin(x + b_1) + a_2 * np.sin(2 * x + b_2) + a_3 * np.sin(3 * x + b_3) + c


                        popt, pcov = curve_fit(fit_func_all, x_data, y_data)
                        fit_all = fit_func_all(x_data, *popt)
                        fit_params = popt

                        A1[energy_index] = fit_params[0]
                        A2[energy_index] = fit_params[2]
                        A3[energy_index] = fit_params[4]

                        for index, key in enumerate(dict_vals.keys()):
                            dict_vals[key].append(fit_params[index])

                        if check_flag:
                            fig = plt.figure(figsize=(12, 6))
                            ax = fig.add_subplot(111)
                            ax.plot(x_data, y_data, label='real')
                            ax.plot(x_data, fit_all, label='fitted')
                            plt.legend()

                            file_name = f'fit_check_{current_energy:.2f}' + '.png'
                            save_folder = full_file_folder + 'fourier/' + 'fit_check/'
                            main_service.save_figure(fig, save_folder, file_name)

                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_subplot(111)
                    ax.plot(config.energy_arr, np.abs(A1), label='A1')
                    ax.plot(config.energy_arr, np.abs(A2), label='A2')
                    ax.plot(config.energy_arr, np.abs(A3), label='A3')
                    plt.legend()

                    file_name = 'params_value' + '.png'
                    save_folder = full_file_folder + 'fourier/'
                    main_service.save_figure(fig, save_folder, file_name)

                    save_folder = config.PROJECT_DIR + 'figs/' + 'fourier/'
                    file_name = f'i_angle={i_angle} betta_mu={betta_mu} mc2={mc2} a_portion={a_portion} fi_0={fi_0}' + '.png'
                    main_service.save_figure(fig, save_folder, file_name)

                    file_name = 'save_params.txt'
                    with open(full_file_folder + 'fourier/' + file_name, 'w') as outfile:
                        for key in dict_vals.keys():
                            outfile.write(f'{key}: ' + ' '.join(list(map(str, dict_vals[key]))) + '\n')
