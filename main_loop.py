import numpy as np
import multiprocessing as mp
from itertools import repeat
import time

import geometricTask.matrix as matrix
import config
import accretionColumnService
from accretionColumn import AccretionColumn
import vectors
import main_service

import plot_funcs.plot_from_main
import plot_funcs.plot_luminosity_in_range
import plot_funcs.plot_L_nu
import plot_funcs.plot_nu_L_nu

if __name__ == '__main__':
    def plot_all():
        plot_funcs.plot_from_main.plot_figs()
        plot_funcs.plot_luminosity_in_range.plot_figs()
        plot_funcs.plot_L_nu.plot_figs()
        plot_funcs.plot_nu_L_nu.plot_figs()

    mc2 = [10, 20, 40]
    a_portion = [0.25, 0.65, 1]
    fi_0 = [10, 20, 70, 140]

    for i in range(len(mc2)):
        for j in range(len(a_portion)):
            for k in range(len(fi_0)):

                config.M_rate_c2_Led = mc2[i]
                config.a_portion = a_portion[j]
                config.phi_accretion_begin_deg = fi_0[k]

                config.update()

                print('calculate for mc=%f a=%f di_0=%f' % (mc2[i], a_portion[j], fi_0[k]))
                R_alfven = (config.mu ** 2 / (
                        2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
                R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
                # print('R_e = %f' % (R_e / config.R_ns))
                R_e_outer_surface, R_e_inner_surface = R_e, R_e  # допущение что толщина = 0
                # вектор на наблюдателя в системе координат двойной системы (условимся что omega и e_obs лежат в пл-ти x0z)
                e_obs = np.array([np.sin(config.obs_i_angle), 0, np.cos(config.obs_i_angle)])
                file_name_variables = "betta_omega=%d betta_mu=%d a_portion=%f M_rate_c2_Led=%d" \
                                      % (config.betta_rotate, config.betta_mu, config.a_portion, config.M_rate_c2_Led)
                approx_method = accretionColumnService.approx_method

                # ------------------- создание папки для графиков --------------------------
                full_file_folder = config.full_file_folder
                main_service.create_file_path(full_file_folder)
                # ------------------- создание папки для графиков --------------------------

                # ----------------- начало инициализации верхней колонки ------------------------
                # от поверхности NS - угол при котором радиус = радиусу НЗ
                theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(
                    R_e_outer_surface)
                theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(
                    R_e_inner_surface)

                top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                                             theta_accretion_begin_inner_surface, True)
                # ----------------- конец инициализации верхней колонки ------------------------

                # ----------------- начало инициализации нижней колонки ------------------------
                theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
                theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

                bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                                             theta_accretion_begin_inner_surface, False)
                # ----------------- конец инициализации нижней колонки ------------------------

                # словарь для того чтобы можно было в цикле ходить по поверхностям
                surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface, 2: bot_column.outer_surface,
                            3: bot_column.inner_surface}

                # ----------------------------- график T_eff --------------------------------
                arr_T_eff = []
                for key, surface in surfaces.items():
                    arr_T_eff.append(surface.T_eff)
                file_name = 'surfaces_T_eff.txt'
                main_service.save_arr_as_txt(arr_T_eff, full_file_folder, file_name)
                # ----------------------------- график T_eff --------------------------------

                file_name = "save_phi_range.txt"
                main_service.save_arr_as_txt(top_column.outer_surface.phi_range, full_file_folder, file_name)

                file_name = "save_theta_range.txt"
                main_service.save_arr_as_txt(top_column.outer_surface.theta_range, full_file_folder, file_name)

                print('phi_theta_range saved')
                # ----------------- углы для нахождения пересечений -------------------------
                theta_accretion_begin = top_column.outer_surface.theta_range[0]
                theta_accretion_end = top_column.outer_surface.theta_range[-1]
                # ----------------- углы для нахождения пересечений -------------------------

                file_name = 'save_values.txt'
                f = open(full_file_folder + file_name, 'w')
                f.write('R_e = %f \nksi_shock = %f\n' % (R_e / config.R_ns, top_column.outer_surface.ksi_shock))
                f.close()

                file_name = 'save_values.txt'
                f = open(full_file_folder + file_name, 'a')

                f.write('beta = %f\n' % top_column.outer_surface.beta)

                power_index = 0
                number = top_column.inner_surface.L_x
                while number > 10:
                    number = number / 10
                    power_index += 1
                print('total L_x = %f * 10**%d' % (number, power_index))

                f.write('total L_x = %f * 10**%d \n' % (number, power_index))

                power_index = 0
                number = top_column.inner_surface.calculate_total_luminosity()

                f.write(
                    'difference L_x / L_calc - 1 : %f ' % (
                            (top_column.inner_surface.L_x / (4 * number) - 1) * 100) + '% \n')

                while number > 10:
                    number = number / 10
                    power_index += 1

                f.write('calculated total L_x of single surface = %f * 10**%d \n' % (number, power_index))
                f.close()

                time_start = time.time()

                # ------------------ начало заполнения матриц косинусов ---------------------------
                # for key, surface in surfaces.items():
                #     surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end, top_column.outer_surface.phi_range,
                #                                bot_column.outer_surface.phi_range, e_obs)

                # распараллелил
                for key, surface in surfaces.items():
                    with mp.Pool(mp.cpu_count()) as pool:
                        result_cos_psi_range = pool.starmap(surface.async_fill_cos_psi_range,
                                                            zip(range(config.t_max), repeat(theta_accretion_begin),
                                                                repeat(theta_accretion_end),
                                                                repeat(top_column.outer_surface.phi_range),
                                                                repeat(bot_column.outer_surface.phi_range),
                                                                repeat(e_obs)))

                    cos_psi_range_final = []
                    for cos_psi in result_cos_psi_range:
                        cos_psi_range_final.append(cos_psi)
                    surface.cos_psi_range = cos_psi_range_final
                # ------------------ конец заполнения матриц косинусов ---------------------------

                time_cos = time.time()

                # ------------------ начало заполнения массивов светимости -----------------------
                arr_simps_integrate = [0] * 4
                sum_simps_integrate = 0
                for key, surface in surfaces.items():
                    arr_simps_integrate[key] = surface.calculate_integral_distribution()
                    # file_name = "%s %s %d.txt" % (file_name_variables, approx_method, key)
                    # np.savetxt(file_name, arr_simps_integrate[key])
                    sum_simps_integrate += np.array(arr_simps_integrate[key])
                # ------------------ конец заполнения массивов светимости -----------------------

                file_name = 'save_values.txt'
                f = open(full_file_folder + file_name, 'a')

                avg_L_on_phase = np.mean(sum_simps_integrate)

                power_index = 0
                number = avg_L_on_phase
                while number > 10:
                    number = number / 10
                    power_index += 1

                print('avg_L_on_phase = %f * 10**%d' % (number, power_index))

                f.write('avg_L_on_phase = %f * 10**%d\n' % (number, power_index))

                f.write('avg_L_on_phase / L_x = %.3f\n' % (avg_L_on_phase / top_column.inner_surface.L_x))

                f.close()

                d0 = accretionColumnService.get_delta_distance(top_column.inner_surface.theta_range[0],
                                                               top_column.inner_surface.R_e)
                l0 = accretionColumnService.get_A_normal(top_column.inner_surface.theta_range[0],
                                                         top_column.inner_surface.R_e) / d0

                file_name = 'save_values.txt'
                f = open(full_file_folder + file_name, 'a')

                f.write('width / length = %f \n' % (d0 / l0))
                f.write('A_perp / 2 d0**2 = %f \n' % (l0 / (2 * d0)))
                f.close()

                time_integral_distribution = time.time()

                print('ksi_shock = %f' % bot_column.outer_surface.ksi_shock)

                # --------------------- вывод графика светимости ----------------------------
                # 1 слитый массив - отдельные поверхности + сумма
                arr_total_luminosity_of_surfaces = np.append(arr_simps_integrate, [sum_simps_integrate], 0)
                print('total_luminosity_of_surfaces.txt saved')
                file_name = 'total_luminosity_of_surfaces.txt'
                main_service.save_arr_as_txt(arr_total_luminosity_of_surfaces, full_file_folder, file_name)
                # --------------------- вывод графика светимости ----------------------------

                # --------------------- вывод графика углов наблюдателя ----------------------------
                observer_theta = [0] * config.t_max
                observer_phi = [0] * config.t_max

                for t in range(config.t_max):
                    phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
                    # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
                    A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                                 config.betta_mu)
                    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
                    observer_phi[t], observer_theta[t] = vectors.get_angles_from_vector(e_obs_mu)

                obs_angles = []
                obs_angles.append(observer_phi)
                obs_angles.append(observer_theta)
                file_name = 'observer_angles.txt'
                main_service.save_arr_as_txt(obs_angles, full_file_folder, file_name)
                # --------------------- вывод графика углов наблюдателя ----------------------------

                # ------------------ цикл для диапазона энергий ----------------------
                # лог сетка по энергии: en+1 = en * const
                energy_step = (config.energy_max / config.energy_min) ** (1 / (config.N_energy - 1))
                energy_arr = list(config.energy_min * energy_step ** i for i in range(config.N_energy - 1))
                # чтобы убрать погрешности и закрыть массив точным числом
                energy_arr.append(config.energy_max)

                PF = [0] * (
                        config.N_energy - 1)  # -1 т.к. берем в диапазоне энергий, следовательно, на 1 меньше чем точек
                data_array = [0] * (config.N_energy - 1)
                folder = 'luminosity_in_range/'
                for energy_index in range(config.N_energy - 1):
                    current_energy_min = energy_arr[energy_index]
                    current_energy_max = energy_arr[energy_index + 1]
                    # ------------------ начало заполнения массивов светимости -----------------------
                    arr_simps_integrate = [0] * 4
                    sum_simps_integrate = 0
                    for key, surface in surfaces.items():
                        arr_simps_integrate[key] = surface.calculate_integral_distribution_in_range(current_energy_min,
                                                                                                    current_energy_max)
                        sum_simps_integrate += np.array(arr_simps_integrate[key])
                    # ------------------ конец заполнения массивов светимости -----------------------
                    PF[energy_index] = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

                    file_name = "sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (
                        current_energy_min, current_energy_max)
                    main_service.save_arr_as_txt(sum_simps_integrate, full_file_folder + folder + 'txt/', file_name)

                    data_array[energy_index] = sum_simps_integrate

                # ------------------ цикл для диапазона энергий ----------------------

                file_name = "PF.txt"
                main_service.save_arr_as_txt(PF, full_file_folder + folder, file_name)

                file_name = "luminosity_in_range.txt"
                main_service.save_arr_as_txt(data_array, full_file_folder + folder, file_name)

                file_name = "energy.txt"
                main_service.save_arr_as_txt(energy_arr, full_file_folder, file_name)

                time_integral_distribution_in_range = time.time()

                # print('спектр')
                print('L_nu')
                PF = [0] * config.N_energy
                data_array = [0] * config.N_energy
                folder = 'L_nu/'
                # ------------------ цикл для Spectrum Distribution ------------------
                for energy_index in range(config.N_energy):
                    current_energy = energy_arr[energy_index]
                    # ------------------ начало заполнения массивов Spectral Energy -----------------------
                    arr_simps_integrate = [0] * 4
                    sum_simps_integrate = 0
                    for key, surface in surfaces.items():
                        arr_simps_integrate[key] = surface.calculate_L_nu_on_energy(current_energy)
                        sum_simps_integrate += np.array(arr_simps_integrate[key])
                    # ------------------ конец заполнения массивов Spectral Energy -----------------------

                    PF[energy_index] = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

                    file_name = "L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % current_energy
                    main_service.save_arr_as_txt(sum_simps_integrate, full_file_folder + folder + 'txt/', file_name)

                    data_array[energy_index] = sum_simps_integrate
                # ------------------ цикл для Spectrum Distribution ------------------

                file_name = "PF.txt"
                main_service.save_arr_as_txt(PF, full_file_folder + folder, file_name)

                file_name = "L_nu.txt"
                main_service.save_arr_as_txt(data_array, full_file_folder + folder, file_name)

                time_calculate_L_nu_on_energy = time.time()

                # print('Spectral Energy Distribution')
                print('nu_L_nu')
                PF = [0] * config.N_energy
                data_array = [0] * config.N_energy
                folder = 'nu_L_nu/'
                # ------------------ цикл для Spectral Energy Distribution ------------------
                for energy_index in range(config.N_energy):
                    current_energy = energy_arr[energy_index]
                    # ------------------ начало заполнения массивов Spectral Energy -----------------------
                    arr_simps_integrate = [0] * 4
                    sum_simps_integrate = 0
                    for key, surface in surfaces.items():
                        arr_simps_integrate[key] = surface.calculate_nu_L_nu_on_energy(current_energy)
                        sum_simps_integrate += np.array(arr_simps_integrate[key])
                    # ------------------ конец заполнения массивов Spectral Energy -----------------------

                    PF[energy_index] = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

                    file_name = "nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % current_energy
                    main_service.save_arr_as_txt(sum_simps_integrate, full_file_folder + folder + 'txt/', file_name)

                    data_array[energy_index] = sum_simps_integrate
                # ------------------ цикл для Spectral Energy Distribution ------------------

                file_name = "PF.txt"
                main_service.save_arr_as_txt(PF, full_file_folder + folder, file_name)

                file_name = "nu_L_nu.txt"
                main_service.save_arr_as_txt(data_array, full_file_folder + folder, file_name)

                time_calculate_nu_L_nu_on_energy = time.time()

                print("execution time of program: %f s" % (
                        time_calculate_nu_L_nu_on_energy - time_start))

                # print("execution time of intersections: %f s" % (time_cos - time_start))
                # print(
                #     "execution time of calculate_integral_distribution: %f s" % (time_integral_distribution - time_cos))
                # print("execution time of calculate_integral_distribution_in_range: %f s" % (
                #         time_integral_distribution_in_range - time_integral_distribution))
                # print("execution time of calculate_L_nu_on_energy: %f s" % (
                #         time_calculate_L_nu_on_energy - time_integral_distribution_in_range))
                # print("execution time of calculate_nu_L_nu_on_energy: %f s" % (
                #         time_calculate_nu_L_nu_on_energy - time_calculate_L_nu_on_energy))

                plot_all()

                if (a_portion[j] == 1):
                    break
