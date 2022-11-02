import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from itertools import repeat
import pathlib

import geometricTask.matrix as matrix
import config
import accretionColumnService
from accretionColumn import AccretionColumn
import vectors
import main_service

if __name__ == '__main__':

    R_alfven = (config.mu ** 2 / (2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
    R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
    print('R_e = %f' % (R_e / config.R_ns))
    R_e_outer_surface, R_e_inner_surface = R_e, R_e  # допущение что толщина = 0
    # вектор на наблюдателя в системе координат двойной системы
    e_obs = np.array([0, np.sin(config.i_angle), np.cos(config.i_angle)])
    file_name_variables = "betta_omega=%d betta_mu=%d a_portion=%f M_rate_c2_Led=%d" \
                          % (config.betta_rotate, config.betta_mu, config.a_portion, config.M_rate_c2_Led)
    approx_method = accretionColumnService.approx_method

    # ------------------- создание папки для графиков --------------------------
    file_folder = 'figs/'
    file_folder_args = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)
    full_file_folder = file_folder + file_folder_args
    full_file_folder = config.full_file_folder
    main_service.create_file_path(full_file_folder)
    # ------------------- создание папки для графиков --------------------------

    # ----------------- начало инициализации верхней колонки ------------------------
    # от поверхности NS - угол при котором радиус = радиусу НЗ
    theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(R_e_outer_surface)
    theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(R_e_inner_surface)

    top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                                 theta_accretion_begin_inner_surface, True)
    # ----------------- конец инициализации верхней колонки ------------------------

    # ----------------- начало инициализации нижней колонки ------------------------
    theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
    theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

    bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                                 theta_accretion_begin_inner_surface, False)
    # ----------------- конец инициализации нижней колонки ------------------------

    # ----------------------------- график T_eff --------------------------------
    fig = main_service.create_figure(top_column.outer_surface.theta_range, top_column.outer_surface.T_eff,
                                     x_axis_label='phase',
                                     y_axis_label=r'$T_{eff}, K$', figure_title=r'$T_{eff}$', is_y_2d=False)
    file_name = 'T_eff.png'
    main_service.save_figure(fig, full_file_folder, file_name)
    # ----------------------------- график T_eff --------------------------------

    file_name = "save_phi_range.txt"
    main_service.save_arr_as_txt(top_column.outer_surface.phi_range, full_file_folder, file_name)
    main_service.save_arr_as_txt(top_column.outer_surface.phi_range, '', file_name)

    file_name = "save_theta_range.txt"
    main_service.save_arr_as_txt(top_column.outer_surface.theta_range, full_file_folder, file_name)
    main_service.save_arr_as_txt(top_column.outer_surface.theta_range, '', file_name)

    print('phi_theta_range saved')
    # ----------------- углы для нахождения пересечений -------------------------
    theta_accretion_begin = top_column.outer_surface.theta_range[0]
    theta_accretion_end = top_column.outer_surface.theta_range[-1]
    # ---------------------------------------------------------------------------

    # словарь для того чтобы можно было в цикле ходить по поверхностям
    surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface, 2: bot_column.outer_surface,
                3: bot_column.inner_surface}

    # ------------------ начало заполнения матриц косинусов ---------------------------
    # распараллелил

    # for key, surface in surfaces.items():
    #     surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end, top_column.outer_surface.phi_range,
    #                                bot_column.outer_surface.phi_range, e_obs)

    for key, surface in surfaces.items():
        with mp.Pool(mp.cpu_count()) as pool:
            result_cos_psi_range = pool.starmap(surface.async_fill_cos_psi_range,
                                                zip(range(config.t_max), repeat(theta_accretion_begin),
                                                    repeat(theta_accretion_end),
                                                    repeat(top_column.outer_surface.phi_range),
                                                    repeat(bot_column.outer_surface.phi_range), repeat(e_obs)))

        cos_psi_range_final = []
        for cos_psi in result_cos_psi_range:
            cos_psi_range_final.append(cos_psi)
        surface.cos_psi_range = cos_psi_range_final
    # ------------------ конец заполнения матриц косинусов ---------------------------

    # ------------------ начало заполнения массивов светимости -----------------------
    arr_simps_integrate = [0] * 4
    sum_simps_integrate = 0
    # for key, surface in surfaces.items():
    #     arr_simps_integrate[key] = surface.calculate_integral_distribution()
    #     # file_name = "%s %s %d.txt" % (file_name_variables, approx_method, key)
    #     # np.savetxt(file_name, arr_simps_integrate[key])
    #     sum_simps_integrate += np.array(arr_simps_integrate[key])

    # распараллелил
    for key, surface in surfaces.items():
        with mp.Pool(mp.cpu_count()) as pool:
            arr_simps_integrate[key] = pool.map(surface.async_calculate_integral_distribution, range(config.t_max))
        sum_simps_integrate += np.array(arr_simps_integrate[key])
    # ------------------ конец заполнения массивов светимости -----------------------

    print('ksi_shock = %f' % bot_column.outer_surface.ksi_shock)

    # --------------------- вывод графика светимости ----------------------------

    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    # 1 слитый массив - отдельные поверхности + сумма
    combined_arrays = np.append(arr_simps_integrate, [sum_simps_integrate], 0)
    arr_to_plt = [0] * len(combined_arrays)
    # нужно расширить массивы, чтобы покрыть фазу [0,2]
    for i in range(len(combined_arrays)):
        arr_to_plt[i] = main_service.extend_arr_for_phase(combined_arrays[i])

    labels_arr = ['top outer', 'top inner', 'bot outer', 'bot inner', 'sum']
    fig_title = 'total luminosity of surfaces'
    fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr, x_axis_label='phase',
                                     y_axis_label='luminosity, erg/s',
                                     figure_title=fig_title)
    file_name = 'total_luminosity_of_surfaces.png'
    main_service.save_figure(fig, full_file_folder, file_name)
    # --------------------- вывод графика светимости ----------------------------

    print('total_luminosity_of_surfaces.png saved')

    file_name = 'total_luminosity_of_surfaces.txt'
    main_service.save_arr_as_txt(arr_to_plt, full_file_folder, file_name)

    # --------------------- вывод графика углов наблюдателя ----------------------------
    observer_theta = [0] * config.t_max_for_plot
    observer_phi = [0] * config.t_max_for_plot

    for t in range(config.t_max_for_plot):
        phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
        # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
        A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                     config.betta_mu)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
        observer_phi[t], observer_theta[t] = vectors.get_angles_from_vector(e_obs_mu)

    labels_arr = [r'$\theta_{observer}$', r'$\phi_{observer}$']
    fig_title = 'Observer angles'
    combined_arrays_for_plot = np.append([observer_theta], [observer_phi], 0)
    fig = main_service.create_figure(phi_for_plot, combined_arrays_for_plot, labels_arr, x_axis_label='phase',
                                     figure_title=fig_title)
    # --------------------- вывод графика углов наблюдателя ----------------------------

    file_name = 'Observer_angles.png'
    main_service.save_figure(fig, full_file_folder, file_name)

    # ------------------ цикл для диапазона энергий ----------------------
    energy_i = 0
    N_energy = 10
    PF = [0] * N_energy
    folder = 'luminosity_in_range/'
    while energy_i < N_energy:
        if energy_i == 0:
            energy_min = 1
            energy_max = 4
        else:
            energy_min = energy_max
            energy_max = energy_max + 4

        # try:
        #     energy_bot = float(input('введите нижний предел в КэВ: '))
        #     energy_top = float(input('введите верхний предел в КэВ: '))
        # except ValueError:
        #     break
        # ------------------ начало заполнения массивов светимости -----------------------
        arr_simps_integrate = [0] * 4
        sum_simps_integrate = 0
        for key, surface in surfaces.items():
            arr_simps_integrate[key] = surface.calculate_integral_distribution_in_range(energy_min, energy_max)
            sum_simps_integrate += np.array(arr_simps_integrate[key])
        # ------------------ конец заполнения массивов светимости -----------------------

        PF[energy_i] = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)
        # --------------------- вывод графика светимости ----------------------------
        phi_for_plot = list(
            config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

        arr_to_plt = main_service.extend_arr_for_phase(sum_simps_integrate)
        fig_title = 'luminosity in range %0.2f - %0.2f KeV of surfaces, PF = %0.3f' % (
            energy_min, energy_max, PF[energy_i])
        fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr='sum', figure_title=fig_title,
                                         is_y_2d=False)

        # --------------------- вывод графика светимости ----------------------------
        file_name = 'luminosity_in_range%0.2f_-_%0.2f_KeV_of_surfaces.png' % (energy_min, energy_max)
        main_service.save_figure(fig, full_file_folder + folder, file_name)

        file_name = "sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (energy_min, energy_max)
        main_service.save_arr_as_txt(arr_to_plt, full_file_folder + folder + 'txt/', file_name)

        energy_i += 1

    # ------------------ цикл для диапазона энергий ----------------------
    file_name = "PF.txt"
    # full_file_name = full_file_folder + folder + file_name
    # np.savetxt(full_file_name, PF)
    main_service.save_arr_as_txt(PF, full_file_folder + folder, file_name)

    print('Spectral Energy Distribution')
    energy_i = 0
    N_energy = 10
    PF = [0] * N_energy
    folder = 'nu_L_nu/'
    # ------------------ цикл для Spectral Energy Distribution ------------------
    while energy_i < N_energy:
        if energy_i == 0:
            energy = 1
        else:
            energy = energy_i * 4
        # try:
        #     energy = float(input('введите энергию в КэВ: '))
        # except ValueError:
        #     break
        # ------------------ начало заполнения массивов Spectral Energy -----------------------
        arr_simps_integrate = [0] * 4
        sum_simps_integrate = 0
        for key, surface in surfaces.items():
            arr_simps_integrate[key] = main_service.fill_arr_with_func(
                AccretionColumn.Surface.calculate_nu_L_nu_on_energy, surface,
                energy)
            sum_simps_integrate += np.array(arr_simps_integrate[key])
        # ------------------ конец заполнения массивов Spectral Energy -----------------------

        PF[energy_i] = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

        # --------------------- вывод графика светимости ----------------------------

        phi_for_plot = list(
            config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
        arr_to_plt = main_service.extend_arr_for_phase(sum_simps_integrate)
        fig_title = 'Spectral Energy Distribution of energy %0.2f KeV of surfaces, PF = %0.3f' % (energy, PF[energy_i])
        fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr=r'$\nu \cdot L_{\nu}(\nu)$',
                                         x_axis_label='phase',
                                         y_axis_label='Spectral energy, erg/s', figure_title=fig_title, is_y_2d=False)

        # --------------------- вывод графика светимости ----------------------------

        file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy
        main_service.save_figure(fig, full_file_folder + folder, file_name)

        file_name = "nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energy
        main_service.save_arr_as_txt(arr_to_plt, full_file_folder + folder + 'txt/', file_name)

        energy_i += 1
        # ------------------ цикл для Spectral Energy Distribution ------------------

    file_name = "PF.txt"
    main_service.save_arr_as_txt(PF, full_file_folder + folder, file_name)

    print('спектр')
    energy_i = 0
    N_energy = 10
    PF = [0] * N_energy
    folder = 'L_nu/'
    # ------------------ цикл для Spectrum Distribution ------------------
    while energy_i < N_energy:
        if energy_i == 0:
            energy = 1
        else:
            energy = energy_i * 4
        # try:
        #     energy = float(input('введите энергию в КэВ: '))
        # except ValueError:
        #     break

        # ------------------ начало заполнения массивов Spectral Energy -----------------------
        arr_simps_integrate = [0] * 4
        sum_simps_integrate = 0
        for key, surface in surfaces.items():
            arr_simps_integrate[key] = main_service.fill_arr_with_func(AccretionColumn.Surface.calculate_L_nu_on_energy,
                                                                       surface,
                                                                       energy)
            sum_simps_integrate += np.array(arr_simps_integrate[key])
        # ------------------ конец заполнения массивов Spectral Energy -----------------------

        PF[energy_i] = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

        # --------------------- вывод графика светимости ----------------------------
        phi_for_plot = list(
            config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

        arr_to_plt = main_service.extend_arr_for_phase(sum_simps_integrate)

        fig_title = 'Spectrum of energy %0.2f KeV of surfaces, PF = %0.3f' % (energy, PF[energy_i])
        fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr=r'$\nu \cdot L_{\nu}(phase)$',
                                         x_axis_label='phase',
                                         y_axis_label=r'$Spectrum, erg \, s^{-1} \, hz^{-1}$', figure_title=fig_title,
                                         is_y_2d=False)

        # --------------------- вывод графика светимости ----------------------------
        file_name = 'L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy
        main_service.save_figure(fig, full_file_folder + folder, file_name)

        file_name = "L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energy
        main_service.save_arr_as_txt(arr_to_plt, full_file_folder + folder + 'txt/', file_name)

        energy_i += 1
        # ------------------ цикл для Spectrum Distribution ------------------

    file_name = "PF.txt"
    main_service.save_arr_as_txt(PF, full_file_folder + folder, file_name)
