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

if __name__ == '__main__':
    R_alfven = (config.mu ** 2 / (2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
    R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
    print('R_e = %f' % (R_e / config.R_ns))
    R_e_outer_surface, R_e_inner_surface = R_e, R_e
    # вектор на наблюдателя в системе координат двойной системы
    e_obs = np.array([0, np.sin(config.i_angle), np.cos(config.i_angle)])
    file_name_variables = "betta_omega=%d betta_mu=%d a_portion=%f M_rate_c2_Led=%d" \
                          % (config.betta_rotate, config.betta_mu, config.a_portion, config.M_rate_c2_Led)
    approx_method = accretionColumnService.approx_method

    file_folder = 'figs/'
    args_folder = 'a=%0.2f fi_0=%d/' % (config.a_portion, config.phi_accretion_begin_deg)
    file_folder = file_folder + args_folder
    pathlib.Path(file_folder).mkdir(parents=True, exist_ok=True)

    # от поверхности NS - угол при котором радиус = радиусу НЗ
    # ----------------- начало инициализации верхней колонки ------------------------
    theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(R_e_outer_surface)
    theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(R_e_inner_surface)

    top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                                 theta_accretion_begin_inner_surface, True)
    # ----------------- конец инициализации верхней колонки ------------------------

    # ----------------- начало инициализации нижней колонки ------------------------
    theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
    theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(top_column.outer_surface.theta_range, top_column.outer_surface.T_eff)
    plt.show()
    plt.close()

    bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                                 theta_accretion_begin_inner_surface, False)
    # ----------------- конец инициализации нижней колонки ------------------------

    surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface, 2: bot_column.outer_surface,
                3: bot_column.inner_surface}

    print('phi_theta_range saved')
    file_name = "save_phi_range.txt"
    full_file_name = file_folder + file_name
    np.savetxt(full_file_name, top_column.outer_surface.phi_range)
    np.savetxt(file_name, top_column.outer_surface.phi_range)
    file_name = "save_theta_range.txt"
    full_file_name = file_folder + file_name
    np.savetxt(full_file_name, top_column.outer_surface.theta_range)
    np.savetxt(file_name, top_column.outer_surface.theta_range)

    # print('T_eff:')
    # print(top_column.outer_surface.T_eff)

    # ----------------- углы для нахождения пересечений -------------------------
    theta_accretion_begin = top_column.outer_surface.theta_range[0]
    theta_accretion_end = top_column.outer_surface.theta_range[-1]
    # ---------------------------------------------------------------------------

    # ------------------ начало заполнения матриц косинусов ---------------------------
    # распараллелил
    for key, surface in surfaces.items():
        with mp.Pool(mp.cpu_count()) as pool:
            result = pool.starmap(surface.async_fill_cos_psi_range,
                                  zip(range(config.t_max), repeat(theta_accretion_begin), repeat(theta_accretion_end),
                                      repeat(top_column.outer_surface.phi_range),
                                      repeat(bot_column.outer_surface.phi_range), repeat(e_obs)))

        cos_psi_range_final = []
        for num in result:
            cos_psi_range_final.append(num)
        surface.cos_psi_range = cos_psi_range_final
    # ------------------ конец заполнения матриц косинусов ---------------------------

    # ------------------ начало заполнения массивов светимости -----------------------
    arr_simps_integrate = [0] * 4
    sum_simps_integrate = 0
    for key, surface in surfaces.items():
        arr_simps_integrate[key] = surface.calculate_integral_distribution()
        # file_name = "%s %s %d.txt" % (file_name_variables, approx_method, key)
        # np.savetxt(file_name, arr_simps_integrate[key])
        sum_simps_integrate += np.array(arr_simps_integrate[key])
    # ------------------ конец заполнения массивов светимости -----------------------

    print('ksi_shock = %f' % bot_column.outer_surface.ksi_shock)

    # --------------------- вывод графика светимости ----------------------------
    fig = plt.figure(figsize=(8, 8))
    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    append_index = config.t_max_for_plot - config.t_max
    for i in range(4):
        arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
    sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

    ax = fig.add_subplot(111)
    ax.plot(phi_for_plot, arr_simps_integrate[0], label='top outer')
    ax.plot(phi_for_plot, arr_simps_integrate[1], label='top inner')
    ax.plot(phi_for_plot, arr_simps_integrate[2], label='bot outer', marker='*')
    ax.plot(phi_for_plot, arr_simps_integrate[3], label='bot inner')
    ax.plot(phi_for_plot, sum_simps_integrate, label='sum')
    ax.legend()
    fig.suptitle('total luminosity of surfaces', fontsize=14)
    # plt.yscale('log')
    plt.show()
    plt.close()
    # --------------------- вывод графика светимости ----------------------------

    file_name = 'total_luminosity_of_surfaces.png'
    full_file_name = file_folder + file_name
    fig.savefig(full_file_name, dpi=fig.dpi)

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

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(phi_for_plot, observer_theta, label=r'$\theta_{observer}$')
    ax.plot(phi_for_plot, observer_phi, label=r'$\phi_{observer}$')
    ax.legend()
    fig.suptitle('Observer angles', fontsize=14)
    plt.show()
    plt.close()
    # --------------------- вывод графика углов наблюдателя ----------------------------

    file_name = 'Observer_angles.png'
    full_file_name = file_folder + file_name
    fig.savefig(full_file_name, dpi=fig.dpi)

    # ------------------ цикл для диапазона энергий ----------------------
    energy_i = 0
    while energy_i < 10:
        if energy_i == 0:
            energy_bot = 1
            energy_top = 4
        else:
            energy_bot = energy_top
            energy_top = energy_top + 4
        energy_i += 1
        # try:
        #     energy_bot = float(input('введите нижний предел в КэВ: '))
        #     energy_top = float(input('введите верхний предел в КэВ: '))
        # except ValueError:
        #     break
        # ------------------ начало заполнения массивов светимости -----------------------
        arr_simps_integrate = [0] * 4
        sum_simps_integrate = 0
        for key, surface in surfaces.items():
            arr_simps_integrate[key] = surface.calculate_integral_distribution_in_range(energy_bot, energy_top)
            sum_simps_integrate += np.array(arr_simps_integrate[key])
        # ------------------ конец заполнения массивов светимости -----------------------

        PF = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

        # --------------------- вывод графика светимости ----------------------------
        fig = plt.figure(figsize=(8, 8))
        phi_for_plot = list(
            config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

        append_index = config.t_max_for_plot - config.t_max
        for i in range(4):
            arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
        sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

        ax = fig.add_subplot(111)
        # ax.plot(phi_for_plot, arr_simps_integrate[0], label='top outer')
        # ax.plot(phi_for_plot, arr_simps_integrate[1], label='top inner')
        # ax.plot(phi_for_plot, arr_simps_integrate[2], label='bot outer', marker='*')
        # ax.plot(phi_for_plot, arr_simps_integrate[3], label='bot inner')
        ax.plot(phi_for_plot, sum_simps_integrate, label='sum')
        ax.legend()
        fig_title = 'luminosity in range %0.2f - %0.2f KeV of surfaces, PF = %0.3f' % (energy_bot, energy_top, PF)
        fig.suptitle(fig_title, fontsize=14)
        # plt.yscale('log')
        # plt.show()
        # --------------------- вывод графика светимости ----------------------------

        folder = 'luminosity_in_range/'
        pathlib.Path(file_folder + folder).mkdir(parents=True, exist_ok=True)

        pathlib.Path(file_folder + folder + 'txt/').mkdir(parents=True, exist_ok=True)
        file_name = "txt/sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (energy_bot, energy_top)
        full_file_name = file_folder + folder + file_name
        np.savetxt(full_file_name, sum_simps_integrate)

        file_name = 'luminosity_in_range%0.2f_-_%0.2f_KeV_of_surfaces.png' % (energy_bot, energy_top)
        full_file_name = file_folder + folder + file_name
        fig.savefig(full_file_name, dpi=fig.dpi)
        plt.close()
    # ------------------ цикл для диапазона энергий ----------------------

    print('спектр')
    energy_i = 1
    while energy_i < 32:
        # try:
        #     energy = float(input('введите энергию в КэВ: '))
        # except ValueError:
        #     break
        energy = energy_i
        energy_i += 1
        # ------------------ начало заполнения массивов светимости -----------------------
        arr_simps_integrate = [0] * 4
        sum_simps_integrate = 0
        for key, surface in surfaces.items():
            arr_simps_integrate[key] = surface.calculate_luminosity_on_energy(energy)
            sum_simps_integrate += np.array(arr_simps_integrate[key])
        # ------------------ конец заполнения массивов светимости -----------------------

        PF = accretionColumnService.get_pulsed_fraction(sum_simps_integrate)

        # --------------------- вывод графика светимости ----------------------------
        fig = plt.figure(figsize=(8, 8))
        phi_for_plot = list(
            config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

        append_index = config.t_max_for_plot - config.t_max
        for i in range(4):
            arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
        sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

        ax = fig.add_subplot(111)
        ax.plot(phi_for_plot, sum_simps_integrate, label=r'$\nu * L_{\nu}(\nu)$')
        ax.legend()
        fig_title = 'luminosity of energy %0.2f KeV of surfaces, PF = %0.3f' % (energy, PF)
        fig.suptitle(fig_title, fontsize=14)
        # plt.yscale('log')
        # plt.show()
        # --------------------- вывод графика светимости ----------------------------

        folder = 'nu_L_nu/'
        pathlib.Path(file_folder + folder).mkdir(parents=True, exist_ok=True)

        pathlib.Path(file_folder + folder + 'txt/').mkdir(parents=True, exist_ok=True)
        file_name = "txt/nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.txt" % energy
        full_file_name = file_folder + folder + file_name
        np.savetxt(full_file_name, sum_simps_integrate)

        file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy
        full_file_name = file_folder + folder + file_name
        fig.savefig(full_file_name, dpi=fig.dpi)
        plt.close()
