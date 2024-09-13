import config
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import matplotlib as mpl

import main_service

plt.style.use(['science', 'notebook', 'grid'])  # для красивых графиков

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'


def plot_save_sky_map(obs_i_angle_arr):
    file_folder = config.PROJECT_DIR + 'figs/sky_map/betta_mu=%d/' % config.betta_mu_deg
    working_folder = config.full_file_folder

    # print('beta_mu=%d' % config.betta_mu_deg)

    # title = 'beta_mu=%d mc2=%d a=%0.2f fi_0=%d' % (
    #     config.betta_mu_deg, config.M_rate_c2_Led, config.a_portion, config.phi_accretion_begin_deg)

    title = (r'$\beta_{\mu}$') + ('=%d' % config.betta_mu_deg) + (r'$^\circ$') + ' ' + (r'$\dot{m}$') + \
            ('=%d ' % config.M_rate_c2_Led) + ('a=%0.2f ' % config.a_portion) + (r'$\phi_0$') + \
            ('=%d' % config.phi_accretion_begin_deg) + (r'$^\circ$')

    with open(working_folder + 'save_values.txt') as f:
        lines = f.readlines()
        L_x = float(lines[3][12:20]) * 10 ** float(lines[3][27:29])
        # total L_x = 4.396383 * 10**38 - 12 это индекс начала числа, 27-29 это степень 10

    file_name = 'total_luminosity_of_surfaces.txt'
    data_array = [0] * len(obs_i_angle_arr)

    for i in range(len(obs_i_angle_arr)):
        config.set_e_obs(obs_i_angle_arr[i], 0)
        data_array[i] = main_service.load_arr_from_txt(config.full_file_folder, file_name)[4]

    phase = np.linspace(0, 1, config.t_max)

    # fig, ax = plt.subplots()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # нормировка на L_nu_avg
    data_to_plot = []
    for arr in data_array:
        data_to_plot.append(arr / L_x)

    im = ax.pcolormesh(phase, obs_i_angle_arr, data_to_plot)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{obs}$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    fig_title = r'$L_{iso}/L_{x}$ ' + title
    fig.suptitle(fig_title, fontsize=14)

    plt.colorbar(im)

    working_folder = file_folder + config.file_folder_accretion_args
    file_name = 'map' + '.png'
    main_service.save_figure(fig, working_folder, file_name)
    plt.close('all')


def plot_save_sky_map_contour(obs_i_angle_arr):
    file_folder = config.PROJECT_DIR + 'figs/sky_map/betta_mu=%d/' % config.betta_mu_deg

    working_folder = config.full_file_folder

    # title = 'beta_mu=%d mc2=%d a=%0.2f fi_0=%d' % (
    #     config.betta_mu_deg, config.M_rate_c2_Led, config.a_portion, config.phi_accretion_begin_deg)

    title = (r'$\beta_{\mu}$') + ('=%d' % config.betta_mu_deg) + (r'$^\circ$') + ' ' + (r'$\dot{m}$') + \
            ('=%d ' % config.M_rate_c2_Led) + ('a=%0.2f ' % config.a_portion) + (r'$\phi_0$') + \
            ('=%d' % config.phi_accretion_begin_deg) + (r'$^\circ$')

    with open(working_folder + 'save_values.txt') as f:
        lines = f.readlines()
        L_x = float(lines[3][12:20]) * 10 ** float(lines[3][27:29])
        # total L_x = 4.396383 * 10**38 - 12 это индекс начала числа, 27-29 это степень 10

    # file_name = "L_x.txt"
    # L_x = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = 'total_luminosity_of_surfaces.txt'
    data_array = [0] * len(obs_i_angle_arr)

    for i in range(len(obs_i_angle_arr)):
        config.set_e_obs(obs_i_angle_arr[i], 0)
        data_array[i] = main_service.load_arr_from_txt(config.full_file_folder, file_name)[4]

    check_index = 2
    # print(np.roll(data_array[check_index], len(data_array[check_index]) // 2))
    # print(data_array[-check_index - 3])

    phase = np.linspace(0, 1, config.t_max)

    # fig, ax = plt.subplots()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # нормировка на L_nu_avg
    data_to_plot = []
    for arr in data_array:
        data_to_plot.append(arr / L_x)

    im = ax.contourf(phase, obs_i_angle_arr, data_to_plot, levels=30)
    # cmap = 'cividis', cmap='Spectral', cmap='Spectral_r'

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{obs} \, [^{\circ}]$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    # fig_title = 'L_iso/L_x'
    # fig.suptitle(fig_title, fontsize=14)

    fig_title = r'$L_{iso}/L_{x}$ ' + title
    # fig.suptitle(fig_title, fontsize=14)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)
    # clb.ax.set_title(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)

    working_folder = file_folder + config.file_folder_accretion_args
    file_name = 'map_contour' + '.png'
    main_service.save_figure(fig, working_folder, file_name)
    plt.close('all')


def try_sky_map(obs_i_angle_arr):
    file_folder = config.PROJECT_DIR + 'figs/sky_map/betta_mu=%d/' % config.betta_mu_deg
    working_folder = config.full_file_folder

    # title = 'beta_mu=%d mc2=%d a=%0.2f fi_0=%d' % (
    #     config.betta_mu_deg, config.M_rate_c2_Led, config.a_portion, config.phi_accretion_begin_deg)

    title = (r'$\beta_{\mu}$') + ('=%d' % config.betta_mu_deg) + (r'$^\circ$') + ' ' + (r'$\dot{m}$') + \
            ('=%d ' % config.M_rate_c2_Led) + ('a=%0.2f ' % config.a_portion) + (r'$\phi_0$') + \
            ('=%d' % config.phi_accretion_begin_deg) + (r'$^\circ$')

    with open(working_folder + 'save_values.txt') as f:
        lines = f.readlines()
        L_x = float(lines[3][12:20]) * 10 ** float(lines[3][27:29])
        # total L_x = 4.396383 * 10**38 - 12 это индекс начала числа, 27-29 это степень 10

    file_name = 'total_luminosity_of_surfaces.txt'
    data_array = [0] * 17

    # roll -- циклическая перестановка - делаем так как симметрическая задача и для угла 180 - theta будет симметрично
    # для сдвига на полфазы -- можно расчитать только до 90 а потом для других переставить и получить для до 180
    for i in range(len(obs_i_angle_arr)):
        config.set_e_obs(obs_i_angle_arr[i], 0)
        file_name = 'total_luminosity_of_surfaces.txt'
        data_array[i] = main_service.load_arr_from_txt(config.full_file_folder, file_name)[4]

        file_name = 'scattered_energy_bot.txt'
        buf = main_service.load_arr_from_txt(config.full_file_folder + 'scattered_on_magnet_lines/', file_name)
        if not np.isnan(buf).any():
            data_array[i] += buf

        file_name = 'scattered_energy_top.txt'
        buf = main_service.load_arr_from_txt(config.full_file_folder + 'scattered_on_magnet_lines/', file_name)
        if not np.isnan(buf).any():
            data_array[i] += buf

        if i != 8:
            data_array[-i - 1] = np.roll(data_array[i], len(data_array[i]) // 2)

    phase = np.linspace(0, 2, 2 * config.t_max)

    # fig, ax = plt.subplots()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # нормировка на L_nu_avg
    data_to_plot = []
    for arr in data_array:
        data_to_plot.append(arr / L_x)

    for i in range(len(data_to_plot)):
        data_to_plot[i] = main_service.extend_arr_for_phase(data_to_plot[i])

    obs_i_angle_arr_to_plot = np.linspace(0, 180, 17)

    im = ax.pcolormesh(phase, obs_i_angle_arr_to_plot, data_to_plot)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{obs}$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    fig_title = r'$L_{iso}/L_{x}$ ' + title
    fig.suptitle(fig_title, fontsize=14)

    plt.colorbar(im)

    working_folder = file_folder + config.file_folder_accretion_args
    file_name = 'try_map' + '.png'
    main_service.save_figure(fig, working_folder, file_name)
    plt.close('all')

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    im = ax.contourf(phase, obs_i_angle_arr_to_plot, data_to_plot, levels=30)
    # cmap = 'cividis', cmap='Spectral', cmap='Spectral_r'

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{obs} \, [^{\circ}]$'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    # fig_title = 'L_iso/L_x'
    # fig.suptitle(fig_title, fontsize=14)

    fig_title = r'$L_{iso}/L_{x}$ ' + title
    # fig.suptitle(fig_title, fontsize=14)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)
    # clb.ax.set_title(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)

    # ax.axhline(config.betta_mu_deg, c='r', linestyle="--")

    ax.scatter([0, 1, 2], [config.betta_mu_deg, config.betta_mu_deg, config.betta_mu_deg], c='white', marker='*', s=300)
    ax.scatter([0.5, 1.5], [180 - config.betta_mu_deg, 180 - config.betta_mu_deg], c='white', marker='*', s=300)

    ax.scatter([0, 1, 2], [config.betta_mu_deg, config.betta_mu_deg, config.betta_mu_deg], c='black', marker='*', s=100)
    ax.scatter([0.5, 1.5], [180 - config.betta_mu_deg, 180 - config.betta_mu_deg], c='black', marker='*', s=100)

    working_folder = file_folder + config.file_folder_accretion_args
    file_name = 'try_map_contour' + '.png'
    main_service.save_figure(fig, working_folder, file_name)
    plt.close('all')


if __name__ == '__main__':

    mc2 = [30, 100]
    a_portion = [0.25, 0.65]
    fi_0 = [60, 120, 180, 220, 240, 300]
    betta_mu = [10 * i for i in range(1, 10)]
    # betta_mu = [60]
    # obs_i_angle_arr = np.linspace(10, 180, 18)
    obs_i_angle_arr = np.linspace(10, 90, 9)

    # for betta_mu_index in range(len(betta_mu)):
    #     for mc_index in range(len(mc2)):
    #         for j in range(len(a_portion)):
    #             for k in range(len(fi_0)):
    #                 config.set_betta_mu(betta_mu[betta_mu_index])
    #
    #                 config.M_rate_c2_Led = mc2[mc_index]
    #                 config.a_portion = a_portion[j]
    #                 config.phi_accretion_begin_deg = fi_0[k]
    #
    #                 config.update()
    #
    #                 plot_save_sky_map(obs_i_angle_arr)
    #                 plot_save_sky_map_contour(obs_i_angle_arr)

    obs_i_angle_arr = np.linspace(10, 90, 9)
    mc2 = [30]
    a_portion = [0.66]
    # fi_0 = [config.fi_0_dict[a_portion[0]]]
    fi_0 = [((config.fi_0_dict[a_portion[i]] + 0) % 360) for i in range(len(a_portion))]
    betta_mu = [20]
    for betta_mu_index in range(len(betta_mu)):
        for mc_index in range(len(mc2)):
            for j in range(len(a_portion)):
                for k in range(len(fi_0)):
                    config.set_betta_mu(betta_mu[betta_mu_index])

                    config.M_rate_c2_Led = mc2[mc_index]
                    config.a_portion = a_portion[j]
                    config.phi_accretion_begin_deg = fi_0[k]

                    config.update()

                    try_sky_map(obs_i_angle_arr)
