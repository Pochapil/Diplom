import config
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

import main_service
import accretionColumnService


def plot_figs():
    plt.style.use(['science', 'notebook', 'grid'])
    working_folder = config.full_file_folder

    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    # -------------------------------------------------------------------------------------------
    file_name = 'total_luminosity_of_surfaces.txt'
    data_array = main_service.load_arr_from_txt(working_folder, file_name)

    arr_to_plt = [0] * len(data_array)
    # нужно расширить массивы, чтобы покрыть фазу [0,2]
    for i in range(len(data_array)):
        arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

    labels_arr = [r'$top_{pol}$', r'$top_{eq}$', r'$bottom_{pol}$', r'$bottom_{eq}$', 'sum']
    # fig_title = 'total luminosity of surfaces'
    fig_title = r'$\theta{obs}$' + ('=%d' % config.obs_i_angle_deg) + (r'$^\circ$') + ' ' + (r'$\beta_{\mu}$') + \
                ('=%d' % config.betta_mu_deg) + (r'$^\circ$') + ' ' + ('a=%0.2f ' % config.a_portion) + \
                (r'$\dot{m}$') + ('=%d ' % config.M_rate_c2_Led) + (r'$\phi_0$') + \
                ('=%d' % config.phi_accretion_begin_deg) + (r'$^\circ$')
    # fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr, x_axis_label='Phase',
    #                                  y_axis_label=r'$Luminosity \: [erg/s]$', figure_title=fig_title)

    folder = 'scattered_on_magnet_lines/'
    working_folder = config.full_file_folder + folder

    file_name = 'scattered_energy_bot.txt'
    bot_scatter = main_service.load_arr_from_txt(working_folder, file_name)
    bot_scatter = main_service.extend_arr_for_phase(bot_scatter)

    file_name = 'scattered_energy_top.txt'
    top_scatter = main_service.load_arr_from_txt(working_folder, file_name)
    top_scatter = main_service.extend_arr_for_phase(top_scatter)

    arr_to_plt[-1] += bot_scatter + top_scatter

    fig = plt.figure(figsize=(21, 10))
    ax = fig.add_subplot(111)
    line = [0] * 5
    i = 0
    ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='red', marker='.', linestyle=':', markersize=12)
    i += 1
    ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='red', marker='*', markersize=12)
    i += 1
    ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='green', marker='+', linestyle=':', markersize=12)
    i += 1
    ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='green', marker='^', markersize=12)
    i += 1
    ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i], color='black')

    ax.plot(phi_for_plot, top_scatter, label='top_scatter', color='purple', marker='*', markersize=12)
    ax.plot(phi_for_plot, bot_scatter, label='bot_scatter', color='blue', marker='^', markersize=12)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$L_{\rm{iso}}$' + ' [erg/s]'
    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)

    # for i in range(len(data_array)):
    #     line[i], = ax.plot(phi_for_plot, arr_to_plt[i], label=labels_arr[i])
    # ax.set_xlabel('Phase', fontsize=24)
    # ax.set_ylabel('L [erg/s]', fontsize=24)

    # fig.suptitle(fig_title, fontsize=24)
    ax.legend()

    file_name = 'total_luminosity_of_surfaces_with_scatter.png'
    main_service.save_figure(fig, working_folder, file_name)

    # ------------------------------
    folder = 'L_nu/'
    working_folder = config.full_file_folder + folder
    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    file_name = "energy.txt"
    energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    file_name = "L_nu.txt"
    data_array = main_service.load_arr_from_txt(working_folder, file_name)

    folder = 'scattered_on_magnet_lines/' + 'L_nu/'
    working_folder = config.full_file_folder + folder

    file_name = "top_column_scatter_L_nu.txt"
    top_column_scattered_data_array = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = "bot_column_scatter_L_nu.txt"
    bot_column_scattered_data_array = main_service.load_arr_from_txt(working_folder, file_name)

    sum_scattered_array = top_column_scattered_data_array + bot_column_scattered_data_array + data_array

    PF = [0] * config.N_energy
    for energy_index in range(config.N_energy):
        PF[energy_index] = accretionColumnService.get_pulsed_fraction(sum_scattered_array[energy_index])

    x_axis_label = r'$h \nu$' + ' [keV]'
    fig = main_service.create_figure(energy_arr, PF, x_axis_label=x_axis_label, y_axis_label='PF', is_y_2d=False)

    file_name = "PF.png"
    main_service.save_figure(fig, working_folder, file_name)

    data_array_to_plt = [0] * len(data_array)
    top_column_scattered_data_array_to_plt = [0] * len(data_array)
    bot_column_scattered_data_array_to_plt = [0] * len(data_array)
    sum_scattered_array_to_plt = [0] * len(data_array)

    for i in range(len(data_array)):
        data_array_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])
        top_column_scattered_data_array_to_plt[i] = main_service.extend_arr_for_phase(
            top_column_scattered_data_array[i])
        bot_column_scattered_data_array_to_plt[i] = main_service.extend_arr_for_phase(
            bot_column_scattered_data_array[i])
        sum_scattered_array_to_plt[i] = main_service.extend_arr_for_phase(sum_scattered_array[i])

    for energy_i in range(config.N_energy):
        fig_title = 'Spectrum of energy %0.2f keV of surfaces, PF = %0.3f' % (energy_arr[energy_i], PF[energy_i])

        y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
        x_axis_label = r'$\Phi$'

        fig = plt.figure(figsize=(21, 10))
        ax = fig.add_subplot(111)

        ax.plot(phi_for_plot, bot_column_scattered_data_array_to_plt[energy_i], label='bot_scatter', color='red')
        ax.plot(phi_for_plot, top_column_scattered_data_array_to_plt[energy_i], label='top_scatter', color='blue')
        ax.plot(phi_for_plot, data_array_to_plt[energy_i], label='without_scatter', color='green')
        ax.plot(phi_for_plot, sum_scattered_array_to_plt[energy_i], label='sum', color='black')
        ax.legend()
        fig.suptitle(fig_title, fontsize=14)

        ax.set_xlabel(x_axis_label, fontsize=24)
        ax.set_ylabel(y_axis_label, fontsize=24)
        # fig = main_service.create_figure(phi_for_plot, sum_scattered_array_to_plt[energy_i], labels_arr=r'$L_{\nu}(\Phi)$',
        #                                  x_axis_label=x_axis_label, y_axis_label=y_axis_label, figure_title=fig_title,
        #                                  is_y_2d=False)

        file_name = 'L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy_arr[energy_i]
        main_service.save_figure(fig, working_folder, file_name)


if __name__ == '__main__':
    i_angle = 60
    betta_mu = 40
    mc2 = 30
    a_portion = 0.65
    fi_0 = 0

    config.set_e_obs(i_angle, 0)
    config.set_betta_mu(betta_mu)
    config.M_rate_c2_Led = mc2
    config.a_portion = a_portion
    config.phi_accretion_begin_deg = fi_0

    config.update()

    plot_figs()
