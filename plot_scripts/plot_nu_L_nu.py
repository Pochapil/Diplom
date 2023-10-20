import numpy as np
import matplotlib.pyplot as plt
import scienceplots

import config
import main_service
import accretionColumnService


def plot_figs():
    plt.style.use(['science', 'notebook', 'grid'])
    # plt.style.use(['default'])

    folder = 'nu_L_nu/'
    working_folder = config.full_file_folder + folder

    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\nu L_{\nu} \: [erg/s]$'
    # -------------------------------------------------------------------------------------------
    file_name = "PF.txt"
    PF = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = "energy.txt"
    energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)
    N_energy = config.N_energy

    fig = main_service.create_figure(energy_arr, PF, x_axis_label=r'$\Phi$', y_axis_label='PF', is_y_2d=False)

    file_name = "PF.png"
    main_service.save_figure(fig, working_folder, file_name)

    # -------------------------------read nu L_nu data---------------------------------------------------
    file_name = "nu_L_nu.txt"
    data_array = main_service.load_arr_from_txt(working_folder, file_name)

    arr_to_plt = [0] * len(data_array)
    for i in range(len(data_array)):
        arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

    # -------------------------------------nu L_nu on dif energy--------------------------------------------
    for energy_i in range(N_energy):
        # fig_title = 'Spectral Energy Distribution of energy %0.2f KeV of surfaces, PF = %0.3f' % (
        #     energy_arr[energy_i], PF[energy_i])
        fig = main_service.create_figure(phi_for_plot, arr_to_plt[energy_i], x_axis_label=x_axis_label,
                                         y_axis_label=y_axis_label, is_y_2d=False)

        file_name = 'nu_L_nu_of_energy_%0.2f_KeV_of_surfaces.png' % energy_arr[energy_i]
        main_service.save_figure(fig, working_folder, file_name)

    # -----------------------------------nu L_nu all energy in one--------------------------------------------
    # по идее переписать!!
    labels_arr = [''] * N_energy
    for i in range(N_energy):
        labels_arr[i] = '%0.2f keV' % energy_arr[i]

    fig_title = r'$\nu L_{\nu}$'
    fig = main_service.create_figure(phi_for_plot, arr_to_plt, labels_arr=labels_arr, x_axis_label=x_axis_label,
                                     y_axis_label=y_axis_label, figure_title=fig_title)

    # ------------------------------------nu L_nu(nu) on phase------------------------------------------------
    x_axis_label = r'$h \nu$' + ' [keV]'
    phase_index = 0  # индекс фазы для nu_L_nu(nu)

    nu_L_nu = [0] * N_energy
    for i in range(N_energy):
        nu_L_nu[i] = arr_to_plt[i][phase_index]

    fig_title = r'$\nu L_{\nu}$'
    fig = main_service.create_figure(energy_arr, nu_L_nu, x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                                     figure_title=fig_title, is_y_2d=False)

    file_name = 'nu_L_nu(nu)' + '.png'
    main_service.save_figure(fig, working_folder, file_name)

    # ------------------------------------------nu L_nu(nu)_avg------------------------------------------------
    nu_L_nu_avg_on_phase = [0] * N_energy
    for i in range(N_energy):
        nu_L_nu_avg_on_phase[i] = np.mean(arr_to_plt[i])

    fig_title = r'$\nu L_{\nu}$'
    fig = main_service.create_figure(energy_arr, nu_L_nu_avg_on_phase, x_axis_label=x_axis_label,
                                     y_axis_label=y_axis_label,
                                     figure_title=fig_title, is_y_2d=False)

    file_name = 'nu_L_nu(nu)_avg' + '.png'
    main_service.save_figure(fig, working_folder, file_name)

    fig_title = r'$\nu L_{\nu}$'
    fig = main_service.create_figure(energy_arr, nu_L_nu_avg_on_phase, x_axis_label=x_axis_label,
                                     y_axis_label=y_axis_label,
                                     figure_title=fig_title, is_y_2d=False, is_x_log_scale=True, is_y_log_scale=True)

    file_name = 'nu_L_nu(nu)_avg_log_log' + '.png'
    main_service.save_figure(fig, working_folder, file_name)

    # -----------------------------------pretty fig-------------------------------------------------
    plt.style.use(['science', 'notebook', 'grid'])

    N_column_plot = config.N_column_plot
    energy_indexes = config.energy_indexes
    plt.rc('font', size=18)
    fig, axes = plt.subplots(N_column_plot, 1, figsize=(12, 3 * N_column_plot), sharex=True)
    for i in range(N_column_plot):
        ax = axes[i]
        label = "%0.1f keV\n PF=%0.3f" % (energy_arr[energy_indexes[i]], PF[energy_indexes[i]])
        ax.tick_params(axis='both', labelsize=12)
        ax.plot(phi_for_plot, arr_to_plt[energy_indexes[i]], color='black', lw=0.8)
        # ax.plot(phi_for_plot, arr_to_plt[i], color='black', lw=0.8, label=label)
        ax.text(0.98, 0.87, label, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black', alpha=1),
                ha='right', va='top')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        # ax.set_xlabel('10 kev PF=0.3', fontsize=16)

        # ax.set_title("%0.1f KeV PF=%0.3f" % (energy_arr[energy_indexes[i]], PF[energy_indexes[i]]), fontsize=22)

        # ax.legend(loc='upper right')

    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel('phase')
    # plt.ylabel('luminosity [erg/s]')
    plt.rc('font', size=24)
    fig.text(0.5, 0.07, r'$\Phi$', ha='center')
    fig.text(0.06, 0.5, r'$\nu \cdot L_{\nu} \: [erg/s]$', va='center', rotation='vertical')
    file_name = 'pretty_fig.png'
    main_service.save_figure(fig, working_folder, file_name)
    plt.rc('font', size=10)

    # ------------------------------log with black body Tmin Tmax--------------------------------------------------
    freq_arr = accretionColumnService.get_frequency_from_energy(energy_arr)

    file_name = "surfaces_T_eff.txt"
    T_eff = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    black_body_max = accretionColumnService.plank_energy_on_frequency(freq_arr, max(T_eff[0]))
    black_body_min = accretionColumnService.plank_energy_on_frequency(freq_arr, min(T_eff[0]))

    coeff_max = max(nu_L_nu_avg_on_phase) / max(black_body_max)
    coeff_min = max(nu_L_nu_avg_on_phase) / max(black_body_min)

    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(111)

    ax.plot(energy_arr, coeff_max * np.array(black_body_max), label='black body max T=%0.1f' % max(T_eff[0]))
    ax.plot(energy_arr, coeff_min * np.array(black_body_min), label='black body min T=%0.1f' % min(T_eff[0]))
    ax.plot(energy_arr, nu_L_nu_avg_on_phase, label=r'$\nu \cdot L_{\nu} \, avg$', color='black')

    plt.xscale('log')
    plt.yscale('log')
    ax.legend()

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)
    plt.ylim(min(nu_L_nu_avg_on_phase))

    file_name = 'nu_L_nu(nu)_avg_and_black_body' + '.png'
    main_service.save_figure(fig, working_folder, file_name)

    # ------------------------------colormesh----------------------------------------------

    phase = np.linspace(0, 1, config.t_max)
    fig, ax = plt.subplots()

    # нормировка на nu_L_nu_avg
    data_to_plot = []
    for arr in data_array:
        data_to_plot.append(arr * max(nu_L_nu_avg_on_phase) / max(arr))

    im = ax.pcolormesh(phase, energy_arr, data_to_plot)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$h \nu$' + ' [keV]'

    ax.set_xlabel(x_axis_label, fontsize=24)
    ax.set_ylabel(y_axis_label, fontsize=24)
    plt.colorbar(im)
    file_name = 'nu_L_nu_colormesh' + '.png'
    main_service.save_figure(fig, working_folder, file_name)


if __name__ == '__main__':
    i_angle = 60
    betta_mu = 40
    mc2 = 200
    a_portion = 0.25
    fi_0 = 0

    config.set_e_obs(i_angle, 0)
    config.set_betta_mu(betta_mu)
    config.M_rate_c2_Led = mc2
    config.a_portion = a_portion
    config.phi_accretion_begin_deg = fi_0

    config.update()

    plot_figs()
