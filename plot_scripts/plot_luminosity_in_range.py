import config
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scienceplots

import main_service


def plot_figs():
    import matplotlib as mpl
    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['font.family'] = 'STIXGeneral'

    plt.style.use(['science', 'notebook', 'grid'])

    folder = 'luminosity_in_range/'
    working_folder = config.full_file_folder + folder

    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))
    # N_energy = 10

    # -------------------------------------------------------------------------------------------
    file_name = "PF.txt"
    PF = main_service.load_arr_from_txt(working_folder, file_name)

    file_name = "energy.txt"
    energy_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)
    N_energy = config.N_energy
    #
    # energy_arr = [1] * N_energy
    # for i in range(1, N_energy):
    #     energy_arr[i] = i * 4

    # fig = main_service.create_figure(energy_arr, PF, is_y_2d=False)

    fig = main_service.create_figure(energy_arr[:-1], PF, x_axis_label=r'$h \nu$ [keV]', y_axis_label='PF',
                                     is_y_2d=False)

    file_name = "PF.png"
    main_service.save_figure(fig, working_folder, file_name)

    # -------------------------------------------------------------------------------------------
    file_name = "luminosity_in_range.txt"
    data_array = main_service.load_arr_from_txt(working_folder, file_name)

    arr_to_plt = [0] * len(data_array)
    for i in range(len(data_array)):
        arr_to_plt[i] = main_service.extend_arr_for_phase(data_array[i])

    # -------------------------------------------------------------------------------------------
    # energy_min = [0] * N_energy
    # energy_max = [0] * N_energy
    # energy_min[0] = 1
    # energy_max[0] = 4
    # for i in range(1, N_energy):
    #     energy_min[i] = 4 * i
    #     energy_max[i] = 4 * (i + 1)

    energy_min = [0] * (N_energy - 1)
    energy_max = [0] * (N_energy - 1)
    for i in range(N_energy - 1):
        energy_min[i] = energy_arr[i]
        energy_max[i] = energy_arr[i + 1]

    for energy_i in range(N_energy - 1):
        fig_title = 'luminosity in range %0.2f - %0.2f keV of surfaces, PF = %0.3f' % (
            energy_min[energy_i], energy_max[energy_i], PF[energy_i])
        fig = main_service.create_figure(phi_for_plot, arr_to_plt[energy_i], x_axis_label=r'$\Phi$',
                                         y_axis_label='L [erg/s]', figure_title=fig_title, is_y_2d=False)

        file_name = 'luminosity_in_range%0.2f_-_%0.2f_KeV_of_surfaces.png' % (
            energy_min[energy_i], energy_max[energy_i])
        main_service.save_figure(fig, working_folder, file_name)

    # -------------------------------------------------------------------------------------------
    N_energy = 6
    labels_arr = [''] * N_energy

    for i in range(N_energy):
        labels_arr[i] = "%0.2f - %0.2f keV" % (energy_min[i], energy_max[i])

    fig = main_service.create_figure(phi_for_plot, arr_to_plt[:N_energy], labels_arr=labels_arr, x_axis_label=r'$\Phi$',
                                     y_axis_label='L [erg/s]')

    file_name = 'sum_of_luminosity_in_range.png'
    main_service.save_figure(fig, working_folder, file_name)

    # -------------------------------------------------------------------------------------------
    plt.style.use(['science', 'notebook', 'grid'])

    N_column_plot = config.N_column_plot
    energy_indexes = config.energy_indexes_luminosity
    fig, axes = plt.subplots(N_column_plot, 1, figsize=(12, 3 * N_column_plot), sharex=True)
    for i in range(N_column_plot):
        ax = axes[i]
        label = "%0.1f - %0.1f keV\n PF=%0.3f" % (
            energy_min[energy_indexes[i]], energy_max[energy_indexes[i]], PF[energy_indexes[i]])
        ax.tick_params(axis='both', labelsize=12)
        ax.plot(phi_for_plot, arr_to_plt[energy_indexes[i]], color='black', lw=0.8)
        # ax.plot(phi_for_plot, arr_to_plt[i], color='black', lw=0.8, label=label)
        ax.text(0.98, 0.87, label, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black'), ha='right',
                va='top')
        # plt.subplots_adjust(hspace=0.15)
        # ax.legend(loc='upper right')

    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel(r'$\Phi$')
    # plt.ylabel('luminosity [erg/s]')
    plt.rc('font', size=24)
    fig.text(0.5, 0.07, r'$\Phi$', ha='center')
    fig.text(0.06, 0.5, 'L [erg/s]', va='center', rotation='vertical')
    file_name = 'pretty_fig.png'
    main_service.save_figure(fig, working_folder, file_name)

    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_subplot(111)
    ax.step(phi_for_plot, arr_to_plt[16], color='black', lw=0.9)
    ax.set_xlabel(r'$\Phi$', fontsize=24)
    ax.set_ylabel('L [erg/s]', fontsize=24)
    file_name = 'step.png'
    main_service.save_figure(fig, working_folder, file_name)
    plt.rc('font', size=10)


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
