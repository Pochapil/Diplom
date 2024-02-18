import matplotlib.pyplot as plt
import numpy as np

import config
import main_service

fi_0 = [0]

mc2 = [30]
a_portion = [0.65]

i_angle = [60]
betta_mu = [40]

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        for i in range(len(mc2)):
            for j in range(len(a_portion)):
                for k in range(len(fi_0)):
                    config.set_e_obs(i_angle[i_angle_index], 0)
                    config.set_betta_mu(betta_mu[betta_mu_index])

                    config.M_rate_c2_Led = mc2[i]
                    config.a_portion = a_portion[j]
                    config.phi_accretion_begin_deg = fi_0[k]

                    config.update()

                    file_name = "save_tau_range.txt"
                    tau_array = main_service.load_arr_from_txt(config.full_file_folder, file_name)
                    min_tau = min(tau_array)

                    file_name = "save_theta_range.txt"
                    theta_range = main_service.load_arr_from_txt(config.full_file_folder, file_name)

                    step_theta_accretion = (np.pi / 2 - theta_range[-1]) / (config.N_theta_accretion - 1)
                    theta_range_for_tau = np.array(
                        [theta_range[-1] + step_theta_accretion * j for j in range(config.N_theta_accretion)])

                    R_alfven = (config.mu ** 2 / (
                            2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
                    R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
                    # print('R_e = %f' % (R_e / config.R_ns))
                    R_e_outer_surface, R_e_inner_surface = R_e, R_e  # допущение что толщина = 0
                    x = R_e / config.R_ns * np.sin(theta_range_for_tau) ** 2

                    fig = plt.figure(figsize=(12, 6))

                    ax = fig.add_subplot(111)
                    #ax.set_xscale('log')
                    #ax.set_yscale('log')
                    ax.plot(x, tau_array, label=f'min_val = {min_tau:.2f}')
                    plt.legend()
                    plt.show()
