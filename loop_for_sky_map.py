import numpy as np
import multiprocessing as mp
from itertools import repeat
import time

import config
import accretionColumnService
from accretionColumn import AccretionColumn
import main_service
import plot_sky_map

if __name__ == '__main__':
    # mc2 = [10, 30, 100]
    mc2 = [30]
    a_portion = [0.25]
    fi_0 = [0]
    #mc2 = [10, 30, 100]
    #a_portion = [0.1, 0.25, 0.65, 1]
    # fi_0 = [20 * i for i in range(1, 18)]
    # betta_mu = [30, 60, 90]
    betta_mu = [20]
    obs_i_angle = np.linspace(0, 180, 19)

    time_start = time.time()

    for betta_mu_index in range(len(betta_mu)):
        for mc_index in range(len(mc2)):
            for j in range(len(a_portion)):
                for k in range(len(fi_0)):
                    config.set_betta_mu(betta_mu[betta_mu_index])

                    config.M_rate_c2_Led = mc2[mc_index]
                    config.a_portion = a_portion[j]
                    config.phi_accretion_begin_deg = fi_0[k]

                    config.update()

                    file_folder = config.PROJECT_DIR + 'figs/sky_map/betta_mu=%d/' % config.betta_mu_deg
                    working_folder = file_folder + config.file_folder_accretion_args

                    print('beta_mu=%d' % config.betta_mu_deg)

                    # ------------------- создание папки для графиков --------------------------
                    main_service.create_file_path(working_folder)
                    # ------------------- создание папки для графиков --------------------------

                    print('calculate for beta_mu=%d mc=%0.2f a=%0.2f fi_0=%0.2f' % (
                        betta_mu[betta_mu_index], mc2[mc_index], a_portion[j], fi_0[k]))

                    for i in range(len(obs_i_angle)):
                        config.set_e_obs(obs_i_angle[i], 0)

                        R_alfven = (config.mu ** 2 / (
                                2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (
                                           2 / 7)
                        R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
                        # print('R_e = %f' % (R_e / config.R_ns))
                        print('%f' % obs_i_angle[i])
                        R_e_outer_surface, R_e_inner_surface = R_e, R_e  # допущение что толщина = 0
                        # вектор на наблюдателя в системе координат двойной системы (условимся что omega и e_obs лежат в пл-ти x0z)
                        # e_obs = np.array([np.sin(config.obs_i_angle), 0, np.cos(config.obs_i_angle)])
                        e_obs = config.e_obs

                        approx_method = accretionColumnService.approx_method

                        # ----------------- начало инициализации верхней колонки ------------------------
                        # от поверхности NS - угол при котором радиус = радиусу НЗ
                        theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(
                            R_e_outer_surface)
                        theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(
                            R_e_inner_surface)

                        top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface,
                                                     R_e_inner_surface,
                                                     theta_accretion_begin_inner_surface, True)
                        # ----------------- конец инициализации верхней колонки ------------------------

                        # ----------------- начало инициализации нижней колонки ------------------------
                        theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
                        theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

                        bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface,
                                                     R_e_inner_surface,
                                                     theta_accretion_begin_inner_surface, False)
                        # ----------------- конец инициализации нижней колонки ------------------------

                        # словарь для того чтобы можно было в цикле ходить по поверхностям
                        surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface,
                                    2: bot_column.outer_surface,
                                    3: bot_column.inner_surface}

                        # ----------------- углы для нахождения пересечений -------------------------
                        theta_accretion_begin = top_column.outer_surface.theta_range[0]
                        theta_accretion_end = top_column.outer_surface.theta_range[-1]
                        # ----------------- углы для нахождения пересечений -------------------------

                        # ------------------ начало заполнения матриц косинусов ---------------------------
                        # for key, surface in surfaces.items():
                        #     surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end, top_column.outer_surface.phi_range,
                        #                                bot_column.outer_surface.phi_range, e_obs)

                        # распараллелил
                        for key, surface in surfaces.items():
                            with mp.Pool(mp.cpu_count()) as pool:
                                result_cos_psi_range = pool.starmap(surface.async_fill_cos_psi_range,
                                                                    zip(range(config.t_max),
                                                                        repeat(theta_accretion_begin),
                                                                        repeat(theta_accretion_end),
                                                                        repeat(top_column.outer_surface.phi_range),
                                                                        repeat(bot_column.outer_surface.phi_range),
                                                                        repeat(e_obs),
                                                                        repeat(config.betta_mu)))

                            cos_psi_range_final = []
                            for cos_psi in result_cos_psi_range:
                                cos_psi_range_final.append(cos_psi)
                            surface.cos_psi_range = cos_psi_range_final
                        # ------------------ конец заполнения матриц косинусов ---------------------------

                        # file_name_for_cos_of_surfaces = {0: 'top_outer', 1: 'top_inner',
                        #                                  2: 'bot_outer', 3: 'bot_inner'}
                        #
                        # cos_file_folder = 'data/cos/' + config.file_folder_angle + config.file_folder_args
                        # for key, surface_name in file_name_for_cos_of_surfaces.items():
                        #     full_cos_file_folder = cos_file_folder + file_name_for_cos_of_surfaces[key] + '/'
                        #     for cos_index in range(config.t_max):
                        #         file_name = 'save_cos_' + surface_name + ('_%d_phase' % cos_index) + '.txt'
                        #         main_service.save_arr_as_txt(surfaces[key].cos_psi_range[cos_index],
                        #                                      full_cos_file_folder, file_name)

                        # ------------------ начало заполнения массивов светимости -----------------------
                        arr_simps_integrate = [0] * 4
                        sum_simps_integrate = 0
                        for key, surface in surfaces.items():
                            arr_simps_integrate[key] = surface.calculate_integral_distribution()
                            # file_name = "%s %s %d.txt" % (file_name_variables, approx_method, key)
                            # np.savetxt(file_name, arr_simps_integrate[key])
                            sum_simps_integrate += np.array(arr_simps_integrate[key])
                        # ------------------ конец заполнения массивов светимости -----------------------

                        file_name = "i=%d.txt" % obs_i_angle[i]
                        main_service.save_arr_as_txt(sum_simps_integrate, working_folder, file_name)

                    file_name = "L_x.txt"
                    f = open(working_folder + file_name, 'w')
                    f.write('%f' % top_column.outer_surface.L_x)
                    f.close()

                    time_calculate_nu_L_nu_on_energy = time.time()

                    print("execution time of program: %f s" % (
                            time_calculate_nu_L_nu_on_energy - time_start))

                    plot_sky_map.plot_save_sky_map(obs_i_angle)
                    plot_sky_map.plot_save_sky_map_contour(obs_i_angle)

                    if (a_portion[j] == 1):
                        break
