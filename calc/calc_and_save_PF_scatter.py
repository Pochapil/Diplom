import geometricTask.matrix as matrix
import config
import accretionColumnService
from accretionColumn import AccretionColumn
import vectors
import main_service

i_angle = [10 * i for i in range(1, 10)]
betta_mu = [10 * i for i in range(1, 10)]
mc2 = [30, 100]
a_portion = [0.25, 0.65]
fi_0 = [20 * i for i in range(18)]

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        for i in range(len(mc2)):
            for j in range(len(a_portion)):
                for k in range(len(fi_0)):
                    # fi_0[k] = 360 - 180 * a_portion[j]

                    print('calculate for i=%d beta_mu=%d mc=%0.2f a=%0.2f fi_0=%0.2f' % (
                        i_angle[i_angle_index], betta_mu[betta_mu_index], mc2[i], a_portion[j], fi_0[k]))

                    config.set_e_obs(i_angle[i_angle_index], 0)
                    config.set_betta_mu(betta_mu[betta_mu_index])

                    config.M_rate_c2_Led = mc2[i]
                    config.a_portion = a_portion[j]
                    config.phi_accretion_begin_deg = fi_0[k]

                    config.update()

                    folder = 'L_nu/'
                    working_folder = config.full_file_folder + folder

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

                    file_name = "PF.txt"
                    main_service.save_arr_as_txt(PF, working_folder, file_name)
