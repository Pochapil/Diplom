import numpy as np
import time

import main_service
import config
import accretionColumnService
from accretionColumn import AccretionColumn

# ------------------------------------------
time_start = time.time()

i_angle = 30
betta_mu = 0

mc2 = 10
a_portion = 0.65
fi_0 = 0

config.set_e_obs(i_angle, 0)
config.set_betta_mu(betta_mu)

config.M_rate_c2_Led = mc2
config.a_portion = a_portion
config.phi_accretion_begin_deg = fi_0

config.update()
# ------------------------------------------
R_alfven = (config.mu ** 2 / (2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
R_e_outer_surface, R_e_inner_surface = R_e, R_e  # допущение что толщина = 0

e_obs = config.e_obs

theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(R_e_outer_surface)
theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(R_e_inner_surface)

top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface,
                             R_e_inner_surface,
                             theta_accretion_begin_inner_surface, True)

theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface,
                             R_e_inner_surface,
                             theta_accretion_begin_inner_surface, False)

surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface, 2: bot_column.outer_surface,
            3: bot_column.inner_surface}

file_name_for_cos_of_surfaces = {0: 'top_outer', 1: 'top_inner', 2: 'bot_outer', 3: 'bot_inner'}
cos_file_folder = 'figs/cos/' + config.file_folder_angle + config.file_folder_args
for key, surface_name in file_name_for_cos_of_surfaces.items():
    cos_psi_range_final = []
    full_cos_file_folder = cos_file_folder + file_name_for_cos_of_surfaces[key] + '/'
    for cos_index in range(config.t_max):
        file_name = 'save_cos_' + surface_name + ('_%d_phase' % cos_index) + '.txt'
        cos_psi_range_on_phase = main_service.load_arr_from_txt(full_cos_file_folder, file_name)
        cos_psi_range_final.append(cos_psi_range_on_phase)
    surfaces[key].cos_psi_range = cos_psi_range_final

working_folder = config.full_file_folder
file_name = 'surfaces_T_eff.txt'

# arr_T_eff = main_service.load_arr_from_txt(working_folder, file_name)
# for i in range(4):
#     surfaces[i].T_eff = arr_T_eff[i]
#     print(surfaces[i].T_eff)

time_calculate = time.time()
print("execution time of program: %f s" % (time_calculate - time_start))
