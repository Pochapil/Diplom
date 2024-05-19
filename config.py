from numpy import pi
import numpy as np
import os


def set_config_params(mu_arg, opacity_above_shock_arg, i_angle_arg, betta_mu_arg, M_rate_c2_Led_arg, a_portion_arg,
                      fi_0_arg):
    global mu
    mu = mu_arg

    global opacity_above_shock
    opacity_above_shock = opacity_above_shock_arg

    global betta_mu, betta_mu_deg
    betta_mu_deg = betta_mu_arg
    betta_mu = betta_mu_deg * grad_to_rad

    global obs_i_angle_deg, obs_i_angle, obs_phi_angle, e_obs
    obs_i_angle_deg = i_angle_arg
    obs_i_angle = obs_i_angle_deg * grad_to_rad
    e_obs = np.array([np.sin(obs_i_angle) * np.cos(obs_phi_angle),
                      np.sin(obs_i_angle) * np.sin(obs_phi_angle),
                      np.cos(obs_i_angle)])

    global M_accretion_rate, M_rate_c2_Led
    M_rate_c2_Led = M_rate_c2_Led_arg
    M_accretion_rate = M_rate_c2_Led * L_edd / c ** 2  # таблица 1

    global lim_phi_accretion, a_portion
    a_portion = a_portion_arg
    lim_phi_accretion = 2 * pi * a_portion

    global phi_accretion_begin, phi_accretion_begin_deg
    phi_accretion_begin_deg = fi_0_arg
    phi_accretion_begin = phi_accretion_begin_deg * grad_to_rad

    update_project_dir()
    update_folder()
    # print(full_file_folder)


def update_folder():
    global file_folder, file_folder_accretion_args, file_folder_angle_args, full_file_folder
    file_folder_angle_args = 'i=%d betta_mu=%d/' % (obs_i_angle_deg, betta_mu_deg)
    file_folder_accretion_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg)
    full_file_folder = PROJECT_DIR + file_folder + file_folder_angle_args + file_folder_accretion_args


def update_project_dir():
    global PROJECT_DIR
    global PROJECT_DIR_ORIGIN

    if new_magnet_lines_flag:
        PROJECT_DIR += 'new_magnet_lines/'
    else:
        PROJECT_DIR += 'old_magnet_lines/'

    PROJECT_DIR += 'new_condition/'
    PROJECT_DIR_ORIGIN = PROJECT_DIR

    buf = mu
    count = 1
    while buf > 1:
        count += 1
        buf //= 10
    # PROJECT_DIR += 'mu=%d/opacity/%0.2f/' % (count, opacity_above_shock)
    if tau_flag:
        if tau_cutoff > 0:
            PROJECT_DIR += 'new_data/' + f'mu=0.1e{count}/' + 'tau_with_cutoff/'
        else:
            PROJECT_DIR += 'new_data/' + f'mu=0.1e{count}/' + 'tau/'
    else:
        PROJECT_DIR += 'new_data/' + f'mu=0.1e{count}/opacity={opacity_above_shock:.2f}/'

    if not NS_shadow_flag:
        PROJECT_DIR += 'NS_shadow_off/'


def update():
    global M_accretion_rate
    M_accretion_rate = M_rate_c2_Led * L_edd / c ** 2  # таблица 1

    global lim_phi_accretion
    lim_phi_accretion = 2 * pi * a_portion

    global phi_accretion_begin
    phi_accretion_begin = phi_accretion_begin_deg * grad_to_rad

    global file_folder, file_folder_accretion_args, file_folder_angle_args, full_file_folder
    file_folder_angle_args = 'i=%d betta_mu=%d/' % (obs_i_angle_deg, betta_mu_deg)
    file_folder_accretion_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg)
    full_file_folder = PROJECT_DIR + file_folder + file_folder_angle_args + file_folder_accretion_args


def get_folder_with_args(i_angle, betta_mu, M_rate_c2_Led, a_portion, fi_0):
    '''получаю имя папки откуда достать данные'''
    file_folder = 'figs/loop/'
    file_folder_angle_args = 'i=%d betta_mu=%d/' % (i_angle, betta_mu)
    file_folder_accretion_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, fi_0)
    full_file_folder = file_folder + file_folder_angle_args + file_folder_accretion_args

    return PROJECT_DIR + full_file_folder


def get_e_obs(i_angle, phi_angle):
    e_obs = np.array([np.sin(i_angle) * np.cos(phi_angle),
                      np.sin(i_angle) * np.sin(phi_angle),
                      np.cos(i_angle)])
    return e_obs


def set_e_obs(i_angle, phi_angle):
    global obs_i_angle_deg, obs_i_angle, obs_phi_angle, e_obs
    obs_i_angle_deg = i_angle
    obs_i_angle = obs_i_angle_deg * grad_to_rad
    obs_phi_angle = phi_angle * grad_to_rad
    e_obs = np.array([np.sin(obs_i_angle) * np.cos(obs_phi_angle),
                      np.sin(obs_i_angle) * np.sin(obs_phi_angle),
                      np.cos(obs_i_angle)])

    update()


def set_betta_mu(betta_mu_deg_arg):
    global betta_mu, betta_mu_deg
    betta_mu_deg = betta_mu_deg_arg
    betta_mu = betta_mu_deg * grad_to_rad

    update()


# Parameters

# глобальные постоянные
M_Sun = 1.9891e33  # масса молнца [г]
G = 6.67e-8  # гравитационная постоянная [см3·с−2·г−1]
c = 2.99792458e10  # скорость света [см/с]
sigm_Stf_Bolc = 5.67e-5  # постоянная Стефана Больцмана в сгс
a_rad_const = 7.5657e-15  # радиационная константа p=aT**4 [эрг см-3 К-4] постоянная излучения - Плотность энергии и давление равновесного излучения
sigma_T = 6.652e-25  # сечение томсона [см-2]
mass_P = 1.67e-24  # масса протона [г]
h_plank_ergs = 6.62607015e-27  # постоянная Планка в [эрг * с]
h_plank_evs = 4.135667669e-15  # постоянная Планка в [эв * с]
k_bolc = 1.380649e-16  # постоянная Больцмана [эрг/К]

grad_to_rad = pi / 180

# параметры НЗ
M_ns = 1.4 * M_Sun  # масса нз [г]
R_ns = 1e6  # радиус нз [см]
# H = 2 * 10 ** 13  # магнитное поле стр 19 над формулой 37
mu = 0.1e31  # магнитный момент [Гаусс * см3]
H = 2 * mu / R_ns ** 3
# p_spin = 3.62  # период вращения, [с]

# параметры аккреционного потока
dRe_div_Re = 0.25  # взял просто число
# M_accretion_rate = 10 ** 38 * R_ns / G / MSun  # темп аккреции
ksi_rad = 3 / 2
ksi_param = 0.5  # между 1 и 2 формулой в статье - размер магнитосферы
k = 0.35  # opacity непрозрачность [см**2 / г]
new_magnet_lines_flag = True
tau_flag = True
tau_cutoff = 0
opacity_above_shock = 0  # непрозрачность вещества над ударной волной: 0 - полностью прозрачное, 1 - непрозрачное
# L_ed = M_ns / MSun * 10 ** 38
L_edd = 4 * pi * G * M_ns * c / k

M_rate_c2_Led = 30
M_accretion_rate = M_rate_c2_Led * L_edd / c ** 2  # таблица 1

a_portion = 0.65  # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ

lim_phi_accretion = 2 * pi * a_portion  # верхний предел по phi
phi_accretion_begin_deg = 0  # нижний предел по phi - phi_0 !!!!
phi_accretion_begin = phi_accretion_begin_deg * grad_to_rad  # нижний предел по phi

# цикл для поворотов, сколько точек на графике интегралов - для фазы от 0 до 2 - с перекрытием чтобы форму макс
omega_ns = 8  # скорость вращения НЗ - будет меняться только угол phi_mu!
max_phase_angle_for_plot = 720  # сколько точек на графике интегралов - для фазы от 0 до 2 - с перекрытием чтобы форму макс
t_max_for_plot = (max_phase_angle_for_plot // omega_ns) + (1 if max_phase_angle_for_plot % omega_ns > 0 else 0)
max_phase = 360
t_max = (max_phase // omega_ns) + (1 if max_phase % omega_ns > 0 else 0)
# цикл для поворотов, сколько точек для фазы от 0 до 1 (полного поворота)

# количество шагов
N_phi_accretion = 100
N_theta_accretion = 100
N_wavelength_range = 10
N_frequency_range = 100

N_energy = 20
energy_min = 1  # [КэВ]
energy_max = 40  # [КэВ]

# лог сетка по энергии: e_n+1 = e_n * const
energy_step = (energy_max / energy_min) ** (1 / (N_energy - 1))
energy_arr = list(energy_min * energy_step ** i for i in range(N_energy - 1))
# чтобы убрать погрешности и закрыть массив точным числом
energy_arr.append(energy_max)

obs_i_angle_deg = 30
obs_i_angle = obs_i_angle_deg * grad_to_rad  # угол между нормалью к двойной системе и наблюдателем
obs_phi_angle = 0 * grad_to_rad
e_obs = np.array([np.sin(obs_i_angle) * np.cos(obs_phi_angle),
                  np.sin(obs_i_angle) * np.sin(obs_phi_angle),
                  np.cos(obs_i_angle)])

# угол между осью вращения системы и собственным вращением НЗ (берем ось z сонаправленно с осью собств. вращения omega)
betta_rotate = 0 * grad_to_rad
phi_rotate = 0 * grad_to_rad
# угол между собственным вращением НЗ и магнитной осью
betta_mu_deg = 30
betta_mu = betta_mu_deg * grad_to_rad
phi_mu_0 = 0 * grad_to_rad

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = PROJECT_DIR.replace('\\', '/') + '/'
PROJECT_DIR_ORIGIN = PROJECT_DIR
NS_shadow_flag = True
update_project_dir()
# print(PROJECT_DIR)

file_folder = 'figs/loop/'
file_folder_angle_args = 'i=%d betta_mu=%d/' % (obs_i_angle_deg, betta_mu_deg)
file_folder_accretion_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg)
full_file_folder = file_folder + file_folder_angle_args + file_folder_accretion_args

# print(PROJECT_DIR + full_file_folder)

# для pretty графиков - индекс по энергии и сколько subfigures
N_column_plot = 5
energy_indexes = [0, 12, 15, 17, 19]
# меньше на 1 т.к. в диапазоне энергий => меньше на 1 чем точек разбиения
energy_indexes_luminosity = [0, 12, 14, 16, 18]

fi_0_dict = {0.11: 340, 0.165: 330, 0.22: 320, 0.275: 310, 0.33: 300, 0.385: 290, 0.44: 280, 0.5: 270, 0.55: 260,
             0.605: 250, 0.66: 240, 0.715: 230, 0.77: 220, 0.825: 210, 0.25: 320, 0.65: 240}

# val_ksi = 100
# val_n = 1e15
# # val = 2 ** (1 / 2) * (G * M_ns) ** (3 / 2) / (val_ksi * R_ns ** (7 / 2) * M_accretion_rate)
# val = G * M_ns * M_accretion_rate / ((val_ksi * R_ns) ** 3 * val_n)
# print(val)
# print(H/10**11)
