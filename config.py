from numpy import pi
import numpy as np

# Parameters

# глобальные постоянные
MSun = 1.9891e33  # масса молнца [г]
G = 6.67e-8  # гравитационная постоянная [см3·с−2·г−1]
c = 2.99792458e10  # скорость света [см/с]
sigmStfBolc = 5.67e-5  # постоянная Стефана Больцмана в сгс
a_rad_const = 7.5657e-15  # радиационная константа p=aT**4 [эрг см-3 К-4]
sigmaT = 6.652e-25  # сечение томсона [см-2]
massP = 1.67e-24  # масса протона [г]
h_plank_ergs = 6.62607015e-27  # постоянная Планка в [эрг * с]
h_plank_evs = 4.135667669e-15  # постоянная Планка в [эв * с]
k_bolc = 1.380649e-16  # постоянная Больцмана [эрг/К]

grad_to_rad = pi / 180

# параметры НЗ
M_ns = 1.4 * MSun  # масса нз [г]
R_ns = 1e6  # радиус нз [см]
# H = 2 * 10 ** 13  # магнитное поле стр 19 над формулой 37
mu = 0.1e30  # магнитный момент [Гаусс * см3]
H = 2 * mu / R_ns ** 3
# p_spin = 3.62  # период вращения, [с]

# параметры аккреционного потока
dRe_div_Re = 0.25  # взял просто число
# M_accretion_rate = 10 ** 38 * R_ns / G / MSun  # темп аккреции
ksi_rad = 3 / 2
a_portion = 0.65  # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
k = 0.35  # opacity непрозрачность
# L_ed = M_ns / MSun * 10 ** 38
L_edd = 4 * pi * G * M_ns * c / k

M_rate_c2_Led = 10
M_accretion_rate = M_rate_c2_Led * L_edd / c ** 2  # таблица 1

ksi_param = 0.5  # между 1 и 2 формулой в статье

lim_phi_accretion = 2 * pi * a_portion  # верхний предел по phi
phi_accretion_begin_deg = 0  # нижний предел по phi
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

obs_i_angle = 0 * grad_to_rad  # угол между нормалью к двойной системе и наблюдателем
obs_phi_angle = 0 * grad_to_rad
e_obs = np.array([np.sin(obs_i_angle) * np.cos(obs_phi_angle),
                  np.sin(obs_i_angle) * np.sin(obs_phi_angle),
                  np.cos(obs_i_angle)])

# угол между осью вращения системы и собственным вращением НЗ
betta_rotate = 30 * grad_to_rad
phi_rotate = 0 * grad_to_rad
# угол между собственным вращением НЗ и магнитной осью
betta_mu = 70 * grad_to_rad
phi_mu_0 = 0 * grad_to_rad

file_folder = 'figs/new_energy/'
file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (M_rate_c2_Led, a_portion, phi_accretion_begin_deg)
full_file_folder = file_folder + file_folder_args

# для pretty графиков - индекс по энергии и сколько subfigures
N_column_plot = 5
energy_indexes = [0, 12, 15, 17, 19]
# меньше на 1 т.к. в диапазоне энергий => меньше на 1 чем точек разбиения
energy_indexes_luminosity = [0, 12, 14, 16, 18]
