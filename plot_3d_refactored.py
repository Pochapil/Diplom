import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox, Button
from matplotlib.animation import PillowWriter
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import scienceplots

import config
import main_service
import vectors

lim_coeff_for_axis = 0.15


def get_projection_on_vector(input_vector, project_vector):
    # cos * project_vector
    return np.dot(input_vector, project_vector) * project_vector / (np.linalg.norm(project_vector)) ** 2


def get_projection_on_surface(input_vector, surface_norm_vector):
    # на нормаль к поверхности и вычитаем из вектора
    projection_on_norm = get_projection_on_vector(input_vector, surface_norm_vector)
    return input_vector - projection_on_norm


def get_roll_angle(i_angle, betta_mu, first_vector, second_vector, phase):
    if phase % 180 == 0:
        roll_angle = 0
        if i_angle < betta_mu and phase % 360 == 0:
            roll_angle = 180
    else:
        roll_angle = -np.arccos(np.dot(first_vector, second_vector) / (
                np.linalg.norm(first_vector) * np.linalg.norm(second_vector))) / config.grad_to_rad

    if betta_mu == 0:
        roll_angle = 0

    if phase > 180 and phase < 360 or phase > 540:
        roll_angle = - roll_angle

    return roll_angle


def calculate_roll_angle(i_angle, betta_mu, e_obs_mu, phase):
    # phase in deg!
    '''надо достать угол между проекциями mu и omega на картинную плоскость - чтобы повернуть картинную плоскость
    на этот угол и перевести в СК связанную с omega'''
    view_plane_normal = [e_obs_mu[0, 0], e_obs_mu[0, 1], e_obs_mu[0, 2]]  # x,y,z
    omega_vector = [np.sin(-betta_mu * config.grad_to_rad) * np.cos(0),
                    np.sin(-betta_mu * config.grad_to_rad) * np.sin(0),
                    np.cos(-betta_mu * config.grad_to_rad)]

    view_plane_normal = np.array(view_plane_normal)
    omega_vector = np.array(omega_vector)

    # беру проекцию на картинную плоскость, оттуда достаю угол.
    omega_projection_on_view_plane = get_projection_on_surface(omega_vector, view_plane_normal)
    mu_projection_on_view_plane = get_projection_on_surface(np.array([0, 0, 1]), view_plane_normal)

    roll_angle = get_roll_angle(i_angle, betta_mu, omega_projection_on_view_plane, mu_projection_on_view_plane, phase)

    if i_angle == 0:
        roll_angle = 0

    return roll_angle


def lim_axes(ax, lim_value):
    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])
    ax.set_aspect("equal")


def plot_NS(ax, phi_range_column, theta_range_column):
    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / config.N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', shade=False, alpha=1)


def plot_accr_columns(ax, phi_range_column, theta_range_column):
    # рисуем силовые линии

    # верх
    theta_range = theta_range_column
    phi_range = phi_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color="r", alpha=0.2)
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    # низ
    ax.plot_wireframe(-x, -y, -z, rstride=4, cstride=4, color="green", alpha=0.2)


def plot_mu_omega_vector(ax, betta_mu, lim_value):
    # вектора
    origin = [0, 0, 0]
    mu_vector = [0, 0, 1]
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    omega_vector = [-np.sin(betta_mu * config.grad_to_rad) * np.cos(0),
                    np.sin(betta_mu * config.grad_to_rad) * np.sin(0),
                    np.cos(betta_mu * config.grad_to_rad)]

    add_vector(ax, origin, omega_vector, 'black', lim_value)


def add_vector(ax, origin, vector, color, lim_value):
    # if sum(vector) != 0:
    # vector = np.array(vector) * lim_value / np.linalg.norm(vector)
    ax.quiver(origin[0], origin[1], origin[2], vector[0], vector[1], vector[2], length=lim_value, color=color)


def hide_axes_and_background(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    ax.set_axis_off()

# vectors.get_angles_from_vector(vector)

def create_gif(i_angle, betta_mu, phi_range_column, theta_range_column):
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    # рисуем звезду
    plot_NS(ax, phi_range_column, theta_range_column)

    plot_accr_columns(ax, phi_range_column, theta_range_column)

    # вектора
    plot_mu_omega_vector(ax, betta_mu, lim_value)

    lim_axes(ax, lim_value)

    phase = 0
    # по сути меняю e_obs
    e_obs = config.e_obs

    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad)

    # Hide axes ticks
    hide_axes_and_background(ax)

    figure_title = ''
    fig.suptitle(figure_title, fontsize=14)

    def animate(i):
        phase = i / 60 * 360
        e_obs = config.e_obs
        A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                     config.betta_mu)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

        azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)
        # if val == 0.5:
        #     observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
        #                           np.sin(elevation) * np.sin(azimuth),
        #                           np.cos(elevation)]
        #
        #     add_vector(ax, origin, observer_mu_vector, 'black', lim_value)

        roll_angle = calculate_roll_angle(i_angle, betta_mu, e_obs_mu, phase)
        # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
        ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad, roll=roll_angle)

        figure_title = r'$\Phi$' + ' = %.2f' % (phase / 360)
        fig.suptitle(figure_title, fontsize=14)

        # Hide axes ticks
        hide_axes_and_background(ax)

    ani = animation.FuncAnimation(fig, animate, frames=60, interval=10)

    file_folder = 'figs/gifs/'
    file_name = 'i=%d betta_mu=%d a=%0.2f m=%d fi0=%d ani.gif' % (
        config.obs_i_angle_deg, config.betta_mu_deg, config.a_portion, config.M_rate_c2_Led,
        config.phi_accretion_begin_deg)

    ani.save(file_folder + file_name, writer='pillow', fps=50, dpi=180)


def visualise_3d_configuration(i_angle, betta_mu, phi_range_column, theta_range_column):
    # fig, ax = plt.subplots()
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    plot_NS(ax, phi_range_column, theta_range_column)

    plot_accr_columns(ax, phi_range_column, theta_range_column)

    # вектора
    plot_mu_omega_vector(ax, betta_mu, lim_value)

    axSlider1 = fig.add_axes([0.25, 0.1, 0.65, 0.05])
    slider1 = Slider(axSlider1, 'phase', valmin=0, valmax=2, valinit=0)

    lim_axes(ax, lim_value)

    phase = 0
    e_obs = config.e_obs
    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

    add_vector(ax, [0, 0, 0], observer_mu_vector, 'purple', lim_value)

    # omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(0 * grad_to_rad),
    #                 np.sin(-betta_mu * grad_to_rad) * np.sin(0 * grad_to_rad),
    #                 np.cos(-betta_mu * grad_to_rad)]
    #
    # add_vector(ax, origin, omega_vector, 'blue', lim_value)
    # config.betta_rotate / grad_to_rad + config.betta_mu / grad_to_rad

    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
    ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad)

    # Hide axes ticks
    hide_axes_and_background(ax)

    def rotate(val):
        phase = slider1.val  # slider1.val
        phase = phase * 360
        e_obs = config.e_obs
        A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                     config.betta_mu)

        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

        azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

        # беру проекцию на картинную плоскость, оттуда достаю угол.
        # не делю на норму т.к. нормированные векторы np.linalg.norm()
        roll_angle = calculate_roll_angle(i_angle, betta_mu, e_obs_mu, phase)

        # if val == 0.5:
        # observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
        #                       np.sin(elevation) * np.sin(azimuth),
        #                       np.cos(elevation)]
        # add_vector(ax, origin, observer_mu_vector, 'black', lim_value)
        #
        # omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(phase * grad_to_rad),
        #                 np.sin(-betta_mu * grad_to_rad) * np.sin(phase * grad_to_rad),
        #                 np.cos(-betta_mu * grad_to_rad)]
        # add_vector(ax, origin, omega_vector, 'blue', lim_value)

        # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
        # ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

        # roll angle добавил
        ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad, roll=roll_angle)
        # ax.view_init(0, phase, roll=config.betta_mu_deg)

        # Hide axes ticks
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.set_zticks([])

    slider1.on_changed(rotate)

    ax_box = fig.add_axes([0.25, 0.03, 0.05, 0.05])
    text_box = TextBox(ax_box, 'phase to set')

    def func_for_text_box(text):
        phase1 = float(text)
        slider1.set_val(phase1)

    text_box.on_submit(func_for_text_box)

    plt.show()


def visualise_3d_configuration_on_phase(phi_range_column, theta_range_column, phase):
    # fig, ax = plt.subplots()

    vector_flag = True

    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    plot_NS(ax, phi_range_column, theta_range_column)

    # рисуем силовые линии
    plot_accr_columns(ax, phi_range_column, theta_range_column)

    # вектора
    if vector_flag:
        plot_mu_omega_vector(ax, betta_mu, lim_value)

    lim_axes(ax, lim_value)

    e_obs = config.e_obs

    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * 2 * np.pi, config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    # omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(0 * grad_to_rad),
    #                 np.sin(-betta_mu * grad_to_rad) * np.sin(0 * grad_to_rad),
    #                 np.cos(-betta_mu * grad_to_rad)]
    #
    # add_vector(ax, origin, omega_vector, 'blue', lim_value)
    # config.betta_rotate / grad_to_rad + config.betta_mu / grad_to_rad
    roll_angle = calculate_roll_angle(i_angle, betta_mu, e_obs_mu, phase * 360)
    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
    ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad, roll=roll_angle)

    # Hide axes ticks
    hide_axes_and_background(ax)

    # plt.show()
    return fig


if __name__ == "__main__":
    gif_flag = False
    config_vectors_flag = False
    phase_flag = True

    lim_coeff_for_axis = 0.1

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

    working_folder = config.full_file_folder

    file_name = "save_phi_range.txt"
    phi_range_column = main_service.load_arr_from_txt(working_folder, file_name)

    # file_name = "save_phi_range.txt"
    # phi_range_column = np.loadtxt(file_name)

    file_name = "save_theta_range.txt"
    theta_range_column = main_service.load_arr_from_txt(working_folder, file_name)
    # theta_range_column = np.loadtxt(file_name)

    if gif_flag:
        create_gif(i_angle, betta_mu, phi_range_column, theta_range_column)

    visualise_3d_configuration(i_angle, betta_mu, phi_range_column, theta_range_column)

    file_folder_angle = 'i=%d betta_mu=%d/' % (config.obs_i_angle_deg, config.betta_mu_deg)
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (
        config.M_rate_c2_Led, config.a_portion, config.phi_accretion_begin_deg)
    save_folder = 'figs/phases/' + file_folder_angle + file_folder_args

    if phase_flag:
        phase = 0.75
        fig = visualise_3d_configuration_on_phase(phi_range_column, theta_range_column, phase)
        file_name = "phase = %0.2f.png" % phase
        main_service.save_figure(fig, save_folder, file_name)

    # if config_vectors_flag:
    #     visualise_3d_angles()
