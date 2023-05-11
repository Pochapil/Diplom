import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox, Button
from matplotlib.animation import PillowWriter
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

import config
import main_service
import vectors

lim_coeff_for_axis = 0.15


def plot_NS(ax):
    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / config.N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', alpha=1)


def get_angles_from_vector(vector):
    x = vector[0, 0]
    y = vector[0, 1]
    z = vector[0, 2]
    theta = np.arccos(z)  # np.arccos(z/r)
    if x > 0:
        if y >= 0:
            phi = np.arctan(y / x)
        else:
            phi = np.arctan(y / x) + 2 * np.pi
    else:
        phi = np.arctan(y / x) + np.pi
    return phi, theta


def add_vector(ax, origin, vector, color, lim_value):
    # if sum(vector) != 0:
    # vector = np.array(vector) * lim_value / np.linalg.norm(vector)
    ax.quiver(origin[0], origin[1], origin[2], vector[0], vector[1], vector[2], length=lim_value, color=color)


def plot_3d_configuration(phi_range_column, theta_range_column, betta_rotate, betta_mu, phase):
    phase = phase * 360

    lim_value = lim_coeff_for_axis
    grad_to_rad = np.pi / 180
    i_angle = 0 * grad_to_rad
    e_obs = config.e_obs
    A_matrix_analytic = matrix.newMatrixAnalytic(0, betta_rotate * grad_to_rad, phase * grad_to_rad,
                                                 betta_mu * grad_to_rad)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = get_angles_from_vector(e_obs_mu)

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    plot_NS(ax)

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
    ax.plot_wireframe(-x, -y, -z, rstride=4, cstride=4, color="r", alpha=0.2)

    # вектора
    origin = [0, 0, 0]

    # откладываем от магнитного вектора
    omega_vector = [np.sin(betta_mu * grad_to_rad) * np.cos(phase * grad_to_rad),
                    np.sin(betta_mu * grad_to_rad) * np.sin(phase * grad_to_rad),
                    np.cos(betta_mu * grad_to_rad)]
    mu_vector = [0, 0, 1]
    observer_vector = [np.sin((betta_rotate + betta_mu) * grad_to_rad) * np.cos(phase * grad_to_rad),
                       np.sin((betta_rotate + betta_mu) * grad_to_rad) * np.sin(phase * grad_to_rad),
                       np.cos((betta_rotate + betta_mu) * grad_to_rad)]

    observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

    add_vector(ax, origin, observer_mu_vector, 'black', lim_value)
    # add_vector(ax, origin, omega_vector, 'green')
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    # ax.view_init(90 - betta_rotate - betta_mu, 0)  # поворот в градусах
    # ax.view_init(90 - betta_rotate - betta_mu, phase)
    ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

    print('elevation = %f' % elevation)
    print('azimuth = %f' % azimuth)

    plt.show()


def animate_3d_configuration(phi_range_column, theta_range_column, betta_rotate, betta_mu):
    lim_value = lim_coeff_for_axis
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    plot_NS(ax)

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
    ax.plot_wireframe(-x, -y, -z, rstride=4, cstride=4, color="r", alpha=0.2)

    # вектора
    origin = [0, 0, 0]

    # откладываем от магнитного вектора
    mu_vector = [0, 0, 1]

    # add_vector(ax, origin, omega_vector, 'green')
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    # ax.view_init(90 - betta_rotate - betta_mu, 0)  # поворот в градусах
    # ax.view_init(90 - betta_rotate - betta_mu, phase)

    i_angle = 0 * grad_to_rad
    e_obs = config.e_obs

    for phase in range(0, 720):

        A_matrix_analytic = matrix.newMatrixAnalytic(0, betta_rotate * grad_to_rad, phase * grad_to_rad,
                                                     betta_mu * grad_to_rad)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
        azimuth, elevation = get_angles_from_vector(e_obs_mu)

        if phase == 0:
            observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
                                  np.sin(elevation) * np.sin(azimuth),
                                  np.cos(elevation)]

            add_vector(ax, origin, observer_mu_vector, 'black', lim_value)

        ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

        # print('elevation = %f' % elevation)
        # print('azimuth = %f' % azimuth)

        # ax.text(9, 0, 0, "red", color='red')
        plt.draw()
        plt.pause(.001)


def create_gif(phi_range_column, theta_range_column):
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    plot_NS(ax)

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

    # вектора
    origin = [0, 0, 0]
    mu_vector = [0, 0, 1]
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    omega_vector = [-np.sin(betta_mu * grad_to_rad) * np.cos(0),
                    np.sin(betta_mu * grad_to_rad) * np.sin(0),
                    np.cos(betta_mu * grad_to_rad)]

    add_vector(ax, origin, omega_vector, 'black', lim_value)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    phase = 0
    e_obs = config.e_obs
    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * grad_to_rad,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    figure_title = ''
    fig.suptitle(figure_title, fontsize=14)

    def animate(i):
        phase = i / 60 * 360
        e_obs = config.e_obs
        A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * grad_to_rad,
                                                     config.betta_mu)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

        azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)
        # if val == 0.5:
        #     observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
        #                           np.sin(elevation) * np.sin(azimuth),
        #                           np.cos(elevation)]
        #
        #     add_vector(ax, origin, observer_mu_vector, 'black', lim_value)

        # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
        ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

        figure_title = r'$\Phi$' + ' = %.2f' % (phase / 360)
        fig.suptitle(figure_title, fontsize=14)

        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    ani = animation.FuncAnimation(fig, animate, frames=60, interval=10)

    file_folder = 'figs/gifs/'
    file_name = 'i=%d betta_mu=%d a=%0.2f m=%d fi0=%d ani.gif' % (
        config.obs_i_angle_deg, config.betta_mu_deg, config.a_portion, config.M_rate_c2_Led,
        config.phi_accretion_begin_deg)

    ani.save(file_folder + file_name, writer='pillow', fps=50, dpi=180)

    # plt.show()


def visualise_3d_configuration(phi_range_column, theta_range_column):
    # fig, ax = plt.subplots()
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    plot_NS(ax)

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

    # вектора
    origin = [0, 0, 0]
    mu_vector = [0, 0, 1]
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(0),
                    np.sin(-betta_mu * grad_to_rad) * np.sin(0),
                    np.cos(-betta_mu * grad_to_rad)]

    add_vector(ax, origin, omega_vector, 'black', lim_value)

    axSlider1 = fig.add_axes([0.25, 0.1, 0.65, 0.05])
    slider1 = Slider(axSlider1, 'phase', valmin=0, valmax=2, valinit=0)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    ax.set_aspect("equal")

    phase = 0
    e_obs = config.e_obs
    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * grad_to_rad,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

    add_vector(ax, origin, observer_mu_vector, 'purple', lim_value)

    # omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(0 * grad_to_rad),
    #                 np.sin(-betta_mu * grad_to_rad) * np.sin(0 * grad_to_rad),
    #                 np.cos(-betta_mu * grad_to_rad)]
    #
    # add_vector(ax, origin, omega_vector, 'blue', lim_value)
    # config.betta_rotate / grad_to_rad + config.betta_mu / grad_to_rad

    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
    ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    def rotate(val):
        phase = slider1.val  # slider1.val
        phase = phase * 360
        e_obs = config.e_obs
        A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * grad_to_rad,
                                                     config.betta_mu)

        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

        azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

        i_vector_in_mu_sk = np.array([np.sin(-config.betta_mu_deg) * np.cos(0 * grad_to_rad),
                                      np.sin(-config.betta_mu_deg) * np.sin(0 * grad_to_rad),
                                      np.cos(-config.betta_mu_deg)])

        A_inverse = np.linalg.inv(A_matrix_analytic)
        omega = np.dot(A_inverse, i_vector_in_mu_sk)  # переход обратно в omega СК

        A_inverse_matrix = matrix.newRy(-config.betta_mu)

        psi = np.arctan(A_matrix_analytic[2, 1] * A_matrix_analytic[2, 2])
        psi_1 = np.arctan(A_matrix_analytic[1, 0] * A_matrix_analytic[0, 0])
        psi_2 = np.arctan(-A_matrix_analytic[2, 0] * A_matrix_analytic[0, 0] / np.cos(psi_1))

        omega = np.dot(A_matrix_analytic, [0, 0, 1])
        omega_vector = [omega[0, 0], omega[0, 1], omega[0, 2]]
        azimuth_1, elevation_1 = vectors.get_angles_from_vector(omega)

        roll_angle = np.arctan(omega[0, 1] / omega[0, 2]) / grad_to_rad

        # беру проекцию на картинную плоскость, оттуда достаю угол.
        view_plane_normal = [e_obs_mu[0, 0], e_obs_mu[0, 1], e_obs_mu[0, 2]]  # x,y,z
        omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(0),
                        np.sin(-betta_mu * grad_to_rad) * np.sin(0),
                        np.cos(-betta_mu * grad_to_rad)]

        view_plane_normal = np.array(view_plane_normal)
        omega_vector = np.array(omega_vector)

        # не делю на норму т.к. нормированные векторы np.linalg.norm()
        omega_projection_on_norm = np.dot(view_plane_normal, omega_vector) * view_plane_normal / (
            np.linalg.norm(view_plane_normal)) ** 2
        omega_projection_on_view_plane = omega_vector - omega_projection_on_norm

        # print(np.linalg.norm(omega_projection_on_norm))

        mu_projection_on_norm = np.dot(view_plane_normal, [0, 0, 1]) * view_plane_normal / (
            np.linalg.norm(view_plane_normal)) ** 2
        mu_projection_on_view_plane = np.array([0, 0, 1]) - mu_projection_on_norm

        # print(np.linalg.norm(mu_projection_on_view_plane))

        # print(np.dot(omega_projection_on_view_plane, mu_projection_on_view_plane))
        if phase % 180 == 0:
            roll_angle = 0
            if i_angle < betta_mu and phase % 360== 0:
                roll_angle = 180
        else:
            roll_angle = -np.arccos(np.dot(omega_projection_on_view_plane, mu_projection_on_view_plane) / (
                    np.linalg.norm(omega_projection_on_view_plane) * np.linalg.norm(
                mu_projection_on_view_plane))) / grad_to_rad

        if phase > 180 and phase < 360 or phase > 540:
            roll_angle = - roll_angle

        x1_vec = [np.sin(90 * grad_to_rad) * np.cos(0),
                  np.sin(90 * grad_to_rad) * np.sin(0),
                  np.cos(90 * grad_to_rad)]
        x2_vec = [0, 0, 1]
        x1_vec_projection_on_norm = np.dot(view_plane_normal, x1_vec) * view_plane_normal / (
            np.linalg.norm(view_plane_normal)) ** 2
        x1_vec_projection_on_view_plane = np.array(x1_vec) - x1_vec_projection_on_norm

        x2_vec_projection_on_norm = np.dot(view_plane_normal, x2_vec) * view_plane_normal / (
            np.linalg.norm(view_plane_normal)) ** 2
        x2_vec_projection_on_view_plane = np.array(x2_vec) - x2_vec_projection_on_norm

        angle = np.arccos(np.dot(x1_vec_projection_on_view_plane, x2_vec_projection_on_view_plane) / (
                np.linalg.norm(x1_vec_projection_on_view_plane) * np.linalg.norm(
            x2_vec_projection_on_view_plane))) / grad_to_rad

        # print(angle)

        # roll_angle = np.arccos(np.dot(np.array([0,0,1]), mu_projection_on_view_plane)/np.linalg.norm(mu_projection_on_view_plane)) / grad_to_rad
        # print(roll_angle)

        azimuth_1, elevation_1 = vectors.get_angles_from_vector_one_dimension(omega_projection_on_view_plane)

        # roll_angle = config.betta_mu_deg - config.obs_i_angle_deg

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
        ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad, roll=roll_angle)
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
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    plot_NS(ax)

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

    # вектора
    origin = [0, 0, 0]
    mu_vector = [0, 0, 1]
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    # omega_vector = [-np.sin(betta_mu * grad_to_rad) * np.cos(0),
    #                 np.sin(betta_mu * grad_to_rad) * np.sin(0),
    #                 np.cos(betta_mu * grad_to_rad)]
    #
    # add_vector(ax, origin, omega_vector, 'black', lim_value)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

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

    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
    ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    # plt.show()
    return fig


def visualise_3d_configuration_angles():
    # fig, ax = plt.subplots()
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    # рисуем звезду
    plot_NS(ax)

    # вектора
    origin = [0, 0, 0]
    mu_vector = [0, 0, 1]
    add_vector(ax, origin, mu_vector, 'red', lim_value)

    omega_vector = [np.sin(-config.betta_mu) * np.cos(0),
                    np.sin(-config.betta_mu) * np.sin(0),
                    np.cos(-config.betta_mu)]
    add_vector(ax, origin, omega_vector, 'black', lim_value)
    omega_vector = [np.sin(np.pi - config.betta_mu) * np.cos(0),
                    np.sin(np.pi - config.betta_mu) * np.sin(0),
                    np.cos(np.pi - config.betta_mu)]
    add_vector(ax, origin, omega_vector, 'black', lim_value)

    # if not (config.betta_rotate + config.betta_mu < np.pi or config.betta_rotate + config.betta_mu > 2 * np.pi):
    #     omega_vector = [np.sin(np.pi - betta_mu * grad_to_rad) * np.cos(0),
    #                     np.sin(np.pi - betta_mu * grad_to_rad) * np.sin(0),
    #                     np.cos(np.pi - betta_mu * grad_to_rad)]
    # add_vector(ax, origin, omega_vector, 'black', lim_value)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    phase = 0
    e_obs = config.e_obs
    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * grad_to_rad,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

    add_vector(ax, origin, observer_mu_vector, 'purple', lim_value)

    # рисуем arc
    theta_range = np.arange(0, -config.betta_mu, -config.betta_mu / config.N_theta_accretion)

    y = [0] * config.N_theta_accretion
    x = lim_value * np.sin(theta_range) * 0.9
    z = lim_value * np.cos(theta_range) * 0.9

    ax.plot(x, y, z, color='black', alpha=1, label=r'$ \beta_mu $')

    # theta_range = np.arange(-betta_mu * grad_to_rad, -(betta_mu + betta_rotate) * grad_to_rad,
    #                         -(betta_rotate) * grad_to_rad / (config.N_theta_accretion - 1))

    theta_range = np.array(
        [-config.betta_mu + -(config.betta_rotate) / (config.N_theta_accretion - 1) * i for i in
         range(config.N_theta_accretion)])

    y = [0] * config.N_theta_accretion
    x = lim_value * np.sin(theta_range) * 0.8
    z = lim_value * np.cos(theta_range) * 0.8

    ax.plot(x, y, z, color='green', alpha=1, label=r'$ \beta_* $')

    y = [0] * config.N_theta_accretion
    x = lim_value * np.sin(theta_range) * 0.7
    z = lim_value * np.cos(theta_range) * 0.7

    ax.plot(x, y, z, color='green', alpha=1)
    # omega_vector = [-np.sin(betta_mu * grad_to_rad) * np.cos(0 * grad_to_rad),
    #                 np.sin(betta_mu * grad_to_rad) * np.sin(0 * grad_to_rad),
    #                 np.cos(betta_mu * grad_to_rad)]
    #
    # add_vector(ax, origin, omega_vector, 'blue', lim_value)
    # config.betta_rotate / grad_to_rad + config.betta_mu / grad_to_rad

    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
    ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)
    ax.legend()
    plt.show()


def visualise_3d_angles():
    # fig, ax = plt.subplots()
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    def add_text_to_vector(ax, vector, label):
        vector = np.array(vector) / np.linalg.norm(vector) * lim_value
        ax.text(vector[0] * 1.1, vector[1] * 1.1, vector[2] * 1.1, s=label, fontsize=14)

    def plot_all(phase):
        # рисуем звезду
        plot_NS(ax)

        # вектора
        origin = [0, 0, 0]
        z_vector = [0, 0, 1]
        add_vector(ax, origin, z_vector, 'blue', lim_value)

        omega_vector = [np.sin(config.betta_rotate) * np.cos(0),
                        np.sin(config.betta_rotate) * np.sin(0),
                        np.cos(config.betta_rotate)]
        add_vector(ax, origin, omega_vector, 'black', lim_value)

        add_text_to_vector(ax, omega_vector, r'$ \omega $')

        omega_vector = [np.sin(np.pi + config.betta_rotate) * np.cos(0),
                        np.sin(np.pi + config.betta_rotate) * np.sin(0),
                        np.cos(np.pi + config.betta_rotate)]
        # add_vector(ax, origin, omega_vector, 'black', lim_value)
        ax.quiver(omega_vector[0] * lim_value, omega_vector[1] * lim_value, omega_vector[2] * lim_value,
                  -omega_vector[0], -omega_vector[1], -omega_vector[2], length=lim_value, color='black')

        # переходим из СК НЗ в OZ так как углы mu в СК НЗ !!!
        A_matrix_analytic = matrix.newMatrixAnalytic(0, -config.betta_rotate, 0, 0)

        mu_vector = [np.sin(config.betta_mu) * np.cos(config.phi_mu_0 + phase),
                     np.sin(config.betta_mu) * np.sin(config.phi_mu_0 + phase),
                     np.cos(config.betta_mu)]

        mu_vector = np.dot(A_matrix_analytic, mu_vector)  # переходим из СК НЗ в OZ так как углы mu в СК НЗ !!!

        azimuth, elevation = vectors.get_angles_from_vector(mu_vector)

        mu_vector = [mu_vector[0, 0], mu_vector[0, 1], mu_vector[0, 2]]
        add_vector(ax, origin, mu_vector, 'red', lim_value)

        ax.quiver(-mu_vector[0] * lim_value, -mu_vector[1] * lim_value, -mu_vector[2] * lim_value,
                  mu_vector[0], mu_vector[1], mu_vector[2], length=lim_value, color='red')

        add_text_to_vector(ax, mu_vector, r'$ \mu $')

        e_obs = config.e_obs
        add_vector(ax, origin, e_obs, 'purple', lim_value)

        # add_text_to_vector(ax, e_obs, 'observer')

        theta_range = np.linspace(-2 * np.pi / 7, 2 * np.pi / 7, config.N_theta_accretion)

        # x = lim_value * np.sin(theta_range) * 0.1
        # y = [0] * config.N_theta_accretion
        # z = -lim_value * np.cos(theta_range) * 0.1 + lim_value * 1.5
        #
        # ax.plot(x, y, z, color='black')
        #
        # x = [x[0], 0, x[-1]]
        # y = [0] * 3
        # z = [z[15], lim_value * 1.7, z[-16]]
        #
        # ax.plot(x, y, z, color='black')

        # рисуем arc
        # theta_range = np.arange(0, config.betta_rotate, config.betta_rotate / config.N_theta_accretion)
        # theta_range = np.arange(0, config.obs_i_angle, config.obs_i_angle / config.N_theta_accretion)

        x = lim_value * np.sin(theta_range) * 0.9
        y = [0] * config.N_theta_accretion
        z = lim_value * np.cos(theta_range) * 0.9

        ax.plot(x, y, z, color='black', alpha=1, label=r'$ i_{obs} $')
        direction = (1, 1, 1)
        ax.text(x[config.N_theta_accretion // 2], y[config.N_theta_accretion // 2],
                z[config.N_theta_accretion // 2] * 1.1, s=r'$ i_{obs} $', fontsize=14)

        # theta_range = np.arange(-betta_mu * grad_to_rad, -(betta_mu + betta_rotate) * grad_to_rad,
        #                         -(betta_rotate) * grad_to_rad / (config.N_theta_accretion - 1))

        # if config.betta_mu + config.betta_rotate < np.pi:
        # elevation - (config.betta_mu)
        # elevation - theta mu в СК OZ, betta_rotate - theta om в СК OZ, рисуем от om до mu в 0Z поэтому:
        theta_range = np.array(
            [config.betta_rotate + (elevation - config.betta_rotate) / (config.N_theta_accretion - 1) * i for i in
             range(config.N_theta_accretion)])
        # else:
        # theta_range = np.array(
        #         [elevation + (config.betta_mu) / (config.N_theta_accretion - 1) * i for i in
        #          range(config.N_theta_accretion)])

        phase_range = np.array([0 + (azimuth) / (config.N_theta_accretion - 1) * i for i in
                                range(config.N_theta_accretion)])
        # if azimuth > np.pi:
        #     phase_range = np.array([azimuth + (2 * np.pi - azimuth) / (config.N_theta_accretion - 1) * i for i in
        #                             range(config.N_theta_accretion)])
        if config.betta_mu + config.betta_rotate > np.pi:
            phase_range = np.array([0 + (azimuth - np.pi) / (config.N_theta_accretion - 1) * i for i in
                                    range(config.N_theta_accretion)])

        x = lim_value * np.sin(theta_range) * np.cos(phase_range) * 0.7
        y = lim_value * np.sin(theta_range) * np.sin(phase_range) * 0.7
        z = lim_value * np.cos(theta_range) * 0.7

        ax.plot(x, y, z, color='green', alpha=1)

        x = lim_value * np.sin(theta_range) * np.cos(phase_range) * 0.8
        y = lim_value * np.sin(theta_range) * np.sin(phase_range) * 0.8
        z = lim_value * np.cos(theta_range) * 0.8

        ax.plot(x, y, z, color='green', alpha=1, label=r'$ \beta_\mu $')

        ax.text(x[config.N_theta_accretion // 2] * 1.2, y[config.N_theta_accretion // 2] * 1.2,
                z[config.N_theta_accretion // 2] * 1.2, s=r'$ \beta_\mu $', color='green', fontsize=14)

        # ax.legend()

        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90

    plot_all(0)

    axSlider1 = fig.add_axes([0.25, 0.1, 0.65, 0.05])
    slider1 = Slider(axSlider1, 'phase', valmin=0, valmax=2, valinit=0)

    def rotate(val):
        ax.cla()

        ax.set_xlim([-lim_value, lim_value])
        ax.set_ylim([-lim_value, lim_value])
        ax.set_zlim([-lim_value, lim_value])

        phase = slider1.val * 2 * np.pi  # slider1.val

        plot_all(phase)

        # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
        # ax.view_init(90 - elevation / grad_to_rad, azimuth / grad_to_rad)

    slider1.on_changed(rotate)

    ax_box = fig.add_axes([0.25, 0.03, 0.05, 0.05])
    text_box = TextBox(ax_box, 'phase to set')

    def func_for_text_box(text):
        phase1 = float(text)
        slider1.set_val(phase1)

    text_box.on_submit(func_for_text_box)

    plt.show()


def plot_sphere():
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    theta_range = np.arange(0, 2 * np.pi, 2 * np.pi / (2 * config.N_theta_accretion))
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(theta_range) * np.cos(0)
    y = r1 * np.sin(theta_range) * np.sin(0)
    z = r1 * np.cos(theta_range)

    ax.plot(x, y, z, color='b', alpha=1)

    theta_range = np.arange(0, np.pi, np.pi / (2 * config.N_theta_accretion))
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(theta_range) * 1
    y = r1 * np.sin(theta_range) * 0 + 0.1
    z = r1 * np.cos(theta_range)

    ax.plot(x, y, z, color='b', alpha=1)

    plt.show()


def visualise_3d_star(phi_range_column, theta_range_column):
    # fig, ax = plt.subplots()
    lim_value = lim_coeff_for_axis * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    plot_NS(ax)

    plt.show()


if __name__ == "__main__":

    gif_flag = False
    config_vectors_flag = False
    phase_flag = False

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
        create_gif(phi_range_column, theta_range_column)

    # visualise_3d_star(phi_range_column, theta_range_column)
    # plot_3d_configuration(phi_range_column, theta_range_column, 40, 60, 0.8)
    visualise_3d_configuration(phi_range_column, theta_range_column)

    file_folder_angle = 'i=%d betta_mu=%d/' % (config.obs_i_angle_deg, config.betta_mu_deg)
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (
        config.M_rate_c2_Led, config.a_portion, config.phi_accretion_begin_deg)
    save_folder = 'figs/phases/' + file_folder_angle + file_folder_args

    if phase_flag:
        phase = 0.5
        fig = visualise_3d_configuration_on_phase(phi_range_column, theta_range_column, phase)
        file_name = "phase = %0.2f.png" % phase
        main_service.save_figure(fig, save_folder, file_name)

    if config_vectors_flag:
        visualise_3d_angles()
    # animate_3d_configuration(phi_range_column, theta_range_column, 40, 30)

    # plot_sphere()
