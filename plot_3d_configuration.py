import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox, Button

import config
import main_service
import vectors


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

    lim_value = 0.2
    grad_to_rad = np.pi / 180
    i_angle = 0 * grad_to_rad
    e_obs = np.array([0, np.sin(i_angle), np.cos(i_angle)])
    A_matrix_analytic = matrix.newMatrixAnalytic(0, betta_rotate * grad_to_rad, phase * grad_to_rad,
                                                 betta_mu * grad_to_rad)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

    azimuth, elevation = get_angles_from_vector(e_obs_mu)

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi, 2 * np.pi / N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', alpha=1)

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
    lim_value = 0.2
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi, 2 * np.pi / N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', alpha=1)

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
    e_obs = np.array([0, np.sin(i_angle), np.cos(i_angle)])

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


def visualise_3d_configuration(phi_range_column, theta_range_column):
    # fig, ax = plt.subplots()
    lim_value = 0.2 * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    N_phi_accretion = len(phi_range_column)
    N_theta_accretion = len(theta_range_column)

    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', alpha=1)

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

    axSlider1 = fig.add_axes([0.25, 0.1, 0.65, 0.05])
    slider1 = Slider(axSlider1, 'phase', valmin=0, valmax=2, valinit=0)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    phase = 0
    e_obs = np.array([0, np.sin(config.obs_i_angle), np.cos(config.obs_i_angle)])
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

    def rotate(val):
        phase = slider1.val  # slider1.val
        phase = phase * 360
        e_obs = np.array([0, np.sin(config.obs_i_angle), np.cos(config.obs_i_angle)])
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

    slider1.on_changed(rotate)

    ax_box = fig.add_axes([0.25, 0.03, 0.05, 0.05])
    text_box = TextBox(ax_box, 'phase to set')

    def func_for_text_box(text):
        phase1 = float(text)
        slider1.set_val(phase1)

    text_box.on_submit(func_for_text_box)

    plt.show()


def visualise_3d_configuration_angles():
    # fig, ax = plt.subplots()
    lim_value = 0.2 * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / config.N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', alpha=1)

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
    e_obs = np.array([0, np.sin(config.obs_i_angle), np.cos(config.obs_i_angle)])
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
    lim_value = 0.2 * config.M_rate_c2_Led / 10
    grad_to_rad = np.pi / 180

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])

    def plot_all(phase):
        # рисуем звезду
        theta_range = np.arange(0, np.pi, np.pi / config.N_theta_accretion)
        phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

        u, v = np.meshgrid(phi_range, theta_range)
        r1 = np.sin(theta_range_column[0]) ** 2
        x = r1 * np.sin(v) * np.cos(u)
        y = r1 * np.sin(v) * np.sin(u)
        z = r1 * np.cos(v)

        ax.plot_surface(x, y, z, color='b', alpha=1)

        # вектора
        origin = [0, 0, 0]
        z_vector = [0, 0, 1]
        add_vector(ax, origin, z_vector, 'blue', lim_value)

        omega_vector = [np.sin(config.betta_rotate) * np.cos(0),
                        np.sin(config.betta_rotate) * np.sin(0),
                        np.cos(config.betta_rotate)]
        add_vector(ax, origin, omega_vector, 'black', lim_value)

        omega_vector = [np.sin(np.pi + config.betta_rotate) * np.cos(0),
                        np.sin(np.pi + config.betta_rotate) * np.sin(0),
                        np.cos(np.pi + config.betta_rotate)]
        add_vector(ax, origin, omega_vector, 'black', lim_value)

        # переходим из СК НЗ в OZ так как углы mu в СК НЗ !!!
        A_matrix_analytic = matrix.newMatrixAnalytic(0, -config.betta_rotate, 0, 0)

        mu_vector = [np.sin(config.betta_mu) * np.cos(config.phi_mu_0 + phase),
                     np.sin(config.betta_mu) * np.sin(config.phi_mu_0 + phase),
                     np.cos(config.betta_mu)]

        mu_vector = np.dot(A_matrix_analytic, mu_vector)  # переходим из СК НЗ в OZ так как углы mu в СК НЗ !!!

        azimuth, elevation = vectors.get_angles_from_vector(mu_vector)

        mu_vector = [mu_vector[0, 0], mu_vector[0, 1], mu_vector[0, 2]]
        add_vector(ax, origin, mu_vector, 'red', lim_value)

        e_obs = config.e_obs
        add_vector(ax, origin, e_obs, 'purple', lim_value)

        # рисуем arc
        theta_range = np.arange(0, config.betta_rotate, config.betta_rotate / config.N_theta_accretion)

        x = lim_value * np.sin(theta_range) * 0.9
        y = [0] * config.N_theta_accretion
        z = lim_value * np.cos(theta_range) * 0.9

        ax.plot(x, y, z, color='black', alpha=1, label=r'$ \beta_* $')

        # theta_range = np.arange(-betta_mu * grad_to_rad, -(betta_mu + betta_rotate) * grad_to_rad,
        #                         -(betta_rotate) * grad_to_rad / (config.N_theta_accretion - 1))

        if config.betta_mu + config.betta_rotate < np.pi:
            # elevation - (config.betta_mu)
            theta_range = np.array(
                [config.betta_rotate + (elevation - config.betta_rotate) / (config.N_theta_accretion - 1) * i for i in
                 range(config.N_theta_accretion)])
        # else:
        #     theta_range = np.array(
        #         [elevation + (config.betta_mu) / (config.N_theta_accretion - 1) * i for i in
        #          range(config.N_theta_accretion)])

        phase_range = np.array([0 + (azimuth) / (config.N_theta_accretion - 1) * i for i in
                                range(config.N_theta_accretion)])

        x = lim_value * np.sin(theta_range) * np.cos(phase_range) * 0.8
        y = lim_value * np.sin(theta_range) * np.sin(phase_range) * 0.8
        z = lim_value * np.cos(theta_range) * 0.8

        ax.plot(x, y, z, color='green', alpha=1, label=r'$ \beta_mu $')

        x = lim_value * np.sin(theta_range) * np.cos(phase_range) * 0.7
        y = lim_value * np.sin(theta_range) * np.sin(phase_range) * 0.7
        z = lim_value * np.cos(theta_range) * 0.7

        ax.plot(x, y, z, color='green', alpha=1)

        ax.legend()

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


if __name__ == "__main__":
    working_folder = config.full_file_folder

    file_name = "save_phi_range.txt"
    phi_range_column = main_service.load_arr_from_txt(working_folder, file_name)

    # file_name = "save_phi_range.txt"
    # phi_range_column = np.loadtxt(file_name)

    file_name = "save_theta_range.txt"
    theta_range_column = main_service.load_arr_from_txt(working_folder, file_name)
    # theta_range_column = np.loadtxt(file_name)

    # plot_3d_configuration(phi_range_column, theta_range_column, 40, 60, 0.8)
    visualise_3d_configuration(phi_range_column, theta_range_column)

    visualise_3d_angles()
    # animate_3d_configuration(phi_range_column, theta_range_column, 40, 30)
