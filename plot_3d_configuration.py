import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt


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
    vector = np.array(vector) * lim_value / np.linalg.norm(vector)
    ax.quiver(origin[0], origin[1], origin[2], vector[0], vector[1], vector[2], color=color)


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


if __name__ == "__main__":
    file_name = "save_phi_range.txt"
    phi_range_column = np.loadtxt(file_name)

    file_name = "save_theta_range.txt"
    theta_range_column = np.loadtxt(file_name)

    # plot_3d_configuration(phi_range_column, theta_range_column, 40, 60, 0.8)

    animate_3d_configuration(phi_range_column, theta_range_column, 40, 30)
