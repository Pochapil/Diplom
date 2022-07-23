import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt


def plot_3d_configuration(phi_range_column, theta_range_column):
    lim_value = 0.2

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

    origin = [0, 0, 0]

    def add_vector(ax, origin, vector):
        vector = np.array(vector) * lim_value / np.linalg.norm(vector)
        ax.quiver(origin[0], origin[1], origin[2], vector[0], vector[1], vector[2])

    omega_vector = [1, 2, 1]
    mu_vector = [0, 0, 1]

    add_vector(ax, origin, omega_vector)
    add_vector(ax, origin, mu_vector)

    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])
    plt.show()


if __name__ == "__main__":
    file_name = "save_phi_range.txt"
    phi_range_column = np.loadtxt(file_name)

    file_name = "save_theta_range.txt"
    theta_range_column = np.loadtxt(file_name)

    plot_3d_configuration(phi_range_column, theta_range_column)
