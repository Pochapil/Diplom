import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

from traits.api import HasTraits, Range, Button, Instance, on_trait_change
from traitsui.api import View, Item, HGroup, VGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

import config
import main_service
import vectors
import geometricTask.matrix as matrix


def get_projection_on_vector(input_vector, project_vector):
    # cos * project_vector
    return np.dot(input_vector, project_vector) * project_vector / (np.linalg.norm(project_vector)) ** 2


def get_projection_on_surface(input_vector, surface_norm_vector):
    # на нормаль к плоскости и вычитаем из вектора
    projection_on_norm = get_projection_on_vector(input_vector, surface_norm_vector)
    return input_vector - projection_on_norm


def get_roll_angle(first_vector, second_vector, phase):
    # first_angle = np.arctan2(first_vector[1], first_vector[0])
    # second_angle = np.arctan2(second_vector[1], second_vector[0])
    # roll_angle = second_angle - first_angle
    #
    # # print(first_vector)
    # # print(second_vector)
    # # print()
    # print('first angle = %f' % (first_angle / config.grad_to_rad))
    # print(second_angle / config.grad_to_rad)
    # print()
    #
    # if (first_vector[0] < 0 and second_vector[0] < 0 or first_vector[0] > 0 and second_vector[0] > 0) \
    #         and first_vector[1] * second_vector[1] < 0:
    #     roll_angle = np.pi - roll_angle
    #
    # return roll_angle / config.grad_to_rad

    if phase % 180 == 0:
        roll_angle = 0
        if i_angle < betta_mu and phase % 360 == 0:
            roll_angle = 180
    else:
        r_first_vector = np.linalg.norm(first_vector)
        r_second_vector = np.linalg.norm(second_vector)
        if r_first_vector == 0 or r_second_vector == 0:
            roll_angle = 0
        else:
            roll_angle = -np.arccos(np.dot(first_vector, second_vector) / (r_first_vector * r_second_vector))
            roll_angle = roll_angle / config.grad_to_rad

    if betta_mu == 0:
        roll_angle = 0

    if phase > 180 and phase < 360 or phase > 540:
        roll_angle = - roll_angle

    return roll_angle


def calculate_roll_angle(i_angle, betta_mu, e_obs_mu, phase):
    # phase in deg!
    '''Метод view init задает координаты картинной плоскости. Поэтому для того, чтобы перейти из магнитной СК обратно
    в СК, связанную с omega, надо достать угол между проекциями mu и omega на картинную плоскость и повернуть картинную плоскость
    на этот угол. Для этого нахожу проекцию на плоскость векторов mu и omega  и нахожу угол между ними через косинус'''
    view_plane_normal = [e_obs_mu[0, 0], e_obs_mu[0, 1], e_obs_mu[0, 2]]  # x,y,z
    omega_vector = [np.sin(-betta_mu * config.grad_to_rad) * np.cos(0),
                    np.sin(-betta_mu * config.grad_to_rad) * np.sin(0),
                    np.cos(-betta_mu * config.grad_to_rad)]

    view_plane_normal = np.array(view_plane_normal)
    omega_vector = np.array(omega_vector)

    # беру проекцию на картинную плоскость, оттуда достаю угол.
    omega_projection_on_view_plane = get_projection_on_surface(omega_vector, view_plane_normal)
    mu_projection_on_view_plane = get_projection_on_surface(np.array([0, 0, 1]), view_plane_normal)

    roll_angle = get_roll_angle(omega_projection_on_view_plane, mu_projection_on_view_plane, phase)
    # print(roll_angle)
    return roll_angle


def lim_axes(ax, lim_value):
    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])
    ax.set_aspect("equal")


def hide_axes_and_background(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    ax.set_axis_off()


def plot_NS_mayavi(phi_range_column, theta_range_column):
    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / config.N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    mlab.mesh(x, y, z, color=(0, 0, 1))


def get_data_for_NS():
    theta_range = np.arange(0, np.pi, np.pi / config.N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi * 1.01, 2 * np.pi * 1.01 / config.N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    return x, y, z


def get_data_for_accretion_disc_side_surface(accretion_disc_con_val):
    # cylinder
    radius = 4
    z = np.linspace(-accretion_disc_con_val * radius, accretion_disc_con_val * radius, config.N_theta_accretion)
    phi = np.linspace(0, 2 * np.pi, config.N_phi_accretion)
    phi_grid, z = np.meshgrid(phi, z)
    x = radius * np.cos(phi_grid)
    y = radius * np.sin(phi_grid)

    return x, y, z


def get_data_for_accretion_disc(accretion_disc_con_val):
    r, phi = np.mgrid[1:4:100j, 0:2 * np.pi:100j]
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = accretion_disc_con_val * (x ** 2 + y ** 2) ** (1 / 2)
    # z, z1 = np.mgrid[-0.003:0.003:100j, -0.003:0.003:100j]
    return x, y, z


def plot_accr_columns_mayavi(phi_range_column, theta_range_column):
    # рисуем силовые линии

    # верх
    theta_range = theta_range_column
    phi_range = phi_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    mlab.mesh(x, y, z, color=(1, 0, 0))
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    # низ
    mlab.mesh(-x, -y, -z, color=(0, 1, 0))


def get_data_for_accretion_columns(theta_range_column, phi_range_column, fi_0):
    step_phi_accretion = config.lim_phi_accretion / (config.N_phi_accretion - 1)
    phi_range = np.array([fi_0 * config.grad_to_rad + step_phi_accretion * i for i in range(config.N_phi_accretion)])
    theta_range = theta_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    return x, y, z


def get_data_for_magnet_lines_with_mask(theta_range_column, phi_range_column, fi_0):
    '''смог реализовать только с помощью маски'''
    theta_array_end = np.pi / 2 + config.betta_mu
    # ограничиваю колонкой
    theta_array_end = min((np.pi - theta_range_column[-1]), theta_array_end)
    theta_array_begin = theta_range_column[-1]

    step_theta_accretion = (theta_array_end - theta_array_begin) / (config.N_theta_accretion - 1)
    theta_range = np.array([theta_array_begin + step_theta_accretion * j for j in range(config.N_theta_accretion)])

    step_phi_accretion = config.lim_phi_accretion / (config.N_phi_accretion - 1)
    phi_range = np.array([fi_0 * config.grad_to_rad + step_phi_accretion * i for i in range(config.N_phi_accretion)])

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    mask = np.zeros_like(x).astype(bool)
    for i in range(len(phi_range)):
        for j in range(len(theta_range)):
            theta_end = np.pi / 2 - config.betta_mu * np.cos(phi_range[i])
            if theta_range[j] > theta_end:
                # pass
                mask[i][j] = True

    return x, y, z, mask


def get_data_for_magnet_lines(theta_range_column, phi_range_column, fi_0, cut_flag=False):
    '''старый метод - уже неправильный'''
    if cut_flag:
        theta_array_end = abs(np.pi / 2 - config.betta_mu)
    else:
        theta_array_end = np.pi / 2

    theta_array_begin = theta_range_column[-1]

    step_theta_accretion = (theta_array_end - theta_array_begin) / (config.N_theta_accretion - 1)
    theta_range = np.array([theta_array_begin + step_theta_accretion * j for j in range(config.N_theta_accretion)])

    step_phi_accretion = config.lim_phi_accretion / (config.N_phi_accretion - 1)
    phi_range = np.array([fi_0 * config.grad_to_rad + step_phi_accretion * i for i in range(config.N_phi_accretion)])

    # верх
    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    return x, y, z

    # ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color="blue", alpha=0.06)
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    # низ
    # ax.plot_wireframe(-x, -y, -z, rstride=4, cstride=4, color="blue", alpha=0.06)


def plot_magnet_lines(ax, phi_range_column):
    step_theta_accretion = (np.pi / 2) / (config.N_theta_accretion - 1)
    theta_range = np.array([step_theta_accretion * j for j in range(config.N_theta_accretion)])

    # верх
    phi_range = phi_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color="blue", alpha=0.06)
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    # низ
    ax.plot_wireframe(-x, -y, -z, rstride=4, cstride=4, color="blue", alpha=0.06)


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


def plot_in_mayavi(i_angle, betta_mu, phi_range_column, theta_range_column):
    step_theta_accretion = (np.pi / 2) / (config.N_theta_accretion - 1)
    theta_range = np.array([step_theta_accretion * j for j in range(config.N_theta_accretion)])

    plot_NS_mayavi(phi_range_column, theta_range_column)
    plot_accr_columns_mayavi(phi_range_column, theta_range_column)

    phase = 0
    A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, config.e_obs)  # переход в магнитную СК

    azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

    observer_mu_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

    # omega_vector = [np.sin(-betta_mu * grad_to_rad) * np.cos(0 * grad_to_rad),
    #                 np.sin(-betta_mu * grad_to_rad) * np.sin(0 * grad_to_rad),
    #                 np.cos(-betta_mu * grad_to_rad)]
    #
    # add_vector(ax, origin, omega_vector, 'blue', lim_value)
    # config.betta_rotate / grad_to_rad + config.betta_mu / grad_to_rad

    # 90 - т.к. находим через arccos (в другой СК - theta от 0Z 0 - 180), а рисовать нужно в СК 90 - -90
    mlab.view(azimuth=azimuth / config.grad_to_rad, elevation=90 - elevation / config.grad_to_rad)

    mlab.show()

    step_phi_accretion = config.lim_phi_accretion / (config.N_phi_accretion - 1)
    phi_range = np.array([config.phi_accretion_begin + step_phi_accretion * i for i in range(config.N_phi_accretion)])

    # fig = mlab.figure(size=(480, 340))
    # mlab.mesh(x, y, z, colormap='cool', figure=fig)


def plot_m(i_angle, betta_mu, phi_range_column, theta_range_column):
    class Visualization(HasTraits):
        slider_fi_0 = Range(0, 360, config.phi_accretion_begin_deg)
        slider_i_angle = Range(0, 90, i_angle)
        slider_betta_mu = Range(0, 90, betta_mu)
        slider_phase = Range(0., 2., 0.)
        slider_distance = Range(0.1, 10, 1)

        button_magnet_line = Button('draw_magnet_lines')
        button_accr_disc = Button('accr_disc_omega_mu')
        button_cut_magnet_lines = Button('cut_magnet_lines')

        scene = Instance(MlabSceneModel, ())

        def __init__(self):
            # Do not forget to call the parent's __init__
            HasTraits.__init__(self)

            # NS
            x, y, z = get_data_for_NS()
            self.NS = self.scene.mlab.mesh(x, y, z, color=(0, 0, 0))

            # рисуем колонки
            x, y, z = get_data_for_accretion_columns(theta_range_column, phi_range_column,
                                                     config.phi_accretion_begin_deg)
            # верх
            self.accretion_column_top = self.scene.mlab.mesh(x, y, z, color=(1, 0, 0))
            # низ
            self.accretion_column_bot = self.scene.mlab.mesh(-x, -y, -z, color=(0, 1, 0))

            self.flag_draw_magnet_lines = True
            self.flag_cut_magnet_lines = False
            # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, config.phi_accretion_begin_deg)
            x, y, z, mask = get_data_for_magnet_lines_with_mask(theta_range_column, phi_range_column,
                                                                config.phi_accretion_begin_deg)

            opacity_for_magnet_line = 0.1
            self.magnet_lines_top = self.scene.mlab.mesh(x, y, z, color=(0, 0, 1), opacity=opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)
            # self.magnet_lines_top = self.scene.mlab.surf(x, y, z, color=(1, 0, 0), warp_scale=0.3,
            #                                              representation='wireframe', line_width=0.5)
            self.magnet_lines_bot = self.scene.mlab.mesh(-x, -y, -z, color=(0, 0, 1), opacity=opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)

            mlab.plot3d([0, 0], [0, 0], [0, 1.0], color=(1, 0, 0), tube_radius=0.005, tube_sides=4)

            omega_vector = [-np.sin(betta_mu * config.grad_to_rad) * np.cos(0),
                            np.sin(betta_mu * config.grad_to_rad) * np.sin(0),
                            np.cos(betta_mu * config.grad_to_rad)]

            self.omega_vector = mlab.plot3d([0, omega_vector[0]], [0, omega_vector[1]], [0, omega_vector[2]],
                                            color=(0, 0, 0), tube_radius=0.005, tube_sides=4)

            # рисуем аккреционный диск
            disc_color = (160, 82, 45)
            disc_color = tuple(rgb / 255 for rgb in disc_color)

            accretion_disc_con_val = 0.01
            x, y, z = get_data_for_accretion_disc(accretion_disc_con_val)

            # сначала отрисовали в СК betta_mu
            self.accretion_disc_top = mlab.mesh(x, y, z, color=disc_color)
            self.accretion_disc_bot = mlab.mesh(x, y, -z, color=disc_color)

            self.flag_accretion_disc_omega_mu = True  # True = omega

            x, y, z = get_data_for_accretion_disc_side_surface(accretion_disc_con_val)

            # боковая поверхность диска (как боковая поверхность цилиндра)
            self.accretion_disc_side_surface = mlab.mesh(x, y, z, color=disc_color)

            # поворачиваем диск на -betta чтобы было перпендикулярно omega -> переходим в СК omega
            self.accretion_disc_rotate_angle = config.betta_mu / config.grad_to_rad
            self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)

        def rotate_accretion_disc(self, rotate_angle):
            self.accretion_disc_top.actor.actor.rotate_y(rotate_angle)
            self.accretion_disc_bot.actor.actor.rotate_y(rotate_angle)
            self.accretion_disc_side_surface.actor.actor.rotate_y(rotate_angle)

        def update_accretion_disc_rotate_angle(self):
            if self.flag_accretion_disc_omega_mu:
                self.rotate_accretion_disc(self.accretion_disc_rotate_angle)
                self.accretion_disc_rotate_angle = config.betta_mu / config.grad_to_rad
                self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)
            else:
                self.accretion_disc_rotate_angle = config.betta_mu / config.grad_to_rad

        def update_magnet_lines(self):
            opacity_for_magnet_line = 0.1
            x, y, z, mask = get_data_for_magnet_lines_with_mask(theta_range_column, phi_range_column,
                                                                self.slider_fi_0)
            # clear
            # self.magnet_lines_top.remove()
            self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0])
            self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0])
            # new
            self.magnet_lines_top = self.scene.mlab.mesh(x, y, z, color=(0, 0, 1), opacity=opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)
            self.magnet_lines_bot = self.scene.mlab.mesh(-x, -y, -z, color=(0, 0, 1), opacity=opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)

        def view_phase(self, phase=0):
            e_obs = config.e_obs
            A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                         config.betta_mu)
            e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

            azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)
            roll_angle = calculate_roll_angle(config.obs_i_angle_deg, config.betta_mu_deg, e_obs_mu, phase)
            # print(roll_angle)

            # ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad, roll=roll_angle)
            distance = self.slider_distance
            self.scene.mlab.view(azimuth=azimuth / config.grad_to_rad, elevation=elevation / config.grad_to_rad,
                                 distance=distance, focalpoint=[0, 0, 0])
            # roll angle считаем для плоскости камеры - поэтому roll там
            # только здесь нашел про камеру https://docs.enthought.com/mayavi/mayavi/mlab_figures_decorations.html
            camera = self.scene.camera
            camera.roll(roll_angle)
            # self.scene.background = (1, 1, 1)  # white background

        @on_trait_change('slider_i_angle')
        def func_change_slider_i_angle_slider(self):
            global betta_mu
            i_angle = self.slider_i_angle

            config.set_e_obs(self.slider_i_angle, 0)

            omega_vector = [-np.sin(config.betta_mu_deg * config.grad_to_rad) * np.cos(0),
                            np.sin(config.betta_mu_deg * config.grad_to_rad) * np.sin(0),
                            np.cos(config.betta_mu_deg * config.grad_to_rad)]

            self.omega_vector.mlab_source.trait_set(x=[0, omega_vector[0]], y=[0, omega_vector[1]],
                                                    z=[0, omega_vector[2]])

            self.update_accretion_disc_rotate_angle()

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('slider_betta_mu')
        def func_change_slider_betta_mu(self):
            global betta_mu
            betta_mu = self.slider_betta_mu

            config.set_betta_mu(self.slider_betta_mu)

            omega_vector = [-np.sin(config.betta_mu_deg * config.grad_to_rad) * np.cos(0),
                            np.sin(config.betta_mu_deg * config.grad_to_rad) * np.sin(0),
                            np.cos(config.betta_mu_deg * config.grad_to_rad)]

            self.omega_vector.mlab_source.trait_set(x=[0, omega_vector[0]], y=[0, omega_vector[1]],
                                                    z=[0, omega_vector[2]])

            self.update_accretion_disc_rotate_angle()

            if self.flag_draw_magnet_lines:
                # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_fi_0,
                #                                     self.flag_cut_magnet_lines)

                # обновляем т.к. зависит от betta_mu
                self.update_magnet_lines()

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('slider_fi_0')
        def func_change_slider_fi_0(self):
            x, y, z = get_data_for_accretion_columns(theta_range_column, phi_range_column, self.slider_fi_0)
            self.accretion_column_top.mlab_source.trait_set(x=x, y=y, z=z, color=(1, 0, 0))
            self.accretion_column_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=(0, 1, 0))

            if self.flag_draw_magnet_lines:
                # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_fi_0,
                #                                     self.flag_cut_magnet_lines)

                self.update_magnet_lines()
                # self.magnet_lines_top.mlab_source.trait_set(x=x, y=y, z=z, mask=mask)
                # self.magnet_lines_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, mask=mask)

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('slider_phase, slider_distance')
        def func_change_slider_phase_slider_distance(self):
            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('button_magnet_line')
        def func_change_button_magnet_line(self):
            self.flag_draw_magnet_lines = not self.flag_draw_magnet_lines
            if self.flag_draw_magnet_lines:
                x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_fi_0,
                                                    self.flag_cut_magnet_lines)
                self.magnet_lines_top.mlab_source.trait_set(x=x, y=y, z=z, color=(0, 0, 1))
                self.magnet_lines_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=(0, 0, 1))
            else:
                self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0], color=(0, 0, 1))
                self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0], color=(0, 0, 1))

        @on_trait_change('button_cut_magnet_lines')
        def func_change_button_cut_magnet_lines(self):
            self.flag_cut_magnet_lines = not self.flag_cut_magnet_lines
            # if self.flag_draw_magnet_lines:
            #     x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_fi_0,
            #                                         self.flag_cut_magnet_lines)
            #     self.magnet_lines_top.mlab_source.trait_set(x=x, y=y, z=z, color=(0, 0, 1))
            #     self.magnet_lines_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=(0, 0, 1))
            # else:
            #     self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0], color=(0, 0, 1))
            #     self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0], color=(0, 0, 1))

            theta = np.pi / 2 - config.betta_mu

            r, phi = np.mgrid[np.sin(theta) ** 2:4:100j, 0:2 * np.pi:100j]
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            z, z1 = np.mgrid[-0.003:0.003:100j, -0.003:0.003:100j]

            # сначала отрисовали в betta_mu
            self.accretion_disc_top.mlab_source.trait_set(x=x, y=y)

        @on_trait_change('button_accr_disc')
        def func_change_button_accr_disc(self):

            if self.flag_accretion_disc_omega_mu:
                self.rotate_accretion_disc(self.accretion_disc_rotate_angle)
            else:
                self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)
            self.flag_accretion_disc_omega_mu = not self.flag_accretion_disc_omega_mu

        # the layout of the dialog created

        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                         height=250, width=300, show_label=False),
                    VGroup(
                        HGroup(
                            '_', 'slider_i_angle', 'slider_betta_mu'
                        ),
                        HGroup(
                            '_', 'slider_fi_0', 'slider_phase'
                        ),
                        HGroup('button_magnet_line', 'button_accr_disc', 'button_cut_magnet_lines', 'slider_distance')
                    )
                    )

    visualization = Visualization()
    visualization.configure_traits()


if __name__ == "__main__":
    gif_flag = False
    config_vectors_flag = False
    phase_flag = False

    plot_magnet_lines_flag = True

    lim_coeff_for_axis = 0.1
    if plot_magnet_lines_flag:
        lim_coeff_for_axis = 0.14

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

    # plot_in_mayavi(i_angle, betta_mu, phi_range_column, theta_range_column)
    plot_m(i_angle, betta_mu, phi_range_column, theta_range_column)
