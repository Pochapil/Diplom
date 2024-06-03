import numpy as np
import math
import matplotlib.pyplot as plt
from mayavi import mlab

from traits.api import HasTraits, Range, Button, Instance, on_trait_change
from traitsui.api import View, Item, HGroup, VGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
import matplotlib.cm as cm

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


def plot_main(i_angle, betta_mu, fi_0):
    class Visualization(HasTraits):
        slider_fi_0 = Range(0, 360, config.phi_accretion_begin_deg)
        slider_i_angle = Range(0, 90, i_angle)
        slider_betta_mu = Range(0, 90, betta_mu)
        slider_phase = Range(0., 2., 0.)
        slider_distance = Range(0.1, 10, 1)

        button_magnet_line = Button('draw_magnet_lines')
        button_accr_disc = Button('accr_disc_omega_mu')
        button_cut_magnet_lines = Button('cut_magnet_lines')
        button_hide_accr_disc = Button('hide_accr_disc')
        button_check_data = Button('check_data')

        scene = Instance(MlabSceneModel, ())

        color_accretion_column_top = (1, 0, 0)
        color_accretion_column_bot = (0, 1, 0)

        color_accretion_column_top_outer = (1, 1, 0)
        color_accretion_column_bot_outer = (0, 1, 1)

        color_magnet_lines_top = (0, 0, 1)
        color_magnet_lines_top_outer = (0, 0.5, 1)
        color_magnet_lines_bot = (0, 0, 1)
        color_magnet_lines_bot_outer = (0, 0, 1)

        mu_vector_tube_radius = 0.005
        omega_vector_tube_radius = 0.005

        mu_vector_color = (1, 0, 0)
        omega_vector_color = (0, 0, 0)

        def get_data_for_indicatrix(self):

            obs_i_angle_arr = np.linspace(10, 90, 9)

            working_folder = config.full_file_folder

            with open(working_folder + 'save_values.txt') as f:
                lines = f.readlines()
                L_x = float(lines[3][12:20]) * 10 ** float(lines[3][27:29])
                # total L_x = 4.396383 * 10**38 - 12 это индекс начала числа, 27-29 это степень 10

            data_array = [0] * 17

            # roll -- циклическая перестановка - делаем так как симметрическая задача и для угла 180 - theta будет симметрично
            # для сдвига на полфазы -- можно расчитать только до 90 а потом для других переставить и получить для до 180
            for i in range(len(obs_i_angle_arr)):
                config.set_e_obs(obs_i_angle_arr[i], 0)
                file_name = 'total_luminosity_of_surfaces.txt'
                data_array[i] = main_service.load_arr_from_txt(config.full_file_folder, file_name)[4]
                file_name = 'scattered_energy_bot.txt'
                data_array[i] += main_service.load_arr_from_txt(config.full_file_folder + 'scattered_on_magnet_lines/',
                                                                file_name)
                file_name = 'scattered_energy_top.txt'
                data_array[i] += main_service.load_arr_from_txt(config.full_file_folder + 'scattered_on_magnet_lines/',
                                                                file_name)
                if i != 8:
                    data_array[-i - 1] = np.roll(data_array[i], len(data_array[i]) // 2)

            config.set_e_obs(i_angle, 0)

            # нормировка на L_nu_avg
            data_to_plot = []
            for arr in data_array:
                data_to_plot.append(arr / L_x)

            # theta_range = np.linspace(0, np.pi, 17)
            theta_range = np.linspace(10 * config.grad_to_rad, np.pi - 10 * config.grad_to_rad, 17)
            """Учесть что вращение против часовой - поэтому предел - 2пи"""
            phi_range = np.linspace(0, -2 * np.pi, config.t_max)
            # theta_range -= config.betta_mu

            u, v = np.meshgrid(phi_range, theta_range)
            r = np.array(data_to_plot)
            # np.repeat(np.array(data_to_plot)[:, np.newaxis], 45, axis=1)
            x = r * np.sin(v) * np.cos(u)
            y = r * np.sin(v) * np.sin(u)
            z = r * np.cos(v)

            return x, y, z, r

        def __init__(self):
            # Do not forget to call the parent's __init__
            HasTraits.__init__(self)
            self.scene.background = (1, 1, 1)

            # NS
            x, y, z, color = self.get_data_for_indicatrix()
            cmaps = cm.jet(color)
            cmaps[:][:][:3] = cmaps[:][:][:3] * 255
            # print(cmaps[0][1][:3])
            self.indicatrix = self.scene.mlab.mesh(x, y, z, scalars=color)

            e_obs = config.e_obs
            A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, 0 * config.grad_to_rad, 0)
            e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
            '''учесть что данные хранятся в омега СК'''
            azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

            obs_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

            x, y, z = obs_vector

            self.obs_point = self.scene.mlab.points3d(x, y, z,
                                                      mode="point",
                                                      scale_mode="none",
                                                      scale_factor=3.0)
            self.obs_point.actor.property.render_points_as_spheres = True
            self.obs_point.actor.property.point_size = 10

            # self.indicatrix1 = self.scene.mlab.mesh(x, y, color)

            # lut = self.indicatrix1.module_manager.scalar_lut_manager.lut.table.to_array()
            # self.indicatrix.module_manager.scalar_lut_manager.lut.table = lut

            # self.indicatrix1.mlab_source.trait_set(x=[0], y=[0], z=[0])
            # self.indicatrix1 = self.scene.mlab.triangular_mesh(x, y, z)

            # self.indicatrix.module_manager.scalar_lut_manager.lut.table = cmaps.reshape(17 * 45,4)
            # print(lut.shape)

            # colormap=

            scale_factor = 2
            '''учесть что данные хранятся в омега СК'''
            self.omega_vector = mlab.quiver3d(0, 0, 1, mode='2ddash', scale_factor=scale_factor,
                                              color=self.omega_vector_color)
            self.omega_vector_1 = mlab.quiver3d(0, 0, -1, mode='2ddash', scale_factor=scale_factor,
                                                color=self.omega_vector_color)

            mu_vector = [np.sin(betta_mu * config.grad_to_rad) * np.cos(0),
                         np.sin(betta_mu * config.grad_to_rad) * np.sin(0),
                         np.cos(betta_mu * config.grad_to_rad)]

            self.mu_vector = mlab.quiver3d(mu_vector[0], mu_vector[1], mu_vector[2], mode='2ddash',
                                              scale_factor=scale_factor, color=self.mu_vector_color)
            self.mu_vector_1 = mlab.quiver3d(-mu_vector[0], -mu_vector[1], -mu_vector[2], mode='2ddash',
                                                scale_factor=scale_factor, color=self.mu_vector_color)

            obs_vector = [
                np.sin(-betta_mu * config.grad_to_rad + config.obs_i_angle_deg * config.grad_to_rad) * np.cos(0),
                np.sin(-betta_mu * config.grad_to_rad + config.obs_i_angle_deg * config.grad_to_rad) * np.sin(0),
                np.cos(-betta_mu * config.grad_to_rad + config.obs_i_angle_deg * config.grad_to_rad)]

            self.obs_vector = mlab.quiver3d(obs_vector[0], obs_vector[1], obs_vector[2], mode='2ddash',
                                            scale_factor=scale_factor, color=(0, 0, 1))

        def view_phase(self, phase=0):
            e_obs = config.e_obs
            A_matrix_analytic = matrix.newMatrixAnalytic(0, 0, phase * config.grad_to_rad, 0)
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
            # camera = self.scene.camera
            # camera.roll(roll_angle)

        @on_trait_change('slider_i_angle')
        def func_change_slider_i_angle_slider(self):
            global betta_mu
            i_angle = self.slider_i_angle

            config.set_e_obs(self.slider_i_angle, 0)

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('slider_phase, slider_distance')
        def func_change_slider_phase_slider_distance(self):
            phase = 360 * self.slider_phase

            e_obs = config.e_obs
            A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, phase * config.grad_to_rad,
                                                         0)
            e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

            azimuth, elevation = vectors.get_angles_from_vector(e_obs_mu)

            obs_vector = [np.sin(elevation) * np.cos(azimuth),
                          np.sin(elevation) * np.sin(azimuth),
                          np.cos(elevation)]

            self.obs_vector.mlab_source.vectors = np.reshape([obs_vector[0], obs_vector[1], obs_vector[2]], (1, 3))

            for i in range(len(obs_vector)):
                obs_vector[i] *= 2

            x, y, z = obs_vector
            self.obs_point.mlab_source.trait_set(x=x, y=y, z=z)

            self.view_phase(phase)

        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                         height=250, width=300, show_label=False),
                    VGroup(
                        HGroup(
                            '_', 'slider_i_angle', 'slider_betta_mu'
                        ),
                        HGroup(
                            '_', 'slider_fi_0', 'slider_phase'
                        ),
                        HGroup('button_magnet_line', 'button_accr_disc', 'button_cut_magnet_lines',
                               'button_hide_accr_disc', 'slider_distance')
                    )
                    )

    visualization = Visualization()
    visualization.configure_traits()


betta_mu = 40
i_angle = 40

mc2 = 100
a_portion = 0.22
fi_0 = 50

config.set_e_obs(i_angle, 0)
config.set_betta_mu(betta_mu)

config.M_rate_c2_Led = mc2
config.a_portion = a_portion
config.phi_accretion_begin_deg = fi_0

config.update()

plot_main(i_angle, betta_mu, fi_0)
