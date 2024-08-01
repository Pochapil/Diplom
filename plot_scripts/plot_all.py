import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import config

from plot_scripts import plot_from_main, plot_luminosity_in_range, plot_nu_L_nu, plot_L_nu, plot_scattered
plt.style.use(['science', 'notebook', 'grid'])  # для красивых графиков

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'


if __name__ == '__main__':
    def plot_all():
        plot_from_main.plot_figs()
        plot_luminosity_in_range.plot_figs()
        plot_L_nu.plot_figs()
        plot_nu_L_nu.plot_figs()
        plot_scattered.plot_figs()
        plt.close('all')  # чтобы очистить память от fig из предыдущих вызовов


    fi_0 = [0]
    i_angle = [10 * i for i in range(2, 10)]
    betta_mu = [10 * i for i in range(1, 10)]
    #  mc2 = [10, 30, 50, 100, 150]
    mc2 = [30, 100]
    a_portion = [0.22, 0.44, 0.66]
    # a_portion = [0.11, 0.33, 0.44]
    # fi_0 = [20 * i for i in range(10, 18)]

    i_angle = [60]
    betta_mu = [40]
    mc2 = [160]
    a_portion = [0.22, 0.44, 0.66]
    # fi_0 = [0]

    # i_angle = [10 * i for i in range(10, 19)]
    #fi_0 = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # a_portion = [0.22, 0.44, 0.66]
    N_big = len(i_angle) * len(betta_mu) * len(mc2) * len(a_portion) * len(fi_0)
    print('to calculate %d loops need about %f hours' % (N_big, 40 * N_big / 3600))

    for i_angle_index in range(len(i_angle)):
        for betta_mu_index in range(len(betta_mu)):
            for i in range(len(mc2)):
                for j in range(len(a_portion)):
                    for k in range(len(fi_0)):
                        # fi_0[k] = 360 - 180 * a_portion[j]

                        fi_0[k] = (config.fi_0_dict[a_portion[j]] + 0) % 360
                        # fi_0[k] = config.fi_0_dict[a_portion[j]]
                        # fi_0[k] = (config.fi_0_dict[a_portion[j]] + 20 * k) % 360

                        config.set_e_obs(i_angle[i_angle_index], 0)
                        config.set_betta_mu(betta_mu[betta_mu_index])

                        config.M_rate_c2_Led = mc2[i]
                        config.a_portion = a_portion[j]
                        config.phi_accretion_begin_deg = fi_0[k]

                        config.update()

                        print('calculate for i=%d beta_mu=%d mc=%0.2f a=%0.2f fi_0=%0.2f' % (
                            i_angle[i_angle_index], betta_mu[betta_mu_index], mc2[i], a_portion[j], fi_0[k]))

                        time_start = time.time()
                        plot_all()
                        time_to_plot = time.time()
                        print("execution time of program: %f s" % (time_to_plot - time_start))
