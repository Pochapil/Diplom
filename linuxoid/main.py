import sys
import numpy as np
import pathlib
import matplotlib.pyplot as plt

import config

def create_file_path(file_path):
    pathlib.Path(file_path).mkdir(parents=True, exist_ok=True)


def save_arr_as_txt(arr, file_folder, file_name):
    create_file_path(file_folder)
    full_file_name = file_folder + file_name
    np.savetxt(full_file_name, arr)


if __name__ == "__main__":

    params_index_dict = {'i_angle': 1, 'betta_mu': 2, 'mc2': 3, 'a_portion': 4, 'fi_0': 5}

    i_angle = int(sys.argv[params_index_dict['i_angle']])
    betta_mu = int(sys.argv[params_index_dict['betta_mu']])

    mc2 = int(sys.argv[params_index_dict['mc2']])
    a_portion = float(sys.argv[params_index_dict['a_portion']])

    fi_0 = int(sys.argv[params_index_dict['fi_0']])

    config.i_angle = i_angle
    config.betta_mu = betta_mu

    config.mc2 = mc2
    config.a_portion = a_portion

    config.fi_0 = fi_0
    config.update()

    file_folder = 'figs/loop/'
    file_folder_angle = 'i=%d betta_mu=%d/' % (i_angle, betta_mu)
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (mc2, a_portion, fi_0)
    full_file_folder = file_folder + file_folder_angle + file_folder_args

    full_file_folder = config.full_file_folder

    arr = [i_angle, betta_mu, mc2, a_portion, fi_0]
    save_arr_as_txt([i_angle, betta_mu, mc2, a_portion, fi_0], full_file_folder, str(arr) + '.txt', )


    # mult = 1
    # for i in range(1, len(sys.argv)):
    #     mult = mult * int(sys.argv[i])
    # print(mult)
    # file_folder = 'try/'
    # save_arr_as_txt([mult], file_folder, str(mult) + '.txt', )

    # for param in sys.argv:
    # print(param)
