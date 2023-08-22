N_cores = 16

i_angle = [10 * i for i in range(10)]
betta_mu = [10 * i for i in range(10)]

mc2 = [10, 14, 20, 27, 38, 50, 74, 100, 143, 200]
# mc2 = [10, 20, 50, 100, 200]
a_portion = [0.1, 0.25, 0.65, 1]

fi_0 = [20 * i for i in range(18)]

file_name = 'run_params.sh'
f = open(file_name, 'w')

f.write('#!/bin/bash\n')
# f.write('source /home/pochapil/PycharmProjects/Diplom/venv/bin/activate\n')
f.write('\n')

count_cores = 0
for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        for mc2_index in range(len(mc2)):
            for a_portion_index in range(len(a_portion)):
                for fi_0_index in range(len(fi_0)):
                    f.write('python3 main_with_args.py ')
                    f.write('%d %d %d %0.2f %d &\n' % (
                        i_angle[i_angle_index], betta_mu[betta_mu_index], mc2[mc2_index], a_portion[a_portion_index],
                        fi_0[fi_0_index]))
                    count_cores += 1

                    if count_cores >= N_cores:
                        count_cores = 0
                        f.write('wait\n')
                        f.write('\n')

                    if (a_portion[a_portion_index] == 1):
                        break
                # f.write('\n')
# f.write('deactivate')
