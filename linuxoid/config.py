i_angle = 10
betta_mu = 40

mc2 = 20
a_portion = 0.65

fi_0 = 20


file_folder = 'figs/loop/'
file_folder_angle_args = 'i=%d betta_mu=%d/' % (i_angle, betta_mu)
file_folder_accretion_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (mc2, a_portion, fi_0)
full_file_folder = file_folder + file_folder_angle_args + file_folder_accretion_args

#print(full_file_folder)

def update():
    global file_folder,file_folder_angle_args,file_folder_accretion_args, full_file_folder
    file_folder = 'figs/loop/'
    file_folder_angle = 'i=%d betta_mu=%d/' % (i_angle, betta_mu)
    file_folder_args = 'mc2=%d/a=%0.2f fi_0=%d/' % (mc2, a_portion, fi_0)
    full_file_folder = file_folder + file_folder_angle + file_folder_args
    #print(full_file_folder)