import numpy as np

# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
G_dat_Bord = np.array([4.9, 5.2, 4.9, 4.6])
G_dat_Breg = np.array([1.8, 2.1, 2.6, 2.7])
G_dat_Btot = np.array([7.3, 7.5, 7.1, 6.3])

#pitch angle data #RM= rotation measure data, M= fitting polarisation angles
RM_dat_po = np.array([30, 29, 26, 27, 27]) * np.pi/180 #pitch angle of ordered field
err_RM_dat_po = np.array([5, 4, 3, 2, 3]) * np.pi/180 #error in po


M_dat_pb = np.array([13, 19, 11, 8]) * np.pi/180 #pitch angle of regular field
err_M_dat_pb = np.array([4, 3, 3, 3]) * np.pi/180
RM_dat_pb = np.array([4, 9, 7, 7, 5]) * np.pi/180 #pitch angle of regular field at another radial range
err_RM_dat_pb = np.array([5, 3, 3, 2, 3]) * np.pi/180

#radial ranges
mrange_endps = np.array([6.8, 9.0, 11.3, 13.6, 15.8]) #lower limit of radial ranges where M is used, given in table 4, Beck+19
mrange = (mrange_endps[1:] + mrange_endps[:-1])/2 #average of each of the intervals given in above array. contains 1 less point than above array
rmrange_endps = np.array([7, 8, 9, 10, 11,12]) #radial ranges where RM is used, given in table 4, Beck+19
rmrange = (rmrange_endps[1:] + rmrange_endps[:-1])/2 #radial ranges where RM is used, given in table 4, Beck+19

#scale height data
kpc_dat_r = np.array([7, 9, 11, 13])
pc_dat_h = np.array([316.4, 371.9, 437.1, 513.7])