import numpy as np

#no po data for ngc 6946

# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
G_dat_Bord = np.array([5.1,5.1])
G_dat_Breg = np.array([19*0.06,13*0.15])
G_dat_Btot = np.array([19,13])

#pitch angle data #RM= rotation measure data, M= fitting polarisation angles
RM_dat_po = np.array([27,21,10,30,32,10]) * np.pi/180 #pitch angle of ordered field
err_RM_dat_po = np.array([2,2,6,2,4,5]) * np.pi/180 #error in po
rmdat_tanpo = np.tan(RM_dat_po)
rm_errdat_tanpo = 1/(np.cos(err_RM_dat_po))**2

# M_dat_pb = np.array(['nan',20,24,22,18,'nan']) * np.pi/180 #pitch angle of regular field
# err_M_dat_pb = np.array(['nan',1,4,4,1,'nan']) * np.pi/180
#RM_dat_pb = np.array([4, 9, 7, 7, 5]) * np.pi/180 #pitch angle of regular field at another radial range
#err_RM_dat_pb = np.array([5, 3, 3, 2, 3]) * np.pi/180
#rmdat_tanpb = np.tan(RM_dat_pb)
#rm_errdat_tanpb = 1/(np.cos(err_RM_dat_pb))**2
# m_errdat_tanpb = 1/(np.cos(err_M_dat_pb))**2

#radial ranges #discontinuous ranges in Beck+19
# mrange_endps = np.array([]) #lower limit of radial ranges where M is used, given in table 4, Beck+19
# mrange = (mrange_endps[1:] + mrange_endps[:-1])/2 #average of each of the intervals given in above array. contains 1 less point than above array
#rmrange = np.arange(7.5, 12, 1) #radial ranges where RM is used, given in table 4, Beck+19

#scale height data
kpc_dat_r = np.array([7, 9, 11, 13])
pc_dat_h = np.array([316.4, 371.9, 437.1, 513.7])