import numpy as np

#no RM data for M51

# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
G_dat_Bord = np.array([8.6,7.6,7.6,7.8])
G_dat_Breg = np.array([17*0.08,16*0.11,15*0.17,13*0.25])
G_dat_Btot = np.array([17,16,15,13])

# pitch angle data #RM= rotation measure data, M= fitting polarisation angles
# RM_dat_po = np.array([20,27,'nan','nan','nan',19]) * np.pi/180 #pitch angle of ordered field
# err_RM_dat_po = np.array([2,2,'nan','nan','nan',5]) * np.pi/180 #error in po

RM_dat_po = np.array([20,27,19]) * np.pi/180 #pitch angle of ordered field
err_RM_dat_po = np.array([2,2,5]) * np.pi/180 #error in po
rmdat_tanpo = np.tan(RM_dat_po)
rm_errdat_tanpo = 1/(np.cos(err_RM_dat_po))**2

# M_dat_pb = np.array(['nan',20,24,22,18,'nan']) * np.pi/180 #pitch angle of regular field
# err_M_dat_pb = np.array(['nan',1,4,4,1,'nan']) * np.pi/180

M_dat_pb = np.array([20,24,22,18]) * np.pi/180 #pitch angle of regular field
err_M_dat_pb = np.array([1,4,4,1,]) * np.pi/180
#RM_dat_pb = np.array([4, 9, 7, 7, 5]) * np.pi/180 #pitch angle of regular field at another radial range
#err_RM_dat_pb = np.array([5, 3, 3, 2, 3]) * np.pi/180
#rmdat_tanpb = np.tan(RM_dat_pb)
#rm_errdat_tanpb = 1/(np.cos(err_RM_dat_pb))**2
m_errdat_tanpb = 1/(np.cos(err_M_dat_pb))**2

#radial ranges
#mrange_endps = np.array([1.2,2.4,3.6,4.8,6.0,7.2]) #lower limit of radial ranges where M is used, given in table 4, Beck+19
mrange_endps = np.array([2.4,3.6,4.8,6.0,7.2]) #lower limit of radial ranges where M is used, given in table 4, Beck+19
mrange = (mrange_endps[1:] + mrange_endps[:-1])/2 #average of each of the intervals given in above array. contains 1 less point than above array
#rmrange = np.arange(7.5, 12, 1) #radial ranges where RM is used, given in table 4, Beck+19

#scale height data
kpc_dat_r = np.array([7, 9, 11, 13])
pc_dat_h = np.array([316.4, 371.9, 437.1, 513.7])

