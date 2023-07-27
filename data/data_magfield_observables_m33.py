import numpy as np

# no RM data for m33

# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
G_dat_Bord = np.array([3.1,3.1])
G_dat_Breg = np.array([1.3,2.4])
G_dat_Btot = np.array([8.7,7.6])

#pitch angle data #RM= rotation measure data, M= fitting polarisation angles
RM_dat_po = np.array([48,40,41,35]) * np.pi/180 #pitch angle of ordered field
err_RM_dat_po = np.array([5,5,5,6]) * np.pi/180 #error in po
rmdat_tanpo = np.tan(RM_dat_po)
rm_errdat_tanpo = 1/(np.cos(err_RM_dat_po))**2

M_dat_pb = np.array([51,41]) * np.pi/180 #pitch angle of regular field
err_M_dat_pb = np.array([2,2]) * np.pi/180
# RM_dat_pb = np.array([4, 9, 7, 7, 5]) * np.pi/180 #pitch angle of regular field at another radial range
# err_RM_dat_pb = np.array([5, 3, 3, 2, 3]) * np.pi/180
# rmdat_tanpb = np.tan(RM_dat_pb)
rmdat_tanpb = np.tan(M_dat_pb) #need to change name from rm to m
# rm_errdat_tanpb = 1/(np.cos(err_RM_dat_pb))**2
rm_errdat_tanpb = 1/(np.cos(err_M_dat_pb))**2 #need to change name from rm to m

#radial ranges
range_po=np.array([1.0,3.0,5.0,7.0,9.0]) #for po
po_mrange = (range_po[1:] + range_po[:-1])/2 #average of each of the intervals given in above array. contains 1 less point than above array

mrange_endps = np.array([1.0,3.0,5.0]) #radial ranges (for B and pB) where M is used, given in table 4, Beck+19
mrange = (mrange_endps[1:] + mrange_endps[:-1])/2 #average of each of the intervals given in above array. contains 1 less point than above array
# rmrange = np.arange(7.5, 12, 1) #radial ranges where RM is used, given in table 4, Beck+19

#scale height data
kpc_dat_r = np.array([2.02,3.99])
pc_dat_h = np.array([268.42,391.90])