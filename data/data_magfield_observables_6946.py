import numpy as np

#no pb data for ngc 6946

# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
G_dat_Bord = np.array([5.1,5.1])
G_dat_Breg = np.array([1.2,1.9])
G_dat_Btot = np.array([19,13])

#pitch angle data #RM= rotation measure data, M= fitting polarisation angles
RM_dat_po_range1 = np.array([27,21,10]) * np.pi/180 #pitch angle of ordered field
RM_dat_po_range2 = np.array([30,32,10]) * np.pi/180 
err_RM_dat_po_range1 = np.array([2,2,6]) * np.pi/180 #error in po
err_RM_dat_po_range2 = np.array([2,4,5]) * np.pi/180 
rmdat_tanpo_range1 = np.tan(RM_dat_po_range1)
rmdat_tanpo_range2 = np.tan(RM_dat_po_range2)
rm_errdat_tanpo_range1 = 1/(np.cos(err_RM_dat_po_range1))**2
rm_errdat_tanpo_range2 = 1/(np.cos(err_RM_dat_po_range2))**2

# M_dat_pb = np.array(['nan',20,24,22,18,'nan']) * np.pi/180 #pitch angle of regular field
# err_M_dat_pb = np.array(['nan',1,4,4,1,'nan']) * np.pi/180
#RM_dat_pb = np.array([4, 9, 7, 7, 5]) * np.pi/180 #pitch angle of regular field at another radial range
#err_RM_dat_pb = np.array([5, 3, 3, 2, 3]) * np.pi/180
#rmdat_tanpb = np.tan(RM_dat_pb)
#rm_errdat_tanpb = 1/(np.cos(err_RM_dat_pb))**2
# m_errdat_tanpb = 1/(np.cos(err_M_dat_pb))**2

#radial ranges #discontinuous ranges in Beck+19
range1_endpoints=np.array([0,6,12,18])
range1=((range1_endpoints[1:] + range1_endpoints[:-1])/2)*(7.72/7)
range2=np.array([1.5,5.5,8.5])*(7.72/7)
mrange_endps = np.array([0,4.7,9.4]) #lower limit of radial ranges where M is used, given in table 4, Beck+19
mrange = ((mrange_endps[1:] + mrange_endps[:-1])/2)*(7.72/7) #average of each of the intervals given in above array. contains 1 less point than above array



#scale height data
kpc_dat_r = np.array([3.01,9.02,13.02])*(7.72/7)
pc_dat_h = np.array([259.2,564.92,923.81])