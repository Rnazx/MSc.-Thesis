import numpy as np

#no RM data for M51
#always correct for distances
#mention inclination angles if not 20 degrees

# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
# *(8.5/7.6) is correction for distance to M51 chosen
G_dat_Bord = np.array([8.6,7.6,7.6,7.8])
G_dat_Breg = np.array([1.3,1.8,2.6,3.2])
G_dat_Btot = np.array([17,16,15,13])

# pitch angle data #RM= rotation measure data, M= fitting polarisation angles

#Surgent+23
#i= 22.5 d=8.58 Mpc
radius_pB_surgent23=np.array([1.4393531,1.913746631,2.377358491,2.851752022,3.326145553,3.811320755,
                              4.285714286,4.749326146,5.223719677,5.687331536,6.161725067,6.64690027,7.121293801,])*(8.5/8.58) #distance correction done
pB_surgent23=np.array([30.75949367,22.78481013,19.36708861,23.92405063,27.72151899,28.86075949,
                       25.44303797,28.10126582,25.06329114,26.20253165,28.86075949,23.16455696,32.27848101])

#Borlaff+23
#i= 22.5 d=8.58 Mpc
radius_pB_borlaff23=np.array([0.4532084,  0.89827749, 1.66880654, 2.3304253,  2.92022675, 3.66114695, 
 4.32883351, 4.91799469 ,5.6388513,  6.34400257, 7.06176041, 7.79785886])*(8.5/8.58) #distance correction done
pB_borlaff23=np.array([80.37383178, 14.95327103, 30.8411215 , 32.71028037 ,28.03738318, 33.64485981,
 45.79439252, 44.85981308, 30.8411215,  30.8411215  ,14.01869159, 33.64485981])

#Borlaff+21
#i= 22
radius_pB_borlaff21=np.array([0.349693006,1.036105179,1.699332621,2.445568607,3.137987186,3.858756006,
4.578323545,5.283355579,6.152602776,6.873011212,7.59197811,8.310704752,9.017058195])*(8.5/8.58) #distance correction done
pB_borlaff21=np.array([70.76775227,36.02349172,24.64388681,22.46129204,23.11478911,30.84570208,31.49706353,26.48585158,
                       29.24933262,34.85638014,31.96796583,27.66364122,30.43993593])

# radius_pB_borlaff23=np.array([  0.89827749, 1.66880654, 2.3304253,  2.92022675, 3.66114695, 
#  4.32883351, 4.91799469 ,5.6388513,  6.34400257, 7.06176041, 7.79785886])
# pB_borlaff23=np.array([ 14.95327103, 30.8411215 , 32.71028037 ,28.03738318, 33.64485981,
#  45.79439252, 44.85981308, 30.8411215,  30.8411215  ,14.01869159, 33.64485981])

RM_dat_po = np.array([20,27,19]) * np.pi/180 #pitch angle of ordered field, not really found using RM
err_RM_dat_po = np.array([2,2,5]) * np.pi/180 #error in po
rmdat_tanpo = np.tan(RM_dat_po)
rm_errdat_tanpo = 1/(np.cos(err_RM_dat_po))**2

M_dat_pb = np.array([20,24,22,18]) * np.pi/180 #pitch angle of regular field
err_M_dat_pb = np.array([1,4,4,1,]) * np.pi/180
m_errdat_tanpb = 1/(np.cos(err_M_dat_pb))**2

#radial ranges
range_po=np.array([1.7,3.0,7.8])*(8.5/7.6) #for po
mrange_endps = np.array([2.4,3.6,4.8,6.0,7.2])*(8.5/7.6) #lower limit of radial ranges where M is used, given in table 4, Beck+19
mrange = (mrange_endps[1:] + mrange_endps[:-1])/2 #average of each of the intervals given in above array. contains 1 less point than above array

#scale height data
kpc_dat_r = np.array([2.96,4.18,5.34,6.56])*(8.5/7.6)
pc_dat_h = np.array([256.54,293.76,336.30,386.81])


