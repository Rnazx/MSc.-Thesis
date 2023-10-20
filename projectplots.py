import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import unicodedata
import math as m

i_m51=m.radians(20)
i_6946=m.radians(38)
plt.rcParams.update(plt.rcParamsDefault)
# plt.rcParams['text.usetex'] = True
# #######################################################################################################################################

# # all plots for m51

# # rotation curve
# RC1_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\RCm51.csv')
# RC2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\RC_sofue+18.csv')

# r1=np.array(RC1_m51['R (kpc)']*(8.5/9.6))
# v1=np.array(RC1_m51['v (km/s)'])

# r2=np.array(RC2_m51['r']*(8.5/9.6))
# v2=np.array(RC2_m51['v'])

# plt.plot(r2,v2,label='Sofue+18')
# plt.plot(r1,v1,label='Sofue+99')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$V_c$ [$kms^{-1}$]',fontsize=13)

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('rotation curve'))
# plt.show()

# #######################################################################################################################################

# # omega
# Omega = unicodedata.lookup('GREEK CAPITAL LETTER OMEGA')

# omega1_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\omega.csv')
# omega2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\omega_sofue+18.csv')

# print(omega1_m51.head())
# print(omega2_m51.head())

# r1=np.array(omega1_m51['r']*(8.5/9.6))
# omega1=np.array(omega1_m51['omega'])

# r2=np.array(omega2_m51['r']*(8.5/9.6))
# omega2=np.array(omega2_m51['omega'])

# plt.plot(r2,omega2,label='Sofue+18')
# plt.plot(r1,omega1,label='Sofue+99')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Omega$ [$kms^{-1} kpc^{-1}$]',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('omega'))
# plt.show()


# #######################################################################################################################################

# # q
# use files n definitions of omega plotting

# q1=-(r1/omega1)*np.gradient(omega1)/np.gradient(r1)
# q2=-(r2/omega2)*np.gradient(omega2)/np.gradient(r2)

# plt.plot(r2,q2,label='Sofue+18')
# plt.plot(r1,q1,label='Sofue+99')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$q=-dln\Omega/dlnR$',fontsize=13)
# # plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('q'))
# plt.show()

########################################################################################################################################

# # sigma_HI, sigma_H2 and sigma_gas 

# HI_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\HI m51...csv')
# HI_m51_b=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\sigma_HI_bigiel.csv')

# H2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\H2 m51..csv')
# H2_m51_b=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\sigma_H2_bigiel.csv')

# sigma_gas=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\HI+H2 m51..csv')
# sigma_gas_b=HI_m51_b+H2_m51_b

# # print(HI_m51.head())
# # print(H2_m51.head())
# # print(sigma_gas.head())

# r1=np.array(HI_m51['r']*(8.5/8.2))
# sigma1=np.array(HI_m51['sigma_HI'])
# r1_b=np.array(HI_m51_b['r']*(8.5/8))
# sigma1_b=np.array(HI_m51_b['sigma_HI'])

# r2=np.array(H2_m51['r']*(8.5/8.2))
# sigma2=np.array(H2_m51['sigma_H2'])
# r2_b=np.array(H2_m51_b['r']*(8.5/8))
# sigma2_b=np.array(H2_m51_b['sigma_H2'])

# r3=np.array(sigma_gas['r']*(8.5/8.2))
# sigma3=np.array(sigma_gas['sigma_gas'])
# # r3_b=np.array(sigma_gas_b['r']*(8.5/8.2))
# # sigma3_b=np.array(sigma_gas_b['sigma_gas'])

# # plt.plot(r1,sigma1,label='$\Sigma_{HI}$ Kumari+20')
# plt.plot(r2,sigma2,label='$\Sigma_{H_2}$ Kumari+20')
# # plt.plot(r3,sigma3,label='$\Sigma_{gas}$',color='blue')
# # plt.plot(r1_b,sigma1_b,label='$\Sigma_{HI}$ Bigiel+08',marker='*')
# plt.plot(r2_b,sigma2_b,label='$\Sigma_{H_2}$ Bigiel+08',marker='*')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Sigma_{H_2}$ $[M_\odot pc^{-2}]$',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('sigma H2'))

# plt.show()


# #######################################################################################################################################

# # velocity dispersion

# vel_disp_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\Tamburro HI velcity dispersion..csv')
# vel_disp2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\CO vel dispersion schuster..csv')

# r1=np.array(vel_disp_6946['r'])*(7.72/6) #boomsma data
# veldisp1=np.array(vel_disp_6946['v'])

# r2=np.array(vel_disp2_6946['r'])*(7.72/5.5) #tamburro data
# veldisp2=np.array(vel_disp2_6946['v disp'])

# plt.plot(r1,veldisp1,label='Tamburro+09')
# plt.plot(r2,veldisp2,label='Boomsma+08')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\sigma_{\star}$ $[km s^{-1}]$',fontsize=13)
# # plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('vel disp'))

# plt.show()

vel_disp_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\not needed\Tamburro HI velcity dispersion..csv')
vel_disp_m51_schuster=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\CO vel dispersion schuster..csv')

vel_disp_m51_schuster['corrected radius']=vel_disp_m51_schuster['r']*(8.5/8.4)
print(vel_disp_m51_schuster.head())

vel_disp_m51['corr radius']=vel_disp_m51['r']*(8.5/7.77)
r1=np.array(vel_disp_m51['corr radius'])
veldisp=np.array(vel_disp_m51['v'])

r2=np.array(vel_disp_m51_schuster['corrected radius'])
veldisp2=np.array(vel_disp_m51_schuster['vel disp'])

plt.plot(r1,veldisp,label='Tamburro+09')
plt.plot(r2,veldisp2,label='Schuster+07')
plt.title('M51 1D velocity dispersion')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('Velocity (km/s)')

plt.show()

# #######################################################################################################################################

# # sigma sfr

# sigma_sfr_FUV_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\SFR_FUV24 m51..csv')
# sigma_sfr_Halpha_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\SFR_Halpha24 m51..csv')
# sigma_sfr_b=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\sigma_sfr_bigiel.csv')

# print(sigma_sfr_FUV_m51.head())
# print(sigma_sfr_Halpha_m51.head())

# r1=np.array(sigma_sfr_FUV_m51['r']*(8.5/8.2))
# sigmasfr1=np.array(sigma_sfr_FUV_m51['sigma_sfr_from_plot'])

# r2=np.array(sigma_sfr_FUV_m51['r']*(8.5/8.2))
# sigmasfr2=np.array(sigma_sfr_Halpha_m51['sigma_sfr_from_plot'])

# r3=np.array(sigma_sfr_b['r']*(8.5/8))
# sigmasfr3=np.array(sigma_sfr_b['sigma_sfr_from_plot'])

# plt.plot(r3,sigmasfr3,label='Bigiel+08')
# plt.plot(r1,sigmasfr1,label='FUV Kumari+20')
# plt.plot(r2,sigmasfr2,label='$H_\\alpha$ Kumari+20')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Sigma_{SFR} \\times 10^{-3}$ $ [M_\odot kpc^{-2} yr^{-1}]$',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('sigma SFR'))
# plt.show()


# #######################################################################################################################################

# # sigma SMD
# SMDf=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\smdf.csv')
# SMDs=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\smds.csv')

# print(SMDf.head())
# print(SMDs.head())

# r1=np.array(SMDf['r']*(8.5/9.6))
# smdf=np.array(SMDf['smdf'])

# r2=np.array(SMDs['r']*(8.5/9.6))
# smds=np.array(SMDs['smds'])

# plt.plot(r1,smdf,label='smd-f')
# plt.plot(r2,smds,label='smd-s')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Sigma_{\star}$ $[M_\odot pc^{-2} Gyr^{-1}]$',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('smd'))

# plt.show()
#######################################################################################################################################
# from scipy import stats

# # temperature #fitting done in curve_fitting.py
# temp_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\csv\temperature.csv')

# r1=np.array(temp_m51['r'])*(8.5/7.6)
# temp=np.array(temp_m51['temp'])
# error_temp=np.array(temp_m51['quant_error'])
# slope, intercept, r_value, p_value, std_err = stats.linregress(r1, temp)

# temp_corr=np.array([r*slope+intercept for r in r1])
# print(error_temp)
# plt.errorbar(r1,temp, yerr=error_temp,fmt=' ',capsize=5,label='Bresolin+04',marker='o')
# plt.plot(r1,temp_corr,label='$T={}r +{}$'.format(np.round(slope,2),np.round(intercept,2)))

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('T [K]',fontsize=13)

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('T'))

# plt.show()

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


# NGC 6946

#######################################################################################################################################

# rotation curve

# RC1_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\RC.csv') #sparc data with error
# RC2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\RC_sofue+18.csv')
# RC3_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\RC_sofue+99.csv')

# r1=np.array(RC1_6946['# Rad']*(7.72/5.52))
# v1=np.array(RC1_6946['Vobs']*(m.cos(i_6946)/m.cos(m.radians(38))))
# error=np.array(RC1_6946['errV']*(m.cos(i_6946)/m.cos(m.radians(38))))

# r2=np.array(RC2_6946['R']*(7.72/5.5))
# v2=np.array(RC2_6946['v'])

# r3=np.array(RC3_6946['r']*(7.72/5.5))
# v3=np.array(RC3_6946['v'])

# plt.plot(r2,v2,label='Sofue+18')
# plt.plot(r3,v3,label='Sofue+99')
# plt.errorbar(r1, v1, yerr=error,label='SPARC data', fmt='.', capsize=1) #to be used when there are rror bars

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$V_c$ [$kms^{-1}$]',fontsize=13)

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('rotation curve'))

# plt.show()

# #######################################################################################################################################

# # omega
# omega1_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\omega.csv')
# omega2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\omega_sofue+18.csv')
# omega3_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\omega_sofue+99.csv')


# r1=np.array(omega1_6946['# Rad']*(7.72/5.52))
# omega1=np.array(omega1_6946['omega']*(m.cos(i_6946)/m.cos(m.radians(38))))
# error=np.array(omega1_6946['error_omega']*(m.cos(i_6946)/m.cos(m.radians(38))))

# r2=np.array(omega2_6946['r']*(7.72/5.5))
# omega2=np.array(omega2_6946['omega'])

# r3=np.array(omega3_6946['r']*(7.72/5.5))
# omega3=np.array(omega3_6946['omega'])

# plt.plot(r2,omega2,label='Sofue+18')
# plt.plot(r3,omega3,label='Sofue+99')
# plt.errorbar(r1, omega1, yerr=error,label='SPARC data', fmt='.', capsize=2) #to be used when there are error bars

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Omega$ [$kms^{-1} kpc^{-1}$]',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('omega'))
# plt.show()

# #######################################################################################################################################

# # q
# use files n definitions of omega plotting

# q1=-(r1/omega1)*np.gradient(omega1)/np.gradient(r1)
# q2=-(r2/omega2)*np.gradient(omega2)/np.gradient(r2)
# q3=-(r3/omega3)*np.gradient(omega3)/np.gradient(r3)

# plt.plot(r2,q2,label='Sofue+18')
# plt.plot(r1,q1,label='Sofue+99')
# plt.plot(r3,q3,label='SPARC data')

# plt.legend(loc='upper center',bbox_to_anchor=(0.3, 1),fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$q=-dln\Omega/dlnR$',fontsize=13)
# # plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('q'))
# plt.show()
# #######################################################################################################################################

# # sigma_HI, sigma_H2 and sigma_gas 

# HI_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\HI 6946.csv')
# HI_m51_b=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\sigma_HI_bigiel.csv')

# H2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\H2 6946.csv')
# H2_m51_b=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\sigma_H2_bigiel.csv')

# sigma_gas=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\HI+H2 6946.csv')

# r1=np.array(HI_m51['r']*(7.72/6.8))
# sigma1=np.array(HI_m51['sigma_HI_from_plot']*(m.cos(i_6946)/m.cos(m.radians(33))))
# r1_b=np.array(HI_m51_b['r']*(7.72/5.9))
# sigma1_b=np.array(HI_m51_b['sigma_HI']*(m.cos(i_6946)/m.cos(m.radians(33))))

# r2=np.array(H2_m51['r']*(7.72/6.8))
# sigma2=np.array(H2_m51['sigma_H2_from_plot']*(m.cos(i_6946)/m.cos(m.radians(33))))
# r2_b=np.array(H2_m51_b['r']*(7.72/5.9))
# sigma2_b=np.array(H2_m51_b['sigma_H2']*(m.cos(i_6946)/m.cos(m.radians(33))))

# r3=np.array(sigma_gas['r']*(7.72/6.8))
# sigma3=np.array(sigma_gas['sigma_from_plot']*(m.cos(i_6946)/m.cos(m.radians(33))))

# # plt.plot(r1_b,sigma1_b,label='$\Sigma_{HI}$ Bigiel+08',marker='*')
# # plt.plot(r1,sigma1,label='$\Sigma_{HI}$ Kumari+20')

# plt.plot(r2_b,sigma2_b,label='$\Sigma_{H_2}$ Bigiel+08',marker='*')
# plt.plot(r2,sigma2,label='$\Sigma_{H_2}$ Kumari+20')

# # plt.plot(r3,sigma3,label='$\Sigma_{gas}$ Kumari+20')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Sigma_{H_2}$ $[M_\odot pc^{-2}]$',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('sigma H2'))

# plt.show()

# #######################################################################################################################################

# # velocity dispersion

# vel_disp_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\HI velocity dispersion-WBD..csv')
# vel_disp2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\velocity disp._tamburro.csv')

# r1=np.array(vel_disp_6946['r'])*(7.72/6) #boomsma data
# veldisp1=np.array(vel_disp_6946['v'])

# r2=np.array(vel_disp2_6946['r'])*(7.72/5.5) #tamburro data
# veldisp2=np.array(vel_disp2_6946['v disp'])

# plt.plot(r1,veldisp1,label='Tamburro+09')
# plt.plot(r2,veldisp2,label='Boomsma+08')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\sigma_{\star}$ $[km s^{-1}]$',fontsize=13)
# # plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('vel disp'))

# plt.show()

# #######################################################################################################################################

# # sigma sfr

# sigma_sfr_FUV_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\SFR_FUV24 6946.csv')
# sigma_sfr_Halpha_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\SFR_Halpha24 6946.csv')
# sigma_sfr_Halpha_6946_b=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\sigma_sfr_bigiel.csv')

# r1=np.array(sigma_sfr_FUV_6946['r']*(7.72/6.8))
# sigmasfr1=np.array(sigma_sfr_FUV_6946['sigma_from_plot']*(m.cos(i_6946)/m.cos(m.radians(33))))

# r2=np.array(sigma_sfr_Halpha_6946['r']*(7.72/6.8))
# sigmasfr2=np.array(sigma_sfr_Halpha_6946['sigma_from_plot']*(m.cos(i_6946)/m.cos(m.radians(33))))

# r3=np.array(sigma_sfr_Halpha_6946_b['r']*(7.72/5.9))
# sigmasfr3=np.array(sigma_sfr_Halpha_6946_b['sigma_sfr_from_plot']*(m.cos(i_6946)/m.cos(m.radians(33))))

# plt.plot(r3,sigmasfr3,label='Bigiel+08')
# plt.plot(r1,sigmasfr1,label='FUV Kumari+20')
# plt.plot(r2,sigmasfr2,label='$H_\\alpha$ Kumari+20')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Sigma_{SFR} \\times 10^{-3}$ $ [M_\odot kpc^{-2} yr^{-1}]$',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('sigma SFR'))

# plt.show()

# #######################################################################################################################################

# # sigma SMD

# SMDf=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\smdf.csv')
# SMDs=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\smd-s.csv')
# JalochaSMDf=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\Jalocha smd..csv')

# r1=np.array(SMDf['R']*(7.72/5.5))
# smdf=np.array(SMDf['smd-f'])*(m.cos(i_6946)/m.cos(m.radians(33)))

# r2=np.array(SMDs['r']*(7.72/5.5))
# smds=np.array(SMDs['smd-s'])*(m.cos(i_6946)/m.cos(m.radians(33)))

# r3=np.array(JalochaSMDf['r']*(7.72/6))
# smdf_jalocha=np.array(JalochaSMDf['smdf']*(m.cos(i_6946)/m.cos(m.radians(30))))

# plt.plot(r1,smdf,label='smd-f Sofue+18')
# plt.plot(r2,smds,label='smd-s Sofue+18')
# plt.plot(r3,smdf_jalocha,label='smd-f Jalocha+10')

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('$\Sigma_{\star}$ $[M_\odot pc^{-2} Gyr^{-1}]$',fontsize=13)
# plt.yscale('log')

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('smd'))

# plt.show()

# #######################################################################################################################################

# from scipy import stats

# # # temperature #fitting done in curve_fitting.py
# temp_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\csv\electron temp.csv')

# r1=np.array(temp_6946['r'])*(7.72/5.9)
# temp=np.array(temp_6946['temp'])
# error_temp=np.array(temp_6946['error_kelvin'])
# slope, intercept, r_value, p_value, std_err = stats.linregress(r1, temp)

# temp_corr=np.array([r*slope+intercept for r in r1])
# print(error_temp)
# plt.errorbar(r1,temp, yerr=error_temp,fmt=' ',capsize=5,label='Bresolin+04',marker='o')
# plt.plot(r1,temp_corr,label='$T={}r +{}$'.format(np.round(slope,2),np.round(intercept,2)))

# plt.legend(fontsize=14)
# plt.minorticks_on()
# plt.tick_params(axis='both',labelsize=12, which='both', width=2)
# plt.xlabel('r [kpc]',fontsize=13)
# plt.ylabel('T [K]',fontsize=13)

# plt.savefig(r'D:\Documents\Gayathri_college\MSc project\plots for supplementary material'+'\{}'.format('T'))

# plt.show()

