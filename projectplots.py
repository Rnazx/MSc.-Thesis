import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import unicodedata
import math as m

i_m51=m.radians(20)
i_6946=m.radians(30)


#######################################################################################################################################

# all plots for m51

# rotation curve
RC1_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\RCm51.csv')
RC2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\RC_sofue+18.csv')

r1=np.array(RC1_m51['R (kpc)']*(8.5/9.6))
v1=np.array(RC1_m51['v (km/s)'])

r2=np.array(RC2_m51['R']*(8.5/9.6))
v2=np.array(RC2_m51['v'])

plt.plot(r1,v1,label='Sofue+99')
plt.plot(r2,v2,label='Sofue+18')
plt.title('M51 rotation curve')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('Velocity (km/s)')


# plt.show()

#######################################################################################################################################

# omega
Omega = unicodedata.lookup('GREEK CAPITAL LETTER OMEGA')

omega1_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\omega.csv')
omega2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\omega_sofue+18.csv')

print(omega1_m51.head())
print(omega2_m51.head())

r1=np.array(omega1_m51['R (kpc)']*(8.5/9.6))
omega1=np.array(omega1_m51['omega'])

r2=np.array(omega2_m51['R']*(8.5/9.6))
omega2=np.array(omega2_m51['omega'])

plt.plot(r1,omega1,label='Sofue+99')
plt.plot(r2,omega2,label='Sofue+18')
plt.title('M51 '+Omega)
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Omega)
plt.yscale('log')

# plt.show()


#######################################################################################################################################

# sihma_HI, sigma_H2 and sigma_gas 

Sigma=unicodedata.lookup('GREEK CAPITAL LETTER SIGMA')

HI_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\HI m51...csv')
H2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\H2 m51..csv')
sigma_gas=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\HI+H2 m51..csv')

print(HI_m51.head())
print(H2_m51.head())
print(sigma_gas.head())

r1=np.array(HI_m51['r']*(8.5/8.2))
sigma1=np.array(HI_m51['sigma_HI'])

r2=np.array(H2_m51['r']*(8.5/8.2))
sigma2=np.array(H2_m51['sigma H2'])

r3=np.array(sigma_gas['r']*(8.5/8.2))
sigma3=np.array(sigma_gas['sigma_HI'])

plt.plot(r1,sigma1,label=Sigma+'HI')
plt.plot(r2,sigma2,label=Sigma+'H2')
plt.plot(r3,sigma3,label=Sigma+'(gas)')
plt.title('M51 '+Sigma+'(gas)')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Sigma+'(gas)')
plt.yscale('log')


# plt.show()


#######################################################################################################################################

# velocity dispersion

vel_disp_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\HI velcity dispersion..csv')
vel_disp_m51_schuster=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\CO vel dispersion schuster..csv')

vel_disp_m51_schuster['corrected radius']=vel_disp_m51_schuster['r']*(8.5/8.4)
print(vel_disp_m51_schuster.head())

r1=np.array(vel_disp_m51_schuster['r'])
veldisp=np.array(vel_disp_m51_schuster['v disp'])

r2=np.array(vel_disp_m51_schuster['corrected radius'])
veldisp2=np.array(vel_disp_m51_schuster['v disp'])

# plt.plot(r1,veldisp,label='Schuster+18')
plt.plot(r2,veldisp2,label='Schuster+18')
plt.title('M51 velocity dispersion')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('Velocity (km/s)')

# plt.show()

#######################################################################################################################################

# sigma sfr

alpha=unicodedata.lookup('GREEK CAPITAL LETTER ALPHA')
Sigma=unicodedata.lookup('GREEK CAPITAL LETTER SIGMA')

# sigma_sfr_FUV_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\SFR_FUV24 m51..csv')
# sigma_sfr_Halpha_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\SFR_Halpha24 m51..csv')

sigma_sfr_FUV_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\SFR_FUV24 6946.csv')
sigma_sfr_Halpha_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\SFR_Halpha24 6946.csv')

print(sigma_sfr_FUV_m51.head())
print(sigma_sfr_Halpha_m51.head())

r1=np.array(sigma_sfr_FUV_m51['r']*(7.72/6.8))
sigmasfr1=np.array(sigma_sfr_FUV_m51['sigma'])

r2=np.array(sigma_sfr_FUV_m51['r']*(7.72/6.8))
sigmasfr2=np.array(sigma_sfr_Halpha_m51['sigma'])

plt.plot(r1,sigmasfr1,label='FUV')
plt.plot(r2,sigmasfr2,label='H_alpha')
plt.title('NGC 6946 '+Sigma+'(SFR)')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Sigma+'(SFR)*10^(-3)')
plt.yscale('log')

# plt.show()


#######################################################################################################################################

# sigma SMD

SMDf=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\smdf.csv')
SMDs=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\smds.csv')

print(SMDf.head())
print(SMDs.head())

r1=np.array(SMDf['r']*(8.5/9.6))
smdf=np.array(SMDf['smdf'])

r2=np.array(SMDs['r']*(8.5/9.6))
smds=np.array(SMDs['smds'])

plt.plot(r1,smdf,label='smd-f')
plt.plot(r2,smds,label='smd-s')
plt.title('M51 Surface mass density')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Sigma+'(total)')
plt.yscale("log")


# plt.show()
#######################################################################################################################################

# temperature #fitting done in curve_fitting.py

temp_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\m51\temperature.csv')

print(temp_m51.head())

temp_m51['corrected radius']=temp_m51['radius']*(8.5/7.6)
r1=np.array(temp_m51['radius'])
temp=np.array(temp_m51['avg temp'])

plt.plot(np.array(temp_m51['corrected radius']),temp,label='corr_bresolin')
plt.plot(r1,temp,label='Bresolin+04')
plt.title('M51 electron temperature')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('temperature (K)')

# plt.show()

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


# NGC 6946

#######################################################################################################################################

# rotation curve

RC1_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\RC.csv') #sparc data with error
RC2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\RC_sofue+18.csv')
RC3_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\RC_sofue+99.csv')

print(RC1_6946.head())
print(RC2_6946.head())
print(RC3_6946.head())

r1=np.array(RC1_6946['# Rad']*(7.72/5.52))
v1=np.array(RC1_6946['Vobs']*(m.cos(i_6946)/m.cos(m.radians(38))))
error=np.array(RC1_6946['errV']*(m.cos(i_6946)/m.cos(m.radians(38))))

r2=np.array(RC2_6946['R']*(7.72/5.5))
v2=np.array(RC2_6946['v'])

r3=np.array(RC3_6946['r']*(7.72/5.5))
v3=np.array(RC3_6946['v'])

# plt.plot(r1,v1,label='SPARC data')
plt.errorbar(r1, v1, yerr=error,label='SPARC data', fmt='.', capsize=1) #to be used when there are rror bars
plt.plot(r2,v2,label='Sofue+18')
plt.plot(r3,v3,label='Sofue+99')
plt.title('NGC 6946 rotation curve')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('Velocity (km/s)')

# plt.show()

#######################################################################################################################################

# omega
Omega = unicodedata.lookup('GREEK CAPITAL LETTER OMEGA')

omega1_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\omega.csv')
omega2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\omega_sofue+18.csv')
omega3_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\omega_sofue+99.csv')

print(omega1_6946.head())
print(omega2_6946.head())
print(omega3_6946.head())

r1=np.array(omega1_6946['# Rad']*(7.72/5.52))
omega1=np.array(omega1_6946['omega']*(m.cos(i_6946)/m.cos(m.radians(38))))
error=np.array(omega1_6946['error_omega']*(m.cos(i_6946)/m.cos(m.radians(38))))

r2=np.array(omega2_6946['R']*(7.72/5.5))
omega2=np.array(omega2_6946['omega'])

r3=np.array(omega3_6946['r']*(7.72/5.5))
omega3=np.array(omega3_6946['omega'])

# plt.plot(r1,omega1,label='SPARC data')
plt.errorbar(r1, omega1, yerr=error,label='SPARC data', fmt='.', capsize=2) #to be used when there are rror bars
plt.plot(r2,omega2,label='Sofue+18')
plt.plot(r3,omega3,label='Sofue+99')
plt.title('NGC 6946 '+Omega)
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Omega)
plt.yscale('log')

# plt.show()

#######################################################################################################################################

# sigma_HI, sigma_H2 and sigma_gas 

HI_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\HI 6946.csv')
H2_m51=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\H2 6946.csv')
sigma_gas=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\HI+H2 6946.csv')

print(HI_m51.head())
print(H2_m51.head())
print(sigma_gas.head())

r1=np.array(HI_m51['r']*(7.72/6.8))
sigma1=np.array(HI_m51['sigma']*(m.cos(i_6946)/m.cos(m.radians(33))))

r2=np.array(H2_m51['r']*(7.72/6.8))
sigma2=np.array(H2_m51['sigma']*(m.cos(i_6946)/m.cos(m.radians(33))))

r3=np.array(sigma_gas['r']*(7.72/6.8))
sigma3=np.array(sigma_gas['sigma']*(m.cos(i_6946)/m.cos(m.radians(33))))

plt.plot(r1,sigma1,label=Sigma+'HI')
plt.plot(r2,sigma2,label=Sigma+'H2')
plt.plot(r3,sigma3,label=Sigma+'(gas)')
plt.title('NGC 6946 '+Sigma+'(gas)')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Sigma+'(gas)')
plt.yscale('log')

# plt.show()

#######################################################################################################################################

# velocity dispersion

vel_disp_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\HI velocity dispersion-WBD..csv')
vel_disp2_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\velocity disp..csv')

vel_disp2_6946['corrected radius']=vel_disp2_6946['r']*(7.72/6)
# print(vel_disp_6946.head())
print(vel_disp2_6946.head())

r1=np.array(vel_disp2_6946['corrected radius'])
# veldisp=np.array(vel_disp_6946['v'])

r2=np.array(vel_disp2_6946['r'])
veldisp2=np.array(vel_disp2_6946['v disp'])

plt.plot(np.array(vel_disp2_6946['corrected radius']),veldisp2,label='Boomsma+08')
#plt.plot(r2,veldisp2,label='Boomsma+08')
plt.title('NGC 6946 velocity dispersion')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('velocity (km/s)')

# plt.show()

#######################################################################################################################################

# sigma sfr

alpha=unicodedata.lookup('GREEK CAPITAL LETTER ALPHA')
Sigma=unicodedata.lookup('GREEK CAPITAL LETTER SIGMA')

sigma_sfr_FUV_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\SFR_FUV24 6946.csv')
sigma_sfr_Halpha_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\SFR_Halpha24 6946.csv')

print(sigma_sfr_FUV_6946.head())
print(sigma_sfr_Halpha_6946.head())

r1=np.array(sigma_sfr_FUV_6946['r']*(7.72/6.8))
sigmasfr1=np.array(sigma_sfr_FUV_6946['sigma']*(m.cos(i_6946)/m.cos(m.radians(33))))

r2=np.array(sigma_sfr_FUV_6946['r']*(7.72/6.8))
sigmasfr2=np.array(sigma_sfr_Halpha_6946['sigma']*(m.cos(i_6946)/m.cos(m.radians(33))))

plt.plot(r1,sigmasfr1,label='FUV')
plt.plot(r2,sigmasfr2,label='H_alpha')
plt.title('NGC 6946 '+Sigma+'(SFR)')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel(Sigma+'(SFR)')
plt.yscale('log')


# plt.show()

#######################################################################################################################################

# sigma SMD

SMDf=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\smdf.csv')
SMDs=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\smd-s.csv')
JalochaSMDf=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\Jalocha smd..csv')

print(SMDf.head())
print(SMDs.head())

r1=np.array(SMDf['R']*(7.72/5.5))
smdf=np.array(SMDf['smd-f'])

r2=np.array(SMDs['r']*(7.72/5.5))
smds=np.array(SMDs['smd-s'])

r3=np.array(JalochaSMDf['r']*(7.72/6))
smdf_jalocha=np.array(JalochaSMDf['smdf']*(m.cos(i_6946)/m.cos(m.radians(38))))

plt.plot(r1,smdf,label='smd-f')
plt.plot(r2,smds,label='smd-s')
plt.plot(r3,smdf_jalocha,label='smd-f Jalocha+10')
plt.title('NGC 6946 SMD')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('SMD')
plt.yscale("log")

# plt.show()

#######################################################################################################################################

# temperature

temp_6946=pd.read_csv(r'D:\Documents\Gayathri_college\MSc project\data\Ngc6946\electron temp.csv')
temp_6946['corrected radius']=temp_6946['r']*(7.72/5.9)

# print(temp_6946.head())


r1=np.array(temp_6946['r'])
r2=np.array(temp_6946['corrected radius'])
temp=np.array(temp_6946['t_NS'])

plt.plot(r1,temp,'o',label='Gusev+12')
plt.plot(r2,temp,'o',label='Gusev+12c')
plt.title('NGC 6946 electron temperature')
plt.legend()
plt.xlabel('Radius (kpc)')
plt.ylabel('Temperature (K)')

# plt.show()