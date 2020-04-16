import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 2 :
	print ("usage: python3 ipole.py ipole.dat")
	quit()

fil = sys.argv[1]

# read in data 
i0, j0, Ia, Is, Qs, Us, Vs, tauF = np.loadtxt(fil, unpack=True)

# set image size
ImRes = int(round(np.sqrt(len(i0))))
print ("Image resolution: ", ImRes)

# size of single pixel in rad: M/pixel . muas/pix . rad/muas
FOV = 40#4.17
print ("assuming FOV = ", FOV, "GM/c^2")
#
# solid angle subtended by pixel
flux = sum(Is)
print ("Flux [Jy]: ", flux)

# recast indices into offset in units of M
i = (np.reshape(i0, (ImRes,ImRes))+1)*FOV/ImRes - FOV/2.
j = (np.reshape(j0, (ImRes,ImRes))+1)*FOV/ImRes - FOV/2.

# LP plot
plt.subplot(2,2,2)
lpfrac = 100.*np.sqrt(Qs*Qs + Us*Us)/Is
z = np.reshape(lpfrac, (ImRes,ImRes))
plt.pcolormesh(i,j,z,cmap='jet', vmin = 0., vmax = 100.)
plt.title('LP [%]')
plt.axis([-FOV/2,FOV/2,-FOV/2,FOV/2])
plt.colorbar()

# EVPA plot
plt.subplot(2,2,3)
# adjusted to IAU conventions
evpa = np.angle(Qs+1j*Us)/2*180/np.pi
z = np.reshape(evpa, (ImRes,ImRes))
plt.pcolormesh(i,j,z,cmap='jet')
plt.title('EVPA [deg]')
plt.axis([-FOV/2,FOV/2,-FOV/2,FOV/2])
plt.colorbar()

# CP plot
plt.subplot(2,2,4)
cpfrac = 100.*Vs/Is
z = np.reshape(cpfrac, (ImRes,ImRes))
plt.pcolormesh(i,j,z,cmap='jet', vmin = -5, vmax = 5.)
plt.title('CP [%]')
plt.axis([-FOV/2,FOV/2,-FOV/2,FOV/2])
plt.colorbar()

# total intensity 
plt.subplot(2,2,1)
z = np.reshape(Is, (ImRes,ImRes))
#plt.pcolormesh(i,j,z,cmap='afmhot', vmin=0., vmax=1.e-4)
plt.pcolormesh(i,j,z,cmap='afmhot', vmin=0., vmax=5.e-4)
plt.colorbar()
plt.title('Stokes I [cgs]')
plt.axis([-FOV/2,FOV/2,-FOV/2,FOV/2])

# superpose EV on total intensity using quiver
Qb = sum(Qs)
Ub = sum(Us)
LP = np.sqrt(Qb*Qb + Ub*Ub)/sum(Is)
print ("LP [%]: ", LP*100.)
CHI = (180./3.14159)*0.5*np.arctan2(Ub,Qb)
print ("EVPA [deg]:",CHI)
CP = sum(Vs)/sum(Is)
print ("CP [%]: ", CP*100.)
amp = np.sqrt(Qs*Qs + Us*Us)
scal = max(amp)


# brand new adjusted to IAU conventions and ehtim
vxp = -np.sqrt(Qs*Qs + Us*Us)*np.sin(evpa*3.14159/180.)/scal
vyp = np.sqrt(Qs*Qs + Us*Us)*np.cos(evpa*3.14159/180.)/scal




vx = np.reshape(vxp, (ImRes,ImRes))
vy = np.reshape(vyp, (ImRes,ImRes))
skip = 4
plt.quiver(i[::skip, ::skip],j[::skip, ::skip],vx[::skip, ::skip],vy[::skip, ::skip], 
	headwidth=1, headlength=1, 
	width=0.005,
	color='green', 
	units='width', 
	scale=16)

plt.subplots_adjust(wspace=0.3,hspace=0.3)

# show, or save
plt.show()
#plt.savefig('tst')

