#!/usr/bin/env python

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.ticker import FuncFormatter, MultipleLocator # for axis multiple of pi

"""
healpix.py: Generator of HEALPix grid

Healpix points are plotted on the sphere or in the lon-lat plane.
Rhomboids are created with its centroids in healpix points.
On the sphere, the edges of the rhomboids are great arcs.
"""

__author__  = """
slavko.brdar@ecmwf.int,
"""

# parameters
ns = 2               # for healpix H<ns>
sph_plot = True      # plot on sphere


pi = np.pi
floor = np.floor
sqrt = np.sqrt
cos = np.cos
sin = np.sin
arccos = np.arccos

def slerp(p1, p2, t):
    omega = arccos( p1.dot(p2) )
    sin_omega = sin(omega)    
    t = t[:, np.newaxis]
    return ( sin( (1-t)*omega )*p1 + sin( t*omega )*p2 )/sin_omega

def linp(p1, p2, t):
    omega = arccos( p1.dot(p2) )
    sin_omega = sin(omega)    
    t = t[:, np.newaxis]
    return (1-t)*p1 + t*p2

npix = 12*ns*ns
nrings = 4*ns-1
rad2deg = 180/pi
deg2rad = 1/rad2deg
ncap_ids = np.arange(2)

# HEALPix dual points
dualp = np.zeros(2*npix).reshape(npix,2)

# north cap
for p in range(0,npix):
    ph = 0.5*(p+1)
    i = floor(sqrt(ph-sqrt(floor(ph)))) + 1
    if( i > ns ):
        ncap_ids[0] = p
        break
    j = p+1-2*i*(i-1)
    dualp[p,0] = (arccos(1-i*i/(3.*ns*ns)))*rad2deg
    dualp[p,1] = (0.25*pi*(2*j-1)/i)*rad2deg

# north belt
for p in range(ncap_ids[0],npix):
    ph = p-2*ns*(ns-1) 
    i = floor(0.25*ph/ns) + ns
    if( i > 2*ns ):
        ncap_ids[1] = p
        break
    j = 1 + ph%(4*ns)
    s = (i-ns+1)%2
    dualp[p,0] = (arccos((4*ns-2*i)/(3.*ns)))*rad2deg
    dualp[p,1] = (0.25*pi*(2*j-s)/ns)*rad2deg

# south belt & south cap
for p in range(ncap_ids[1],npix):
    dualp[p,0] = 180-dualp[npix-p-1,0]
    dualp[p,1] = dualp[npix-p-1,1]

# HEALPix vertex points
vertp = np.zeros(2*(npix+2)).reshape(npix+2,2)

vertp[0,0] = 0
vertp[0,1] = 180
vertp[npix+1,0] = 180
vertp[npix+1,1] = 180

# north cap
for p in range(0,ncap_ids[0]):
    ph = 0.5*(p+1)
    i = floor(sqrt(ph-sqrt(floor(ph)))) + 1
    last_id = 2*i*(i-1)
    j = p + 1 - last_id
    vertp[p+1,0] = dualp[p,0]
    pp = int(p-1 if j>1 else last_id+4*i-1)
    vertp[p+1,1] = dualp[pp,1] + dualp[p,1]
    vertp[p+1,1] = 0.5*vertp[p+1,1] if p>pp else 0.5*(vertp[p+1,1]+360.)

# north belt
for p in range(ncap_ids[0],ncap_ids[1]):
    ph = p-2*ns*(ns-1) 
    i = floor(0.25*ph/ns) + ns
    j = 1 + ph%(4*ns)
    vertp[p+1,0] = dualp[p,0]
    pp = int(p-1 if j>1 else ncap_ids[0]+4*ns*(i-ns)-1)
    vertp[p+1,1] = dualp[pp,1] + dualp[p,1]
    vertp[p+1,1] = 0.5*vertp[p+1,1] if p>pp else 0.5*(vertp[p+1,1]+360.)

# south belt & south cap
for p in range(ncap_ids[1],npix):
    vertp[p+1,0] = 180-vertp[npix-p,0]
    vertp[p+1,1] = vertp[npix-p,1]

# HEALPix cells
el_vid = np.zeros(4*npix).reshape(npix,4)

# north cap for cells
for el in range(0,ncap_ids[0]):
    ph = 0.5*(el+1)
    i = floor(sqrt(ph-sqrt(floor(ph)))) + 1
    j = el+1-2*i*(i-1)
    el_vid[el,0] = ((j-int((j-1)/i) if j<4*i else 1) if i>1 else 0) + ( 2*(i-1)*(i-2) if i>2 else 0 )
    el_vid[el,1] = j+2*i*(i-1)
    el_vid[el,3] = ( j+1+2*i*(i-1) if j<4*i else 1+2*i*(i-1) )
    if( el < 2*(ns-1)*ns ):
        el_vid[el,2] = 2*(i+1)*i + j + 1 + int((j-1)/i)
    else:
        el_vid[el,2] = 2*(i+1)*i + j

# north belt for cells
for el in range(ncap_ids[0],6*ns*ns+2*ns):
    ph = el-2*ns*(ns-1) 
    i = floor(0.25*ph/ns) + ns
    j = 1 + ph%(4*ns)
    s = (i-ns+1)%2
    el_vid[el,0] = 2*ns*(ns-1)+4*ns*(i-ns-1)+1 + (0 if j==4*ns and not s else j) - s
    el_vid[el,1] = el+1
    el_vid[el,3] = el+2 - (4*ns if j==4*ns else 0)
    if( el < 6*ns*ns-2*ns ):
        el_vid[el,2] = 1 + (2*ns*(ns-1)+4*ns*(i-ns+1) if j==4*ns and not s else el+4*ns+1) - s

# equator for cells
for el in range(6*ns*ns-2*ns, 6*ns*ns+2*ns):
    el_vid[el,2] = npix + 1 - el_vid[el,0]

# south belt & south cap for cells
for el in range(6*ns*ns+2*ns, npix):
    el_vid[el,0] = npix + 1 - el_vid[npix-el-1,0]
    el_vid[el,1] = npix + 1 - el_vid[npix-el-1,1]
    el_vid[el,3] = npix + 1 - el_vid[npix-el-1,3]
    if( el_vid[npix-el-1,2] > 6*ns*ns-2*ns and el_vid[npix-el-1,2] <= 6*ns*ns+2*ns ):
        el_vid[el,2] = el_vid[npix-el-1,2]
    else:
        el_vid[el,2] = npix + 1 - el_vid[npix-el-1,2]

# PLOTTING

# to radians
dualp = dualp * deg2rad
vertp = vertp * deg2rad
fig = plt.figure()

if (sph_plot):
    xx = sin(dualp[:,0])*cos(dualp[:,1])
    yy = sin(dualp[:,0])*sin(dualp[:,1])
    zz = cos(dualp[:,0])
    xx_v = sin(vertp[:,0])*cos(vertp[:,1])
    yy_v = sin(vertp[:,0])*sin(vertp[:,1])
    zz_v = cos(vertp[:,0])

    #eind = int(npix/2) # plot just north hemisphere
    eind = npix # plot full sphere

    # plot dual and vertex points
    ax = fig.add_subplot(111, projection='3d')
    ax.set_axis_off() # hide axis
    ax.scatter(xx[0:eind],yy[0:eind],zz[0:eind],color="red",s=20) # dual points
    ax.scatter(xx_v[0:eind],yy_v[0:eind],zz_v[0:eind],color="blue",s=10) # vertex points

    # plot rectanges of each element
    pt = np.zeros(12).reshape(4,3)
    t = np.linspace(0,1,30)
    for el in range(0, eind ):
        for ipt in range(0,4):
            pt[ipt,:] = np.array([xx_v[int(el_vid[el,ipt])],
                yy_v[int(el_vid[el,ipt])],zz_v[int(el_vid[el,ipt])]])
        for ipt in range(0,4):
            arc = slerp(pt[ipt],pt[(ipt+1)%4],t)
            ax.plot(arc[:,0], arc[:,1],arc[:,2],c="k")

    # label HEALPix dual points
    #for p in range(0, eind ):
    #    ax.text(xx[p],yy[p],zz[p],str(p))

    # label HEALPix vertex points
    for p in range(0,eind+2):
        ax.text(xx_v[p],yy_v[p],zz_v[p],str(p))

    # transparent sphere
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    x = sin(phi)*cos(theta)
    y = sin(phi)*sin(theta)
    z = cos(phi)
    #ax.plot_surface(
    #    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.4, linewidth=0)

else: # lonlat plot
    #for p in range(0,npix+2):
    #    vertp[p,1] = (vertp[p,1] if vertp[p,1]<2*pi+1e-10 else vertp[p,1]-2*pi)
    dualp[:,1] = 2.*pi-dualp[:,1]
    dualp[:,0] = pi/2.-dualp[:,0]
    vertp[:,1] = 2.*pi-vertp[:,1]
    vertp[:,0] = pi/2.-vertp[:,0]

    ax = fig.add_subplot(111)
    #ax.set_xticklabels(['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$',
    #      r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$', r'$2\pi$'])
    #ax.set_yticklabels([r'$-\pi/2$', r'$-\pi/4$', '$0$', r'$\pi/4$', r'$\pi/2$'])
    ax.set_axis_off() # hide axis
    ax.scatter(dualp[:,1], dualp[:,0],color="red",s=5) # dual points
    ax.scatter(vertp[:,1], vertp[:,0],color="blue",s=20) # vertex points

    # plot rectanges of each element
    pt = np.zeros(8).reshape(4,2)
    t = np.linspace(0,1,5)
    for el in range(0,npix):
        for ipt in range(0,4):
            indx = int(el_vid[el,ipt])
            pt[ipt,:] = np.array([vertp[indx,1], vertp[indx,0]])
        for ipt in range(0,4):
            #ax.plot( pt[ipt], pt[(ipt+1)%4])
            arc = linp(pt[ipt],pt[(ipt+1)%4],t)
            if( abs(pt[ipt,0]-pt[(ipt+1)%4,0])<1.25*pi ):
                ax.plot(arc[:,0], arc[:,1],c="k")

    # label HEALPix dual points
    for p in range(0,npix):
        ax.text(dualp[p,1], dualp[p,0], str(p))

    # label HEALPix vertex points
    for p in range(0,npix+2):
        ax.text(vertp[p,1]-pi/2.,2.*pi-vertp[p,0],str(p))

plt.tight_layout()
plt.show()
