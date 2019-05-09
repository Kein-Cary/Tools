#!/usr/bin/env python3

from astropy.cosmology import LambdaCDM
import numpy as np

d = open( './c_shape' ).readlines()
temp = open( './sqlcl.in' ).readlines()
cos = LambdaCDM( H0=73, Om0=0.27, Ode0=0.73 )

ra = []
dec = []
z = []
#print( d )
for i in d:
    t = i[:-1].split()
    ra.append( float(t[0]) )
    dec.append( float(t[1]) )
    z.append( float(t[-3]) )

cl = 'dbo.fGetNearbyObjEq(%f, %f, %f) AS GN\n'
R = 5.0

for i in range( len(ra) ):
    f = open( './sql/sql_%i'%(i), 'w' )
    r = R / cos.angular_diameter_distance( z[i] ).value / np.pi * 180 * 60
    for j in temp:
        if ( j[0:3] == 'dbo' ):
            f.write( cl%(ra[i], dec[i], r) )
        else:
            f.write( j )
    f.close()
    #break
