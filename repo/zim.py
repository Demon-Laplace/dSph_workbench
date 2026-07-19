import numpy as np
import math as m

def haveing_circle() :
    _N=50 ; dpi=np.pi/4. ; A = np.linspace(0, 2*np.pi,_N) + dpi ; _circlex = np.cos(A); _circley = np.sin(A)
    return(_circlex, _circley)



def Rx(theta):
  return(np.matrix([[ 1, 0           , 0           ],
                   [ 0, m.cos(theta),-m.sin(theta)],
                   [ 0, m.sin(theta), m.cos(theta)]])  )
  
def Ry(theta):
  return(np.matrix([[ m.cos(theta), 0, m.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-m.sin(theta), 0, m.cos(theta)]]) )
  
def Rz(theta):
  return(np.matrix([[ m.cos(theta), -m.sin(theta), 0 ],
                   [ m.sin(theta), m.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]]))





# ;-------------------------------------------------------------
# ;+
# ; NAME:
# ;       ROT_3D
# ; PURPOSE:
# ;       Rotate 3-d coordinate system.
# ; CATEGORY:
# ; CALLING SEQUENCE:
# ;       rot_3d, axis, x1, y1, z1, ang, x2, y2, z2
# ; INPUTS:
# ;       axis=Axis number to rotate about: 1=X, 2=Y, 3=Z.     in
# ;       x1, y1, z1 = arrays of original x,y,z vector comp.   in
# ;       ang = rotation angle in radians.                     in
# ; KEYWORD PARAMETERS:
# ;       Keywords:
# ;         /DEGREES means angle is in degrees, else radians.
# ; OUTPUTS:
# ;       x2, y2, z2 = arrays of new x,y,z vector components.  out
# ; COMMON BLOCKS:
# ; NOTES:
# ;       Note: Right-hand rule is used: Point thumb along +axis.
# ;         Fingers curl in vector rotation direction (for +ang).
# ;         This is for coordinate system rotation.  To rotate the
# ;         vectors in a fixed coord. system use the left hand rule.
# ; MODIFICATION HISTORY:
# ;       R. Sterner.  28 Jan, 1987.
# ;       6 May, 1988 --- modified to work with any shape arrays.
# ;       R. Sterner, 6 Nov, 1989 --- converted to SUN.
# ;       RES 13 Feb, 1991 --- added /degrees.
# ;       Johns Hopkins University Applied Physics Laboratory.
# ;
# ; Copyright (C) 1987, Johns Hopkins University/Applied Physics Laboratory
# ; This software may be used, copied, or redistributed as long as it is not
# ; sold and this copyright notice is reproduced on each copy made.  This
# ; routine is provided as is without any express or implied warranties
# ; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
# ;-
# ;-------------------------------------------------------------
def rot_3d(axis,x1,y1,z1,ang,degrees=True) :  
# 	pro rot_3d, axis, x1, y1, z1, ang, help=hlp, $
# 	  degrees=degrees
 
# 	if (n_params(0) ne 5) or keyword_set(hlp) then begin
# 	  print,' Rotate 3-d coordinate system.'
# 	  print,' rot_3d, axis, x1, y1, z1, ang, x2, y2, z2'
# 	  print,'   axis=Axis number to rotate about: 1=X, 2=Y, 3=Z.     in'
# 	  print,'   x1, y1, z1 = arrays of original x,y,z vector comp.   in'
# 	  print,'   ang = rotation angle in radians.                     in'
# 	  print,'   x2, y2, z2 = arrays of new x,y,z vector components.  out'
# 	  print,' Keywords:'
#           print,'   /DEGREES means angle is in degrees, else radians.'
# 	  print,' Note: Right-hand rule is used: Point thumb along +axis.'
# 	  print,'   Fingers curl in vector rotation direction (for +ang).'
# 	  print,'   This is for coordinate system rotation.  To rotate the'
# 	  print,'   vectors in a fixed coord. system use the left hand rule.'
# 	  return
# 	endif
 
    if degrees :
        radeg = (180./np.pi)
        c = np.cos(ang/radeg)
        s = np.sin(ang/radeg)
    else :
        c = np.cos(ang)
        s = np.sin(ang)
# 	case axis of			; depending on axis.
    if (axis == 1) : 
        x2 =  x1
        y2 =  c*y1 + s*z1
        z2 = -s*y1 + c*z1
    if (axis == 2) : 
        x2 = c*x1 - s*z1
        y2 = y1
        z2 = s*x1 + c*z1
    if (axis == 3) : 
        x2 =  c*x1 + s*y1
        y2 = -s*x1 + c*y1
        z2 = z1
    return x2,y2,z2





def unitconversion(u0,ne0):
# ;/*
# ;  double GRAVITY, BOLTZMANN, PROTONMASS;
# ;  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
# ;  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
# ;  double G, Xh, HubbleParam;
# ;*/

# ;  int i;
# ;  double MeanWeight, u, gamma, u1;

# ;  /* physical constants in cgs units */
    u0 = np.asarray(u0, dtype=np.float64)
    ne0 = np.asarray(ne0, dtype=np.float64)

    GRAVITY   = 6.672e-8;
    BOLTZMANN = 1.3806e-16;
    PROTONMASS = 1.6726e-24;

# ;  /* internal unit system of the code */
    UnitLength_in_cm= 3.085678e21;  # /*  code length unit in cm/h */
    UnitMass_in_g= 1.989e43;        # /*  code mass unit in g/h */
    UnitVelocity_in_cm_per_s= 1.0e5;

    UnitTime_in_s= UnitLength_in_cm/ UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs= UnitMass_in_g/(UnitLength_in_cm**3);
    UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ (UnitTime_in_s**2);
    UnitEnergy_in_cgs = UnitMass_in_g * (UnitLength_in_cm**2) / (UnitTime_in_s**2);


    # ;  G=GRAVITY/(UnitLength_in_cm^3) * UnitMass_in_g * (UnitTime_in_s^2);

    Xh= 0.76;  #/* mass fraction of hydrogen */
    # ;  HubbleParam= 0.65;

    # ;  for(i=1; i<=NumPart; i++)
    # ;    {
    # ;      if(P[i].Type==0)  /* gas particle */
    # ;	{
    # ;	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;
    MeanWeight = 4.0/(3*Xh+1+4*Xh*ne0) * PROTONMASS
    # ;	  /* convert internal energy to cgs units */

    # ;	  u  = (P[i].U + P[i].Fd) * UnitEnergy_in_cgs/ UnitMass_in_g;
    # ;  u = (u0+fd0) * UnitEnergy_in_cgs/ UnitMass_in_g

    gamma= 5.0/3                  ;

    # ;	  /* get temperature in Kelvin */

    # ;	  P[i].Tempeff= MeanWeight/BOLTZMANN * (gamma-1) * u;
    # ;  tempeff=MeanWeight/BOLTZMANN * (gamma-1) * u

    # ;	  u  = (P[i].U ) * UnitEnergy_in_cgs/ UnitMass_in_g;
    u=u0* UnitEnergy_in_cgs/ UnitMass_in_g

    # ;	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
    temp = MeanWeight/BOLTZMANN * (gamma-1) * u
    return temp



def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


def projcolm(x0,y0,d0,m0,rho0,v0,bd,nsize,vrange,dv,fhi=0.76) :
    GRAVITY   = 6.672e-8;
    BOLTZMANN = 1.3806e-16;
    PROTONMASS = 1.6726e-24;
    #;  /* internal unit system of the code */
    UnitLength_in_cm= 3.085678e21;  # /*  code length unit in cm/h */
    UnitMass_in_g= 1.989e43;       #  /*  code mass unit in g/h */
    UnitVelocity_in_cm_per_s= 1.0e5;

    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs= UnitMass_in_g/(UnitLength_in_cm**3);
    UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ (UnitTime_in_s**2);
    UnitEnergy_in_cgs= UnitMass_in_g * (UnitLength_in_cm**2) / (UnitTime_in_s**2);

    rhofactor = UnitDensity_in_cgs / PROTONMASS 
#     print(4*4)
#   if not keyword_set(cloudsize) then begin
#      cloudsize=30 ;; degree  
#   endif else print,'filter clouds size > ',cloudsize, 'degree'
 
    hrfactor=1
    
    nx=(nsize[0]*hrfactor)
    ny=(nsize[1]*hrfactor)

# ;; determinlocations
    x1=bd[0] ; x2= bd[1] ; y1=bd[2] ;  y2=bd[3]
    dx=(x2-x1)/nx
    dy=(y2-y1)/ny
    print(x1,x2,dx)
#     img=np.array(nx,ny)
    ix=(x0-x1)/dx
    iy=(y0-y1)/dy
    ix=ix.astype(int)
    iy=iy.astype(int)
    
# # ;;; find v index
# #    dv=5                         ;km/s
# #   if keyword_set(dv0) then dv=dv0
    vmax=vrange[1]
    vmin=vrange[0]
    nv=int((vmax-vmin) / dv)
    iv=(v0 - vmin) / dv
    iv=iv.astype(int)
    

# #     # ;; create cube and project data.
    cube = np.zeros( (nx,ny,nv) )
    
    idx = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (iv >= 0) & (iv < nv) 

    
    if sum(idx) == 0 :
        print( "no points in side")
        return(cube)
           
    ix = ix[idx]
    iy = iy[idx]
    iv = iv[idx]
    m  = m0[idx]  / 1e10   # to code unit
    d  = d0[idx]  
    v  = v0[idx]
    rho= rho0[idx] / rhofactor  # to code unit

    

# # ;calculate volme --> scale of cloud
    vol=m/rho # kpc^3    ## also in code unit
    scl=vol**(1./3) # kpc
    ll = (scl/d) /np.pi * 180 * 60 # in arcmin  ;; matching the grid asked.

    beamfactor=6
    
    
    NHI = fhi * rho * rhofactor * scl * UnitLength_in_cm  # atoms / cm^3 * cm = atom / cm^2
    ## add hydrogen fraction 0.76, will affect all clauclation done before 08/09/2021. 

    
    ## calcualte extension of each point
    
    dkll = (ll  / dx / 2).astype(int) 
    djll = (ll  / dy / 2).astype(int) 
    lx1 = ix-dkll 
    lx2 = ix+dkll 
    ly1 = iy-djll 
    ly2 = iy+djll 

#     lx1[lx1<0]=0
#     lx2[lx2>nx]=nx-1
#     ly1[ly1<0]=0
#     ly2[ly2>ny]=ny-1
    

    for i in range(len(ll)):
        if (lx1[i] >=1)  & (lx2[i] < nx-1)  & (ly1[i] >=1) & (ly2[i] < ny-1): 
#             print(lx2[i]-lx1[i]+1,dkll[i]/2.0)
#             break
            clip = makeGaussian (lx2[i]-lx1[i],fwhm=dkll[i]/beamfactor)            
            #plt.imshow(clip)
            clip = clip / np.sum(clip) *  NHI[i] #clip.size *
            #print(clip.shape)
            #print(lx1[i],lx2[i],ly1[i],ly2[i])
            cube[lx1[i]:lx2[i],ly1[i]:ly2[i],iv[i]] += clip #NHI[i]
            
    return(cube)
