import numpy as np
# import robustats


#copy from get_mean_pm.py, Vasiliev
def angular_distance(ra0, dec0, ra1, dec1):
    '''
    Compute the angular distance between two points on a sphere (coordinates expressed in degrees)
    '''
    d2r = np.pi/180  # degrees to radians
    return 2 * np.arcsin( (np.sin( (dec0-dec1)*0.5 * d2r )**2 +
        np.cos(dec0 * d2r) * np.cos(dec1 * d2r) * np.sin( (ra0-ra1)*0.5 * d2r )**2 )**0.5 ) / d2r


# to calcualte area of sector of an ellipse
#
#https://keisan.casio.com/exec/system/1343722259
#
def Fthe(a,b,theta_deg):
    theta = np.deg2rad(theta_deg)
    tt = (b - a) * np.sin(2*theta) / (  b + a + (b-a) * np.cos(2*theta)  )
    y = 0.5 * a * b * (theta - np.arctan(tt))
    return(y)

def Asecc(a,b,the1,the2):
    if the1 >= the2 : 
        print('require the1 < the2')
        return None
    return(Fthe(a,b,the2) - Fthe(a,b,the1) )

def rot_2d(x,y,theta,deg=True):
    the = np.radians(theta)
    c, s = np.cos(the), np.sin(the)
    rx = c*x - s*y
    ry = s*x + c*y
    return(rx,ry)


def rot_3d(axis,x1,y1,z1,ang,deg=True):
    if deg : 
        ang = np.radians(ang)
    c = np.cos(ang)
    s = np.sin(ang)
    
    if axis==1 :
        x2 =  x1
        y2 =  c*y1 + s*z1
        z2 = -s*y1 + c*z1
    if axis==2 :
        x2 = c*x1 - s*z1
        y2 = y1
        z2 = s*x1 + c*z1
    if axis==3 :
        x2 =  c*x1 + s*y1
        y2 = -s*x1 + c*y1
        z2 = z1
    return x2,y2,z2


def locateAM (x,y,z,vx,vy,vz,idx) :
    jx = y[idx] * vz[idx] - z[idx] * vy[idx]
    jy = z[idx] * vx[idx] - x[idx] * vz[idx]
    jz = x[idx] * vy[idx] - y[idx] * vx[idx]

    jtot = np.sqrt(jx**2 + jy**2 + jz**2)

    jx /= jtot
    jy /= jtot
    jz /= jtot

    jxm = np.median(jx)
    jym = np.median(jy)
    jzm = np.median(jz)

    theta = np.arccos(jzm) / np.pi * 180.

    JR = np.sqrt( jxm**2 + jym**2 )

    phi = np.arccos( jxm / JR )  / np.pi * 180.

    if jym < 0 :
        phi = 360 - phi

    jl =  phi
    jb =  90-theta
    
    return (jl,jb)

def calshape(x,y):
    x_1 = np.median(x)
    y_1 = np.median(y)
    dx =  x-x_1
    dy =  y-y_1 
    x_2 = np.mean(dx ** 2) - x_1**2
    y_2 = np.mean(dy ** 2) - y_1**2
    xy  = np.mean(x*y) - x_1*y_1 
    theta0=0.5 * np.arctan(2*(xy/(x_2-y_2))) * 180./np.pi
    if y_2 > x_2 :
        theta0=theta0+90
    a2 = 0.5*(x_2+y_2) +  np.sqrt( 0.25*(x_2-y_2)**2 +xy**2 )
    b2 = 0.5*(x_2+y_2) -  np.sqrt( 0.25*(x_2-y_2)**2 +xy**2 )
    ell = 1 - np.sqrt(b2/a2)
    print,"ELL ELLE :",ell
    return (np.sqrt(a2),np.sqrt(b2),ell,theta0)

def myweighted_average(x,err) :
    w=1./err**2
    wtot=sum(w)
    xm = sum(x*w)/wtot
    std = np.sqrt( sum( w * (x-xm)**2 ) / wtot )
    xs  = std / np.sqrt(len(x))
    return(xm,xs)

# def myweighted_median(x,err) :
#     w=1./err**2
#     wtot=sum(w)
#     robustats.weighted_median(x,w)
#     xm = robustats.weighted_median(x,w)
#     std = np.sqrt( sum( w * (x-xm)**2 ) / wtot )
#     xs  = std / np.sqrt(len(x))
#     return(xm,xs)


# def kinemaps(x,err,median=False) :
#     w=1./err**2
#     wtot=sum(w)
#     if median :
#         robustats.weighted_median(x,w)
#         xm = robustats.weighted_median(x,w)
#     else :
#         xm = sum(x*w)/wtot
#     std = np.sqrt( sum( w * (x-xm)**2 ) / wtot )
#     xs  = std / np.sqrt(len(x))
#     return(xm,xs,std)

def haveing_circle() :
    _N=50 ; dpi=np.pi/4. ; A = np.linspace(0, 2*np.pi,_N) + dpi ; _circlex = np.cos(A); _circley = np.sin(A)
    return(_circlex, _circley)


def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)





def mysigma_clip(x0,clipat=3.0) :
    x = x0
#    print(c)
    while True:
        n0 = len(x)
        ms = np.median(x)
        ss = np.std(x)
        print(ms, clipat * ms)
        print(np.abs(x-ms))
#        idx = np.where(np.abs(x0-ms) < clipat * ss) 
        idx = np.where( (np.abs(x0-ms) < clipat * ss) )
        x = x0[idx]
        nx = len(x)
#         print(idx)
#         print(n0,nx)
        if n0 == nx : 
            break
    return(ms,ss)


