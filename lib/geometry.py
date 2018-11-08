#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ingress Maxfield - geometry

GNU Public License
http://www.gnu.org/licenses/
Copyright(C) 2016 by
Jonathan Baker; babamots@gmail.com
Trey Wenger; tvwenger@gmail.com

Triangles and the like

original version by jpeterbaker
29 Sept 2014 - tvw V2.0 major updates
26 Feb 2016 - tvw v3.0
              merged some new stuff from jpeterbaker's new version
30 Mai 2018 - Goyet Christopher
              - convex hull with scipy
              - point in Polygon in a sphere (+sides)
              - spherePolygonArea and area of triangle in a plane
              - distance**2
              - isIntersect(A, B, S, T): True if [AB] intersect [ST]
              # TODO : replace with spheric geometry
"""

import math
import numpy as np
from scipy.spatial import ConvexHull
from itertools import combinations

def LLtoRads(pts):
    pts = pts.astype(float)
    pts *= np.pi / (180.*1.e6)
    return pts

def radstoxyz(pts,R=1):
    # Converts degree latitude/longitude to xyz coords
    # Returns corresponding n x 3 array

    pts = pts.reshape([-1,2])

    lat = pts[:,0]
    lng = pts[:,1]

    # The radius of the latitude line
    r = np.cos(lat)

    x = np.cos(lng)*r
    y = np.sin(lng)*r
    z = np.sin(lat)
    
    xyz = np.column_stack([x,y,z])
    xyz *= R
    return xyz

def xyztorads(pts,R=1):
    pts = pts.reshape([-1,3])
    pts = pts/R
    x = pts[:,0]
    y = pts[:,1]
    z = pts[:,2]

    lat = np.arcsin(z)
    lng = np.arctan2(y,x)

    return np.column_stack([lat,lng])

def greatArcAng(x,y):
    '''
    x,y should be nx2 arrays expressing latitude,longitude (in radians)
    Great arc angle between x and y (in radians)
    '''

    # If either is a single point (not in a list) return a 1-d array
    flatten = y.ndim==1 or x.ndim==1

    # Formula taken from Wikipedia, accurate for distances great and small
    x = x.reshape([-1,2])
    y = y.reshape([-1,2])
    
    nx = x.shape[0]
    ny = y.shape[0]

    # After reshaping, arithmetic operators produce distance-style matrices
    latx = np.tile(x[:,0],[ny,1])
    lngx = np.tile(x[:,1],[ny,1])
    laty = np.tile(y[:,0],[nx,1]).T
    lngy = np.tile(y[:,1],[nx,1]).T

    dlng = np.abs(lngx-lngy)

    sinx = np.sin(latx)
    cosx = np.cos(latx)

    siny = np.sin(laty)
    cosy = np.cos(laty)
    
    sind = np.sin(dlng)
    cosd = np.cos(dlng)

    numer = np.sqrt( (cosx*sind)**2 + (cosy*sinx-siny*cosx*cosd)**2 )
    denom = siny*sinx + cosy*cosx*cosd

    # great arc angle containing x and y
    angles = np.arctan2(numer,denom)

    if flatten:
        angles.shape = -1

    return angles

def sphereDist(x,y,R=6371000):
    '''
    x,y are n x 2 arrays with lattitude, longitude in radians
    '''
    sigma = greatArcAng(x,y)
    return R*sigma
    
def sphereTriAngles(vertsLL):
    '''
    vertsLL is 3 x 2 array with lattitude, longitude in radians
    ''' 
    n = len(vertsLL) # hope n = 3
    sides = [greatArcAng(vertsLL[i], vertsLL[(i+1)%n])[0] for i in range(n)]
    #print sides
    s = 0.5*np.sum(sides)
    #print "s = ",s
    sin_s = np.sin(s)
    sin_sa = np.sin(s-sides)
    sin_sb = np.roll(sin_sa, -1)
    sin_sc = np.roll(sin_sb, -1)
    numer = np.sqrt(sin_sa*sin_sb)
    denom = np.sqrt(sin_s*sin_sc)
    angles = 2*np.arctan2(numer,denom)
    #print "angles=",angles
    return angles

def sphereTriArea(vertsLL, R=6371000):
    '''
    vertsLL is 3 x 2 array with lattitude, longitude in radians
    '''
    n = len(vertsLL) # hope n = 3
    sides = [greatArcAng(vertsLL[i], vertsLL[(i+1)%n])[0] for i in range(n)]
    #print sides
    s = 0.5*np.sum(sides)
    tan_s = np.tan(0.5*s)
    tan_ds = np.tan(0.5*(s-sides))
    if (tan_s*np.prod(tan_ds)<0):
        print "tan_s*np.prod(tan_ds)=",tan_s*np.prod(tan_ds)
    # e = 4*np.arctan( np.sqrt(np.absolute(tan_s*np.prod(tan_ds))) )
    e = 4*np.arctan( np.sqrt(tan_s*np.prod(tan_ds)) )
    # print "Aire =", e*(R**2)
    # n = len(vertsLL)
    # e = np.sum(sphereTriAngles(vertsLL)) - (n-2)*np.pi
    # print "Aire =", e*(R**2)
    return int(e*(R**2))
        

def sphereTriContains(pts,x):
    '''
    pts is a 3 x 3 array representing vertices of a triangle
        pts[i] contains the x,y,z coords of vertex i
    x is a 2-array representing the test point

    points should be represented in xyz format

    returns True iff x is inside the triangle
        yes, three points make two triangles, but we assume the small one

    behavior in border cases ont guaranteed
    '''
    x = x.reshape([-1,3])
    # Find vectors orthogonal to the planes through origin and triangle sides
    crosses = np.cross( pts[[1,2,0]] , pts[[2,0,1]] )
    # crosses = [cross(1, 2), cross(2, 0), cross(0, 1) ]
    xsign = np.dot( crosses,x.T )
    psign = np.sum(crosses*pts,1).reshape([3,1])
    # Check whether opposite vertex is always on same side of plane as x
    return np.all( xsign*psign > 0,0)
    
def spherePolygonContains(pts, x, close=True):
    '''
    pts is a n x 3 array representing vertices of a polygon
        pts[i] contains the x,y,z coords of vertex i
    x contains the x,y,z coords of the test point

    points should be represented in xyz format

    returns True iff x is inside the polygon

    behavior in border cases ont guaranteed
    A points on border seems exclude ?
    '''
    
    x = x.reshape([-1,3]) #x.T
    # Vectors orthogonal to planes through origin and sides (i, i+1)
    crosses = np.cross( np.roll(pts, 0) , np.roll(pts,-1, 0) )
    # crosses = [cross(0, 1), cross(1, 2), cross(2, 3) ... cross(n, 0)]
    if not close:
        # in case of demi-plan for example
        crosses = np.delete(crosses, -1, 0) #remove last
    xsign = np.dot( crosses,x.T )
    if close:
        # we also accepte clockwise 
        return np.all( xsign > 0 , 0) or np.all( xsign < 0, 0)  
    else:
        # be careful to points sequence
        return np.all( xsign > 0 , 0)

def makeLace(n):
    # sequence of perimeter nodes to hit for a lacing-style triangulation
    # begins with the edge 1,-1
    lace = np.arange(1,n//2)
    lace = np.vstack([lace,(-lace)%n])
    lace = lace.T.reshape(-1)
    lace = list(lace)
    if n%2==1:
        lace.append(n//2)
    return lace

def rotate(x):
    # rotate the vector(s) in x by one quarter turn counter-clockwise
    if x.ndim == 1:
        x[[0,1]] = [-x[1],x[0]]
    else:
        x[:,[0,1]] = x[:,[1,0]]
        x[:,0] *= -1
        
        

def gnomonicProj(pts,ptsxyz=None):
    '''
    pts should be in lat/lng
    Uses the centroid of pts as the center, North Pole as positive y-direction
    This is only guaranteed to work if no two points are more than 90 degrees apart (great arcwise)
    This is about 9700 km across the surface of Earth
    '''
    if ptsxyz is None:
        ptsxyz = radstoxyz(pts)

    # We'll project onto the plane tangent at base
    basexyz = ptsxyz.mean(0)
    basexyz /= np.linalg.norm(basexyz)

    base = xyztorads(basexyz).reshape(-1)

    # We'll us the triangle base - portal - North Pole
    # The angles at these vertices are, respectively A - B - C
    # The corresponding lowercase letter is arc-angle of the opposite edge
    
    a = np.pi/2-pts[:,0]
    b = np.pi/2-base[0]
    c = greatArcAng(base,pts)
    C = base[1] - pts[:,1]

    # http://en.wikipedia.org/wiki/Spherical_trigonometry#Identities
    # A = arcsin[   sin(a)*sin(C)          /   sin(c)          ]
    # A = arccos[ { cos(a)-cos(c)*cos(b) } / { sin(c)*sin(b) } ]
    sinA = np.sin(a)*np.sin(C) / np.sin(c)
    cosA= (np.cos(a)-np.cos(c)*np.cos(b))/(np.sin(c)*np.sin(b))
    
    # arcsin can only fall in [-pi/2,pi/2]
    # we can find obtuse angles this way
    theta = np.arctan2(sinA,cosA)

    # Distance from base
    r = np.tan(c).reshape([-1,1])
    
    # theta measures counter-clockwise from north
    xy = np.column_stack([ -np.sin(theta) , np.cos(theta) ])*r
    return xy



def arc(a,b,c):
    '''
    Finds the arc through three points in a plane
    returns z,r,ta,tb,tc

    z = [x,y] is the center of the arc
    r is the radius of the arc
    a = z+r*[cos(ta),sin(ta)]
    b = z+r*[cos(tb),sin(tb)]
    c = z+r*[cos(tc),sin(tc)]
    '''
    
    # center points on b
    ab = a-b
    cb = c-b
    ac = a-c

    # squared lengths
    slab =  ab[0]**2+ab[1]**2
    slcb =  cb[0]**2+cb[1]**2
    # length
    lac = (ac[0]**2+ac[1]**2)**.5

    # this is from wikipedia http://en.wikipedia.org/wiki/Circumscribed_circle
    D = 2*(ab[0]*cb[1]-ab[1]*cb[0])
    z = np.array([ cb[1]*slab - ab[1]*slcb ,\
                   ab[0]*slcb - cb[0]*slab ])/D + b

    # the angle a,b,c
    t = np.abs( np.arctan(ab[1]/ab[0]) - np.arctan(cb[1]/cb[0]) )

    # the angle a,z,c is 2*t
    # the angles a,c,z and c,a,z are equal (isosolecsescscs triangle)
    # a,c,z + c,a,z + a,z,c = 180
    acz = np.pi/2-t

    # d is the midpoint of ac
    lad = lac/2 # the length of ad

    # d,c,z is a right triangle with hypoteneuse az
    # and since a,c,z = a,d,z
    r = lad/np.cos(acz)

    az = a-z
    bz = b-z
    cz = c-z
    ta = np.arctan2(az[1],az[0])
    tb = np.arctan2(bz[1],bz[0])
    tc = np.arctan2(cz[1],cz[0])
    
    return z,r,ta,tb,tc

def orthplane(xyz):
    '''
    xyz should be a 3 x 3 numpy array
    returns the vector orthogonal to the plane passing through the rows of xyz such that all( np.dot(xyz,p) > 0 )
    '''
    # treat xyz[0] as origin
    a,b,c = tuple(xyz[1]-xyz[0])
    d,e,f = tuple(xyz[2]-xyz[0])

    # cross product of other shifted vectors
    p = np.array( [b*f-c*e,
                   c*d-a*f,    
                   a*e-b*d])   
           
    return p/np.linalg.norm(p)

def commonHemisphere(xyz,getDisproof=False):
    '''
    xyz should be an n x 3 numpy array with point coordinates
    if it exists, returns (p,None)
        p is a 3-vector such that all( np.dot(xyz,p) > 0 )
        p is orthogonal to the plane through the points indexed by inds
    otherwise, returns (None,inds)
        inds are the indices of 4 non-co-hemispherical points
        inds are None if getDisproof is False (since these take extra time to compute with this implementation)

    the plane through the origin and orthogonal to p has all points of xyz on the same side
    this defines a hemisphere appropriate for gnomic projection
    '''
    n = xyz.shape[0]
    if n < 4:
        if n == 0:
            return (np.array([1,0,0]),None)
        if n == 1:
            return (xyz,None)
        if n == 2:
            return (np.mean(xyz,0),None)
        if n == 3:
            return (orthplane(xyz),None)

    for tri in combinations(xyz,3):
        p = orthplane(tri)
        if np.all(np.dot(xyz,p) > 0):
            # print np.dot(xyz,p)
            return (p,None)

    if not getDisproof:
        return (None,None)

    range1_4 = range(1,4)
    range4 = range(4)

    for quad in combinations(range(n),4):
        for j in range4:
            noj = [ quad[j-i] for i in range1_4 ]
            p = orthplane(xyz[noj])
            if np.dot(xyz[quad[j]],p) > 0:
                # The convex hull of these four don't contain the origin
                break
        else:
            # The loop exited normally
            # The current quad is a counter example
            return (None,quad)

    print xyz
    print "We shouldn't be here"

'''
********************************************************************************
    
    Plane Geometry
    DONT USE 
'''

def planeDist(x,y=None):
    x = x.reshape([-1,2])
    if y == None:
        y = x
    else:
        y = y.reshape([-1,2])

    #TODO this is not a clever way of makeing the matrix
    return np.sqrt(np.array([ [sum( (a-b)**2 ) for a in y] for b in x ]))

    
def norms(x):
    'Norm per row of x, a matrix of vectors ?'
    return np.sqrt(np.sum(x**2,1))
    
def norms2(x):
    res = 0
    try:
        res = np.sum(x**2,1)
    except Exception:
        'Norm of vector x'
        res = np.sum(x**2)
    return res
    
def distance2(a,b):
    return norms2(b-a)
    
def distance(a,b):
    return np.sqrt(distance2(a,b))
    
def planeArea2(a, b, c):
    #with determinant/2 ? 
    l1 = distance(a, b)
    l2 = distance(a, c)
    l3 = distance(b, c)
    dp = 0.5*(l1+l2+l3)
    return dp*(dp-l1)*(dp-l2)*(dp-l3)
    
def planeArea(a, b, c):
    #with determinant/2 ?
    return np.sqrt(planeArea2(a, b, c))
    
def planeTriArea(pts):
    return 0.5*abs(determinant(pts[1]-pts[0], pts[2]-pts[0]))
    
def determinant(vectA, vectB):
    return (vectA[0]*vectB[1]-vectA[1]*vectB[0])
        

# def areParalleles(A, B, S, T):
#     xA, yA = A[0], A[1]
#     xB, yB = B[0], B[1]
#     xS, yS = S[0], S[1]
#     xT, yT = T[0], T[1]
#     vectST = (xT-xS, yT-yS)
#     vectAB = (xB-xA, yB-yA)
#     res = abs(determinant(vectAB, vectST))
#     return res==0
#     #print res==0
#     #return res<1.e-8

def positif_N_over_D(n, d):#sign of n/d > 0 -> True
    if n == 0:
        return False
    if (n>0):
        return (d>0)
    else:
        return (d<0)
        
def isIntersect(A, B, S, T):
    '''
    Returns True if [AB] intersect [ST]"
    A, B, S, T are coords of points in plane
    
    This is for planar points (spherical points should be get Gnomonic projection first)
    '''
    #print A, B, S, T
    xA, yA = A[0], A[1]
    xB, yB = B[0], B[1]
    xS, yS = S[0], S[1]
    xT, yT = T[0], T[1]
    vectST = (xT-xS, yT-yS) # s
    vectAB = (xB-xA, yB-yA) # r
    d = determinant(vectAB, vectST) # r x s
    if abs(d) == 0: #if areParalleles( A, B, S, T):
        return False #TODO they may be colineaire or parallele and intersect
    vectAS = (xS-xA, yS-yA) # q-p
    # with x the 2-dimensional vector cross product
    n = determinant(vectAS, vectST) # (q-p) x s
    m = determinant(vectAS, vectAB) # (q-p) x r
    #print "n="+str(n)+"   ;    d="+str(d)+"    ;    m="+str(m)
    #print "n/d="+str(n/d)+"   ;    m/d="+str(m/d)
    if positif_N_over_D(n, d) and positif_N_over_D(m, d) and (n/d<1) and (m/d<1):
    # t = n/d so, if 0 <= t <= 1  and
    # u = m/opd, if 0 <= u <= 1
        return True
    return False
    
'''
********************************************************************************
    
    convex hull and Spliting Perimeter functions 
    
'''
# def between(a,b,pts):
#     # For use with gerperim
#     # Returns the index of a point in pts "left" of the ray a-b
#     # diff will be orthogonal to the line through a,b
#     diff = pts[a]-pts[b]
#     rotate(diff)
#     
#     # maximum inner product with diff
#     c = np.argmax(np.dot(pts,diff))
# 
#     if c == a or c == b:
#         return None
#     else:
#         return c
# 
# def getPerim(pts): #enveloppe convexe, algorithme peu efficace en n^2 (meilleur en nlog(n) )
#     '''
#     Returns a list of indices of the points on the "outside" (in the boundary of the convex hull)
#    
#     This is for planar points (spherical points should be get Gnomonic projection first)
#     '''
#     # Point with the greatest x-coordinate is certainly outside
#     hix = np.argmax(pts[:,0])
#     # Same goes for the point with the least x-coordinate
#     lox = np.argmin(pts[:,0])
# 
#     perim = {hix:lox , lox:hix}
#     perimlist = []
# 
#     a = hix
#     b = lox
#     
#     aNeverChanged = True
# 
#     while a != hix or aNeverChanged:
#         c = between(a,b,pts) 
#         #print c
#         if c is None:
#             # there is no perimeter point between a and b
#             # proceed to the next adjacent pair
#             perimlist.append(a)
#             a = b
#             b = perim[b]
# 
#             aNeverChanged = False
#         else:
#             # c is on the perimeter after a
#             # we will next look for another point between a,c
#             perim[a] = c
#             perim[c] = b
#             b = c
#             
#     hull = ConvexHull(pts)
#     print hull.vertices
#     print perimlist
#     return perimlist
def getPerim(pts):
    '''
    Returns a list of indices of the points on the "outside" (in the boundary of the convex hull)
    
    This is for planar points (spherical points should be get Gnomonic projection first)
    '''
    hull = ConvexHull(pts)
    return hull.vertices
    

    
def split_to_partition(demi_plan, inter, max=None):
    # list all possibles partitions of the perim
    if max == 0:
        return []
    if np.all([len(i) < 2 for i in inter]):
        return [demi_plan]
    part = []
    for s in range( len(demi_plan)-1 ):       
        n = len(inter[s]) if max == None else min(len(inter[s]), max)
        if n>0:
            for i in range(n):
                dp = [demi_plan[j] for j in range(s)]
                if i == 0:
                    dp += [demi_plan[s]]                
                else:
                    dp += [demi_plan[s][:-i]]
                dp += [np.append(demi_plan[s+1][:2],demi_plan[s+1][1+n-i:])]
                dp += [demi_plan[j] for j in range(s+2, len(demi_plan))]
                it = [inter[j] for j in range(s)]
                it.append([])
                it += [inter[j] for j in range(s+1, len(demi_plan))]
                part = part + split_to_partition(dp, it, max=(max-1 if max != None else None))
            return part
    s = len(demi_plan)-1
    n = len(inter[s]) if max == None else min(len(inter[s]), max)
    for i in range(n):
        dp = [ np.append(demi_plan[0][:2],demi_plan[0][1+n-i:])]
        dp += [demi_plan[j] for j in range(1, s)]
        if i == 0:
            dp += [demi_plan[s]]
        else:
            dp += [demi_plan[s][:-i]]
        part.append(dp)
    return part
    # if len(inter[0])>0:
    #     n = len(inter[0])
    #     for i in range(n):
    #         if i == 0:
    #             dp = [ demi_plan[0], np.append(demi_plan[1][:2],demi_plan[1][1+n-i:]), demi_plan[2] ]
    #         else:
    #             dp = [ demi_plan[0][:-i], np.append(demi_plan[1][:2],demi_plan[1][1+n-i:]), demi_plan[2] ]
    #         it = [[], inter[1], inter[2]]
    #         part = part + split_to_partition(dp, it)
    # elif len(inter[1])>0:
    #     n = len(inter[1])
    #     for i in range(n):
    #         if i == 0:
    #             dp = [ demi_plan[0], demi_plan[1], np.append(demi_plan[2][:2],demi_plan[2][1+n-i:]) ]
    #         else:
    #             dp = [ demi_plan[0], demi_plan[1][:-i],  np.append(demi_plan[2][:2],demi_plan[2][1+n-i:]) ]
    #         it = [inter[0], [], inter[2] ]
    #         part = part + split_to_partition(dp, it)
    # elif len(inter[2])>0:
    #     n = len(inter[2])
    #     for i in range(n):
    #         if i == 0:
    #             dp = [ np.append(demi_plan[0][:2],demi_plan[0][1+n-i:]), demi_plan[1], demi_plan[2] ]
    #         else:
    #             dp = [ np.append(demi_plan[0][:2],demi_plan[0][1+n-i:]), demi_plan[1], demi_plan[2][:-i] ]
    #         #it = [inter[0], inter[1], [] ]
    #         part.append(dp)
    # return part
    
# def split_perim(a, perim, triangle):
#     demi_plan = []
#     inter = []
#     for i in range(3):
#         new_perim = [triangle[(i+1)%3], triangle[i]]
#         new_coord = [a.node[i]['xyz'] for i in new_perim]
#         # points in new perim must be in same side than their edges
#         part_perim = np.array([p for p in perim if (not p in new_perim) and \
#                 geometry.spherePolygonContains(new_coord, a.node[p]['xyz'], close=False)], np.int32)
#         demi_plan.append( np.append(new_perim, part_perim) )
#     #  points in multiples demi-plan -> multiples split
#     for i in range(3):
#         inter.append(np.intersect1d(demi_plan[i][2:], demi_plan[(i+1)%3][2:], assume_unique=True))
#     #print "triangle, perim, inter =",triangle, demi_plan, inter
#     return split_to_partition(demi_plan, inter)
    
def split_perim(a, perim, polygon, max=None):
    #print "perim = ", perim
    rotperim = perim[:]
    n = len(polygon)
    demi_plan = []
    inter = []
    for i in range(n):
        new_perim = [polygon[(i+1)%n], polygon[i]]
        new_coord = [a.node[i]['xyz'] for i in new_perim]
        # points in new perim must be in same side than their edges
        bool_perim = []
        part_perim = []
        for p in rotperim :
            if (not p in new_perim) and \
                spherePolygonContains(new_coord, a.node[p]['xyz'], close=False):
                bool_perim.append(1)
            else:
                bool_perim.append(0)
        while not np.all([bool_perim[i] >= bool_perim[i+1] for i in range(len(bool_perim)-1)]):
            bool_perim = np.roll(bool_perim, -1)
            rotperim = np.roll(rotperim, -1)
        # part_perim = np.array([p for p in perim if (not p in new_perim) and \
        #         spherePolygonContains(new_coord, a.node[p]['xyz'], close=False)], np.int32)
        for i in range(len(bool_perim)):
            if bool_perim[i] == 0:
                break
            part_perim.append(rotperim[i])
        demi_plan.append( np.append(new_perim, np.array(part_perim,np.int32) ) )
    #  points in multiples demi-plan -> multiples split
    for i in range(n):
        inter.append(np.intersect1d(demi_plan[i][2:], demi_plan[(i+1)%n][2:], assume_unique=True))
    #print "demi_plan = ", demi_plan
    #print "inter = ", inter
    return split_to_partition(demi_plan, inter, max=max)
 

if __name__ == '__main__':
    # Test common hemisphere finder

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # xyz = np.random.randn(7,3)
    # xyz = (xyz.T/norms(xyz)).T
    xyz = np.array([[-0.30581918,-0.46686818,-0.82976426],
                    [ 0.59465481, 0.19030562, 0.78113342],
                    [ 0.8265863 ,-0.56278406,-0.00540285],
                    [-0.50141151, 0.78501969, 0.36377271],
                    [ 0.23231895, 0.90232697,-0.36308944],
                    [-0.33705904,-0.56767828, 0.75108759],
                    [-0.32538217, 0.94383169, 0.05751689]])
    # p,pts = commonHemisphere(xyz,True)

    ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],'bo')
    """
    if p == None:
        print 'disproof found'
        ax.plot([0,xyz[pts[0],0]],[0,xyz[pts[0],1]],[0,xyz[pts[0],2]],'bo-')
        ax.plot([0,xyz[pts[1],0]],[0,xyz[pts[1],1]],[0,xyz[pts[1],2]],'ko-')
        ax.plot([0,xyz[pts[2],0]],[0,xyz[pts[2],1]],[0,xyz[pts[2],2]],'go-')
        ax.plot([0,xyz[pts[3],0]],[0,xyz[pts[3],1]],[0,xyz[pts[3],2]],'ro-')
    else:
        print 'plane found'
        print np.dot(xyz,p)
        ax.plot([0,p[0]],[0,p[1]],[0,p[2]],'ko-')
        ax.plot([0],[0],[0],'r*')
    ax.plot(xyz[pts,0],xyz[pts,1],xyz[pts,2],'r*')
   """
        
    plt.show()
    plt.close('all')

