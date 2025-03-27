# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:52:19 2025

Numerical methods for Pore Calculation.

@author: samda
"""
import math
import numpy as np
from matplotlib import path

from AtomicPoint import AtomicPoint


def threeAtomCircumcenter(p1, p2, p3):
    '''
    Function for calculating the center of a circle described by how it touches three non-linear atoms.
    Prof. Grimm's version of the equation, explinations in paper/file on three atom circumcenter

    Mathmactical formula for the coordinates of a circle center based on https://math.stackexchange.com/questions/3100828/calculate-the-circle-that-touches-three-other-circles    

    Parameters
    ----------
    p1 : AtomicPoint
        Point 1.
    p2 : AtomicPoint
        Point 2.
    p3 : AtomicPoint
        Point 3.

    Returns
    -------
    AtomicPoint
        Potential candidate sphere. The AtomicPoint describing the center and 
        radius of the circle defined by the 3 surrounding atoms.
    '''
    r1 = p1.VanDerWaalsRadius
    r2 = p2.VanDerWaalsRadius
    r3 = p3.VanDerWaalsRadius
    
    a = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y)
    b = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z)
    c = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
    d =  a * p1.x + b * p1.y + c * p1.z
    
    # Check added to prevent devide by 0 errors occuring due to small c --> Grimm will want print out of case
    if abs(c) < 0.000001: return AtomicPoint(0, 0, 0, '', vdwr=-2)
    
    l = p2.x - p1.x + (a/c) * p1.z - (a/c) * p2.z
    m = p2.y - p1.y + (b/c) * p1.z - (b/c) * p2.z
    n = r1 - r2
    p_ = r1**2 - r2**2 + p2.x**2 - p1.x**2 + p2.y**2 - p1.y**2 + p2.z**2 - p1.z**2
    p = 0.5*p_ + (d/c) * p1.z - (d/c) * p2.z

    q = p3.x - p1.x + (a/c) * p1.z - (a/c) * p3.z
    s = p3.y - p1.y + (b/c) * p1.z - (b/c) * p3.z
    t = r1 - r3
    w_ = r1**2 - r3**2 + p3.x**2 - p1.x**2 + p3.y**2 - p1.y**2 + p3.z**2 - p1.z**2
    w = 0.5*w_ + (d/c) * p1.z - (d/c) * p3.z
    
    A = (s*n - m*t) / (s*l - m*q)
    B_ = (s*p - m*w) / (s*l - m*q)
    B = B_ - p1.x
    C = (q*n - l*t) / (q*m - l*s)
    D_ = (q*p - l*w) / (q*m - l*s)
    D = D_ - p1.y
    E = (-b/c) * C - (a/c) * A
    F = (d/c) - (b/c) * D_ - (a/c) * B_ - p1.z
    
    if ((A*B + C*D + E*F - r1)**2 - (A**2 + C**2 + E**2 -1) * (B**2 + D**2 + F**2 - r1**2)) < 0:
        return AtomicPoint(0, 0, 0, '', vdwr=-1)
    
    #Getting a negative of a sqrt
    value = -(A*B + C*D + E*F - r1) / (A**2 + C**2 + E**2 - 1)
    uncertainty = math.sqrt((A*B + C*D + E*F - r1)**2 - (A**2 + C**2 + E**2 -1) * (B**2 + D**2 + F**2 - r1**2)) / (A**2 + C**2 + E**2 - 1)

    r_ans1 = value + uncertainty
    r_ans2 = value - uncertainty

    x_1 = A*r_ans1 + B_
    x_2 = A*r_ans2 + B_

    y_1 = C*r_ans1 + D_
    y_2 = C*r_ans2 + D_

    z_1 = (d/c) - (b/c) * y_1 - (a/c) * x_1
    z_2 = (d/c) - (b/c) * y_2 - (a/c) * x_2
    
    # Solution 1
    sol_1 = AtomicPoint(x_1, y_1, z_1) 
    sol_1.VanDerWaalsRadius = abs(r_ans1)
    rCheck1 = threeAtomRadiusCheck(p1, p2, p3, sol_1) #Checks if the circle touches the surface of each atom
    
    # Solution 2
    sol_2 = AtomicPoint(x_2, y_2, z_2) 
    sol_2.VanDerWaalsRadius = abs(r_ans2)
    rCheck2 = threeAtomRadiusCheck(p1, p2, p3, sol_2) #Checks if the circle touches the surface of each atom

    if rCheck1 and rCheck2: 
        if abs(r_ans1) < abs(r_ans2):
            return sol_1
        else:
            return sol_2
    elif rCheck1 and not rCheck2:
        return sol_1
    elif not rCheck1 and rCheck2:
        return sol_2
    else:
        if abs(r_ans1) < abs(r_ans2):
            return sol_1
        else:
            return sol_2
        

def threeAtomRadiusCheck(a, b, c, cs):
    '''
    Function takes three AtomicPoints and a candidate sphere and checks if the 
    candidate sphere touches the surface of each of the AtomicPoints.

    Parameters
    ----------
    a : AtomicPoint
        Point 1.
    b : AtomicPoint
        Point 2.
    c : AtomicPoint
        Point 3.
    cs : AtomicPoint
        Candidate sphere.

    Returns
    -------
    bool
        If the candidate sphere is touching the surface (within 0.001 Angstroms),
        return True. Else returns False.

    '''
    if (abs(cs.distBetween(a, True)) < 0.001 and
        abs(cs.distBetween(b, True)) < 0.001 and
        abs(cs.distBetween(c, True)) < 0.001):
        return True
    return False


def circleSurrounded(a, b, c, cs):
    '''
    Function checks if the center of candidate sphere cs is within the triangle
    created between AtomicPoints a, b, and c. 

    Parameters
    ----------
    a : AtomicPoint
        Point 1.
    b : AtomicPoint
        Point 2.
    c : AtomicPoint
        Point 3.
    cs : AtomicPoint
        Candidate sphere.

    Returns
    -------
    bool
        If the candidate sphere center is within the trangle return True. Else 
        returns False.
    ''' 
    p = path.Path([(a.x,a.y),(b.x,b.y),(c.x,c.y)])
    point = [cs.x,cs.y]
    
    if p.contains_point(point): return True
    else: return False


def getPlanarValues(a, b, c):
    '''
    Function takes 3 AtomicPoints and finds the 2d plane that intersects them 
    all. Returns the principle values to convert (x, y, z) points into points 
    in the new plane.

    Based on this: https://www.quora.com/What-is-the-fastest-way-to-find-the-equation-of-a-plane-given-three-points

    Parameters
    ----------
    a : AtomicPoint
        Point 1. Becomes origin of new plane.
    b : AtomicPoint
        Point 2.
    c : AtomicPoint
        Point 3.

    Returns
    -------
    list consisting of [origin, surface normal, ab, in plane normal]
        Origin is the AtomicPoint representing the origin of the new plane. The
        Surface Normal is the unit vector surface normal of the new plane. ab
        is the unit vector from a toward be (one of the in plane axes). The In
        Plane Normal is the second in plane axis.
    '''
    vab = [b.x-a.x,b.y-a.y,b.z-a.z]  #vector from a --> b
    vac = [c.x-a.x,c.y-a.y,c.z-a.z]  #vector from a --> c
    
    ab = normalize(vab) #normaized ab vector
    ac = normalize(vac) #normaized ac vector
    
    #Cross product vector between ab and ac
    surface_normal = [ab[1]*ac[2] - ab[2]*ac[1], -1*(ab[0]*ac[2] - ab[2]*ac[0]), ab[0]*ac[1] - ab[1]*ac[0]]
    sn = normalize(surface_normal)
    
    #Second planar axis along with ab, cross product of surface normal and ab
    in_plane_normal = [sn[1]*ab[2] - sn[2]*ab[1], -1*(sn[0]*ab[2] - sn[2]*ab[0]), sn[0]*ab[1] - sn[1]*ab[0]]
    ipn = normalize(in_plane_normal)
    
    #Declares planar origin = point a
    origin = AtomicPoint(a.x, a.y, a.z)
    
    return [origin, sn, ab, ipn]


def getPlanarValuesV2(a, n, u):
    '''
    Function takes origin, plane normal, and an in-plane vector to calculate 
    the principle values to convert (x, y, z) points into points in the new plane.

    Parameters
    ----------
    a : AtomicPoint
        Coordinate Origin.
    n : list of floats [x,y,z]
        Plane normal.
    u : list of floats [x,y,z]
        In-plane vector to act as one in-plane axis

    Returns
    -------
    List consisting of [origin, surface normal, ab, in plane normal]
        Origin is the AtomicPoint representing the origin of the new plane. The
        Surface Normal is the unit vector surface normal of the new plane. ab
        is the unit vector from a toward be (one of the in plane axes). The In
        Plane Normal is the second in plane axis.

    '''
    #ensure input vectors are normalized - i.e. unit vectors
    u_n = normalize(n)
    u_u = normalize(u)
    
    #calculate 3rd axis vector
    w = cross(u_n, u_u)
    u_w = normalize(w)
    
    #Planar constant --> n1 ∗ x + n2 ∗ y + n3 ∗ z + K = 0
    k = -1 * (a.x*u_n[0] + a.y*u_n[1] + a.z*u_n[2])
    
    #Declares planar origin = point a
    origin = AtomicPoint(a.x, a.y, a.z, "")

    return [origin, u_n, u_u, u_w, k]


def getReversePlanarValues(pv):
    '''
    Function constructs list of values necessary to return points converted 
    using this list of planar values to the standard basis.

    Parameters
    ----------
    pv : [Plane Origin as AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, Othoganal Axis Unit Vector 2,
                 Constant for cartesian planar equation with regards to surface normal k].
        Principle planar values of plane to convert back from.

    Returns
    -------
    rppv : [Plane Origin as AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, Othoganal Axis Unit Vector 2,
                 Constant for cartesian planar equation with regards to surface normal k].
            List of values necessary to return points converted using this list of planar values to the standard basis.
    '''

    r_origin = convertPoint(AtomicPoint(0,0,0,''), pv)
    r_normal = convertPoint(AtomicPoint(0,0,1,''), pv)
    r_Xaxis = convertPoint(AtomicPoint(1,0,0,''), pv)
    r_Yaxis = convertPoint(AtomicPoint(0,1,0,''), pv)

    r_n_vec = [r_normal.x - r_origin.x, r_normal.y - r_origin.y, r_normal.z - r_origin.z]
    u_rn = normalize(r_n_vec)

    r_Xvec = [r_Xaxis.x - r_origin.x, r_Xaxis.y - r_origin.y, r_Xaxis.z - r_origin.z]
    u_rX = normalize(r_Xvec)

    r_Yvec = [r_Yaxis.x - r_origin.x, r_Yaxis.y - r_origin.y, r_Yaxis.z - r_origin.z]
    u_rY = normalize(r_Yvec)

    return [r_origin, u_rn, u_rX, u_rY]


def convertPoint(a, pv):
    '''
    Function takes an AtomicPoint a and a set of values produced by 
    getPlanarValues, and returns the location of point a in the basis of the 
    new plane. 

    Based on this: https://stackoverflow.com/questions/23472048/projecting-3d-points-to-2d-plane

    Parameters
    ----------
    a : AtomicPoint
        Point to be converted.
    pv : List
        Set of values produced by getPlanarValues --> [Plane Origin as 
        AtomicPoint, Surface Normal Unit Vector, Axis Unit Vector 1, Othoganal 
        Axis Unit Vector 2].

    Returns
    -------
    AtomicPoint
        Location of a in the basis of the new plane.
    '''
    origin = pv[0]
    surface_normal = np.array(pv[1])
    axis1 = np.array(pv[2])
    axis2 = np.array(pv[3])
    
    origin_vector = [a.x - origin.x, a.y - origin.y, a.z - origin.z]
    
    s = np.dot(surface_normal, origin_vector)
    u = np.dot(axis1, origin_vector)
    v = np.dot(axis2, origin_vector)

    return AtomicPoint(u, v, s, vdwr=a.VanDerWaalsRadius) 


def normalize(v):
    '''
    Function normalizes a given vector and returns it. A normalized vector has
    a magnitude of 1. (Norm 2). --> Can change to us numpy function

    Parameters
    ----------
    v : List of Floats
        Vector in list form.

    Returns
    -------
    normalized_vector : List of Floats
        Unit vector in direction v in the form of a list.
    '''
    magnitude = 0
    #calculates magnitude of vector from components, runs for all components
    for component in v: magnitude += (float(component)**2)
    final_magnitude = math.sqrt(magnitude)
    #Divides each component of the vector by the vector's magnitude
    normalized_vector = [x/final_magnitude for x in v]
    return normalized_vector


def cross(a, b):
    '''
    Calculates cross product of two given vectors, each expected to have 3
    dimensions. [x,y,z] --> Can change to us numpy function

    Parameters
    ----------
    a : list of floats
        Vector 1.
    b : list of floats
        Vector 2.

    Returns
    -------
    Unit vector calculated by the cross product.
    '''
    # caluclate cross product values
    i = a[1]*b[2] - (a[2]*b[1])
    j = a[2]*b[0] - (a[0]*b[2])
    k = a[0]*b[1] - (a[1]*b[0])

    # construct new vector
    n = [i, j, k]

    # return normalized vector
    return normalize(n)