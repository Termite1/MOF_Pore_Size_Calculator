# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:57:16 2025

Code used to generate a spherical path from a candidate sphere and see if its 
path can connect the top and bottom of the pore.

@author: samda
"""
import numpy as np
import math

from AtomicPoint import AtomicPoint
import NumericMethods as nm
import PoreCalculatorFunctions as pcf

def pathGen(A, n, sa, epl, **kwargs):
    '''Given a starting point (a) tries to generate a spherical path to the top and bottom of the environment.

    Parameters
    ----------
    A : AtomicPoint
        Starting point of the spherical path.
    n : List of floats [x,y,z]
        Normal to circle defining starting point location. (extendCylinder direction of travel for 1st step, afterwords
        continues in +z or -z direction)
    sa : List of 3 AtomicPoints
        Atoms surrounding candidate sphere initial position from circle generation. Only used in 1st step.
    epl : List of AtomicPoints
        Environment being traversed.
    **kwargs
        plot : if == true, creates plots before return showing locations along spherical path
        
        print : if == true, prints out information
        
        tolerance : extendCylinder tolerance = to given value
        
        pdb : if == true, prints out path_plotting as list of atomic points for addition to VMD file
        
        GeoGebraPath : if == Candidate Sphere, prints out path_plotting as list of GeoGebra coordinates
        
        GeoGebraEnvironment :
            Must have a value for GeoGebraPath. Value float, returns in 
            GeoGebra notation all environmental atoms withing given distance of
            GeoGebraPath center.

    Returns [boolean, condition]
    -------
    True if path to top and bottom of environment could be found
    False if no path could be found
    
    Condition = 0 if no further action needs to be taken
    If [False, 1], error has occured, special return given so these answers/frames can be identified
    '''
    # Could generate normal from list of surounding points 
    
    # Make sure initial move is in theright direction
    if n[2] < 0: n = [-1 * n[0], -1 * n[1], -1 * n[2]]
      
    up = [0,0,1] # direction of travel upwards
    down = [0,0,-1] # direction of travel downwards
    dt_list = [up,down]
    
    # Used to generate list of path locations
    path_plotting = []
    path_plotting.append(A)
    
    # Tries to move to top then bottom of environment
    for dt in dt_list:
        #print('dt')
        running = True
        
        temp_n = [0,0,0] # Used to get correct initial travel direction
        if dt[2] == 1: temp_n = n
        elif dt[2] == -1: temp_n = [num * -1 for num in n]
        
        B = None
        if 'tolerance' in kwargs:
            B = extendCylinder(A, temp_n, epl, sa, tolerance=kwargs['tolerance'])
        else: B = extendCylinder(A, temp_n, epl, sa)
    
        # If extendCylinder initially is intersecting with an atom, terminate this run
        if B[0] == -1: 
            # print("Initial collision error")
            return [False,1]
        
        C = [-1,AtomicPoint(0,0,0,'')]

        # Prevent infinite loop
        max_cycles = 200 # could make a kwargs argument
        current_cycles = 0
        
        while running:
            #print('cycle')
            if current_cycles > max_cycles: 
                # print('Ran out of pathGen cycles')
                return [False,0]
            current_cycles = current_cycles + 1
            
            if B[0] == 0: # cylinder extended without issue, end cycle move to next dt or signal end
                #print('Valid run')
                running = False
            elif B[0] == 1: # cylinder interrupted
                push_point(B[1], path_plotting, dt)

                cs = B[1] #candidate sphere location
                ri = B[2] #recent interrupt
                
                # Helper function to deal with recursive orbits
                if 'print' in kwargs: C = orbiting(cs, ri, dt, epl, path_plotting, print = True)
                elif 'GeoGebraOut' in kwargs: C = orbiting(cs, ri, dt, epl, path_plotting, GeoGebraOut=kwargs['GeoGebraOut'])
                else: C = orbiting(cs, ri, dt, epl, path_plotting)
                
                if C[0] == -1:  # Error return -->
                    #if 'plot' in kwargs: test.plotSphericalPath(epl, path_plotting)
                    # print("pathGen error return")
                    return [False,1]
                elif C[0] == -2: # Termination return
                    #if 'plot' in kwargs: test.plotSphericalPath(epl, path_plotting)
                    #print('pathGen termination return')
                    return [False,0]
                elif C[0] == -3: # Ran out of orbiting cycles
                    return [False,0]
            
                i_list = [] 
                
                # During extend cylinder, ignores orbit center of unobstructed oneAtomOrbits. 
                # 90 degrees to angle of travel, only can collide if too close to starting location, shouldn't be.
                if C[0] == 1: i_list = [C[2]]
                
                # twoAtomOrbit linear set return
                if C[0] == 3: 
                    i_list = [C[3]]
                    B = extendCylinder(C[1], C[2], epl, i_list)
                elif 'tolerance' in kwargs:
                    B = extendCylinder(C[1], dt, epl, i_list, tolerance=kwargs['tolerance'])
                else: B = extendCylinder(C[1], dt, epl, i_list)
    
    if 'plot' in kwargs:
        if kwargs['plot'] == True:
            pass
            #test.plotSphericalPath(epl, path_plotting)
    if 'pdb' in kwargs:
        pass
        #pdb.pdbWrite(path_plotting, 1)
    if 'GeoGebraPath' in kwargs:
        for item in path_plotting: 
            item.printCoords(GeoGebraOffset=kwargs['GeoGebraPath'])
        # Also print environment spheres within given istance
        if 'GeoGebraEnvironment' in kwargs:
            for item in epl: 
                if kwargs['GeoGebraPath'].distBetween(item, False) < kwargs['GeoGebraEnvironment']:
                    item.printCoords(GeoGebraOffset=kwargs['GeoGebraPath'])
        
    return [True,0]
    
    
def orbiting(cs, ri, dt, epl, path_plotting, debug = False, **kwargs):
    ''' Function manages orbiting loop in pathGen, continues until orbit functions return no collision

    Parameters
    ----------
    cs : AtomicPoint
        Where orbiting behavior starts from. Begins touching ri
    ri : AtomicPoint
        Recent interrupt, where cs ran into environmental point ri as it was trying to travel through the environment
        space. Atom cs initially orbits
    dt : List of floats [x,y,z]
        Direction of travel. Will be [0,0,1] or [0,0,-1] for the moment
    epl : List of AtomicPoints
        Environment being traversed. List of environmental points.
    path_plotting : List of AtomicPoints
        Used for plotting candidate sphere travel paths through.
    debug : Boolean, optional
        Prints some debug info if true. Default = False
    **kwargs
        if print=True, prints values upon termination
        if GeoGebraOut=float, print atoms within float distance in GeoGebra Notation

    Returns
    -------
    [int, AtomicPoint, **]
        int = -1 if error
        int = -2 if terminated
        int = 1 if from oneAtomOrbit --> ** present as last collided atom
        int = 2 if from twoAtomOrbit
    '''
    col_list = [ri, None, None] # Used to remember most recent collisions
    t_count = 0 #Termination counter, should cause exit when reaches 3
    
    b = oneAtomOrbit(ri, cs, epl, dt)
    push_point(b[1], path_plotting, dt)  # For kwargs evaluating spherical path trajectory
    
    max_cycles = 200 # could make a kwargs argument
    current_cycles = 0

    while (True): # Runs until return value is generated
        # Prevents calc from getting trapped in an infinite loop
        if current_cycles > max_cycles: 
            if debug: print("ran out of orbiting cycles")
            return [-3, AtomicPoint(0,0,0,'',vdwr=-3)]
        current_cycles = current_cycles + 1
        
        if b[0] == -1: # Error
            if debug: print('orbiting error')
            return [-1, AtomicPoint(0,0,0,'',vdwr=-1)]
        elif b[0] == 0: # No interrupt
            push_point(b[1], path_plotting, dt) # For kwargs evaluating spherical path trajectory
            if b[2] == 1: return [1, b[1], col_list[0]]
            elif b[2] == 2: return [2, b[1]]
            # Location rotated to returned, along with code if from 1 atom orbit or two atom orbit
        elif b[0] == 1: # Interrupt
            new_col = b[2] # Most recent collision
            new_cs = b[1]  # Where candidate sphere ended up after most recent orbit

            # Looking for trap, should be indicated by 3 consecutive increases in t_count, will need to test
            if new_col in col_list: t_count += 1
            else: t_count = 0
            # Terminate current run
            if t_count == 3:
                if 'print' in kwargs:
                    print("Termination Conditions: Point at")
                    new_cs.printCoords()
                    print("Blocking atoms locations and distance:")
                    for atom in col_list:
                        atom.printCoords()
                        print(atom.distBetween(new_cs, True))
                if 'GeoGebraOut' in kwargs:
                    print(dt[2])
                    new_cs.printCoords(GeoGebraOffset=new_cs)
                    temp_val = kwargs['GeoGebraOut']
                    for P in epl:
                        if P.distBetween(b[1], False) < temp_val:
                            P.printCoords(GeoGebraOffset=new_cs)
                return [-2, AtomicPoint(0,0,0,'',vdwr=-2)]
            
            shift(new_col, col_list) # Add collision point to col_list, moves all points in list >> and discards last
            
            c1 = oneAtomOrbit(new_col, new_cs, epl, dt)
            c2 = twoAtomOrbit(col_list[0], col_list[1], new_cs, epl, dt)

            # error checking 
            if c1[0] == -1: 
                if debug: print('oneAtomOrbit no answer error')
                return [-1, AtomicPoint(0,0,0,'',vdwr=-1)]
                #b = c2 # If other is error, caught at top of while loop
            elif c2[0] == -1: 
                if debug: print('twoAtomOrbit no answer error')
                return [-1, AtomicPoint(0,0,0,'',vdwr=-2)]
                #b = c1
            # Check if twoAtomOrbit returns a linear set, requires a special extendCylinder
            elif c2[0] == 2:
                # 3 is an idicated, c2[1] candidate sphere location, c2[2] extendCylinder direction of travel, c2[3] ignore list for extendCylinder
                return [3, c2[1], c2[2], c2[3]]
            # Regular height checking, value that moved furthest in desired z direction kept
            elif dt[2] == 1: # Traveling upward +z
                if c1[1].z < c2[1].z: b = c2
                else: b = c1
            elif dt[2] == -1: # Traveling downward -z
                if c1[1].z > c2[1].z: b = c2
                else: b = c1

            push_point(b[1], path_plotting, dt) # For kwargs evaluating spherical path trajectory


def extendCylinder(A, n, epl, il, **kwargs):
    '''Function calculates how far a cylinder could be extended from the 
    given starting point (a) in a given direction of travel (n) before coming 
    into contact with a sphere in the environment. Returns given sphere's new 
    location where it contacts the interrupting environmental sphere.
    
    Parameters
    ----------
    A : AtomicPoint
        Starting location for calculation
    n : List of floats [x,y,z]
        Unit vector direction of travel
    epl : List of AtomicPoints
        Environmental point list. All the spheres that could interrupt the 
        extention of the cylinder from a in n direction.
    il : List of AtomicPoints
        List of AtomicPoints to ignore during the cylinder extention calculation.
    **kwargs
        tolerance : float
            Set tolerance factor, used to increase minimum distance between point need to register as colliding with
            extended cylinder, default = 0
    
    Returns
    -------
    [-1] : Error
    [0] : No interrupt
        Cylinder extended without interruption in given direction
        of travel.
    [1, AtomicPoint closest_point, AtomicPoint closest_interrupt] : Interrupt
        Cylinder extention interrupted, AtomicPoint closest_point
        where given sphere a contacts interrupting sphere. AtomicPoint
        closest_interrupt interrupting sphere.
    '''
    #print('extendCylinder')
    # Check for collision at starting location - should be redundant, may be able to remove later
    # Keeping for now as sanity check. If collide, return error
    if atomCollide(A, epl): 
        #print('extendCylinder initial collision error')
        return [-1]
    
    # Sets tolerance factor - Currently unused
    tolerance = 0
    if 'tolerance' in kwargs: tolerance = kwargs['tolerance']
    
    closest_point = AtomicPoint(0, 0, 0, '')            # Tracks closest collision
    closest_interrupt = AtomicPoint(0, 0, 0, '')        # Tracks closest interrupting atom
    closest_dist = 1000                                 # Tracks distance to closest collision
    
    u_n = nm.normalize(n) # Makes sure directional vector is unit vector
    
    # Check every environmental atom for potential collisions with cylindrical path
    for P in epl:
        # Ignore values in ignor list
        skip_check = False
        for ignore in il: #May want to make ignore list optional. Can you run for loop on empty list?
            if ignore == P: 
                skip_check = True
        if skip_check: continue
        
        vec_AP = [P.x - A.x, P.y - A.y, P.z - A.z]      # Vector from A to P
        d = np.dot(np.array(u_n), np.array(vec_AP))                              # Distance from A to B
        
        # Point on the center line of the cylinder closest to P
        B = AtomicPoint(A.x + d*u_n[0], A.y + d*u_n[1], A.z + d*u_n[2], '', vdwr=A.VanDerWaalsRadius)   
        
        # Ensure B is in the direction of travel n relative to A. Otherwise discard and move on
        vec_AB = [B.x - A.x, B.y - A.y, B.z - A.z]      # Vector from A to B
        counter = 0
        for i in range(len(vec_AB)): # Check each component of AB against n
            if sign(u_n[i]) != sign(vec_AB[i]): counter += 1
            elif u_n[i] == 0 and vec_AB[i] == 0: counter += 1
        if counter == 3: continue
    
        vec_BP = [P.x - B.x, P.y - B.y, P.z - B.z]      # Vector from B to P
        mag_BP = magnitude(vec_BP)                      # Magnitude of vector BP
        r = A.VanDerWaalsRadius + P.VanDerWaalsRadius   # If P is closer to B than this distance, it intersects the cylinder
    
        # Check for intersection. If not, continue
        if mag_BP >= r - tolerance: continue
    
        bq = math.sqrt(r**2 - mag_BP**2)                # Distance from B to Q
        # Where candidate sphere touches B
        Q = AtomicPoint(B.x - bq*u_n[0], B.y - bq*u_n[1], B.z - bq*u_n[2], '', vdwr=A.VanDerWaalsRadius)
        mag_AQ = A.distBetween(Q, False)                # Distance from A to Q
    
        # Checks if Q is the closest collision. Keeps best values
        if mag_AQ < closest_dist:
            closest_point = Q
            closest_interrupt = P
            closest_dist = mag_AQ

    if closest_dist == 1000: return[0]                  # No collisions occured
    else: return[1, closest_point, closest_interrupt]   # Collision occured
    
    
def oneAtomOrbit(A, B, epl, dt):
    '''Function orbits candidate sphere B around obstruction A until B lies
    at the same z as sphere A, or until the first obstruction B touches as
    travels along the orbit. Returns [-1] if there is an error, [0, B2] if
    there is no obstruction, [1, B2, co] if there is an obstruction. B2 is
    location B was calculated to orbit to and co is the obstructing sphere 
    the orbiting B came into contact with. 
    
    Direction of travel
    
    Parameters
    ----------
    A : AtomicPoint
        Sphere B is being rotated around
    B : AtomicPoint
        Sphere B's initial position
    epl : List of AtomicPoints
        Environmental point list. All the spheres that could 
        interrupt the orbit of B about A.
    dt : List of floats [x,y,z]
        Unit vector or desired direction of travel. Currently limited to
        [0,0,1] or [0,0,-1], the +z and -z directions of travel.
    
    Returns (1 at end of list indicates oneAtomOrbit
    -------
    [-1, 1] : Error
    [0, B2, 1] : No interrupt
        B orbited a without interruption. B2 is AtomicPoint of location
        B orbited to.
    [1, B2, co, 1] : Interrupt
        B orbited with interruption. B2 is AtomicPoint of location B orbited
        to when it contacted with closest obstruction co. co AtomicPoint
        the interrupted B's orbit.
    '''    
    #print("oneAtomOrbit")
    nrl = -0.175 # Negative radian limit. Returns more negative than this ignorned. Accepts returns with up to ~-10 degrees. Want to catch angle pushback, but not excessivly
    
    r_tot = A.VanDerWaalsRadius + B.VanDerWaalsRadius
    u_AB = nm.normalize([B.x - A.x, B.y - A.y, B.z - A.z])
    
    # Distance modifications so B2 is in dt direction of travel at r_tot distance from A
    x_mod = dt[0] * r_tot
    y_mod = dt[1] * r_tot
    z_mod = dt[2] * r_tot
    
    # Where calculation is atempting to rotate to
    B2 = AtomicPoint(A.x + x_mod, A.y + y_mod, A.z + z_mod, '', vdwr=B.VanDerWaalsRadius)
    
    n = nm.cross(u_AB, dt) # Surface normal
    oppv = nm.getPlanarValuesV2(A, n, u_AB) # planar values
    rpv = nm.getReversePlanarValues(oppv)   # values to reverse planar transformation
    
    o_A = AtomicPoint(0, 0, 0, A.element, vdwr = A.VanDerWaalsRadius)               # origin in orbit bases
    
    o_B = nm.convertPoint(B, oppv) # B and B2 in orbit basis
    o_B2 = nm.convertPoint(B2, oppv)
    
    val = [] # valid angle list
    ial = [] # invalid angle list
    ccl = [] # close collision list --> used in frequently. 
             # Helps avoid errors where there would be no valid answer due to candidate sphere
             # touching environment atoms at initial location
    
    # Add angle to B2 to valid angle list 
    # (only if B2 does not collide with the environment, otherwise add to invalid angle list)
    # used to see if there is a free travel path
    angle_B2 = angle(o_B, o_B2)
    if atomCollide(B2, epl): ial.append(angle_B2)
    else: val.append([angle_B2])
    
    # Check collisions against evironmental atoms
    for P in epl:
        if (P.x == A.x and P.y == A.y and P.z == A.z): continue  # Skips checking atom being orbited
        
        o_P = nm.convertPoint(P, oppv) # Atom in orbit basis
        
        if o_P.x == 0 and o_P.y == 0: # Shouldn't happen, but in theory could occur.
            # Should have a form of error handling   
            print('oneAtomOrbit exact match error')
            return [-1, 1]
        
        u_AP_star = nm.normalize([o_P.x, o_P.y])
        # Closest point on orbit to P
        F = AtomicPoint(r_tot * u_AP_star[0], r_tot * u_AP_star[1], 0,'', vdwr=B.VanDerWaalsRadius) 
        # Should just grab angle to F, will need Ref and B2 value. Otherwise will recalculate in orbitCollisionLocation
        
        # If false, P close enough to interrupt orbit
        if F.distBetween(o_P, True) >= 0: continue
        
        C = orbitCollisionLocation(o_A, o_B, o_P, o_B2)     # result from orbit collision location, location B can rotate to before colliding with P
        
        # Makes sure there is a single system to test orcbitCollisionLocation returns. Either 1 or 2 segments returned in list
        for ts in C:
            # Catches errors from orbitCollision
            if ts[1].VanDerWaalsRadius == -1: 
                print('oneAtomOrbit orbitCollision error')
                return [-1,1]
            
            test_point = nm.convertPoint(ts[1], rpv)         # location B can rotate to before colliding with P in standard basis
    
            # Colliding atoms will also inhibit answers innately, add to list as well
            ial.append(ts[2])
    
            # Check if answer collides with any other atoms
            if atomCollide(test_point, epl) == True:            
                # If intersects environment, note angle as invalid
                ial.append(ts[0])
                if abs(ts[0]) < 0.01:
                    #print('close collision values')
                    # If having a near collision at the starting location, save in close collision list
                    if abs(P.distBetween(B, True)) < 0.001: ccl.append(P)
                continue
        
            # Other, don't add points to list whose angle is larger than the B2 angle, since B2 shoudl be the maximum return
            # or points with angles less than the negative radian limit, see top of function (currently ~ -10 degrees)            
            if ts[0] < nrl: continue
            elif ts[0] > angle_B2: continue
        
            # Makes sure colliding atom included in list, then adds point to valid angle list
            ts.append(P)
            val.append(ts)
        
    # Sorts lists in decending order
    val = pcf.sortList(val)
    ial.sort(reverse=True)

    sia = None          # Smallest invalid angle
    
    # Gets the smallest invalid angle greater than 0 --> Should I let 0 be an invalid angle? Sure, this does that, only breaks if less than 0
    for ia in ial:
        if ia < 0: break
        sia = ia
    
    # If there are no invalid angles, then B2 can be the only return
    if sia == None: return [0, B2, 1]
    
    p_sol = None        # Potential solution
    
    # Selects the answer with the largest angle smaller than the smallest invalid angle
    for va in val:
        # + 0.00001 is there to prevent rounding errors, particularly for 1 location returns. 
        # Small enough that it shouldn't mess with regular solution returns, since answers must be valid to begin with
        if va[0] >= sia + 0.00001: continue      
        p_sol = va
        break
    
    if p_sol == None: # Error, ending calc for now
        if len(ccl) > 0: # If there were close collisions an no other valid anwsers -->
            # should check if environment sphere is in orbit path, but
            # there being no other valid answers should already place at least
            # on environmental sphere in the orbit path
            dt_dist = 0
            es = ccl[0]
            # Checks every environment sphere in ccl
            for P in ccl: 
                # Gets environmnet sphere furthest in direction of travel (only directions of travel right now +z and -z)
                if (P.z - B.z) * dt[2] > dt_dist:
                    dt_dist = (P.z - B.z) * dt[2]
                    es = P
            # returns keeping B location and giving best environment sphere
            return [1, B, es, 1]
        print("oneAtomOrbit : No answer!")
        #A.printCoords(GeoGebraOffset=A)
        #B.printCoords(GeoGebraOffset=A)
        #B2.printCoords(GeoGebraOffset=A)
        #geoGebraOut(A, B, epl, oppv, r_tot, offset = -0.1)
        #print('ial')
        #for item in ial:
        #   if abs(item) < 0.5: print(item)
        return [-1, 1]
            
    # If return is B2 angle, return B2
    if len(p_sol) == 1: # Only 1 entry will have length 1, the angle of B2
        return [0, B2, 1]
    
    # Needs to evaluate angle to colliding atom if there are multiple 0 angle returns. 
    # In theory would be the same if there were multiple returns at any angle,
    # But at 0 degrees is really the only place I think I need to be concerned about it
    if p_sol[0] == 0:
        # Get item index
        ps_index = val.index(p_sol)
        
        angle0_list = [] # Holds all 0 angle returns
        
        # make list of every entry with 0 degree return 
        while val[ps_index][0] == 0:
            angle0_list.append(val[ps_index])
            ps_index += 1
            if ps_index >= len(val): break  # Dont't want to go out of list index range

        # Choose answer with largest angle_F, need to sort by index 2 --> don't need to sort just get largest answer?
        angle0_list = pcf.sortList(angle0_list, sort_index=2)
        # largest index at top of list
        p_sol = angle0_list[0]
    
    r_point = nm.convertPoint(p_sol[1], rpv)    # Returns answer to standard basis
    return [1, r_point, p_sol[3], 1]
       

def twoAtomOrbit(A, B, C1, epl, dt): # Might need to improve this method/logic
    '''Calculates how far a candidate sphere can orbit around pair of environmental spheres it is touching.

    Parameters
    ----------
    A : AtomicPoint
        One of the spheres being rotated around
    B : AtomicPoint
        Second sphere being rotated around
    C1 : AtomicPoint
        Initial position of candidate sphere orbiting the two environmental spheres
    epl : List of AtomicPoints
        Environmental point list. All the spheres that could
        interrupt the orbit of B about A.
    dt : List of floats [x,y,z]
        Unit vector or desired direction of travel. Currently limited to
        [0,0,1] or [0,0,-1], the +z and -z directions of travel.


    Returns
    -------
    [-1, 2] : Error
    [0, C2, 2] : No interrupt
        C1 orbited A and B without interruption. C2 is AtomicPoint of location
        C2 orbited to.
    [1, C2, co, 2] : Interrupt
        C1 orbited with interruption. C2 is AtomicPoint of location C orbited
        to when it contacted with closest obstruction co. co AtomicPoint
        the interrupted C1's orbit.
    [2, C1, extCy_dt, [A,B], 2] : Linear set
        A, B, and C1 situated in a line. extCy_dt unit vector in direction 
        extendCylinder should be applied in.
    '''
    
    #print("twoAtomOrbit")
    nrl = -0.175 # Negative radian limit. Returns more negative than this ignorned. Accepts returns with up to ~-10 degrees. Want to catch angle pushback, but not excessivly
    
    # Heron's formula
    a = A.VanDerWaalsRadius + C1.VanDerWaalsRadius
    b = B.VanDerWaalsRadius + C1.VanDerWaalsRadius
    c = A.distBetween(B, False)
    s = (a + b + c) / 2
    
    # Used to check for linear arangement of atom and orbit center. 
    # Should double check distances to see that path unobscured, then use extendCylinder 
    # in the optimal twoAtomTarget direction
    if (abs(a - b - c) < 0.001 or abs(c - a - b) < 0.001 or abs(b - a - c) < 0.001):
        # Need an in-plane vector to define oppv --> cross product of AB with any non-parallel vector should give an in-plane vector
        vec_AB = [B.x - A.x, B.y - A.y, B.z - A.z]  # Vector from A to B
        u_AB = nm.normalize(vec_AB)                    # Unit vector from A to B
        
        # Just need a vector that is not parallel to u_AB 
        dummy_vec = [u_AB[2], -1 * u_AB[0], u_AB[1]]
        
        # Arbitrary planar values
        oppv = nm.getPlanarValuesV2(C1, u_AB, dummy_vec)
        rpv = nm.getReversePlanarValues(oppv)
        
        # Target direction of travel parallel to u_AB
        target_list = twoAtomTarget(u_AB, dt, 1, C1.VanDerWaalsRadius, oppv, rpv)
        C2 = target_list[0]
        
        # Desired direction of travel for cylinder
        u_CC = nm.normalize([C2.x - C1.x, C2.y - C1.y, C2.z - C1.z])
        
        return [2, C1, u_CC, [A, B], 2]
    
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    orbit_r = (area * 2) / c                    # Orbit radius
    
    d = math.sqrt(a**2 - orbit_r**2)            # Distance from A to orbit center O
    vec_AB = [B.x - A.x, B.y - A.y, B.z - A.z]  # Vector from A to B
    u_AB = nm.normalize(vec_AB)                    # Unit vector from A to B
    O = AtomicPoint(A.x + d * u_AB[0], A.y + d * u_AB[1], A.z + d * u_AB[2], '', vdwr = orbit_r - C1.VanDerWaalsRadius) # Orbit center, radius may be negative right now
    vec_OC1 = [C1.x - O.x, C1.y - O.y, C1.z - O.z] # Vector from O to C1
    
    O_side = False
    
    # Check for which side of A orbit origin O lies on (shouldn't have tolerance issues like dot product, still confirm with dot product)
    if twoAtomOSide(A, B, O, C1, d):
        O = AtomicPoint(A.x - d * u_AB[0], A.y - d * u_AB[1], A.z - d * u_AB[2], '', vdwr = orbit_r - C1.VanDerWaalsRadius) # Orbit center, radius may be negative right now
        vec_OC1 = [C1.x - O.x, C1.y - O.y, C1.z - O.z] # Vector from O to C1
        O_side = True
    
    # a and b not accurate enough, recalculate the correct values and O location
    if 0.0001 < abs(np.dot(np.array(vec_AB), np.array(vec_OC1))) < 0.1 or O.distBetween(C1, False) < 0.1:
        a = A.distBetween(C1, False)
        b = B.distBetween(C1, False)
        s = (a + b + c) / 2
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        orbit_r = (area * 2) / c 
        d = math.sqrt(a**2 - orbit_r**2)            # Distance from A to orbit center O
        if O_side: O = AtomicPoint(A.x - d * u_AB[0], A.y - d * u_AB[1], A.z - d * u_AB[2], '', vdwr = orbit_r - C1.VanDerWaalsRadius)
        else: O = AtomicPoint(A.x + d * u_AB[0], A.y + d * u_AB[1], A.z + d * u_AB[2], '', vdwr = orbit_r - C1.VanDerWaalsRadius) # Orbit center, radius may be negative right now
        vec_OC1 = [C1.x - O.x, C1.y - O.y, C1.z - O.z] # Vector from O to C1
    
    # Check if vec_OC1 perpendicular to vec_AB
    if abs(np.dot(np.array(vec_AB), np.array(vec_OC1))) > 0.1: 
        print('twoAtom dot error')
        #print(dot(vec_AB, vec_OC1))
        #A.printCoords(GeoGebraOffset=O)
        #B.printCoords(GeoGebraOffset=O)
        #O.printCoords(GeoGebraOffset=O)
        #C1.printCoords(GeoGebraOffset=O)
        return [-1, 2]
    
    oppv = nm.getPlanarValuesV2(O, u_AB, nm.normalize(vec_OC1)) # orbit principle planar values, used to convert points to orbit basis
    rpv = nm.getReversePlanarValues(oppv)                    # values to return points from orbit basis to standard basis
    
    # [target location in standard basis, target location in orbit basis]
    target_list = twoAtomTarget(u_AB, dt, orbit_r, C1.VanDerWaalsRadius, oppv, rpv)
    C2 = target_list[0]
    
    o_O = AtomicPoint(0, 0, 0, '', vdwr=O.VanDerWaalsRadius)               # origin in orbit bases
    o_C2 = target_list[1]
    o_C1 = nm.convertPoint(C1, oppv)    # Starting location in orbit basis
    
    val = [] # valid angle list
    ial = [] # invalid angle list
    ccl = [] # close collision list --> used in frequently. 
             # Helps avoid errors where there would be no valid answer due to candidate sphere
             # touching environment atoms at initial location
    
    # Add angle to C2 to valid angle list 
    # (only if C2 does not collide with the environment, otherwise add to invalid angle list)
    # used to see if there is a free travel path
    angle_C2 = angle(o_C1, o_C2)
    if atomCollide(C2, epl): 
        ial.append(angle_C2)
    else: val.append([angle_C2])

    # Check collisions against evironmental atoms
    for P in epl:
        if (P.x == A.x and P.y == A.y and P.z == A.z) or (P.x == B.x and P.y == B.y and P.z == B.z): continue  # Skips checking atom being orbited
        
        o_P = nm.convertPoint(P, oppv) # Atom in orbit basis
        
        if o_P.x == 0 and o_P.y == 0: # Shouldn't happen, but in theory could occur.
            # Should have a form of error handling    
            print('twoAtomOrbit exact match error')
            return [-1, 2]
        
        u_AP_star = nm.normalize([o_P.x, o_P.y])
        # Closest point on orbit to P
        F = AtomicPoint(orbit_r * u_AP_star[0], orbit_r * u_AP_star[1], 0,'', vdwr=C1.VanDerWaalsRadius) 
        
        # If false, P close enough to interrupt orbit
        if F.distBetween(o_P, True) >= 0: continue
        
        D = orbitCollisionLocation(o_O, o_C1, o_P, o_C2)     # result from orbit collision location, location B can rotate to before colliding with P
        #test_point = proc.planarizePoint(D[1], rpv)         # location B can rotate to before colliding with P in standard basis
    
        for ts in D:
            # Catches errors from orbitCollision
            if ts[1].VanDerWaalsRadius == -1: 
                print('twoAtomOrbit orbitCollision error')
                return [-1,2]
            
            test_point = nm.convertPoint(ts[1], rpv)         # location C can rotate to before colliding with P in standard basis            
            
            # Colliding atoms will also inhibit answers innately, add to list as well
            ial.append(ts[2])
            
            # Check if answer collides with any other atoms
            if atomCollide(test_point, epl) == True:            
                # If intersects environment, note angle as invalid    
                ial.append(ts[0])
                if abs(ts[0]) < 0.01:
                    #print('close collision values')
                    # If having a near collision at the starting location, save in close collision list
                    if abs(P.distBetween(C1, True)) < 0.001: ccl.append(P)
                continue
            
            # Other, don't add points to list whose angle is larger than the C2 angle, since C2 should be the maximum return
            # or points with angles less than the negative radian limit, see top of function (currently ~ -10 degrees)
            if ts[0] < nrl: continue
            elif ts[0] > angle_C2: continue
        
            # Makes sure colliding atom included in list, then adds point to valid angle list
            ts.append(P)
            val.append(ts)  
        
    # Sorts lists in decending order
    val = pcf.sortList(val)
    ial.sort(reverse=True)

    sia = None          # Smallest invalid angle
    
    # Gets the smallest invalid angle greater than 0 --> Should I let 0 be an invalid angle? Sure, this does that, only breaks if less than 0
    for ia in ial:
        if ia < 0: break
        sia = ia
    
    # If there are no invalid angles, then C2 can be the only return
    if sia == None: 
        #print("C2 angle return - 1")  
        return [0, C2, 2]
    
    p_sol = None        # Potential solution
    
    # Selects the answer with the largest angle smaller than the smallest invalid angle
    for va in val:
        # + 0.00001 is there to prevent rounding errors, particularly for 1 location returns. 
        # Small enough that it shouldn't mess with regular solution returns, since answers must be valid to begin with
        if va[0] >= sia + 0.00001: continue
        p_sol = va
        break

    if p_sol == None: # Error, ending calc for now
        if len(ccl) > 0: # If there were close collisions an no other valid anwsers -->
            # should check if environment sphere is in orbit path, but
            # there being no other valid answers should already place at least
            # on environmental sphere in the orbit path
            dt_dist = 0
            es = ccl[0]
            # Checks every environment sphere in ccl
            for P in ccl: 
                # Gets environmnet sphere furthest in direction of travel (only directions of travel right now +z and -z)
                if (P.z - C1.z) * dt[2] > dt_dist:
                    dt_dist = (P.z - C1.z) * dt[2]
                    es = P
            # returns keeping C1 location and giving best environment sphere
            return [1, C1, es, 2]
    
        print("twoAtomOrbit : No answer!")
        #print(atomCollide(O, epl))
        #print('ial')
        #for item in ial: 
        #    if abs(item) < 0.5: print(item)
        return [-1, 2]
    
    # If return is C2 angle, return C2
    if len(p_sol) == 1: # Only 1 entry will have length 1, the angle of C2
        #print("C2 angle return - 2")    
        return [0, C2, 2]
    
    # Needs to evaluate angle to colliding atom if there are multiple 0 angle returns. 
    # In theory would be the same if there were multiple returns at any angle,
    # But at 0 degrees is really the only place I think I need to be concerned about it
    if p_sol[0] == 0:
        # Get item index
        ps_index = val.index(p_sol)
        
        angle0_list = [] # Holds all 0 angle returns
        
        # make list of every entry with 0 degree return 
        while val[ps_index][0] == 0:
            angle0_list.append(val[ps_index])
            ps_index += 1
            if ps_index >= len(val): break  # Dont't want to go out of list index range
        
        # Choose answer with largest angle_F, need to sort by index 2 --> don't need to sort just get largest answer?
        angle0_list = pcf.sortList(angle0_list, sort_index=2)
        # largest index at top of list
        p_sol = angle0_list[0]
    
    r_point = nm.convertPoint(p_sol[1], rpv)  # Returns answer to standard basis
    
    return [1, r_point, p_sol[3], 2]
    

def orbitCollisionLocation(A, B1, C, B2, **kwargs):
    '''Function finds location an orbiting candidate sphere
    contacts an interrupting environmental sphere.

    Parameters
    ----------
    A : AtomicPoint
        Orbit origin.
    B1 : AtomicPoint
        Orbiting candidate sphere starting location, generally referred to as o_B.
    C : AtomicPoint
        Sphere interrupting orbit. Projection into orbit plane.
    B2 : AtomicPoint
        Location attempting to be rotated to.
    kwargs
        small_return : if present, garentees return with smallest angle with respect to starting location
        print : prints some values

    Returns
    -------
    [[angle, inter, angle_C]] or [[angle, inter, angle_C] x 2]
        angle : float
            Angle to interrupting sphere from B1
        inter : AtomicPoint
            AtomicPoint on orbit in orbit arc contacting interrupting sphere C.
        angle_C : float
            Angle to colliding atom from B1
    '''    
    # May be worth just putting values here and put determination/selection entirely in oneAtom or twoAtom orbit
    P = AtomicPoint(C.x, C.y, 0, '', vdwr = C.VanDerWaalsRadius) # Projection of C onto orbit plane
    
    # B1 coords (+, 0, 0). Ref allows me to determine which side points are on 
    # so I can label angles as + or -. Side with B2 is +. Check relative to ref 
    # by seeing if > or < than pi / 2
    Ref = AtomicPoint(0, 1, 0,'') 
    q = 0       # Used to keep track of if angles are to be labled as + or -
    if angle(Ref, B2) < math.pi / 2: q = 1
    elif angle(Ref, B2) > math.pi / 2: q = -1
    else:
        print("180 degree orbit")
        #print(angle(Ref, B2))
        #Ref.printCoords()
        #B1.printCoords()
        #B2.printCoords()
        return [[0, AtomicPoint(0, 0, 0, '', vdwr=-1), 0]]
        #quit() Handle this!
        # Should have some sort of error handling here in theory
    
    a = P.distBetween2D(A, False)           # Distance between A and P
    b = abs(C.z)                            # Distance from C to orbit plane
    
    if A.VanDerWaalsRadius < 0: c = B1.distBetween(A, False)             # c is the orbit radius
    else: 
        #c = B1.distBetween(A, False)
        c = A.VanDerWaalsRadius + B1.VanDerWaalsRadius

    # B1.distBetween(C, False)
    d = B1.VanDerWaalsRadius + C.VanDerWaalsRadius          # What will be the distance between answer Bi and C, distance between centers if spheres touching
    e = math.sqrt(d**2 - b**2)                              # Distance between Bi and P
    
    # Find F, closest point on orbit to C. Angle to collision atom used by oneAtomOrbit and twoAtomOrbit

    u_AP = nm.normalize([P.x, P.y])
    F = AtomicPoint(c * u_AP[0], c * u_AP[1], 0, '')
    angle_F = 0
    if angle(Ref, F) < math.pi / 2: angle_F = q * angle(B1, F)
    elif angle(Ref, F) > math.pi / 2: angle_F = -1 * q * angle(B1, F)
    else: angle_F = angle(B1, F) # Only when angle is 0 or pi, don't really need sign then
    
    
    
    # Heron's Formula - Finding two possible solutions from this point on
    s = (a + c + e) / 2
    
    # Check if A, Bi, and P in a line. If so only 1 point can be calculated
    # if : 
    if s * (s - a) * (s - c) * (s - e) < 0 or (abs(a - c - e) < 0.001 or abs(c - a - e) < 0.001 or abs(e - a - c) < 0.001):
        u_AP = nm.normalize([P.x, P.y])  # Unit vector in direction from A to P
        # Location of orbit 'collision'
        I = AtomicPoint(c * u_AP[0], c * u_AP[1], 0, '', vdwr=B1.VanDerWaalsRadius)
        
        # Check if distances align, if not, correct position
        if abs(I.distBetween(C, True)) > 0.001:
            I = AtomicPoint(P.x - e * u_AP[0], P.y - e * u_AP[1], 0, '', vdwr=B1.VanDerWaalsRadius)
        
        angle_B1I = 0
        # Gets angle to point, makes sure sign is correct as well relative to B2's positon
        if angle(Ref, I) < math.pi / 2: angle_B1I = q * angle(B1, I)
        elif angle(Ref, I) > math.pi / 2: angle_B1I = -1 * q * angle(B1, I)
        else: angle_B1I = angle(B1, I)  # Only when angle is 0 or pi, don't really need sign then
        
        # [[Angle to collision, orbit collision location, Angle to colliding atom]]
        return [[angle_B1I, I, angle_F]]
        
    area = math.sqrt(s * (s - a) * (s - c) * (s - e))
    h = (2 * area) / a              # Distance between solutions Bi and T
    
    # Initialize these
    Bi_1 = None
    Bi_2 = None
    
    # Check needed. Occationaly angle between c and a will be near 90 degrees, must use different method to find Bi
    if h > c: 
        h = (2 * area) / e
        f = math.sqrt(c**2 - h**2)
        f_ = e - f
        s_ = (h + f_ + a) / 2
        area_ = math.sqrt(s_ * (s_ - h) * (s_ - f_) * (s_ - a))
        h_ = (area_ * 2) / a
        g = math.sqrt(h**2 - h_**2)
        
        u_AP = nm.normalize([P.x, P.y])
        T_ = AtomicPoint(g * u_AP[0], g * u_AP[1], 0, '')
        p_AP = nm.normalize([P.y, -1 * P.x]) # Vector perpendicular to u_AP
        T1 = AtomicPoint(T_.x + h_ * p_AP[0], T_.y + h_ * p_AP[1], 0)
        T2 = AtomicPoint(T_.x - h_ * p_AP[0], T_.y - h_ * p_AP[1], 0)
        
        p_AT = nm.normalize([T_.y, -1 * T_.x]) # Unit vector perpendicular to AT
        
        Bi_11 = AtomicPoint(T1.x + f * p_AT[0], T1.y + f * p_AT[1], 0)
        Bi_12 = AtomicPoint(T1.x - f * p_AT[0], T1.y - f * p_AT[1], 0)
        Bi_21 = AtomicPoint(T2.x + f * p_AT[0], T2.y + f * p_AT[1], 0)
        Bi_22 = AtomicPoint(T2.x - f * p_AT[0], T2.y - f * p_AT[1], 0)
        
        if abs(Bi_11.distBetween(C, True)) < 0.001: Bi_1 = Bi_11
        else: Bi_1 = Bi_12
        if abs(Bi_21.distBetween(C, True)) < 0.001: Bi_2 = Bi_21
        else: Bi_2 = Bi_22
        
    else:
        mag_a = P.distBetween(A, False)
        vec_I = [(h / mag_a) * P.y, -1 * (h / mag_a) * P.x]     # Vector from T to one Bi
        
        f = math.sqrt(c**2 - h**2)
        u_AP = nm.normalize([P.x, P.y])
        T = AtomicPoint(f * u_AP[0], f * u_AP[1], 0, '')        # Point from which Bi's are found
        
        # Check to make sure T is located in correct direction, if not correct position - Confirm this?
        if (abs(T.distBetween(C, False) - math.sqrt(d**2 - h**2)) > 0.001):
            T = AtomicPoint(-1 * f * u_AP[0], -1 * f * u_AP[1], 0, '')
        
        # The two potential solutions
        Bi_1 = AtomicPoint(T.x + vec_I[0], T.y + vec_I[1], 0, '', vdwr=B1.VanDerWaalsRadius)
        Bi_2 = AtomicPoint(T.x - vec_I[0], T.y - vec_I[1], 0, '', vdwr=B1.VanDerWaalsRadius)
    
    angle_B1I1 = 0
    # Gets angle to point, makes sure sign is correct as well relative to B2's positon
    if angle(Ref, Bi_1) < math.pi / 2: angle_B1I1 = q * angle(B1, Bi_1)
    elif angle(Ref, Bi_1) > math.pi / 2: angle_B1I1 = -1 * q * angle(B1, Bi_1)
    else: angle_B1I1 = angle(B1, Bi_1)  # Only when angle is 0 or pi, don't really need sign then
    
    angle_B1I2 = 0
    # Gets angle to point, makes sure sign is correct as well relative to B2's positon
    if angle(Ref, Bi_2) < math.pi / 2: angle_B1I2 = q * angle(B1, Bi_2)
    elif angle(Ref, Bi_2) > math.pi / 2: angle_B1I2 = -1 * q * angle(B1, Bi_2)
    else: angle_B1I2 = angle(B1, Bi_2)  # Only when angle is 0 or pi, don't really need sign then
    
    if 'print' in kwargs:
        print('set')
    
    # Perform Angle checks here then decide
    # Or return everything and decide in orbit functions?
    # [[Angle to collision, orbit collision location, Angle to colliding atom] x 2]
    return [[angle_B1I1, Bi_1, angle_F], [angle_B1I2, Bi_2, angle_F]]
        

def twoAtomTarget(n, dt, o_r, cs_r, oppv, rpv, **kwargs):
    '''Helper function that calculates the location twoAtomOrbit wants to 
    rotate to.
    
    Parameters
    ----------
    n : List of floats
        Orbit normal
    dt : List of floats [0,0,1] or [0,0,-1] for now
        Desired direction of travel
    o_r : float
        Orbit radius
    cs_r : float
        Candidate sphere radius
    oppv : List
        Principle planar values from lcfProcessing. Used to convert points to
        orbit basis.
    rpv : List
        Priciple planar values to revert points from orbit basis to standard
        basis.
    kwargs 
        for_testing : list for testing
            [A, B]
    
    Returns
    -------
    [AtomicPoint in standard basis, AtomicPoint in orbit basis]
    '''
    DT = AtomicPoint(oppv[0].x + dt[0], oppv[0].y + dt[1], oppv[0].z + dt[2], '') # Direction of travel as a unit vector from the orbit origin
    o_dt = nm.convertPoint(DT, oppv) # direction of travel as a unit vector from the origin in the orbit basis

    theta = math.atan(o_dt.y / o_dt.x)   # Angle to C2 from 0 degrees (+x direction) in orbit basis
    
    x1 = o_r * math.cos(theta)
    y1 = o_r * math.sin(theta)

    # atan only goes up to pi, need posibilites on both sides
    x2 = -1 * o_r * math.cos(theta)
    y2 = -1 * o_r * math.sin(theta)
    
    # Probably don't need to check having pi, but want to cover my bases
    x3 = o_r * math.cos(theta + math.pi)  
    y3 = o_r * math.sin(theta + math.pi)

    x4 = -1 * o_r * math.cos(theta + math.pi)
    y4 = -1 * o_r * math.sin(theta + math.pi)
    
    # Test point in orbit basis
    o_T1 = AtomicPoint(x1, y1, 0, '', vdwr=cs_r)
    o_T2 = AtomicPoint(x2, y2, 0, '', vdwr=cs_r)
    o_T3 = AtomicPoint(x3, y3, 0, '', vdwr=cs_r)
    o_T4 = AtomicPoint(x4, y4, 0, '', vdwr=cs_r)
    
    o_T_list = [o_T1, o_T2, o_T3, o_T4]
    
    # Convert points to standard basis and get z height
    T_list = [nm.convertPoint(o_T, rpv) for o_T in o_T_list]
    
    Tz_list = [T.z for T in T_list]
    
    if 'for_testing' in kwargs:
        print('2at for_testing')
        for item in T_list: 
            print(item.distBetween(kwargs['for_testing'][0], True))
            print(item.distBetween(kwargs['for_testing'][1], True))
        print("orbit dist")
        for item in o_T_list:
            print(item.distBetween(AtomicPoint(0,0,0,''), False))
    
    
    # Take one furthest in desired direction of travel
    ret_index = None    # Return index
    if dt[2] == 1: ret_index = Tz_list.index(max(Tz_list))
    elif dt[2] == -1: ret_index = Tz_list.index(min(Tz_list))
    else: print("Error! Invalid dt provided to twoAtomTarget : SphericalPathV2")
    
    # Returns answer in orbit and standard basis
    return [T_list[ret_index], o_T_list[ret_index]]
    

def twoAtomOSide(A, B, O, C1, d):
    '''Helper function to test which side of A O lies on.

    Parameters
    ----------
    A : AtomicPoint
        First atom defining orbit
    B : AtomicPoint
        Second atom defining orbit
    O : AtomicPoint
        Potential Orbit Origin location
    C1 : AtomicPoint
        Candidate Sphere starting location
    d : float
        Distance from A to C1
        
    Returns
    -------
    True if O on opposite side of A
    False otherwise
    '''
    #vec_BA = [A.x - B.x, A.y - B.y, A.z - B.z]      # Vector from B to A
    #u_BA = nm.normalize(vec_BA)                        # Unit vector from B toward A
    #vec_BC1 = [C1.x - B.x, C1.y - B.y, C1.z - B.z]  # Vector from B to C1
    vec_OC1 = [C1.x - O.x, C1.y - O.y, C1.z - O.z]
    m = magnitude(vec_OC1)
    test_C1 = AtomicPoint(O.x + m * vec_OC1[0], O.y + m * vec_OC1[1], O.z + m * vec_OC1[2])

    # Calcuated location should align with C1, otherwise O in wrong direction with respect to A
    if test_C1.distBetween(C1, False) > 0.001: return True
    else: return False

    #d = np.dot(np.array(vec_BC1), np.array(u_BA))                          # Distance from B to theoretical origin location
    
    # If d larger than the magnitude of vec_BA, O lies on the far side of A rather than between A and B
    #if d > magnitude(vec_BA): return True
    #return False
    

def magnitude(a):
    '''
    Calculates magnitude of given vector and returns it.
    
    Parameters
    ----------
    a : List of floats
        Vector whose magnitude is to be calculated.
    
    Returns
    -------
    mag_tot : float
        Magnitude of given vector.
    '''
    mag_tot = 0
    for n in a: mag_tot += n**2 # Runs through each number in vector a
    return math.sqrt(mag_tot)


def angle(a, b):
    '''
    Calculates angle between given vectors (or AtomicPoints which are converted into vectors).

    Parameters
    ----------
    a : list of floats [x,y,z]
        Vector 1.
    b : list of floats [x,y,z]
        Vector 2.

    Returns
    -------
    Float (0 to pi) : theta angle between the vectors.
    '''
    if type(a) == AtomicPoint: a = [a.x,a.y,a.z]
    if type(b) == AtomicPoint: b = [b.x,b.y,b.z]
    
    a = nm.normalize(a)
    b = nm.normalize(b)
    
    equal = True;
    opposite = True;
    # Checks if vectors are the same or opposite
    for i in range(len(a)):
        if(abs(a[i] - b[i]) > 0.0001): equal = False
        if(abs(a[i] + b[i]) > 0.0001): opposite = False
    
    # If vectors are the same
    if (equal): return 0
    # If vectors are opposite
    if (opposite): return math.pi
    
    return math.acos(np.dot(np.array(a), np.array(b))/(magnitude(a) * magnitude(b)))


def sign(val):
    '''
    Returns sign of given value.

    Parameters
    ----------
    val : float
        Given value

    Returns
    -------
    -1 is val < 0, 1 if val > 0, 0 if val == 0
    '''
    if val == 0: return 0
    elif val > 0: return 1
    else: return -1


def shift(val, l):
    ''' Function shifts each value in the list down one, kicking off the last value of the list

    Parameters
    ----------
    val : Any
        Value to be added to front of list.
    l : list of any
        List being modified.
    '''
    place_val = val
    temp_val = None

    # Runs through each item in list
    for i in range(len(l)):
        temp_val = l[i]
        l[i] = place_val
        place_val = temp_val


def push_point(val, l, dt):
    '''Helps with assembling path_plotting list. Adds to front of list if dt positive, adds to end if dt negative.

    val : AtomicPoint
        AtomicPoint to add to path_lotting
    l : list of AtomicPoints
        List val is to be added to
    dt : list of floats [x,y,z]
        Desired direction of travel
    '''
    if dt[2] > 0: l.insert(0, val)
    else: l.append(val)
    
    
def atomCollide(a, pointList, **kwargs):
    '''
    Function checks if given AtomicPoint intersects with any other atoms in the
    pointList. Returns True if it does so, false otherwise.
    
    @peram a AtomicPoint being checked
    @peram pointList List of AtomicPoints that p is being chekd against
    @return True if p itersects with any of the atoms in the pointList, returns 
    false otherwise
    '''
    for p in pointList:
        if a.distBetween(p, True) < -0.001: 
            if 'pv' in kwargs: nm.convertPoint(p, kwargs['pv']).printCoords()
            if 'pd' in kwargs: print(a.distBetween(p, True))
            return True
    return False


def geoGebraOut(a, cs, epl, pv, r, **kwargs):
    '''
    Function prints list of atoms within a certain distance of the orbit in
    GeoGebra Sphere notation

    Parameters
    ----------
    a : AtomicPoint
        Orbit orgin.
    cs : AtomicPoint
        Candidate sphere.
    epl : List of AtomicPoints
        Environmental point list.
    pv : List
        Principle planar values.
    r : Float
        Orbit radius.
    kwargs
        offset : Float
            Modifies collision calc distance

    Returns
    -------
    None.

    '''
    #Tests each atom in epl
    for p in epl:
        # Converts plane into orbit coordinate system
        o_p = nm.convertPoint(p, pv)
        
        a_vec = [o_p.x, o_p.y]
        a_nvec = nm.normalize(a_vec)
        
        #Find nearest point on orbit to o_p 
        proj = AtomicPoint(r * a_nvec[0], r * a_nvec[1], 0, '', vdwr=cs.VanDerWaalsRadius)
        
        dist = 0
        # If a collision would occur
        if 'offset' in kwargs: dist = o_p.distBetween(proj, True) + kwargs['offset']
        else: dist = o_p.distBetween(proj, True)
        
        # Print atom in GeoGebra form
        if dist < 0:
            p.printCoords(GeoGebraOffset = a)