# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:15:18 2025

This file contains the code the principally handles reading the .pdb file and
returning calculated answers as to the largest solid sphere that can pass
through the given pore each frame.

@author: samda
"""
import pandas as pd

from AtomicPoint import AtomicPoint
import PoreCalculatorFunctions as pcf
import NumericMethods as nm
import SphericalPath as sp


def PoreCalculator(file, i=1, n=150, name="PoreCalculator", noReturn=False, start = 1, **kwargs):
    '''
    Function calculates the largest solid sphere that can pass through the pore
    decribed by the .pdb file for each frame of the molecular dynamics simulation.

    Parameters
    ----------
    file : String
        File path to .pdb file being used.
    i : Int, optional
        Frame index to start on. Default = 1.
    n : Int, optional
        Number of atoms used to generate candidate spheres. Default = 150.
    name : String, optional
        Name of file to create (do not include file extension). Should allow 
        to save to a file path.
    noReturn : Boolean, optional
        If True, does not create a return excel file. Default = False
    start : Int, optional
        .pdb frame to start on. Indexed from 1. Default = i
    --- **kwargs ---    
    working_frame : Int
        Only runs calculation on given frame. Useful for debugging.
    testing : Float
        working_frame must be present. Returns path and environment atoms
        within given distance of the candidate sphere in GeoGebra notation.

    Returns
    -------
    None. Creates excel file containing information on the largest solid sphere
    that could pass through the given pore each frame. Information represented 
    as: Frame #, Radius, X, Y, Z

    '''
    frame = [] # Where atoms decribed by .pdb are held
    
    df = pd.DataFrame({"Frame" : [], "Radius" : [], "x" : [], "y" : [], "z" : []})
    
    first_line = True # Ensures first line is skipped
    with open(file) as f:
        # Continues running until all lines of file have been read. Stops upon
        # reaching "END" characters, knows it has a full .pdb list to run calculations on
        for line in f: # Gets line
            if first_line:
                first_line = False
                continue
        
            t_line = line.split()       # splits line based on white space
            
            # Reached 'END' indictor, frame fully built, time to run algorithm on
            if len(t_line) == 11:
                #5 is x coord, 6 y coord, 7 z coord, 10 element, 2 varient of element by bond type
                frame.append(AtomicPoint(float(t_line[5]), #x coord
                                               float(t_line[6]), #y coord
                                               float(t_line[7]), #z coord
                                               t_line[10],       #element
                                               bondType=t_line[2])) #bond type
            else:
                # allows me to test specific frames if I want
                if 'working_frame' in kwargs and kwargs['working_frame'] != i:
                    frame = []
                    i += 1
                    continue
                # Skips calculating frames up to indicated frame
                if i < start - 1:
                    frame = []
                    i += 1
                    continue
                
                # Get Center
                bounds = pcf.boundsFinder(frame)
                center = AtomicPoint((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2, "")
                
                # transform frame so center is located at (0,0,z)
                frame = pcf.centering(center, frame)
                
                # Gets list of candidate spheres --> generateCandidateSpheres(frame, n)
                cs_list = pcf.generateCandidateSpheres(frame)
                
                j = 0 # used for backtracking. Previous index
                
                sol = [i,0,0,0,0]
                
                # For testing
                testing_cs = None 
                substitute = [False]
                
                # Tests each candidate sphere with sphericalPath
                for cs in cs_list:
                    oppv = nm.getPlanarValues(cs[2][0], cs[2][1], cs[2][2])
                    sa = [cs[2][0], cs[2][1], cs[2][2]]
                    
                    sp_result = sp.pathGen(cs[1], oppv[1], sa, frame)
                    if sp_result[0]: 
                        sol = [i, cs[0], cs[1].x, cs[1].y, cs[1].z]
                        testing_cs = cs
                        break
                    else: # Try a mild backtrack to avoid corners
                        sub_cs = cs_list[j-1]
                        oppv = nm.getPlanarValues(sub_cs[2][0], sub_cs[2][1], sub_cs[2][2])
                        sa = [sub_cs[2][0], sub_cs[2][1], sub_cs[2][2]]
                        temp_cs = AtomicPoint(sub_cs[1].x, sub_cs[1].y, sub_cs[1].z, '', vdwr = cs[1].VanDerWaalsRadius)
                        
                        sp_result = sp.pathGen(temp_cs, oppv[1], sa, frame)
                        if sp_result[0]: 
                            sol = [i, cs[0], cs[1].x, cs[1].y, cs[1].z]
                            testing_cs = cs
                            substitute = [True, sub_cs]
                            break
                    j += 1
                
                # Revert to original positions relative to center
                temp_sol = [AtomicPoint(sol[2],sol[3],sol[4],'',vdwr=sol[1])]
                rev_sol = pcf.centering(center, temp_sol, revert=True)[0]
                sol = [sol[0], sol[1], rev_sol.x, rev_sol.y, rev_sol.z]
                
                if noReturn:
                    print(sol)
                    print("Answer is")
                    print(sol[1])
                # Immediately writes to sheet 1 and saves to prevent data lose
                else:
                    new_row = pd.DataFrame({"Frame" : [sol[0]], "Radius" : [sol[1]], "x" : [sol[2]], "y" : [sol[3]], "z" : [sol[4]]})
                    df = pd.concat([df, new_row], ignore_index=True)
                    # Write the DataFrame to an Excel file
                    df.to_excel(name + '.xlsx', index=False)
                
                # Cuts of caluclation after reaching working frame if using working frame
                if 'working_frame' in kwargs and kwargs['working_frame'] == i:
                    # For testing
                    if 'testing' in kwargs:
                        if substitute[0]:
                            oppv = nm.getPlanarValues(substitute[1][2][0], substitute[1][2][1], substitute[1][2][2])
                            sa = [substitute[1][2][0], substitute[1][2][1], substitute[1][2][2]]
                            temp_cs = AtomicPoint(substitute[1][1].x, substitute[1][1].y, substitute[1][1].z, '', vdwr = testing_cs[1].VanDerWaalsRadius)
                            sp_result = sp.pathGen(temp_cs, oppv[1], sa, frame, GeoGebraPath=temp_cs, GeoGebraEnvironment=kwargs['testing'])
                        else:
                            oppv = nm.getPlanarValues(testing_cs[2][0], testing_cs[2][1], testing_cs[2][2])
                            sa = [testing_cs[2][0], testing_cs[2][1], testing_cs[2][2]]
                            sp_result = sp.pathGen(cs[1], oppv[1], sa, frame, GeoGebraPath=cs[1], GeoGebraEnvironment=kwargs['testing'])
                    if 'cs_count' in kwargs:
                        print('cs_count information')
                        print(len(cs_list))
                        for item in cs_list: print(item[0])
                    break
                
                frame = []
                i += 1
            
        # Closes file to clean up
        f.close()