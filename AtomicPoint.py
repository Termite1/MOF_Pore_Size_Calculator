# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 19:10:37 2023

@author: samda
"""

import math
import VanDerWaalsRadius as vdwrf

class AtomicPoint:    
    
    
    """
    
    
    @param x 
    @param y 
    @param z 
    @param element 
    """
    def __init__(self, x, y, z, element='', **kwargs):
        '''
        Constructor for the AtomicPoint class.

        Parameters
        ----------
        x : float
            X coordinate of atomic point.
        y : float
            Y coordinate of atomic point.
        z : float
            Z coordinate of atomic point.
        element : String
            String atomic symbol of element, given in all caps. Default is '' 
            which represents no element (default VanDerWaalsRadius = 0).
        --- **kwargs --- 
        bondType : String
            Additional information about relative location of atom within 
            molecular structure, represented in .pdb file by column 3.
        vdwr : Float
            Allows for manual setting of VanDerWaalsRadius 

        Returns
        -------
        Constructed AtomicPoint.
        '''
        self.x = x
        self.y = y
        self.z = z
        self.element = element
        self.VanDerWaalsRadius = vdwrf.findRadius(element)
        
        if 'bondType' in kwargs:
            self.bondType = kwargs['bondType']
            if self.bondType == '':
                self.bondType = 'H_'
        else: 
            self.bondType = 'H_'
            
        if 'vdwr' in kwargs:
            self.VanDerWaalsRadius = kwargs['vdwr'] 


    def distBetween(self, p2, includeVanDerWaalRadius):
        '''
        Calculates the euclidean distance between two given AtomicPoints.

        Parameters
        ----------
        p2 : AtomicPoint
            Second point distance is to be caluclated to.
        includeVanDerWaalRadius : Boolean
            If True, subtracts Van Der Waals Radius of both points from 
            distance. If False, returns distance as is.

        Returns
        -------
        Float euclidean distance between the two AtomicPoints.
        '''
        temp_x = self.x - p2.x 
        temp_y = self.y - p2.y
        temp_z = self.z - p2.z
        if(includeVanDerWaalRadius): return (math.hypot(temp_x,temp_y,temp_z) - (self.VanDerWaalsRadius + p2.VanDerWaalsRadius))
        else: return math.hypot(temp_x,temp_y,temp_z)
        

    def distBetween2D(self, p2, includeVanDerWaalRadius):
        '''
        Calculates the euclidean distance in the (x,y) plane between two given 
        AtomicPoints.

        Parameters
        ----------
        p2 : AtomicPoint
            Second point distance is to be caluclated to.
        includeVanDerWaalRadius : Boolean
            If True, subtracts Van Der Waals Radius of both points from 
            distance. If False, returns distance as is.

        Returns
        -------
        Float euclidean distance in (x,y) plane between the two AtomicPoints.
        '''
        temp_x = self.x - p2.x 
        temp_y = self.y - p2.y
        if(includeVanDerWaalRadius):
            return (math.hypot(temp_x,temp_y) - (self.VanDerWaalsRadius + p2.VanDerWaalsRadius))
        else:
            return math.hypot(temp_x,temp_y)
        
    """
    Prints coordinates of an AtomicPoint.
    
    @param self The AtomicPoint
    """    
    def printCoords(self, **kwargs):
        '''
        Prints coordinates and Van Der Waals Radius of AtomicPoint to 3 decimal 
        places.

        Parameters
        ----------
        --- **kwargs ---
        full : Boolean
            If in **kwargs, does not limit decimal places of printed value.
        GeoGeobra : Boolean
            If in **kwargs, prints data in comma seperated format with no 
            decimal limitation. For use in copying data from console to a 
            GeoGebra spreadsheet for easy visualization. Use E1=(A1,B1,C1) to 
            generate points, then F1=Sphere(E1,D1) to generate visualization in
            GeoGebra spreadsheet.
        GeoGebraOffset : AtomicPoint
            Same as above, but adjusts the printed coordinates to be with 
            respect to the given AtomicPoint, which will be located at (0,0,0).

        Returns
        -------
        None.

        '''
        if 'full' in kwargs: print(f'({self.x},{self.y},{self.z}) with vdwr = {self.VanDerWaalsRadius}')
        elif 'GeoGebra' in kwargs: print(f'{self.x},{self.y},{self.z},{self.VanDerWaalsRadius}')
        elif 'GeoGebraOffset' in kwargs: 
            x = self.x - kwargs['GeoGebraOffset'].x
            y = self.y - kwargs['GeoGebraOffset'].y
            z = self.z - kwargs['GeoGebraOffset'].z
            print(f'{x},{y},{z},{self.VanDerWaalsRadius}')
        else: print(f'({self.x:.3f},{self.y:.3f},{self.z:.3f}) with vdwr = {self.VanDerWaalsRadius:.3f}')