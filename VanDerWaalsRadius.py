# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 19:13:30 2023

Unless otherwise specified, values were aquired from "Consistent van der Waals 
Radii for the Whole Main Group" at https://doi.org/10.1021/jp8111556

@author: samda
"""

VanDerWaalsRadius_dict = {
    "":0,
    "H":1.10,
    "HE":1.40,
    "LI":1.81,
    "BE":1.53,
    "B":1.92,
    "C":1.70,
    "N":1.55,
    "O":1.52,
    "F":1.47,
    "NE":1.54,
    "NA":2.27,
    "MG":1.73,
    "AL":1.84,
    "SI":2.10,
    "P":1.80,
    "S":1.80,
    "CL":1.75,
    "AR":1.88,
    "K":2.75,
    "CA":2.31,
    "ZN":1.39, #Aquired from https://en.wikipedia.org/wiki/Van_der_Waals_radius
    "GA":1.87,
    "GE":2.11,
    "AS":1.85,
    "SE":1.90,
    "BR":1.83,
    "KR":2.02,
    "RB":3.03,
    "SR":2.49,
    "ZR":2.52, #Aquired from https://webelements.com/periodicity/van_der_waals_rad/
    "IN":1.93,
    "SN":2.17,
    "SB":2.06,
    "TE":2.06,
    "I":1.98,
    "XE":2.16,
    "CS":3.43,
    "BA":2.68,
    "TI":1.96,
    "PB":2.02,
    "BI":2.07,
    "PO":1.97,
    "AT":2.02,
    "RN":2.20,
    "FR":3.48,
    "RA":2.83
    }


"""
Returns the Van der Waals radius of the given element in angstroms

@param element The element as a string, all caps
@return float The Van der Waals radius of the given element in angstroms
"""
def findRadius(element):
    return VanDerWaalsRadius_dict[element]