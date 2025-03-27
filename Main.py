# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:04:24 2025

Below is a sample implementation of the MOF Pore Size Calculator. The principle function is
PoreCalculator.PoreCalculator(...), which can be imported as needed for your pore size calculation needs. The attached
README in the GitHub explains the main parameters of the function and what should be done if debugging is necessary.
For a full explanation of the math behind the code, please see the paper at ...

@author: samda
"""
import PoreCalculator
import time

if __name__ == "__main__":
    extension = "..." # Extension for folder .pdb files have been placed in
    fl = ["..."] # File list. Name of .pdb files to open
    nl = [f[:-4] for f in fl] # Name List. Removes the .pdb extension from file names

    for i in range(len(fl)):
        tStart = time.perf_counter()

        print('loading ' + fl[i])
        PoreCalculator.PoreCalculator(extension + fl[i], name=nl[i])

        tEnd = time.perf_counter()
        time_diff = tEnd - tStart

        print(f'The program took {time_diff:.3f} seconds to run' + fl[i])