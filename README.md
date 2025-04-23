# MOF_Pore_Size_Calculator

The purpose of this algorithm is to calculate the largest solid sphere that can pass through a given pore for each frame of a molecular dynamics simulation (MDS) of the pore. The principle function is PoreCalculator(), which is contained in the PoreCalculator.py file. The description of that function is given below:

    PoreCalculator(file, i=1, n=150, name="PoreCalculator", noReturn=False, start = 1, **kwargs):
    
    Function calculates the largest solid sphere that can pass through the pore
    described by the .pdb file for each frame of the molecular dynamics simulation.

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
        Useful for debugging.

    Returns
    -------
    None. Creates excel file containing information on the largest solid sphere
    that could pass through the given pore each frame. Information represented 
    as: Frame #, Radius, X, Y, Z


A more extensive tutorial on how to use this code in the form of a Word Document titled 'How to Use MOF Pore Size Calculator' is located in Auxiliary Files, along with a sample .pdb file to test the code with titled 'Sample_MOF-5_TPAA.pdb'. The radius produced by applying the Pore Size Calculator algorithm to this file should be 3.45958416925372.
