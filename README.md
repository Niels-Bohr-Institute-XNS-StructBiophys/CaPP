# CaPP 1.0
Calculating Pair distance distribution functions for Proteins.  
The program calculates the PDDF from a high-resolution protein structure in PDB format,  
and the scattering intensity can be calculated by Fourier transform of the PDDF.  
A water layer can be added as an option. 

# License
CaPP is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.          
                                                                     
CaPP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).  
                                                                     
# Citing the program  
If you use the CaPP software in your work, please cite the Github download address:                                    
CaPP, available at github.com/Niels-Bohr-Institute-XNS-StructBiophys/CaPP                                                  

# Dependencies for the GUI  
- python 2.7  
- wxpython  
- matplotlib  
python packages can easily be managed by e.g. Enthought Canopy  

# Dependencies for developers and users of Linux
- To recompile the source code, a c-compiler is needed, e.g. gcc (Linux/MacOS) or Pelles C (Windows).  

# Different platforms  
Executables are made for MacOS and Windows  
Users of other OS should:  
1) compile Mainfunction.c and call the executable "capp"  
        >> gcc Mainfunction.c -o capp  
2) place the executable, capp, in the same folder as CaPP.py  
3) Run CaPP  
        >> python CaPP_1.0.py  

CaPP has been tested on MacOS 10.12 (Sierra) and Windows Vista  

# Running the program, GUI mode
To start the GUI, type in the terminal:  
        >> python CaPP_1.0.py  

# Running the program, batch mode
Type in the terminal:  
        >> capp [options] PDBFILE.pdb  
  
Options:  
  
- c [input: Contrast of water layer]  
Add a water layer with (c)ontrast between 0 and 2 times the solvent scattering length.  
Typically 0.1.  
Default: No water layer. 
  
- d [no input]  
Only relevant for membrane proteins.  
Removes water layer from the bilayer region.  
Choose -d if the pdb is from the OPM (d)atabase, that provides the bilayer thickness.  
  
- m [input: Bilayer Thickness]  
Only relevant for membrane proteins.  
Removes water layer from the bilayer region.  
Choose -m to (m)anually provide the bilayer thickness in Aangstrom.  
Typically 30 Aangstrom.  
NB: Remember to place the TMD perpendicular to the xy-plane, in z=0!
  
- s [input: prc D20 in the solvent]  
Choose SANS contrast and enter the D2O-content (between 0 and 1) of the (s)olvent.  
SAXS contrast asssumed if option is not chosen.  
  
- r [input: Resolution of p(r) function]  
Change the (r)esolution, i.e. the binsize (in Aangstrom) of the p(r) function.  
Default: 3.0 Aangstrom.  

# About the calculations
The PDDF is calculated using the positions of each atom in the PDB file.  
X-ray scattering length are simply calculated as the number of electrons times the electron scattering length.  
The atomic form factor of each atom is approximated by the carbon atomic formfactor.  
The form factor of the excluded solvent is given as a Gaussian sphere with volume equal to the atomic volume. The Van der Waals radii are used for most atoms, except the volumes for H,C,D,N, and O, which are found for proteins experimentally by Fraser et al. (J Appl Cryst(1978), 11, p693).  
Neutron scattering length are imported from the ILL neutron data booklet, and the nucleai are assumed to be point-like, i.e. with form factor of unity.  
The hydrogens and deuteriums are included implicitely.  
Water layer is added explicitely and included in the PDDF and thus in the calculated scattering intensity.  
