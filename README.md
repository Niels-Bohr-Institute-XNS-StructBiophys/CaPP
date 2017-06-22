# CaPP 1.0
Calculating Pair distance distribution functions for Proteins.  
The program calculates the PDDF from a high-resolution protein structure in PDB format,  
and the scattering intensity can be calculated by Fourier transform of the PDDF.  
A water layer can be added as an option. 

# License
CaPP is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.          
                                                                     
CaPP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).                        
                                                                     
# Citing the program
If you use the CaPP software in your work, please cite:                                         
CaPP, github.com/Niels-Bohr-Institute-XNS-StructBiophys/CaPP                                                  

# Dependencies
- python 2.7  
- wxpython  
- to recompile the source code c-compiler is needed, e.g. gcc  

# Running the program
To start the GUI, type in the terminal:  
        >> python CaPP_1.0.py  
the GUI is self-explaining  

# Different platforms
Executables are made for MacOS and Windows  
Users of other OS should:  
1) compile Mainfunction.c and call the executable "capp"  
        >> gcc Mainfunction.c -o capp
2) place the executable, capp, in the same folder as CaPP.py  
3) Run CaPP  
        >> python CaPP.py 

# On the calculations
The PDDF is calculated using the positions of each atom in the PDB file.  
X-ray scattering length are simply calculated as the number of electrons times the electron scattering length.  
The atomic form factor of each atom is approximated by the carbon atomic formfactor.  
The form factor of the excluded solvent is given as a Gaussian sphere with volume equal to the atomic volume. The Van der Waals radii are used for most atoms, except the volumes for H,C,D,N, and O, which are found for proteins experimentally by Fraser et al. (J Appl Cryst(1978), 11, p693).  
Neutron scattering length are imported from the ILL neutron data booklet, and the nucleai are assumed to be point-like, i.e. with form factor of unity.  
The hydrogens and deuteriums are included implicitely.  
Water layer is added explicitely and included in the PDDF and thus in the calculated scattering intensity.  
