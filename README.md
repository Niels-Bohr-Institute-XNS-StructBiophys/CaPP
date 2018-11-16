# CaPP 3.9
Calculating Pair distance distribution functions (PDDF) for Proteins.  
The program calculates the PDDF from a high-resolution protein structure in PDB format,  
and the scattering intensity can be calculated by Fourier transform of the PDDF.  
A water layer can be added as an option. 

## Install and run the program

### Installation
Download python files (CaPP_A.B.py, where A.B is the version), the sourcefile and the executable and place it all in a folder.  
Check that the denpendencies are properly installed (see below)  

### Running the program, GUI mode
To start the GUI, type in the terminal

        >> python path-to-folder/CaPP_A.B.py  

### Running the program, batch mode
Type in the terminal  

        >> ./path-to-folder/capp_mac [options] PDBFILE.pdb  (mac OS)

        >> ./path-to-folder/capp_windows [options] PDBFILE.pdb  (windows)

        >> ./path-to-folder/capp [options] PDBFILE.pdb  (linux)

##### Options  

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

- p [input: perdeuteration of protein]
Choose SANS contrast and enter the (p)erdeuteration (between 0 and 1) of the protein. 

- r [input: Resolution of p(r) function]  
Change the (r)esolution, i.e. the binsize (in Aangstrom) of the p(r) function.  
Default: 1.0 Aangstrom.  
Too small binsize will not an oscillating p(r), and too large will give wrong results for P(q).  

##### Example:  

        >> capp -c 0.1 -m 30 -s 0.4 -r 2.0 4u2p.pdb  

## Dependencies

### Dependencies for the GUI  
- python 2.7  
- wxpython 2.9
- matplotlib  
- scipy 0.17 or newer (0.17 have bounds option in the function scipy.optimize.least_squares)

## Platforms

### Tested OS 
CaPP_A.B (latest version) has been tested on   
- MacOS 10.12, compiled with gcc, and with python from Enthought Canopy 2.7 (64-bit)  
- May work on Windows, but without option of plotting

##### Install on Windows (with cygwin terminal) type:  
    $ wget.exe http://peak.telecommunity.com/dist/ez_setup.py  
    $ python ez_setup.py  
    $ /cygdrive/c/Python27/Scripts/easy_install-2.7.exe pip  
    $ /cygdrive/c/Python27/Scripts/pip.exe install wxpython  
    $ /cygdrive/c/Python27/Scripts/pip.exe install matplotlib 
    § /cygdrive/c/Python27/Scripts/pip.exe install scipy

### Developers and other OS
- To recompile the source code, a c-compiler is needed,  
e.g. gcc (Linux/MacOS) or Pelles C (Windows).  
  
Executables have been made exclusively for MacOS  
Users of other OS should:  
1) compile MainFunction.c (with flag -lm) and name the executable "capp"  

        >> gcc MainFunction.c -o capp  -lm

2) place the executable, capp, in the same folder as CaPP.py  
3) Run CaPP  

        >> python CaPP_A.B.py  

## Output
assuming ABC.pdb as input pdb file  

##### ABC_w.pdb  
A new pdb with water beads added as lines to ABC.pdb (after CONNECT, before END), e.g. 

    ATOM  26485  Q   WAT W1497     -24.394   2.153 -28.888  1.00 10.00           Q  
    
##### ABC_w_only.pdb  
A pdb file with the water bead lines only (for visualization, e.g. in PyMOL)  

##### ABC_pr.dat/ABC_w_pr.dat
A data files with the PDDF for the pdb file without/with water layer, columns:  
1. r     = pair distance  
2. P(r)  = PDDF = sum ( dBa(i) * dBa(j) ) for i \neq j  
3. G(r)  = sum ( Ba(i) * Ba(j) )  
4. H(r)  = sum ( Ba(i) * Ba(j) )  
5. J(r)  = sum ( Ba(i) * Ba(j) )  
6. K(r)  = sum ( Ba(i) * Ba(j) )  

where  
Ba[i]: scattering length of atom i  
Bs[i]: scattering length of solvent, excluded by atom i  
dB[i]: excess scattering length of atom i  
columns 3-6 are used to calculate the form factor P(q)  

##### ABC_Pq  
A datafile with the form factor for the protein, columns:  
1. q, scattering vector  
2. P(q), form factor  
3. A00(q)^2, the square of the partial amplitude with l=m=0 (Svergun et al, 1995, J. Appl. Cryst., 28, 768-773)  
4. Beta, factor used in the decoupling approximation, given by Beta = A00(q)^2/P(q)  (Hoiberg-Nielsen2009, 2009, Biophys. J., 97, 1445-1453)

## About the calculations
The PDDF is calculated using the positions of each atom in the PDB file.  

X-ray scattering length are calculated as the number of electrons times the electron scattering length. 

Neutron scattering length are imported from the ILL neutron data booklet, and the nucleai are assumed to be  
point-like, i.e. with form factor of unity.  

The hydrogens and deuteriums are included implicitely. This is done for the 20 natural amino acids and the  
hetero atoms currently included in the HETATM library (MainFunction.c). If a hetero atom in the PDB does not  
exist in the library, no H's or D's will be added. The H's and D's are added according to the Protein Data  
Base (www4.rcsb.org/ligand/).  

### Calculating the form factor, P(q) from the p(r) fucntion
The atomic form factor of each atom is approximated by the carbon atomic formfactor. 

The form factor of the excluded solvent is given as a Gaussian sphere with volume equal to the atomic volume.  
The Van der Waals radii are used for most atoms, except the volumes for H, D, C, N and O, which were found  
experimentally for proteins by Fraser et al. (J. Appl. Cryst.(1978), 11, 693).  

The P(q) is affected by the binsize, r. Try decreasing r until convergence is obtained (Default is 2.0 aa).    

### Water layer
If chosen, a water layer is added explicitely and included in the PDDF and thus in the calculated scattering  
intensity.  
The water layer is included as described in the S.I. of Midtgaard et al. (2017), FEBS, DOI: 10.1111/febs.14345.  
Briefly, dummy beads are created with the chosen scattering contrast at the edge of the protein. If chosen,  
the beads are not included in the transmembrane (TM) part, which is relevant for membrane proteins.  
The water layer is explicitely written to a "pseudo" PDB file  and can thus be visualized, e.g. with PyMOL.  
("pseudo") because the dummy beads are not real atoms.  

### Atomic volumes and scattering lengts (copied from Readme.h)
                 x    y    z       X-ray[cm]        Neutron[cm]                                    V[AA^3]  Mw[Da]  Name
[HYDROGEN]  =  { 0.0, 0.0, 0.0,    1 * 2.82e-13,    3.741e-13,                                       5.15,   1.0,   'A'} ,    
[DEUTERIUM]  = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13*SolventD2O-3.741e-13*(1.0-SolventD2O), 5.15,   2.0,   'A'} ,    
[CARBON]     = { 0.0, 0.0, 0.0,    6 * 2.82e-13,    6.646e-13,                                      16.44,  12.0,   'A'} ,    
[NITROGEN]   = { 0.0, 0.0, 0.0,    7 * 2.82e-13,    9.360e-13,                                       2.49,  14.0,   'A'} ,    
[OXYGEN]     = { 0.0, 0.0, 0.0,    8 * 2.82e-13,    5.803e-13,                                       9.13,  16.0,   'A'} ,    
[FLUORINE]   = { 0.0, 0.0, 0.0,    9 * 2.82e-13,    5.654e-13,                                      12.04,  19.0,   'A'} ,    
[SODIUM]     = { 0.0, 0.0, 0.0,   11 * 2.82e-13,    3.630e-13,                                      59.00,  23.0,   'A'} ,    
[MAGNESIUM]  = { 0.0, 0.0, 0.0,   12 * 2.82e-13,    5.375e-13,                                      45.09,  24.3,   'A'} ,    
[PHOSPHORUS] = { 0.0, 0.0, 0.0,   15 * 2.82e-13,    5.130e-13,                                       5.81,  31.0,   'A'} ,    
[SULFUR]     = { 0.0, 0.0, 0.0,   16 * 2.82e-13,    2.847e-13,                                      25.31,  32.1,   'A'} ,    
[CHLORINE]   = { 0.0, 0.0, 0.0,   17 * 2.82e-13,    9.577e-13,                                      24.49,  35.5,   'A'} ,    
[CALCIUM]    = { 0.0, 0.0, 0.0,   20 * 2.82e-13,    4.700e-13,                                      59.07,  40.1,   'A'} ,    
[MANGANESE]  = { 0.0, 0.0, 0.0,   25 * 2.82e-13,    -3.730-13,                                      37.05,  54.9,   'A'} ,    
[IRON]       = { 0.0, 0.0, 0.0,   26 * 2.82e-13,    9.450e-13,                                      35.88,  55.8,   'A'} ,    
[COPPER]     = { 0.0, 0.0, 0.0,   29 * 2.82e-13,    7.718e-13,                                      35.13,  63.5,   'A'} ,    
[ZINK]       = { 0.0, 0.0, 0.0,   30 * 2.82e-13,    5.680e-13,                                      37.79,  65.4,   'A'} ,    
[GOLD]       = { 0.0, 0.0, 0.0,   79 * 2.82e-13,    7.630e-13,                                      12.41, 196.7,   'A'} ,    
[WATER]      = { 0.0, 0.0, 0.0,n_W*8 * 2.82e-13,n_W*5.803e-13,                          n_W*(a*30-2*5.15),   0.0,   'A'} ,    
[UNKNOWN]    = { 0.0, 0.0, 0.0,    0 * 2.82e-13,    0.000e-13,                                       1.00,   0.0,   'A'}    


## Note on Dmax
The Dmax is the true Dmax, i.e. the furhest distance betweeen any two atoms (or water beads) in the protein.  
However, this value is not detectable, since only very few atom pairs have this distance. Therefore, Dmax  
found with CaPP will be a larger value than the experimentally determined Dmax, and the values are thus not  
directly comparable. A pragmatic solution is to use as Dmax the value of r, where p(r) has decreased to 1% of 
its maximal value. This number is comparable to experimentally determined Dmax values, and mathes well with a 
visual assesment of Dmax for the calculated pair distance distribution functions. The Dmax found using this 1%
threshold is thus good prediction of the experimentally determined Dmax.  

## Fitting
CaPP can automatically fit the theoretical formfactor P(q) to a dataset, if this is provided. 
Two parameters, scaling (S) and background (B), are fitted to obtain a theoretical scattering given as:  

        I(q) = S * P(q) + B

where q = 4pi sin(t)/l and 2t and l are the scattering angle and wavelength respectively. The q-values are  
imported from the data file and assumed to be in units of aa (nm can be chosen in the GUI).  
Data are fitted with chi2-minimization, using the error bars from data (3rd column).  
No resolution effects are included for SAXS, nor SANS data.  
The density of the water layer cannot (yet) be fitted.  

## Versions  
#### CaPP 1  
- Release June 2017  

#### CaPP 2  
- Release January 2018 
- Fitting of scale and bg  
- other minor improvements  

#### CaPP 3
- Release August 2018
- Fitting of WL contrast  
- always calculate P(q) when p(r) is calculated  
- option to exclude the first q-points from the fit  
- improvements of GUI  
- more precise calc of P(q), better form factor approximation
- treat RNA and DNA
- bugfix: not include TER lines when counting #atoms (function CheckNumberOfAtoms()). Fixed!
- robustness: pdb was not read correctly after PyMOL manipulation, so WL was placed wrong. Fixed! 
- stop use of structs for treatment of residues - difficult, unintutive, ineffective, etc... 
- clean-up in source codes - significant simplifications!
- less printing to terminal
- Option to fit with a linear combination of 2 pdb files  
- allows for sucrose SAXS contrast variation  
- include resolution effects (from 4th column in data) in SANS fitting (3.6) 
- check number of columns in data files (3.6)
- handle datafile with header in the buttom of data file (3.6)
- GUI implementation of Rg calculator was buggy - removed (3.8)
- Option for perdeuteration included (3.8)
- include info about excluded WL in header (3.8)
- calculate A00 functions (3.8)
- fix minor bugs, expand HETATM list (3.9)

## License
CaPP is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.          

CaPP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).  

## Citing the program  
If you use the CaPP software in your work, please cite the GitHub addres:                                 

        github.com/Niels-Bohr-Institute-XNS-StructBiophys/CaPP  

## Acknowledgements
Søren Kynde made most of the c subroutines that were edited and gathered to form the program.  
Thanks to Lasse Sander Dreyer, and to Martin Schmiele for valuable contributions:  
debugging, testing, writing smarter subfunctions etc.  
To CoNeXT and University of Copenhagen for co-funding the project.   
To my supervisor, Lise Arleth, for supporting the project.  

## Known bugs (to be fixed)
- the program cannot find the pdb file if the path contains white space.  
- the plotting is not working on Windows (the program crashes) and is therefore disabled.  
- the program is generally not stable under Windows.   

## Future development   
- automatically give Dmax comparable with exp. detectable Dmax (see "Note on Dmax")   
- write a paper     
