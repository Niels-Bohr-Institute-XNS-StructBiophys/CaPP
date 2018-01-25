# CaPP 2.0
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

        >> python path-to-folder/CaPP_2.0.py  

### Running the program, batch mode
Type in the terminal  

        >> ./path-to-folder/capp_mac [options] PDBFILE.pdb  (mac OS)

        >> ./path-to-folder/capp_windows [options] PDBFILE.pdb  (windows)

        >> ./path-to-folder/capp [options] PDBFILE.pdb  (linux)

##### Options:  

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
Default: 2.0 Aangstrom.  
Too small will not give a smooth p(r), too large will give wrong results for P(q) if  
this is calculated from the p(r) function.  

##### Example:  

        >> capp -c 0.1 -m 30 -s 0.4 -r 2.0 4u2p.pdb  

## Dependencies

### Dependencies for the GUI  
- python 2.7  
- wxpython 2.9
- matplotlib  
- scipy
tested on MacOS 10.12 under Enthought Canopy 2.7 (64-bit)   
and with python under cygwin on windows 7 (no build-in plotting available)  

##### Windows (with cygwin terminal) type:  
    $ wget.exe http://peak.telecommunity.com/dist/ez_setup.py  
    $ python ez_setup.py  
    $ /cygdrive/c/Python27/Scripts/easy_install-2.7.exe pip  
    $ /cygdrive/c/Python27/Scripts/pip.exe install wxpython  
    $ /cygdrive/c/Python27/Scripts/pip.exe install matplotlib 
    § /cygdrive/c/Python27/Scripts/pip.exe install scipy

### Dependencies for developers and users of Linux
- To recompile the source code, a c-compiler is needed,  
e.g. gcc (Linux/MacOS) or Pelles C (Windows).  

### Different platforms  
Executables have been made for MacOS and Windows  
Users of other OS should:  
1) compile MainFunction.c (with flag -lm) and name the executable "capp"  

        >> gcc MainFunction.c -o capp  -lm

2) place the executable, capp, in the same folder as CaPP.py  
3) Run CaPP  

        >> python CaPP_A.B.py  

CaPP has been tested on  
- MacOS 10.12 (Sierra), compiled with gcc
- Windows Vista and Windows 7, compiled with Pelles C  
- Ubuntu 16.04 LTS, compiled with gcc  

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
columns 3-6 used only to calculate th form factor P(q)  
(optional, see CaPP_A.B.py for details)    

##### ABC_Pq  
A datafile with the form factor for the protein  

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

## Note on Dmax
The Dmax is the true Dmax, i.e. the furhest distance betweeen any two atoms (or water beads) in the protein.  
However, this value is not detectable, since only very few atom pairs have this distance. Therefore, Dmax  
found with CaPP will be a larger value than the experimentally determined Dmax, and the values are thus not  
directly comparable. A pragmatic solution is to use as Dmax the value of r, where p(r) has decreased to 1% of 
its maximal value. This number is comparable to experimentally determined Dmax values, and mathes well with a 
visual assesment of Dmax for the calculated pair distance distribution functions. The Dmax found using this 1%
threshold is thus good prediction of the experimentally determined Dmax.  

## Fitting (New feature in CaPP 2.0)
CaPP can automatically fit the theoretical formfactor P(q) to a dataset, if this is provided. 
Two parameters, scaling (S) and background (B), are fitted to obtain a theoretical scattering given as:  

        I(q) = S * P(q) + B

where q = 4pi sin(t)/l and 2t and l are the scattering angle and wavelength respectively. The q-values are  
imported from the data file and assumed to be in units of aa (nm can be chosen in the GUI).  
Data are fitted with chi2-minimization, using the error bars from data (3rd column).  
No resolution effects are included for SAXS, nor SANS data.  
The density of the water layer cannot (yet) be fitted.  

## Versions  
#### CaPP_1.0  
- from summer 2017
#### CaPP_2.0  
- January 2018  
- Fitting of scale and bg
- other minor improvements  

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
- the program cannot find the pdb file if the path contain white space.  
- the plotting is not working on Windows (the program crashes) and is therefore disabled.  
- the program is generally not stable under Windows.   

## Future development  
- fitting option for the WL contrast
- include resolution effects (from 4th column in data) in SANS fitting
- automatically give Dmax comparable with exp. detectable Dmax (see "Note on Dmax")  
- expand HETATM library  
- include info about excluded WL in header  
- write paper to present the program to the world  
- add the 20 natural occuring aa to the HETATM library (these can be ligands, thus listed as hetero atoms)  
- possibility to plot p(r) and P(q), if already calculated  
