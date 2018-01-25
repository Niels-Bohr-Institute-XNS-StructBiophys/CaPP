///                           ///
///  CaPP 1.0                 ///
///  Copyright 2015,          ///
///  University of Copenhagen ///
///                           ///
///  Søren Kynde              ///
///  Andreas N Larsen         ///
///                           ///
///  anlarsen@nbi.ku.dk       ///
///                           ///

//Compile with:
//gcc cap_1.0.c -o capp_mac
//gcc cap_1.0.c -o capp_windows
//gcc cap_1.0.c -o capp_linux

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include "Helpfunctions.h"
#include "AddWaterLayerToPDB.h"

int main(int argc, char **argv)
{
    clock(); /*Measure computation time */
    
    
            // DEFAULT PARAMETER VALUES
    
    double SolventD2O = -1.0;
    int PointsInPofr;
    double Delta_r = 2.0; //default "resolution" is 2 Å
    double WaterLayerContrast = 0.0; //contrast of hydration layer
    double HalfBilayerThickness = 0.0;
    int OPTION_d_CHOSEN = 0;
    int OPTION_m_CHOSEN = 0;
    int OPTION_g_CHOSEN = 0;
    int NO_OPTIONS_CHOSEN = 1;
    
    double NonExchNH = 0.1; // cannot be changed (in GUI or in batch mode).
                            // give rise to minor, almost q-independent contribution,
                            // since NH groups are distributed all over the protein.
    
            // INSTRUCTIONS FOR USE
    
/*    printf("\n \
           Welcome to CaPP (Calculate p(r) function from a PDB file)\n \
           -------------------------------------------------------------------- \n \
           This program calculates the pair distance distribution function, p(r),\n \
           of a PDB file. The distances are weighted with the scattering contrast\n \
           for X-rays in an aqueous solution or for Neutrons in a D2O/H2O solvent,\n \
           as specified by the user. \n \
           Optionally, a water layer can be added \n \
           -------------------------------------------------------------------- \n \
           Usage:\n\n \
           capp [options] NAME_of_PDB \n\n \
           Options: \n\n \
            -c [input: Contrast of water layer] \n\
               Add a water layer with (c)ontrast between 0 and 2 times the solvent scattering length.\n\
               Typically 0.1.\n\
               Default: No water layer.\n\n\
            -d [no input] \n\
               Only relevant for membrane proteins.\n\
               Removes water layer from the bilayer region. \n\
               Choose -d if the pdb is from the OPM (d)atabase, that provides the bilayer thickness. \n\n \
            -m [input: Bilayer Thickness] \n\
               Only relevant for membrane proteins. \n\
               Removes water layer from the bilayer region. \n\
               Choose -m to (m)anually provide the bilayer thickness in Aangstrom. \n\
               Typically 30 Aangstrom. \n\
               NB: Remember to place the TMD perpendicular to the xy-plane, in z=0!\n\n \
            -s [input: prc D20 in the solvent] \n\
               Choose SANS contrast and enter the D2O-content (between 0 and 1) of the (s)olvent. \n\
               SAXS contrast asssumed if option is not chosen. \n\n \
            -r [input: Resolution of p(r) function] \n\
               Change the (r)esolution, i.e. the binsize (in Aangstrom) of the p(r) function. \n\
               Default: 2.0 Aangstrom. ");
 */

    if(argc == 1){
        printf("\n\n");
        printf("\n            ******************************************************");
        printf("\n            Please provide a PDB file name.");
        printf("\n            See above for instructions on how to use the program.");
        printf("\n            ******************************************************\n\n\n\n");
        exit(-1);
    }
    
            // IMPORT PDB FILE AND ADD .pdb EXTENSION IF NEEEDED
    
    FILE *PointerToFile;

    char *filename = argv[argc-1];
    if (strcmp(filename,"-help")==0 || strcmp(filename,"-h")==0 || strcmp(filename,"-use")==0 || strcmp(filename,"help")==0 || strcmp(filename,"h")==0 || strcmp(filename,"use")==0 ){
        printf("\n\n");
        printf("\n            ******************************************************");
        printf("\n            See above for instructions on how to use the program.");
        printf("\n            ******************************************************\n\n\n\n");
        exit(-1);
    }
    
    if( (PointerToFile = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename
    }
    
    if( (PointerToFile = fopen(filename,"r"))==0){
        printf("\n\n");
        printf("\n            *************************************************************");
        printf("\n            Cannot find the PDB file: %s, is the name and path correct?",filename);
        printf("\n            - is the name and path correct?");
        printf("\n            The pdb file should be the last input, after the options.");
        printf("\n            See above for instructions on how to use the program.");
        printf("\n            *************************************************************\n\n\n\n");
        exit(-1);
    }
    
            // READING OPTIONS
    
        printf("\n");
        printf("\n            ************* Options chosen: *************\n");

    char ch;
    const char *ValidOpts = "c:dm:s:r:g";
    char *CheckOptArg;
    while((ch=options(argc,argv,ValidOpts))!=-1)
    {
        NO_OPTIONS_CHOSEN = 0;
        switch(ch)
        {
            case 'c':
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -c, a number between 0 and 1.\n");
                    printf("\n            Input argument for option -c was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                WaterLayerContrast = char2double(OptArg);
                printf("\n            (-c) Water Layer Contrast = %6.2f prc of solvent scattering length\n", WaterLayerContrast*100.0);
                if (WaterLayerContrast == 0) {
                    printf("\n\n            !!! NB !!! Since you have chosen 0 as water layer contrast (input -c), no water layer will be added\n\n");
                }
                if (WaterLayerContrast < -2.00 || WaterLayerContrast > 2.00){ printf("\n\n\n            !!!ERROR!!! The water Layer contrast, option -c, should be between -2.0 and 2.0\n\n\n"); exit(-1);}

                break;
            case 'd':
                if (OPTION_m_CHOSEN == 1) {
                    printf("\n\n\n            !!!ERROR!!! Dear user, you cannot choose BOTH option -d and option -m.\n\n\n\n");
                    exit(-1);
                }
                HalfBilayerThickness = CheckHalfBilayerThickness(filename);
                printf("\n            (-d) Bilayer Thickness from OPM = %6.2f\n", 2 * HalfBilayerThickness);
                OPTION_d_CHOSEN = 1;
                break;
            case 'm':
                if (OPTION_d_CHOSEN == 1) {
                    printf("\n\n\n            !!!ERROR!!! Dear user, you cannot choose BOTH option -d and option -m.\n\n\n\n");
                    exit(-1);
                }
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -m, a bilayer thickness in Angstrom.\n");
                    printf("\n            Input argument for option -m was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                HalfBilayerThickness = char2double(OptArg)*0.5;
                printf("\n            (-m) Manual Bilayer Thickness = %6.2f Angstrom.\n            NB: Remember to place the TMD perpendicular to the xy-plane, in z=0!\n", 2 * HalfBilayerThickness);
                OPTION_m_CHOSEN = 1;
                if (HalfBilayerThickness == 0) {
                    printf("\n\n\n            !!! Error !!! Dear user, you have chosen a bilayer thickness of 0 Angstrom or not given any input after option -m. Please try again with an input. \n\n\n");
                    exit(-1);
                }
                break;
            case 's':
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -s, a D20 fraction between 0 and 1 for SANS. SAXS is assumed if option -s is not used.\n");
                    printf("\n            Input argument for option -s was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                SolventD2O = char2double(OptArg);
                if (SolventD2O < 0.0) { printf("\n        (-s) SAXS\n"); }
                else{ printf("\n            (-s) SANS with %4.1f prc D2O in the solvent\n", SolventD2O*100.0); }
                if (SolventD2O > 1.00){ printf("\n\n\n            !!!ERROR!!! The solvent D2O content, option -s, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n"); exit(-1);}
                break;
            case 'r':
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -r, resoluton in Aangstrom.\n");
                    printf("\n            Input argument for option -r was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                Delta_r = char2double(OptArg);
                if (Delta_r == 0) {
                    printf("\n\n\n            !!! Error !!! Dear user, you have chosen 0 Angstrom resolution or not given any resolution after option -r. Please try again with an input. \n\n\n");
                    exit(-1);
                }
                printf("\n            (-r) Resolution = %6.4f Angstrom\n", Delta_r);
                break;
            case 'g':
                OPTION_g_CHOSEN = 1;
        }
    }
    if(NO_OPTIONS_CHOSEN == 1) {printf("\n            No options chosen - default values used\n");}
    
    if( (OPTION_d_CHOSEN || OPTION_m_CHOSEN) && WaterLayerContrast == 0) {
        printf("\n\n");
        printf("\n            *************************************************************");
        printf("\n            You have chosen option -m or option -d (exclude WL at TMD),");
        printf("\n            but the contrast of the water layer is 0, i.e. there is no WL,");
        printf("\n            i.e. there is nothing to exclude.");
        printf("\n            See above for instructions on how to use the program.");
        printf("\n            *************************************************************\n\n\n\n");
        exit(-1);
    }

    if (SolventD2O >= 0.0 && SolventD2O <= 1.0) {
        //printf("\n            SANS contrast with %2.0f prc D2O in the solvent\n", SolventD2O*100);
    } else {
        printf("\n            SAXS contrast (default)\n");
    }
    
        printf("\n            ********************************************\n");

    if(WaterLayerContrast > 0.0 || WaterLayerContrast < 0.0) {
        printf("\n\n");
        printf("\n            ************ Adding Water Layer ************\n");
        filename = AddWaterLayerToPDB(filename, HalfBilayerThickness);
        PointerToFile = fopen(filename,"r");
        printf("\n            Bilayer Thickness = %6.4f\n", 2 * HalfBilayerThickness);
        printf("\n            Water has been added with contrast = %6.4f\n", WaterLayerContrast);
        printf("\n            PDB with protein and water generated: %s \n", filename);
        int PDBFilenameLength = strlen(filename) - 4; // length of first part of pdb filename, e.g. 6 for 1x2y_w.pdb
        char WaterOnlyFileName[PDBFilenameLength];
        memcpy(WaterOnlyFileName, filename, PDBFilenameLength);
        WaterOnlyFileName[PDBFilenameLength] = '\0'; // NULL terminates the new string
        printf("\n            PDB with water only generated: %s_only.pdb \n", WaterOnlyFileName);
        printf("\n            ********************************************\n");
    }
    
            // READ DATA FOM PDB FILE
    
    double a = 1.0/(1.0 + WaterLayerContrast);

    double n_W = 4.133; // number of watermolecules per water residue
    double r, D, Dmax=0.0;
    int index;
    enum {DEUTERIUM, CARBON, NITROGEN, OXYGEN, FLUORINE, SODIUM, MAGNESIUM, PHOSPHORUS, SULFUR, CHLORINE, CALCIUM, MANGANESE, IRON, COPPER, ZINK, WATER, UNKNOWN};
    
    // ATOM LIST
    /*                                      Scattering-length                                                                       */
    /*                   x    y    z       X-ray[cm]        Neutron[cm]                                    V[AA^3] Mw[Da]  Name*/
    struct Atom element[] = {
        [DEUTERIUM]  = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13*SolventD2O-3.741e-13*(1.0-SolventD2O), 5.15,  2.0,   'A'} ,
        [CARBON]     = { 0.0, 0.0, 0.0,    6 * 2.82e-13,    6.646e-13,                                      16.44, 12.0,   'A'} ,
        [NITROGEN]   = { 0.0, 0.0, 0.0,    7 * 2.82e-13,    9.360e-13,                                       2.49, 14.0,   'A'} ,
        [OXYGEN]     = { 0.0, 0.0, 0.0,    8 * 2.82e-13,    5.803e-13,                                       9.13, 16.0,   'A'} ,
        [FLUORINE]   = { 0.0, 0.0, 0.0,    9 * 2.82e-13,    5.654e-13,                                      12.04, 19.0,   'A'} ,
        [SODIUM]     = { 0.0, 0.0, 0.0,   11 * 2.82e-13,    3.630e-13,                                      59.00, 23.0,   'A'} ,
        [MAGNESIUM]  = { 0.0, 0.0, 0.0,   12 * 2.82e-13,    5.375e-13,                                      45.09, 24.3,   'A'} ,
        [PHOSPHORUS] = { 0.0, 0.0, 0.0,   15 * 2.82e-13,    5.130e-13,                                       5.81, 31.0,   'A'} ,
        [SULFUR]     = { 0.0, 0.0, 0.0,   16 * 2.82e-13,    2.847e-13,                                      25.31, 32.1,   'A'} ,
        [CHLORINE]   = { 0.0, 0.0, 0.0,   17 * 2.82e-13,    9.577e-13,                                      24.49, 35.5,   'A'} ,
        [CALCIUM]    = { 0.0, 0.0, 0.0,   20 * 2.82e-13,    4.700e-13,                                      59.07, 40.1,   'A'} ,
        [MANGANESE]  = { 0.0, 0.0, 0.0,   25 * 2.82e-13,    -3.730-13,                                      37.05, 54.9,   'A'} ,
        [IRON]       = { 0.0, 0.0, 0.0,   26 * 2.82e-13,    9.450e-13,                                      35.88, 55.8,   'A'} ,
        [COPPER]     = { 0.0, 0.0, 0.0,   29 * 2.82e-13,    7.718e-13,                                      35.13, 63.5,   'A'} ,
        [ZINK]       = { 0.0, 0.0, 0.0,   30 * 2.82e-13,    5.680e-13,                                      37.79, 65.4,   'A'} ,
        [WATER]      = { 0.0, 0.0, 0.0,n_W*8 * 2.82e-13,n_W*5.803e-13,                          n_W*(a*30-2*5.15),  0.0,   'A'} ,
        [UNKNOWN]    = { 0.0, 0.0, 0.0,    0 * 2.82e-13,    0.000e-13,                                       1.00,  0.0,   'A'}
    };
    // Note: [WATER] is without hydrogen/deuterium... will be added later on (implicit way)
    // Note: weight of [WATER] is set to 0, since the WL should not be included in the calculation of the atomic weight
    // Volumes from:
    //    (H,C,D,N,O):                              Fraser, J Appl Cryst(1978), 11, p693.
    //                                               Experimentally determined values
    //    (P, S, Mn, Ca, F, Mg, Cl, Cu, Fe, Zn):    Batsanov, Inorg Mat(2001), vol37, no9, p871.
    //                                               Van der Waals radii, mean of values from Table 1 and 2 (P from Table 5)
    
    // Scattering length of solvent
    double ScatteringLengthSolvent = 2.0 * element[DEUTERIUM].NeutronScatteringLength + element[OXYGEN].NeutronScatteringLength;
    if(SolventD2O<0.0){
           ScatteringLengthSolvent = 2.0 * element[DEUTERIUM].XRayScatteringLength + element[OXYGEN].XRayScatteringLength;
    }
    double VolumeSolvent = 30.0;
    double rhoSolvent = ScatteringLengthSolvent/VolumeSolvent;
    
    double Sum_dB = 0.0, SumAtomWeight=0.0;
    int R,i=0,j,u=0;
    char buffer[256];

    int GuessNumberOfAtoms = CheckNumberOfAtomsInPDBFile(filename);
    struct Atom *Atoms = (struct Atom *)  calloc(GuessNumberOfAtoms, sizeof(struct Atom));
    struct Atom COM   = {0,0,0};
    double Rg2=0;
    double SqDist, SqDistMax=0, GuessOnDmax;
    double SumOfVolume = 0.0;
    
    char AminoName, AtomName, LongAtomName,AlternativeAtomPosition;
    
    double n_D,n_H; // no of Deuterium and Hydrogen attachd to main atom

    double *Ba = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // Ba = scattering length vector for atom
    double *Bs = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // Bs = scattering length vector for excluded solvent
    double *dB = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // dB = Ba - Bs

    int NumberOfH = 0;
    int NumberOfAlternativeAtoms = 0;
    int NumberOfOtherLines = 0;
    
    while(fgets(buffer,sizeof(buffer),PointerToFile)!=NULL){
        R = 0;
        
        // read atom name.
        if (sscanf(buffer, "ATOM%*73c%s",&Atoms[i].Name) == 1 || sscanf(buffer, "HETATM%*71c%s",&Atoms[i].Name) == 1) {}
        // read alternative atom pos
        if (sscanf(buffer, "ATOM%*12c%c",&AlternativeAtomPosition) == 1 || sscanf(buffer, "HETATM%*10c%c",&AlternativeAtomPosition) == 1) {}
        
        // ignore if H or D and ignore if it is a repeated atom with alternative position, only use position A. Else, read atom.
        if (strcmp(&Atoms[i].Name,"H") == 0 || strcmp(&Atoms[i].Name,"D") == 0) {NumberOfH++;}
        else if (AlternativeAtomPosition == 'B') {NumberOfAlternativeAtoms++;}
        else if (sscanf(buffer,"ATOM%*18c%d%*4c%lf%lf%lf%*22c%s",&R,&Atoms[i].x,&Atoms[i].y,&Atoms[i].z,&AtomName)==5 ||
            sscanf(buffer,"HETATM%*16c%d%*4c%lf%lf%lf%*22c%s",&R,&Atoms[i].x,&Atoms[i].y,&Atoms[i].z,&AtomName)==5 ) {
            n_D = 0;
            n_H = 0;
            
            // Read element and give values from atom list
            if (strcmp(&AtomName,"C")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[CARBON]);}
            else if (strcmp(&AtomName,"N")==0  || strcmp(&AtomName,"N1+")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[NITROGEN]);}
            else if (strcmp(&AtomName,"O")==0 || strcmp(&AtomName,"O1-")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[OXYGEN]);}
            else if (strcmp(&AtomName,"F")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[FLUORINE]);}
            else if (strcmp(&AtomName,"NA")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[SODIUM]);}
            else if (strcmp(&AtomName,"MG")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[MAGNESIUM]);}
            else if (strcmp(&AtomName,"P")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[PHOSPHORUS]);}
            else if (strcmp(&AtomName,"S")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[SULFUR]);}
            else if (strcmp(&AtomName,"CL")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[CHLORINE]);}
            else if (strcmp(&AtomName,"CA")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[CALCIUM]);}
            else if (strcmp(&AtomName,"MN")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[MANGANESE]);}
            else if (strcmp(&AtomName,"FE")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[IRON]);}
            else if (strcmp(&AtomName,"CU")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[COPPER]);}
            else if (strcmp(&AtomName,"ZN")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[ZINK]);}
            else if (strcmp(&AtomName,"Q")==0) {
                CopyPhysicalParameters(&Atoms[i], &element[WATER]); n_D = 2.0 * n_W;}
            else {
                CopyPhysicalParameters(&Atoms[i], &element[UNKNOWN]); u++;}
            
            // Add H/D to amino acids. Ordered alphabetically
            if (sscanf(buffer,"ATOM%*13c%s", &AminoName) == 1)
            {
                if (strcmp(&AminoName,"ALA")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"ASN")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"ND2")==0) {n_D = 2;}
                }
                else if (strcmp(&AminoName,"ARG")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"NE")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"NH1")==0) {n_D = 2;}
                    else if (strcmp(&LongAtomName,"NH2")==0) {n_D = 2; Atoms[i].XRayScatteringLength = Atoms[i].XRayScatteringLength - 2.82e-13; } // kation
                }
                else if (strcmp(&AminoName,"ASP")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"OD2")==0) {Atoms[i].XRayScatteringLength = Atoms[i].XRayScatteringLength + 2.82e-13; } // anion
                }
                else if (strcmp(&AminoName,"CYS")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"SG")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"GLN")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"NE2")==0) {n_D = 2;}
                }
                else if (strcmp(&AminoName,"GLU")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"OE2")==0) {Atoms[i].XRayScatteringLength = Atoms[i].XRayScatteringLength + 2.82e-13; } // anion}
                }
                else if (strcmp(&AminoName,"GLY")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"HIS")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CE1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"NE2")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"ILE")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CG1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG2")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CD1")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"LEU")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CD1")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CD2")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"LYS")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CE")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"NZ")==0) {n_D = 3; Atoms[i].XRayScatteringLength = Atoms[i].XRayScatteringLength - 2.82e-13;} // kation
                }
                else if (strcmp(&AminoName,"MET")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CE")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"PHE")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CD2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CE1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CE2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CZ")==0) {n_H = 1;}
                }
                else if (strcmp(&AminoName,"PRO")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {} // proline is the only aa without an amide NH group
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CG")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"SER")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"OG")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"THR")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CG2")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"OG1")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"TRP")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CE3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CZ2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CZ3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CH2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"NE1")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"TYR")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CD1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CD2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CE1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CE2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"OH")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"VAL")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    //printf("\n\nLongAtomName = %3c\n\n",LongAtomName);
                    //puts(&LongAtomName);
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CG1")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CG2")==0) {n_H = 3;}
                }
            }
            // H/D added to HETATM according to PHENIX
            else if (sscanf(buffer,"HETATM%*11c%s", &AminoName) == 1)
            {
                if (strcmp(&AminoName,"HOH")==0) {
                    n_D = 2; // HOH is a water atom in the crystal
                    Atoms[i].Volume = 30 - 2*5.15; // implicit solvent, therefore the volume of two H/D is subtracted
                }
                else if (strcmp(&AminoName,"ZK1")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"NAP")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"CL")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CAJ")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CAK")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CAL")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CAM")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CAN")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CAO")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"CLA")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C18")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2A")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3A")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C9")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CAB")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CHB")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CHC")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CHD")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C10")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C14")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C15")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C16")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C17")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C19")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C20")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CAA")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CAC")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CBA")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CBB")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CBC")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CED")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CMA")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CMB")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CMC")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CMD")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"LMG")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"O2")==0) {n_D = 2;}
                    else if (strcmp(&LongAtomName,"O3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O5")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C45")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C44")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C43")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C42")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C41")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C40")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C39")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C38")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C37")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C36")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C35")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C34")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C33")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C32")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C30")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C29")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C27")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C25")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C23")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C22")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C20")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C19")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C18")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C17")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C16")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C15")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C14")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C13")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C10")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C9")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                }
                else if (strcmp(&AminoName,"LMT")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"O2")==0) {n_D = 2;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C9")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C10")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O1'")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O6'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O1B")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2B")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3B")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O5B")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O6B")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C1B")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2B")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3B")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4B")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5B")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6B")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"LHG")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"O2")==0) {n_D = 2;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C10")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C9")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"BCR")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C40")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C39")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C38")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C37")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C36")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C35")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C29")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C28")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C27")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C23")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C21")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C20")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C19")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C17")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C16")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C15")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C14")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C34")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C32")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C31")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C33")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C10")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"MES")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N4")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"NAG")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"CO3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"CO4")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"CO6")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"BMA")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"CO3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"CO4")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"TG1")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C34")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C9")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C14")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C15")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C16")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C17")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C18")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C19")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C20")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C23")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C25")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C26")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C31")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C33")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C28")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C29")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C30")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"O6")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"O11")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C25")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"PCW")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"PEG")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"PG0")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OTT")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"GOL")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"CRS")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                }
            }
            // Update volume and scattering length with implicit H and D
            Atoms[i].Volume += (n_D + n_H) * 5.15;
            SumOfVolume += Atoms[i].Volume;
            
            if(SolventD2O < 0.0){
                Atoms[i].XRayScatteringLength += (n_D + n_H) * element[DEUTERIUM].XRayScatteringLength;
                Ba[i] = Atoms[i].XRayScatteringLength;
            }
            else {
                Atoms[i].NeutronScatteringLength += n_D * element[DEUTERIUM].NeutronScatteringLength - n_H * 3.741e-13;
                Ba[i] = Atoms[i].NeutronScatteringLength;
            }
            Bs[i] = rhoSolvent * Atoms[i].Volume;
            dB[i] = Ba[i] - Bs[i];
            Sum_dB += dB[i];
            SumAtomWeight += Atoms[i].Weight;
            COM.x += dB[i] * Atoms[i].x;
            COM.y += dB[i] * Atoms[i].y;
            COM.z += dB[i] * Atoms[i].z;
            
            i++;
        } // end read atom
        else {NumberOfOtherLines++;}
    } // end While loop
    fclose(PointerToFile);
    
    int NumberOfAtomsRead = i;
    double MeanVolume = SumOfVolume/NumberOfAtomsRead;

    COM.x /= Sum_dB;
    COM.y /= Sum_dB;
    COM.z /= Sum_dB;

    
            // CREATE FILE FOR p(r) FUNCTION
    
    FILE *outputFile;
    char *inputPDB = GetCStringBeforeLastDelimiter(filename, '.'); //extracts 1HEJ from 1HEJ.pdb
    char *outputfilename = strcat(inputPDB, "_pr.dat"); // adds pr.dat to 1HEJ
    outputFile = fopen(outputfilename,"w");
    
    
            // PRINT OUTPUT TO TERMINAL AND p(r) FILE
    
    printf("\n\n# Location of PDB-file: %s\n", filename); fprintf(outputFile, "# Location of PDB-file: %s\n", filename);
    printf("# Number of atoms read: %d\n",NumberOfAtomsRead);
    printf("# Unknown atoms encoutered: %d\n",u);
    printf("# Number of explicit Hs and Ds ignored: %d \n", NumberOfH);
    printf("# Number of other lines ignored: %d \n", NumberOfOtherLines);
    printf("# Number of alternative atoms (alternative position B) ignored: %d \n", NumberOfAlternativeAtoms);
    if (SolventD2O >= 0.0 && SolventD2O <= 1.0) {
        printf("# SANS with D2O fraction in solvent: %.2lf\n", SolventD2O); fprintf(outputFile, "# SANS with D2O fraction in solvent: %.2lf\n", SolventD2O);
    } else {
        printf("# SAXS contrast\n"); fprintf(outputFile, "# SAXS contrast\n");
    }
    if (WaterLayerContrast > 0.0 || WaterLayerContrast < 0.0) {
        printf("# Water layer added with contrast: %6.4lf\n", WaterLayerContrast); fprintf(outputFile, "# Water layer added with contrast: %6.4lf\n", WaterLayerContrast);
    } else {
        printf("# No water layer added\n"); fprintf(outputFile, "# No water layer added\n");
    }
    
    printf("# Center of excess scattering (%g, %g, %g)\n",COM.x,COM.y,COM.z);
    printf("# Atomic weight: %g\n", SumAtomWeight); fprintf(outputFile, "# Atomic weight: %g\n", SumAtomWeight);

    double SqDistMaxPrevious = 0.0;
    for(i=0;i<NumberOfAtomsRead;i++){
        SqDist = SquareDistance(Atoms[i],COM);
        Rg2 += dB[i]*SqDist;
        SqDistMax = max(SqDistMax,SqDist);
        SqDistMaxPrevious = SqDistMax;
    }
    GuessOnDmax = 2*sqrt(SqDistMax);
    PointsInPofr=floor(1.07 * GuessOnDmax/Delta_r);
    if (OPTION_g_CHOSEN == 0) { printf("# PointsInPofr = %d\n", PointsInPofr); fprintf(outputFile, "# PointsInPofr = %d\n", PointsInPofr); }
    
    double *Pofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Gofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Hofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Jofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Kofr = (double *) calloc( PointsInPofr, sizeof(double));
    if(OPTION_g_CHOSEN == 1){
        printf("\n\n###########################################\n");
    }
    printf("# Radius of gyration    %8.2lf Angstrom\n", sqrt(Rg2/Sum_dB)); fprintf(outputFile, "# Radius of gyration    %8.2lf\n", sqrt(Rg2/Sum_dB));

    if (OPTION_g_CHOSEN == 1) {
        printf("###########################################\n\n");
        exit(-1);
    }
    printf("# Upper limit for Dmax %8.2lf\n", GuessOnDmax);
    if(!PointsInPofr)
        exit(1);
    printf("# Progression [prc]:\n# ");
    
    int progression = 0;

    for(i=0;i<NumberOfAtomsRead;i++){
        Gofr[0] += Ba[i] * Ba[i];
        Hofr[0] += Ba[i] * Bs[i];
        Jofr[0] += Bs[i] * Ba[i];
        Kofr[0] += Bs[i] * Bs[i];
        for(j=i+1;j<NumberOfAtomsRead;j++){
            D = Distance(Atoms[i],Atoms[j]);
            index = floor(D/Delta_r);
            Pofr[index] +=       dB[i] * dB[j];
            Gofr[index] += 2.0 * Ba[i] * Ba[j];
            Hofr[index] += 2.0 * Ba[i] * Bs[j];
            Jofr[index] += 2.0 * Bs[i] * Ba[j];
            Kofr[index] += 2.0 * Bs[i] * Bs[j];
            if (D > Dmax){
                Dmax = D;
            }
        }
        if(i > 1                  && progression == 0){printf("0 > "); progression++;}
        else if(i > NumberOfAtomsRead/10.0 && progression == 1){printf("10 > "); progression++;}
        else if(i > NumberOfAtomsRead/5.00 && progression == 2){printf("20 > "); progression++;}
        else if(i > NumberOfAtomsRead/3.33 && progression == 3){printf("30 > "); progression++;}
        else if(i > NumberOfAtomsRead/2.50 && progression == 4){printf("40 > "); progression++;}
        else if(i > NumberOfAtomsRead/2.00 && progression == 5){printf("50 > "); progression++;}
        else if(i > NumberOfAtomsRead/1.67 && progression == 6){printf("60 > "); progression++;}
        else if(i > NumberOfAtomsRead/1.43 && progression == 7){printf("70 > "); progression++;}
        else if(i > NumberOfAtomsRead/1.25 && progression == 8){printf("80 > "); progression++;}
        else if(i > NumberOfAtomsRead/1.11 && progression == 9){printf("90 > "); progression++;}
        else if(i > NumberOfAtomsRead-1000 && progression ==10){printf("100!\n"); progression++;}
        fflush(stdout);
    }
    fprintf(outputFile, "# Dmax %8.2lf\n", Dmax);
    fprintf(outputFile, "# MeanVolume %8.2lf\n", MeanVolume);
    printf("# Pair distance distribution function:\n"); fprintf(outputFile, "# Pair distance distribution function:\n");
    printf("#      r [AA]   P(r) [cm^2]\n"); fprintf(outputFile, "#     r [AA]  P(r) [cm^2]   G(r) [cm^2]   H(r) [cm^2]   J(r) [cm^2]   K(r) [cm^2]\n");
    printf("%8.4g  %12.6g\n", 0.0, 0.0);
    fprintf(outputFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for(i=0;i<PointsInPofr;i++)
    {
        r=(i+0.5)*Delta_r;
        printf("%12.4g  %-12.6g\n", r, Pofr[i]);
        fprintf(outputFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", r, Pofr[i], Gofr[i], Hofr[i], Jofr[i], Kofr[i]);
    }
    printf("# Dmax %8.2lf\n\n", Dmax);
    printf("# Done\n# time: %lf sec\n\n",(double) clock()/CLOCKS_PER_SEC); /* print execution time  */
    printf("# MeanVolume %8.2lf\n", MeanVolume);
    printf("# SumOfVolume %8.2lf\n", SumOfVolume);
    
    return 0;
}
