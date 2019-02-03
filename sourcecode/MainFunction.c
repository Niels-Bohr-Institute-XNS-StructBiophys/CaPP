///                           ///
///  CaPP A.B                 ///
///  Copyright 2019           ///
///  University of Copenhagen ///
///  Niels Bohr Institute     ///
///                           ///
///  Andreas Haahr Larsen     ///
///                           ///
///  andreas.larsen@nbi.ku.dk ///
///                           ///

// Compile:
// gcc sourcecode/MainFunction.c -lm -o ./capp
// Run GUI:
// python CaPP_A.B.py (A.B is the version)

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include "Structs.h"
#include "Helpfunctions.h"
#include "GeneratePDB.h"
#include "ReadPDB.h"
#include "CalcCrossTerms.h"

int main(int argc, char **argv)
{
    clock(); /*Measure computation time */
    
            // DEFAULT PARAMETER VALUES
    
    double SolventD2O = -1.0; // D2O content in solvent (SANS)
    double Perdeuteration = 0.0; // perdeuterated fraction of the protein (SANS)
    double Perdeuteration_B = 0.0; // perdeuterated fraction of chain B (SANS)
    double Perdeuteration_C = 0.0;
    double Perdeuteration_D = 0.0;
    double Perdeuteration_E = 0.0;
    double Perdeuteration_F = 0.0;
    double Perdeuteration_G = 0.0;
    
    double PrcSucrose = 0.0; // sucrose in solvent (SAXS) in 0.01*% = 0.01*g/100ml
    double Delta_r = 1.0; //default "resolution" (binsize in r) is 1 A
    double WaterLayerContrast = 0.0; //contrast of hydration layer
    double HalfBilayerThickness = 0.0; // exclude water layer from this region
    int OPTION_d_CHOSEN = 0; // take halfbilayer thickness from the OPM database
    int OPTION_m_CHOSEN = 0; // type halfbilayer thickness manually
    int OPTION_g_CHOSEN = 0; // only calc Rg
    int OPTION_WL_CHOSEN = 0; // include a water layer
    int NO_OPTIONS_CHOSEN = 1;
    
    char *filename = argv[argc-1];
    if (strcmp(filename,"-help")==0 || strcmp(filename,"-h")==0 || strcmp(filename,"-use")==0 || strcmp(filename,"help")==0 || strcmp(filename,"h")==0 || strcmp(filename,"use")==0 ){
        printf("\n \
               *********************************************************************  \n \
               Welcome to CaPP (Calculate p(r) function from a PDB file)              \n \
               --------------------------------------------------------------------   \n \
               This program calculates the pair distance distribution function, p(r), \n \
               of a PDB file. The distances are weighted with the scattering contrast \n \
               for X-rays in an aqueous solution or for Neutrons in a D2O/H2O solvent,\n \
               as specified by the user. Proteins can be deuterated in SANS.          \n \
               A water layer can be added.                                            \n \
               --------------------------------------------------------------------   \n \
               Questions/bug-reports/suggustions/collaborations/etc, contact:         \n \
               Andreas Haahr Larsen, andreas.larsen@nbi.ku.dk                         \n \
               *********************************************************************  \n \
               Usage:\n\n \
               capp [options] NAME_of_PDB \n\n\
               Options: \n\n\
               -c [input: Contrast of water layer] \n\
               Add a water layer with (c)ontrast between 0 and 2 times the solvent scattering length.\n\
               Typically 0.1.\n\
               Default: No water layer.\n\n\
               -d [no input] \n\
               Only relevant for membrane proteins.\n\
               Removes water layer from the bilayer region. \n\
               Choose -d if the pdb is from the OPM (d)atabase, that provides the bilayer thickness. \n\n\
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
               Default: 2.0 Aangstrom. \n\n \
               -A,-B,...,-G [input: prc perdeuteration] \n\
               Perdeuteration of chain (A),...,(G). Enter percent (between 0 and 1). \n\
               For single-chain proteins (no chain labels), use -A, which is then used globally.\n\
               \n\n");
        exit(-1);
    }

    if(argc == 1){
        printf("\n\n");
        printf("\n            ******************************************************");
        printf("\n            Please provide a PDB file name.");
        printf("\n            Type capp -h for help.");
        printf("\n            ******************************************************\n\n\n\n");
        exit(-1);
    }
    
            // IMPORT PDB FILE AND ADD .pdb EXTENSION IF NEEEDED
    
    FILE *PointerToFile;

    if( (PointerToFile = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename if not given
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
    
            // OPTIONS
    
        //printf("\n");
        //printf("\n            ************* Options chosen: *************\n");

    char ch;
    const char *ValidOpts = "c:dm:s:x:r:A:B:C:D:E:F:G:g";
    char *CheckOptArg;
    while((ch=options(argc,argv,ValidOpts))!=-1)
    {
        NO_OPTIONS_CHOSEN = 0;
        switch(ch)
        {
            case 'c':
                OPTION_WL_CHOSEN = 1;
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -c, a number between -2 and 2.\n");
                    printf("\n            Input argument for option -c was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                WaterLayerContrast = char2double(OptArg);
                //printf("\n            (-c) Water Layer Contrast = %6.2f prc of solvent scattering length\n", WaterLayerContrast*100.0);
                if (WaterLayerContrast == 0) {
                    printf("\n\n            !!! NB !!! Since you have chosen 0 as water layer contrast (input -c), effectively no water layer will be added\n\n");
                }
                break;
            case 'd':
                if (OPTION_m_CHOSEN == 1) {
                    printf("\n\n\n            !!!ERROR!!! You cannot choose BOTH option -d and option -m.\n\n\n\n");
                    exit(-1);
                }
                HalfBilayerThickness = CheckHalfBilayerThickness(filename);
                //printf("\n            (-d) Bilayer Thickness from OPM = %6.2f\n", 2 * HalfBilayerThickness);
                OPTION_d_CHOSEN = 1;
                break;
            case 'm':
                if (OPTION_d_CHOSEN == 1) {
                    printf("\n\n\n            !!!ERROR!!! You cannot choose BOTH option -d and option -m.\n\n\n\n");
                    exit(-1);
                }
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -m, a bilayer thickness in Angstrom.\n");
                    printf("\n            Input argument for option -m was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                HalfBilayerThickness = char2double(OptArg)*0.5;
                //printf("\n            (-m) Manual Bilayer Thickness = %6.2f Angstrom.\n            NB: Remember to place the TMD perpendicular to the xy-plane, in z=0!\n", 2 * HalfBilayerThickness);
                OPTION_m_CHOSEN = 1;
                if (HalfBilayerThickness == 0) {
                    printf("\n\n\n            !!! Error !!! You have chosen a bilayer thickness of 0 Angstrom or not given any input after option -m. Please try again with an input. \n\n\n");
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
                //if (SolventD2O < 0.0) { printf("\n        (-s) SAXS\n"); }
                //else{ printf("\n            (-s) SANS with %4.1f prc D2O in the solvent\n", SolventD2O*100.0); }
                if (SolventD2O > 1.00){ printf("\n\n\n            !!!ERROR!!! The solvent D2O content, option -s, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n"); exit(-1);}
                break;
            case 'x':
                CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
                if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                    printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -x, a sucrose content (g/100ml) for SAXS.\n");
                    printf("\n            Input argument for option -x was \"%s\"\n\n\n\n", OptArg);
                    exit(-1);
                }
                PrcSucrose = char2double(OptArg);
                if (PrcSucrose < 0.0){
                    printf("\n\n\n            !!!ERROR!!! The percent sucrose, option -x, should be positive.\n\n\n");
                    exit(-1);}
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
                    printf("\n\n\n            !!! Error !!! You have chosen 0 Angstrom resolution or not given any resolution after option -r. Please try again with an input. \n\n\n");
                    exit(-1);
                }
                //printf("\n            (-r) Resolution = %6.4f Angstrom\n", Delta_r);
                break;
            case 'A':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -A, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -A was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -A, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -A, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'B':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -B, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -B was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -B, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration_B = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -B, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'C':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -C, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -C was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -C, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration_C = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -C, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'D':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -D, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -D was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -D, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration_D = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -D, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'E':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -E, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -E was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -E, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration_E = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -E, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'F':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -F, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -F was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -F, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration_F = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -F, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'G':
            CheckOptArg = strchr("abcdefgahcdefghijklmnopqrstuvxyz-*!%&/<>)(][{}", OptArg[0]);
            if (strcmp(OptArg,filename)==0 || CheckOptArg != NULL){
                printf("\n\n\n            !!!ERROR!!! Please provide a proper input for option -G, a perdeuteration fraction between 0 and 1 for SANS.\n");
                printf("\n            Input argument for option -G was \"%s\"\n\n\n\n", OptArg);
                exit(-1);
            }
            if (SolventD2O < 0.0) {
                printf("\n\n\n            WARNING! You have chosen the perdeuteration, option -G, and SAXS contrast, so the perdeuteration has no effect.\n");
            }
            Perdeuteration_G = char2double(OptArg);
            if (Perdeuteration > 1.00){ printf("\n\n\n            !!!ERROR!!! The perdeuteration, option -G, should be between 0 and 1 (SANS), or omitted (SAXS)\n\n\n");
                exit(-1);
            }
                break;
            case 'g':
                OPTION_g_CHOSEN = 1;
        }
    }
    //if(NO_OPTIONS_CHOSEN == 1) {printf("\n            No options chosen - default values used\n");}
    
    if( (OPTION_d_CHOSEN || OPTION_m_CHOSEN) && WaterLayerContrast == 0) {
        printf("\n\n\
                    *************************************************************\n\
                    You have chosen option -m or option -d (exclude WL at TMD),\n\
                    but the contrast of the water layer is 0, i.e. there is no WL,\n\
                    i.e. there is nothing to exclude.\n\
                    See above for instructions on how to use the program.\n\
                    *************************************************************\n\n\n\n");
        exit(-1);
    }

    //if (SolventD2O >= 0.0 && SolventD2O <= 1.0) {
        //printf("\n            SANS contrast with %2.0f prc D2O in the solvent\n", SolventD2O*100);
    //} else {
        //printf("\n            SAXS contrast (default)\n");
    //}
    
        //printf("\n            ********************************************\n");

            // CREATE PDB WITH WATER LAYER
    char *waterfilename;
    char *totalfilename;
    if(OPTION_WL_CHOSEN == 1) {
        //printf("\n\n");
        //printf("\n            ************ Adding Water Layer ************\n");
        waterfilename = GenerateWaterLayerPDB(filename, HalfBilayerThickness);
        PointerToFile = fopen(filename,"r");
        //printf("\n            Bilayer Thickness = %6.4f\n", 2 * HalfBilayerThickness);
        //printf("\n            Water has been added with contrast = %6.4f\n", WaterLayerContrast);
        //printf("\n            PDB with water only generated: %s\n", waterfilename);
        //printf("\n            ********************************************\n");
        totalfilename = GenerateTotalPDB(filename);
    }
    
            // READ DATA FOM PDB FILE
    
    int GuessNumberOfAtoms = CheckNumberOfAtomsInPDBFile(filename)*2;
    struct Atom *Atoms = (struct Atom *)  calloc(GuessNumberOfAtoms, sizeof(struct Atom));
    double *Da = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // Da = distance from COM to atom
    double *Ba = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // Ba = scattering length vector for atom
    double *Bs = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // Bs = scattering length vector for excluded solvent
    double *dB = (double *) calloc( GuessNumberOfAtoms, sizeof(double)); // dB = Ba - Bs
    int i_start = 0;
    int PointsInPofr = 0;
    double DmaxPDB = 0.0;
    double DmaxWL = 0.0;
    double Dmax = 0.0;
    double SumAtomWeightPDB = 0.0;
    double SumAtomWeightWL = 0.0;
    double RgPDB = 0.0;
    double RgWL = 0.0;
    double MeanVolumePDB = 0.0;
    double MeanVolumeWL = 0.0;
    int TOTAL_pr = 0;
    
    int NumberOfAtomsPDB = ReadPDB(filename, Da, Ba, Bs, dB, Atoms, WaterLayerContrast, SolventD2O, Perdeuteration, Perdeuteration_B, Perdeuteration_C, Perdeuteration_D, Perdeuteration_E, Perdeuteration_F, Perdeuteration_G, PrcSucrose, Delta_r, HalfBilayerThickness, OPTION_g_CHOSEN, OPTION_WL_CHOSEN, i_start, &PointsInPofr, &DmaxPDB, &SumAtomWeightPDB, &RgPDB, &MeanVolumePDB,TOTAL_pr);

    if (OPTION_WL_CHOSEN == 1){
        
        int NumberOfAtomsWL = ReadPDB(waterfilename, Da, Ba, Bs, dB, Atoms, WaterLayerContrast, SolventD2O, Perdeuteration, Perdeuteration_B, Perdeuteration_C, Perdeuteration_D, Perdeuteration_E, Perdeuteration_F, Perdeuteration_G, PrcSucrose, Delta_r, HalfBilayerThickness, OPTION_g_CHOSEN, OPTION_WL_CHOSEN, NumberOfAtomsPDB, &PointsInPofr, &DmaxWL, &SumAtomWeightWL, &RgWL, &MeanVolumeWL,TOTAL_pr);
        
        CalcCrossTerms(filename, NumberOfAtomsPDB, NumberOfAtomsWL, Ba, Bs, dB, Atoms, WaterLayerContrast, Delta_r, PointsInPofr);
        
        
        // generate filenames
        char buf[99];
        char buf_tmp[99];
        if(SolventD2O >= 0.0 || SolventD2O >= 1.0){
            sprintf(buf,"_N%1.0f_P%1.0f",SolventD2O*100,Perdeuteration*100);
            if (Perdeuteration_B > 0.0 && Perdeuteration_B <= 1.0) {
                sprintf(buf_tmp,"_B%1.0f",Perdeuteration_B*100.0);
                strcat(buf,buf_tmp);
            }
            if (Perdeuteration_C > 0.0 && Perdeuteration_C <= 1.0) {
                sprintf(buf_tmp,"_C%1.0f",Perdeuteration_C*100.0);
                strcat(buf,buf_tmp);
            }
            if (Perdeuteration_D > 0.0 && Perdeuteration_D <= 1.0) {
                sprintf(buf_tmp,"_D%1.0f",Perdeuteration_D*100.0);
                strcat(buf,buf_tmp);
            }
            if (Perdeuteration_E > 0.0 && Perdeuteration_E <= 1.0) {
                sprintf(buf_tmp,"_E%1.0f",Perdeuteration_E*100.0);
                strcat(buf,buf_tmp);
            }
            if (Perdeuteration_F > 0.0 && Perdeuteration_F <= 1.0) {
                sprintf(buf_tmp,"_F%1.0f",Perdeuteration_F*100.0);
                strcat(buf,buf_tmp);
            }
            if (Perdeuteration_G > 0.0 && Perdeuteration_G <= 1.0) {
                sprintf(buf_tmp,"_G%1.0f",Perdeuteration_G*100.0);
                strcat(buf,buf_tmp);
            }
        } else {
            sprintf(buf,"_X%1.0f",PrcSucrose);
        }
        char *input1 = GetCStringBeforeLastDelimiter(filename, '.');        // extracts 1HEJ from 1HEJ.pdb
        strcat(input1, buf);                                 //adds e.g. _X10 to 1HEJ -> 1HEJ_X10
        char *protein_pr_filename = strcat(input1, "_pr.dat");               // adds _pr.dat to 1HEJ_X10 -> 1HEJ_X10_pr.dat
        char *input2 = GetCStringBeforeLastDelimiter(filename, '.');         // extracts 1HEJ from 1HEJ.pdb
        char *wl_pr_filename = strcat(input2, "_w_only_pr.dat");             // adds _w_only_pr.dat to 1HEJ -> 1HEJ_w_only_pr.dat
        char *input3 = GetCStringBeforeLastDelimiter(filename, '.');         // extracts 1HEJ from 1HEJ.pdb
        char *cross_pr_filename = strcat(input3, "_cross_pr.dat");           // adds _cross_pr.dat to 1HEJ -> 1HEJ_cross_pr.dat
        char *input4 = GetCStringBeforeLastDelimiter(totalfilename, '.');   // extracts 1HEJ_w from 1HEJ_w.pdb
        strcat(input4, buf);                                 //adds e.g. _X10 to 1HEJ_w -> 1HEJ_w_X10
        char *total_pr_filename = strcat(input4, "_pr.dat");                 // adds _pr.dat to 1HEJ_w_X10 -> 1HEJ_w_X10_pr.dat
        
        // import from water, protein and cross term files. Sum up the terms to get the total pr.
        FILE * ProteinFile;
        FILE * WLFile;
        FILE * CrossFile;
        FILE * TotalFile;
        
        ProteinFile = fopen(protein_pr_filename,"r"); // PDB file
        WLFile = fopen(wl_pr_filename,"r"); // WL file
        CrossFile = fopen(cross_pr_filename,"r"); // cross file
        
        char buffer[256];
        double r;
        double *pr = (double *) calloc( PointsInPofr, sizeof(double));
        double *gr = (double *) calloc( PointsInPofr, sizeof(double));
        double *hr = (double *) calloc( PointsInPofr, sizeof(double));
        double *jr = (double *) calloc( PointsInPofr, sizeof(double));
        double *kr = (double *) calloc( PointsInPofr, sizeof(double));
        
        double c1, c2, c3, c4, c5, c6;
        int i = 0;
        while(fgets(buffer,sizeof(buffer),ProteinFile)!=NULL){
            if(sscanf(buffer, "      %lf %lf %lf %lf %lf %lf",&c1, &c2, &c3, &c4, &c5, &c6 )==6) {
                pr[i] += c2;
                gr[i] += c3;
                hr[i] += c4;
                jr[i] += c5;
                kr[i] += c6;
                i++;}
        }
        i = 0;
        while(fgets(buffer,sizeof(buffer),WLFile)!=NULL){
            if(sscanf(buffer, "%lf %lf %lf %lf %lf %lf",&c1, &c2, &c3, &c4, &c5, &c6 )==6) {
                pr[i] += c2;
                gr[i] += c3;
                hr[i] += c4;
                jr[i] += c5;
                kr[i] += c6;
                i++;}
        }
        i = 0;
        while(fgets(buffer,sizeof(buffer),CrossFile)!=NULL){
            if(sscanf(buffer, "%lf %lf %lf %lf %lf %lf",&c1, &c2, &c3, &c4, &c5, &c6 )==6) {
                pr[i] += c2;
                gr[i] += c3;
                hr[i] += c4;
                jr[i] += c5;
                kr[i] += c6;
                i++;}
        }
        
        double MeanVolumeTot  = 0.0;
        double RgTot = 0.0;
        double SumAtomWeightTot = 0.0;
        double DmaxTot = 0.0;
        TOTAL_pr = 1;
        
        // Generate pr for protein + water and print headerlines until (and inclusive) Rg//
        NumberOfAtomsPDB = ReadPDB(totalfilename, Da, Ba, Bs, dB, Atoms, WaterLayerContrast, SolventD2O, Perdeuteration, Perdeuteration_B, Perdeuteration_C, Perdeuteration_D, Perdeuteration_E, Perdeuteration_F, Perdeuteration_G, PrcSucrose, Delta_r, HalfBilayerThickness, OPTION_g_CHOSEN, OPTION_WL_CHOSEN, i_start, &PointsInPofr, &DmaxTot, &SumAtomWeightTot, &RgTot, &MeanVolumeTot,TOTAL_pr);
        
        // Open file with pr for protein + water
        TotalFile = fopen(total_pr_filename,"a"); // total file (append)
        
        // Find largest Dmax from either WL file or from PDB file, and print it to total pr file //
        Dmax = DmaxWL;
        if (DmaxPDB > DmaxWL){Dmax = DmaxPDB;}
        fprintf(TotalFile, "# Dmax: %8.2lf\n", Dmax);
        
        // print mean volume for protein alone (not much point in printing volume for water beads)
        //fprintf(TotalFile, "# Mean volume (without water layer): %8.2lf A**3. Total volume (without water layer): %8.2lf A**3\n", MeanVolumePDB,MeanVolumePDB*NumberOfAtomsPDB);
        fprintf(TotalFile, "# Mean volume (without water layer): %8.2lf\n", MeanVolumePDB);
        
        // print total pr function = sum of water + pure protein + cross terms
        fprintf(TotalFile, "# Pair distance distribution function:\n");
        fprintf(TotalFile, "#     r [AA]  P(r) [cm^2]   G(r) [cm^2]   H(r) [cm^2]   J(r) [cm^2]   K(r) [cm^2]\n");
        fprintf(TotalFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        
        for(i=1;i<PointsInPofr+1;i++)
        {
            r=(i-0.5)*Delta_r;
            if (i >= PointsInPofr) {
                fprintf(TotalFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", r, 0.0, 0.0, 0.0, 0.0, 0.0);
            } else {
                fprintf(TotalFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", r, pr[i], gr[i], hr[i], jr[i], kr[i]);
            }
        }
        fclose(ProteinFile);
        fclose(WLFile);
        fclose(CrossFile);
        fclose(TotalFile);
        
    } // end if WL
    
    // Write file with dB and Da, used to calculate A00 and beta for decoupling approximation
    char *inputPDB = GetCStringBeforeLastDelimiter(filename, '.'); //extracts 1HEJ from 1HEJ.pdb
    char *betafilename = strcat(inputPDB, "_beta.list");
    FILE * BetaFile;
    BetaFile = fopen(betafilename,"w"); // total file (write)
    fprintf(BetaFile, "# List of distances and excess scattering length for each atom.\n");
    fprintf(BetaFile, "# Used to calculate beta(q) for the decoupling approximation:\n");
    fprintf(BetaFile, "# r [A]          Delta b [cm]\n");
    for(int i=0;i<NumberOfAtomsPDB;i++)
    {
        fprintf(BetaFile, "%12.4g %12.4g\n", Da[i], dB[i]);
    }
    fclose(BetaFile);
    return 0;
}
