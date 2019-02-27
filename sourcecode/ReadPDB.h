int ReadPDB(char *filename, double *Da, double *Ba, double *Bs, double *dB, struct Atom *Atoms, double WaterLayerContrast, double SolventD2O, double Perdeuteration, double Perdeuteration_B, double Perdeuteration_C, double Perdeuteration_D, double Perdeuteration_E, double Perdeuteration_F, double Perdeuteration_G, double PrcSucrose, double Delta_r, double HalfBilayerThickness, int OPTION_g_CHOSEN, int OPTION_WL_CHOSEN, int OPTION_EXPLICIT_H_CHOSEN, int i_start, int * PointsInPofrPointer, double * DmaxPointer, double * SumAtomWeightPointer, double * RgPointer, double *MeanVolumePointer, int TOTAL_pr)
{
    double NonExchNH = 0.1; // cannot be changed (in GUI or in batch mode).
                            // give rise to minor, almost q-independent contribution,
                            // since NH groups are distributed all over the protein.
        
    int PointsInPofr = * PointsInPofrPointer;
    double Dmax = * DmaxPointer;
    double Rg = * RgPointer;
    double Rg2 = pow(Rg,2.0);
    FILE *PointerToFile;
    PointerToFile = fopen(filename,"r");

    // READ DATA FOM PDB FILE
    double Sum_dB = 0.0, Sum_B = 0.0;
    double SumAtomWeight = * SumAtomWeightPointer;
    double a = 1.0/(1.0 + WaterLayerContrast);
    double n_W = 4.133; // number of watermolecules per water residue
    double r, D;
    int index;
    enum {HYDROGEN, DEUTERIUM, HD, CARBON, NITROGEN, OXYGEN, FLUORINE, SODIUM, MAGNESIUM, PHOSPHORUS, SULFUR, CHLORINE, CALCIUM, MANGANESE, IRON, COPPER, ZINK, GOLD, WATER, WATERBULK, UNKNOWN};
    
    // SUCROSE CONTRAST VARIATION
    double VolSolution = 0.1; // liter
    double MW_Suc = 342.3; // Da = g/mol, sucrose is C12 H22 O11
    double MoleSuc = PrcSucrose/MW_Suc; // mol
    //printf("Molar sucrose = %f\n", MoleSuc/VolSolution);
    double VolSuc = 358.16; // aa**3, from Mia's WillItFit code (constraints.h, variable VolumeOfContrast)
    double VolH20 = 30.0; // aa**3
    double VolSucLiter = VolSuc * 1e-27;
    double NA = 6.022e23; // Avogadros number
    double VolH2OinSolution = VolSolution - VolSucLiter * MoleSuc * NA;
    double VolSucinSolution = VolSolution - VolH2OinSolution;
    double SucPerSolv = (VolSucinSolution/VolSuc)/(VolH2OinSolution/VolH20);
    //printf("Vol H20 in solution = %f, vol suc in sol = %f\n",VolH2OinSolution, VolSucinSolution);
    double ScatLenSuc = 182 * 2.82e-13;
    
    // SOLVENT VOLUME (with sucrose).
    double VolumeSolvent = 30.0;
    if(SolventD2O<0.0){VolumeSolvent = VolumeSolvent + SucPerSolv * VolSuc;}
    
    // ATOM LIST
    /*                                      Scattering-length                                                                       */
    /*                   x    y    z       X-ray[cm]        Neutron[cm]                                    V[AA^3] Mw[Da]  Name*/
    struct Atom element[] = {
        [HYDROGEN]   = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    3.741e-13,                                       5.15,   1.0,   'A'} ,
        [DEUTERIUM]  = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13,                                       5.15,   2.0,   'A'} ,
        [HD]         = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13*SolventD2O-3.741e-13*(1.0-SolventD2O), 5.15,   1.5,   'A'} ,
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
        [WATER]      = { 0.0, 0.0, 0.0,n_W*8 * 2.82e-13,n_W*5.803e-13,                        n_W*(a*30.00-2*5.15),  0.0,   'A'} ,
        [WATERBULK]  = { 0.0, 0.0, 0.0,   10 * 2.82e-13,    5.803e-13+2.0*(6.671e-13*SolventD2O-3.741e-13*(1.0-SolventD2O)),    30.00,   0.0,   'A'} ,
        [UNKNOWN]    = { 0.0, 0.0, 0.0,    0 * 2.82e-13,    0.000e-13,                                       1.00,   0.0,   'A'}
    };
    // Note: [WATER] is without hydrogen/deuterium... will be added later on (implicit way)
    // Note: weight of [WATER] is set to 0, since the WL should not be included in the calculation of the atomic weight
    // Volumes from:
    //    (H,C,D,N,O):                              Fraser, J Appl Cryst(1978), 11, p693.
    //                                               Experimentally determined values
    //    (P, S, Mn, Ca, F, Mg, Cl, Cu, Fe, Zn):    Batsanov, Inorg Mat(2001), vol37, no9, p871.
    //                                               Van der Waals radii, mean of values from Table 1 and 2 (P from Table 5)
    // HD is a libile Hydrogen, that can be exchanged with Deuterium
    //
    
    // SOLVENT SCATT LENGTH (with sucrose)
    double ScatteringLengthSolvent = element[WATERBULK].NeutronScatteringLength;
    if(SolventD2O<0.0){
        ScatteringLengthSolvent = element[WATERBULK].XRayScatteringLength;
        ScatteringLengthSolvent = ScatteringLengthSolvent + SucPerSolv * ScatLenSuc;
    }
    double rhoSolvent = ScatteringLengthSolvent/VolumeSolvent;
    
    int R,j,u=0; // R = residue number, j = index, u = unknown atom
    char buffer[256];
    
    //int GuessNumberOfAtoms = CheckNumberOfAtomsInPDBFile(filename);
    
    struct Atom COM   = {0,0,0};
    double SqDist, SqDistMax=0, GuessOnDmax;
    double SumOfVolume = 0.0;
    
    char AminoName, LongAtomName, AlternativeAtomPosition, Chain;
    
    double n_H; // number of non-exchangable (strongly bound) Hydrogen
    double n_D; // number of exchangeable Hydrogen
    double n_DD; // number of non-exchangeable (strongly bound) Deuterium (after perdeuteration)
    double n_DB; // number exchangeable Hydrogen attacehd to water bead
    
    // VOLUME CORRECTION FACTORS (to be multiplied on volume of excluded water by atom, relative to the atom in a protein)
    int DNA,RNA,LIP,SUC; // is it a lip, a DNA, RNA or sugar atom?
    double LipCorrFactor = 1.06; // Calculated with Matlab script: CalcLipVolume.m
    double DNACorrFactor = 1.00; // not determined
    double RNACorrFactor = 1.00; // not determined
    double SucCorrFactor = 1.00; // not determined
    //double SucCorrFactor = 1.00/1.49; // Biophysical Journal, 97, 2009, 1445–1453. Interrelationship of Steric Stabilization and Self-Crowding of a Glycosylated Protein. R. Høiberg-Nielsen, P. Westh, L. K. Skov, and L. Arleth

    int NumberOfIgnoredHD = 0;
    int NumberOfAlternativeAtoms = 0;
    int NumberOfOtherLines = 0;
    int i = i_start;
    while(fgets(buffer,sizeof(buffer),PointerToFile)!=NULL){
        
        // read atom name:
        if (sscanf(buffer, "ATOM%*73c%s",&Atoms[i].Name) == 1 || sscanf(buffer, "HETATM%*71c%s",&Atoms[i].Name) == 1) {}
        
        // read alternative atom pos:
        if (sscanf(buffer, "ATOM%*12c%c",&AlternativeAtomPosition) == 1 || sscanf(buffer, "HETATM%*10c%c",&AlternativeAtomPosition) == 1) {}
        
        // read chain:
        if (sscanf(buffer, "ATOM%*17c%s",&Chain) == 1 || sscanf(buffer, "HETATM%*15c%s",&Chain) == 1) {} else {sscanf("U","%s",&Chain);} // U for unknown chain
        
        // ignore H and D:
        if ((strcmp(&Atoms[i].Name,"H") == 0 || strcmp(&Atoms[i].Name,"D") == 0) && OPTION_EXPLICIT_H_CHOSEN == 0) {NumberOfIgnoredHD++;}
        // ignore alternative position atoms (only use position A):
        else if (AlternativeAtomPosition == 'B') {NumberOfAlternativeAtoms++;}
        // read all other atoms:
        else if (sscanf(buffer,"ATOM%*18c%d%*4c%lf%lf%lf",&R,&Atoms[i].x,&Atoms[i].y,&Atoms[i].z)==4 ||
            sscanf(buffer,   "HETATM%*16c%d%*4c%lf%lf%lf",&R,&Atoms[i].x,&Atoms[i].y,&Atoms[i].z)==4 ) {
            n_D = 0;
            n_H = 0;
            n_DD = 0;
            n_DB = 0;
            DNA = 0;
            RNA = 0;
            LIP = 0;
            //printf("%s",buffer);
            //printf("    R       ,            : %d\n",R);
            //printf("    x,                   : %f\n",Atoms[i].x);
            //printf("    y,                   : %f\n",Atoms[i].y);
            //printf("    z,                   : %f\n",Atoms[i].z);
            //printf("    Chain,         repeat: %s\n",&Chain);
            //printf("    Buffer[12-15],       : %c%c%c%c\n",buffer[12],buffer[13],buffer[14],buffer[15]);
            
            // Read element and give values from atom list
            if (buffer[13] == 'C') {
                CopyPhysicalParameters(&Atoms[i], &element[CARBON]);
                //printf("Carbon\n");
            } else if (buffer[13] == 'H') {
                CopyPhysicalParameters(&Atoms[i], &element[HYDROGEN]);
                //printf("Hydrogen\n");
            } else if (buffer[13] == 'D') {
                CopyPhysicalParameters(&Atoms[i], &element[DEUTERIUM]);
                //printf("Deuterium\n");
            } else if (buffer[12] != 'M' && buffer[12] != 'Z' && buffer[13] == 'N') {
                CopyPhysicalParameters(&Atoms[i], &element[NITROGEN]);
                //printf("Nitrogen\n");
            } else if (buffer[13] == 'O') {
                CopyPhysicalParameters(&Atoms[i], &element[OXYGEN]);
                //printf("Oxygen\n");
            } else if (buffer[13] == 'F') {
                CopyPhysicalParameters(&Atoms[i], &element[FLUORINE]);
                //printf("Flourine\n");
            } else if (buffer[12] == 'N' && buffer[13] == 'A') {
                CopyPhysicalParameters(&Atoms[i], &element[SODIUM]);
                //printf("Sodium\n");
            } else if (buffer[12] == 'M' && buffer[13] == 'G') {
                CopyPhysicalParameters(&Atoms[i], &element[MAGNESIUM]);}
                //printf("Magnesium\n");
            else if (buffer[13] == 'P') {
                CopyPhysicalParameters(&Atoms[i], &element[PHOSPHORUS]);}
                //printf("Phosporus\n");
            else if (buffer[13] == 'S') {
                CopyPhysicalParameters(&Atoms[i], &element[SULFUR]);}
                //printf("Sulfur\n");
            else if (buffer[12] == 'C' && buffer[13] == 'L') {
                CopyPhysicalParameters(&Atoms[i], &element[CHLORINE]);}
                //printf("Chlorine\n");
            else if (buffer[12] == 'C' && buffer[13] == 'A') {
                CopyPhysicalParameters(&Atoms[i], &element[CALCIUM]);
                //printf("Calcium\n");
            } else if (buffer[12] == 'M' && buffer[13] == 'N') {
                CopyPhysicalParameters(&Atoms[i], &element[MANGANESE]);
                //printf("Manganese\n");
            } else if (buffer[12] == 'F' && buffer[13] == 'E') {
                CopyPhysicalParameters(&Atoms[i], &element[IRON]);
                //printf("Iron\n");
            } else if (buffer[12] == 'C' && buffer[13] == 'U') {
                CopyPhysicalParameters(&Atoms[i], &element[COPPER]);
                //printf("Copper\n");
            } else if (buffer[12] == 'Z' && buffer[13] == 'N') {
                CopyPhysicalParameters(&Atoms[i], &element[ZINK]);
                //printf("Zink\n");
            } else if (buffer[13] == 'Q') {
                CopyPhysicalParameters(&Atoms[i], &element[WATER]); n_D = 2.0 * n_W; n_DB = 2.0 * n_W;
                //printf("WaterLayerBead\n");
            } else {
                CopyPhysicalParameters(&Atoms[i], &element[UNKNOWN]); u++;
                //printf("UnknownAtom\n");
            }
            // Add H/D. According to http://www.rcsb.org/ligand
            if (sscanf(buffer,"ATOM%*13c%s", &AminoName) == 1 || sscanf(buffer,"HETATM%*11c%s", &AminoName) == 1)
            {
                // ADD H/D to AMINO ACIDS
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
                    if (strcmp(&LongAtomName,"N")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CA")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CB")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CG1")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"CG2")==0) {n_H = 3;}
                }
                // ADD H/D to DNA nucleotides
                else if (strcmp(&AminoName,"DA")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C5")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C6")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"N1")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N3")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N6")==0) {n_D = 2*(1.0-NonExchNH); n_H = 2*NonExchNH;}
                    //else if (strcmp(&LongAtomName,"N7")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N9")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    DNA = 1;
                }
                else if (strcmp(&AminoName,"DC")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"C2")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"N1")==0) {n_D = 0}
                    //else if (strcmp(&LongAtomName,"N3")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N4")==0) {n_D = 2*(1.0-NonExchNH); n_H = 2*NonExchNH;}
                    //else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    DNA = 1;
                }
                else if (strcmp(&AminoName,"DG")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"C2")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C5")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C6")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"N1")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"N2")==0) {n_D = 2*(1.0-NonExchNH); n_H = 2*NonExchNH;}
                    //else if (strcmp(&LongAtomName,"N3")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N7")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N9")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O6")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    DNA = 1;
                }
                else if (strcmp(&AminoName,"DT")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 3;}
                    //else if (strcmp(&LongAtomName,"N1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N3")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    //else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O4")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    DNA = 1;
                }
                // ADD H/D to RNA nucleotides
                else if (strcmp(&AminoName,"A")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C5")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C6")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"N1")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N3")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N6")==0) {n_D = 2*(1.0-NonExchNH); n_H = 2*NonExchNH;}
                    //else if (strcmp(&LongAtomName,"N7")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N9")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    RNA = 1;
                }
                else if (strcmp(&AminoName,"C")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"C2")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"N1")==0) {n_D = 0}
                    //else if (strcmp(&LongAtomName,"N3")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N4")==0) {n_D = 2*(1.0-NonExchNH); n_H = 2*NonExchNH;}
                    //else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    RNA = 1;
                }
                else if (strcmp(&AminoName,"G")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"C2")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C5")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C6")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"N1")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"N2")==0) {n_D = 2*(1.0-NonExchNH); n_H = 2*NonExchNH;}
                    //else if (strcmp(&LongAtomName,"N3")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N7")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N9")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O6")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    RNA = 1;
                }
                else if (strcmp(&AminoName,"U")==0) {
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4'")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5'")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"C4")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"N1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N3")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    //else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O4")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O2'")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3'")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O4'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"O5'")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"OP1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"OP2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"OP3")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_D = 0;}
                    RNA = 1;
                }
                
                // ADD H/D to Sugars
                else if (strcmp(&AminoName,"Fuc")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O5")==0) {n_D = 0;}
                }
                else if (strcmp(&AminoName,"Gal")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3")==0) {n_D = 0.5;} //else if (strcmp(&LongAtomName,"O3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O5")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O6")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"Glc")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 0;} //else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O5")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O6")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"N2")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"CN2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"CAN2")==0) {n_H = 3;}
                    //else if (strcmp(&LongAtomName,"OCN2")==0) {n_D = 0;}
                }
                else if (strcmp(&AminoName,"Man")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O2")==0) {n_D = 0.5;} //else if (strcmp(&LongAtomName,"O2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O3")==0) {n_D = 0.5;} //else if (strcmp(&LongAtomName,"O3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 0.5;} //else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O5")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O6")==0) {n_D = 0.5;} //else if (strcmp(&LongAtomName,"O6")==0) {n_D = 1;}
                }
                else if (strcmp(&AminoName,"NeuPDB")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C3")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C5")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C7")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C8")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C10")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C16")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C26")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C28")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"O16")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"O11")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"O21")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O22")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O24")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"N2")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    //else if (strcmp(&LongAtomName,"N4")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N6")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"N9")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N13")==0) {n_D = 2.0*(1.0-NonExchNH); n_H = 2.0*NonExchNH;}
                }
                else if (strcmp(&AminoName,"Neu")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C2")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C9")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"CAN5")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"O4")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O7")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O8")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"O9")==0) {n_D = 1;}
                    //else if (strcmp(&LongAtomName,"1O1")==0) {n_D = 0;}
                    //else if (strcmp(&LongAtomName,"2O1")==0) {n_D = 0;}
                    else if (strcmp(&LongAtomName,"N5")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                }
                else if (strcmp(&AminoName,"UZ9")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 1;}
                }
                else if (strcmp(&AminoName,"CRS")==0) {
                    SUC = 1;
                    sscanf(buffer,"ATOM%*9c%3s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"C3")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C4")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"O1")==0) {n_D = 1;}
                }
                // ADD H/D to RNA to other heteroatomic molecueles
                else if (strcmp(&AminoName,"HOH")==0) {
                    if (OPTION_EXPLICIT_H_CHOSEN) {
                        CopyPhysicalParameters(&Atoms[i], &element[WATERBULK]);
                    } else {
                        n_D = 2; // HOH is a water atom in the crystal
                        Atoms[i].Volume = 30.0 - 2*5.15; // implicit solvent, therefore the volume of two H/D is subtracted
                    }
                }
                else if (strcmp(&AminoName,"ZK1")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"NAP")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
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
                    if (strcmp(&LongAtomName,"N4")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C5")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C6")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C7")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C8")==0) {n_H = 2;}
                }
                else if (strcmp(&AminoName,"NAG")==0) {
                    sscanf(buffer,"HETATM%*7c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
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
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
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
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
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
                    if (strcmp(&LongAtomName,"N2")==0) {n_D = 1.0-NonExchNH; n_H = NonExchNH;}
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
                else if (strcmp(&AminoName,"POPE")==0) {
                    LIP = 1;
                    sscanf(buffer,"HETATM%*6c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"NH3")==0) {n_D = 3.0*(1.0-NonExchNH); n_H = 3.0*NonExchNH;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"P")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"O11")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"O12")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"O13")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"O14")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    //else if (strcmp(&LongAtomName,"O21")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C21")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"O22")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C22")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    //else if (strcmp(&LongAtomName,"O31")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"O32")==0) {n_H = 0;}
                    //else if (strcmp(&LongAtomName,"C31")==0) {n_H = 0;}
                    else if (strcmp(&LongAtomName,"C32")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C23")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C25")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C26")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C27")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C28")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C29")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"0C21")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"1C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"2C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"3C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"4C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"5C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"6C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"7C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"8C21")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C33")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C34")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C35")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C36")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C37")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C38")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C39")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"0C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"1C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"2C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"3C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"4C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"5C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"6C31")==0) {n_H = 3;}
                }
                else if (strcmp(&AminoName,"POPG")==0) {
                    LIP = 1;
                    sscanf(buffer,"HETATM%*6c%s", &LongAtomName);
                    if (strcmp(&LongAtomName,"C1")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C1H")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"OC3")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C12")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"OC2")==0) {n_D = 1;}
                    else if (strcmp(&LongAtomName,"C11")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C22")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C2")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"C22")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C3")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C32")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C23")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C24")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C25")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C26")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C27")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C28")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C29")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"0C21")==0) {n_H = 1;}
                    else if (strcmp(&LongAtomName,"1C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"2C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"3C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"4C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"5C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"6C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"7C21")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"8C21")==0) {n_H = 3;}
                    else if (strcmp(&LongAtomName,"C33")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C34")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C35")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C36")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C37")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C38")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"C39")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"0C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"1C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"2C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"3C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"5C31")==0) {n_H = 2;}
                    else if (strcmp(&LongAtomName,"6C31")==0) {n_H = 3;}
                }
            }
            if (OPTION_EXPLICIT_H_CHOSEN) {n_H=0; n_D=n_DB; n_DD=0;} // use only implicit H for added water layer beads
            // Update volume and scattering length with implicit H and D (adjust if DNA or RNA)
            Atoms[i].Volume += (n_D + n_H) * 5.15;
            
            // Volume corrections factors for DNA, RNA and LIP
            if (DNA == 1) {Atoms[i].Volume *= DNACorrFactor;}
            else if (RNA == 1) {Atoms[i].Volume *= RNACorrFactor;}
            else if (LIP == 1) {Atoms[i].Volume *= LipCorrFactor;}
            else if (SUC == 1) {Atoms[i].Volume *= SucCorrFactor;}
            
            SumOfVolume += Atoms[i].Volume;
            
            if(SolventD2O < 0.0){
                Atoms[i].XRayScatteringLength += (n_D + n_H) * element[HD].XRayScatteringLength;
                Ba[i] = Atoms[i].XRayScatteringLength;
            }
            else {
                if (strcmp(&Chain,"B")==0 && Perdeuteration_B >= 0.0) {n_DD = Perdeuteration_B * n_H;}
                else if (strcmp(&Chain,"C")==0 && Perdeuteration_C >= 0.0) {n_DD = Perdeuteration_C * n_H;}
                else if (strcmp(&Chain,"D")==0 && Perdeuteration_D >= 0.0) {n_DD = Perdeuteration_D * n_H;}
                else if (strcmp(&Chain,"E")==0 && Perdeuteration_E >= 0.0) {n_DD = Perdeuteration_E * n_H;}
                else if (strcmp(&Chain,"F")==0 && Perdeuteration_F >= 0.0) {n_DD = Perdeuteration_F * n_H;}
                else {n_DD = Perdeuteration * n_H;}
                n_H = n_H - n_DD;
                Atoms[i].NeutronScatteringLength += n_D * element[HD].NeutronScatteringLength - n_H * element[HYDROGEN].NeutronScatteringLength + n_DD * element[DEUTERIUM].NeutronScatteringLength;
                Ba[i] = Atoms[i].NeutronScatteringLength;
            }
            Bs[i] = rhoSolvent * Atoms[i].Volume;
            dB[i] = Ba[i] - Bs[i];
            Sum_dB += dB[i];
            Sum_B += Ba[i];
            SumAtomWeight += Atoms[i].Weight;
            COM.x += dB[i] * Atoms[i].x;
            COM.y += dB[i] * Atoms[i].y;
            COM.z += dB[i] * Atoms[i].z;
            
            i++;
        } // end read atom
        else {NumberOfOtherLines++;}
    } // end While loop
    fclose(PointerToFile);
    
    int TotalNumberOfAtoms = i;
    int NumberOfAtomsRead = i-i_start;
    double MeanVolume = SumOfVolume/NumberOfAtomsRead;
    double p_PDB = Sum_B/SumOfVolume;

    if (Sum_dB == 0){
        COM.x = 0.0;
        COM.y = 0.0;
        COM.z = 0.0;
    } else {
    COM.x /= Sum_dB;
    COM.y /= Sum_dB;
    COM.z /= Sum_dB;
    }
    
    
    // CALFULATE RG AND CALCULATE DISTANCE FOR EACH ATOM
    double SqDistMaxPrevious = 0.0;
    for(i=i_start;i<TotalNumberOfAtoms;i++){
        SqDist = SquareDistance(Atoms[i],COM);
        Da[i] = sqrt(SqDist);
        Rg2 += dB[i]*SqDist;
        SqDistMax = max(SqDistMax,SqDist);
        SqDistMaxPrevious = SqDistMax;
    }
    if (Sum_dB == 0.0){
        Rg = 0.0;
    } else {
        Rg = sqrt(Rg2/Sum_dB);
    }

    PointsInPofr = *PointsInPofrPointer;
    GuessOnDmax = 2*sqrt(SqDistMax);
    if (PointsInPofr == 0){
        PointsInPofr = floor(1.2 * GuessOnDmax/Delta_r);
    }
    
    // CHECK IF FILE FOR p(r) FUNCTION EXIST
    FILE *outputFile;
    char *inputPDB = GetCStringBeforeLastDelimiter(filename, '.'); //extracts 1HEJ from 1HEJ.pdb
    char buf[50];
    char buf_tmp[99];
    if(i_start == 0){
        if(SolventD2O >= 0.0 && SolventD2O <= 1.0){
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
            strcat(buf,"_pr.dat");
        } else {
            sprintf(buf,"_X%1.0f_pr.dat",PrcSucrose);
        }
    } else {sprintf(buf,"_pr.dat");}
    char *outputfilename = strcat(inputPDB, buf); // adds string to 1HEJ
    
    if (i_start == 0 && fopen(outputfilename,"r") && TOTAL_pr != 1) {
        if(OPTION_g_CHOSEN == 1){
            printf("\n\n###########################################\n");
            printf("# Radius of gyration    %8.2lf\n", Rg);
            printf("###########################################\n\n");
            exit(-1);
        }
    } else {
        outputFile = fopen(outputfilename,"w");
        
        // PRINT OUTPUT TO TERMINAL AND p(r) FILE
        fprintf(outputFile, "# Location of PDB-file: %s\n", filename);
        if (SolventD2O >= 0.0 && SolventD2O <= 1.0) {
            fprintf(outputFile, "# SANS with %1.0f prc D2O in the solvent\n", SolventD2O*100);
        } else {
            fprintf(outputFile, "# SAXS contrast with %1.0f prc sucrose in the solvent\n",PrcSucrose);
        }
        if (i_start == 0 && TOTAL_pr == 0) {
            fprintf(outputFile, "# p(r) for PDB without water layer\n");
        } else if (OPTION_WL_CHOSEN == 1 && i_start > 0){
            fprintf(outputFile, "# p(r) for water layer with contrast = %1.0lf prc of solvent scattering length, excluded from hydrophobic bilayer with thickness %1.0lf A\n", WaterLayerContrast*100, HalfBilayerThickness*2);
        } else {
            fprintf(outputFile, "# p(r) for PDB and water layer with contrast = %1.0lf prc of solvent scattering length, excluded from hydrophobic bilayer with thickness %1.0lf A\n", WaterLayerContrast*100, HalfBilayerThickness*2);
        }
        fprintf(outputFile, "# Atomic weight: %g\n", SumAtomWeight);
        
        fprintf(outputFile, "# PointsInPofr: %d\n", PointsInPofr+1);
        
        double *Pofr = (double *) calloc( PointsInPofr, sizeof(double));
        double *Gofr = (double *) calloc( PointsInPofr, sizeof(double));
        double *Hofr = (double *) calloc( PointsInPofr, sizeof(double));
        double *Jofr = (double *) calloc( PointsInPofr, sizeof(double));
        double *Kofr = (double *) calloc( PointsInPofr, sizeof(double));
        double *Aofr = (double *) calloc( PointsInPofr, sizeof(double));
        fprintf(outputFile, "# Radius of gyration: %8.2lf\n",Rg);
        if(OPTION_g_CHOSEN == 1 && OPTION_WL_CHOSEN == 0){
            fprintf(outputFile, "%8.2lf\n",Rg);
            printf("\n\n###########################################\n");
            printf("# Radius of gyration    %8.2lf\n",Rg);
            printf("###########################################\n\n");
            fclose(outputFile);
            exit(-1);
        }
        if(OPTION_g_CHOSEN == 1 && OPTION_WL_CHOSEN == 1 && TOTAL_pr == 1){
            fprintf(outputFile, "%8.2lf\n",Rg);
            printf("\n\n###########################################\n");
            printf("# Radius of gyration    %8.2lf\n",Rg);
            printf("###########################################\n\n");
        }
        
        if (TOTAL_pr == 1) {
            fclose(outputFile);
        } else {
            if (i_start == 0) {printf("# Calculating p(r), progression, [%%]: ");}

            int progression = 0;
            for(i=0;i<TotalNumberOfAtoms-i_start;i++){
                Gofr[0] += Ba[i+i_start] * Ba[i+i_start];
                Hofr[0] += Ba[i+i_start] * Bs[i+i_start];
                Jofr[0] += Bs[i+i_start] * Ba[i+i_start];
                Kofr[0] += Bs[i+i_start] * Bs[i+i_start];
                for(j=i+1;j<TotalNumberOfAtoms-i_start;j++){
                    D = Distance(Atoms[i+i_start],Atoms[j+i_start]);
                    index = floor(D/Delta_r);
                    Pofr[index] +=       dB[i+i_start] * dB[j+i_start];
                    Gofr[index] += 2.0 * Ba[i+i_start] * Ba[j+i_start];
                    Hofr[index] += 2.0 * Ba[i+i_start] * Bs[j+i_start];
                    Jofr[index] += 2.0 * Bs[i+i_start] * Ba[j+i_start];
                    Kofr[index] += 2.0 * Bs[i+i_start] * Bs[j+i_start];
                    if (D > Dmax){Dmax = D;}
                }
                if (i_start == 0) {
                    if(i > (TotalNumberOfAtoms-i_start)/10.0 && progression == 0){printf("10 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/5.00 && progression == 1){printf("20 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/3.33 && progression == 2){printf("30 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/2.50 && progression == 3){printf("40 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/2.00 && progression == 4){printf("50 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/1.67 && progression == 5){printf("60 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/1.43 && progression == 6){printf("70 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/1.25 && progression == 7){printf("80 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)/1.11 && progression == 8){printf("90 > "); progression++;}
                    else if(i > (TotalNumberOfAtoms-i_start)-1000 && progression ==9){printf("100!\n"); progression++;}
                    fflush(stdout);
                }
            }
            fprintf(outputFile, "# Dmax: %8.2lf\n", Dmax);
            //fprintf(outputFile, "# Mean volume: %8.2lf A**3. Total Volume = %8.2lf A**3\n", MeanVolume,SumOfVolume);
            fprintf(outputFile, "# Mean volume: %8.2lf\n", MeanVolume);
            fprintf(outputFile, "# Pair distance distribution function:\n");
            fprintf(outputFile, "#     r [AA]  P(r) [cm^2]   G(r) [cm^2]   H(r) [cm^2]   J(r) [cm^2]   K(r) [cm^2]\n");
            fprintf(outputFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for(i=0;i<PointsInPofr;i++)
            {
                r=(i+0.5)*Delta_r;
                fprintf(outputFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", r, Pofr[i], Gofr[i], Hofr[i], Jofr[i], Kofr[i]);
            }
        
            fclose(outputFile);
        } // end the [if TOTAL_pr !=] statement
    } // end the [if p(r) file does not exist] statement
    
    // return values
    *PointsInPofrPointer = PointsInPofr;
    *DmaxPointer = Dmax;
    *SumAtomWeightPointer = SumAtomWeight;
    *RgPointer = Rg;
    *MeanVolumePointer = MeanVolume;
    return NumberOfAtomsRead;
}
