int ReadPDB(char *filename, double *Da, double *Ba, double *Bs, double *dB, struct Atom *Atoms, double WaterLayerContrast, double SolventD2O, double Perdeuteration, double PrcSucrose, double Delta_r, double HalfBilayerThickness, int OPTION_g_CHOSEN, int OPTION_WL_CHOSEN, int i_start, int * PointsInPofrPointer, double * DmaxPointer, double * SumAtomWeightPointer, double * RgPointer, double *MeanVolumePointer, int TOTAL_pr)
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
<<<<<<< HEAD
    enum {HYDROGEN, DEUTERIUM, HD, CARBON, NITROGEN, OXYGEN, FLUORINE, SODIUM, MAGNESIUM, PHOSPHORUS, SULFUR, CHLORINE, CALCIUM, MANGANESE, IRON, COPPER, ZINK, GOLD, WATER, UNKNOWN};
=======
    enum {HYDROGEN, DEUTERIUM, CARBON, NITROGEN, OXYGEN, FLUORINE, SODIUM, MAGNESIUM, PHOSPHORUS, SULFUR, CHLORINE, CALCIUM, MANGANESE, IRON, COPPER, ZINK, GOLD, WATER, UNKNOWN};
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
    
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
<<<<<<< HEAD
        [HYDROGEN]   = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    3.741e-13,                                       5.15,   1.0,   'A'} ,
        [DEUTERIUM]  = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13,                                       5.15,   2.0,   'A'} ,
        [HD]         = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13*SolventD2O-3.741e-13*(1.0-SolventD2O), 5.15,   1.5,   'A'} ,
=======
        [HYDROGEN]  =  { 0.0, 0.0, 0.0,    1 * 2.82e-13,    3.741e-13,                                       5.15,   1.0,   'A'} ,
        [DEUTERIUM]  = { 0.0, 0.0, 0.0,    1 * 2.82e-13,    6.671e-13*SolventD2O-3.741e-13*(1.0-SolventD2O), 5.15,   2.0,   'A'} ,
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
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
    };
    // Note: [WATER] is without hydrogen/deuterium... will be added later on (implicit way)
    // Note: weight of [WATER] is set to 0, since the WL should not be included in the calculation of the atomic weight
    // Volumes from:
    //    (H,C,D,N,O):                              Fraser, J Appl Cryst(1978), 11, p693.
    //                                               Experimentally determined values
    //    (P, S, Mn, Ca, F, Mg, Cl, Cu, Fe, Zn):    Batsanov, Inorg Mat(2001), vol37, no9, p871.
    //                                               Van der Waals radii, mean of values from Table 1 and 2 (P from Table 5)
<<<<<<< HEAD
    // HD is a libile Hydrogen, that can be exchanged with Deuterium
    //
    
    // SOLVENT SCATT LENGTH (with sucrose)
    double ScatteringLengthSolvent = 2.0 * element[HD].NeutronScatteringLength + element[OXYGEN].NeutronScatteringLength;
    if(SolventD2O<0.0){
        ScatteringLengthSolvent = 2.0 * element[HD].XRayScatteringLength + element[OXYGEN].XRayScatteringLength;
=======
    
    // SOLVENT SCATT LENGTH (with sucrose)
    double ScatteringLengthSolvent = 2.0 * element[DEUTERIUM].NeutronScatteringLength + element[OXYGEN].NeutronScatteringLength;
    if(SolventD2O<0.0){
        ScatteringLengthSolvent = 2.0 * element[DEUTERIUM].XRayScatteringLength + element[OXYGEN].XRayScatteringLength;
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
        ScatteringLengthSolvent = ScatteringLengthSolvent + SucPerSolv * ScatLenSuc;
    }
    double rhoSolvent = ScatteringLengthSolvent/VolumeSolvent;
    
    int R,j,u=0; // R = residue number, j = index, u = unknown atom
    char buffer[256];
    
    //int GuessNumberOfAtoms = CheckNumberOfAtomsInPDBFile(filename);
    
    struct Atom COM   = {0,0,0};
    double SqDist, SqDistMax=0, GuessOnDmax;
    double SumOfVolume = 0.0;
    
    char AminoName, AtomName, LongAtomName,AlternativeAtomPosition;
    
    double n_H; // number of non-exchangable (strongly bound) Hydrogen
    double n_D; // number of exchangeable Hydrogen attachd to main atom
<<<<<<< HEAD
    double n_DD; // number of non-exchangeable (strongly bound) Deuterium (after perdeuteration)
    
    // VOLUME CORRECTION FACTORS
    int DNA,RNA,LIP; // is it a lip, a DNA or a RNA atom?
    double LipCorrFactor = 1.06; // from Matlab: CalcLipVol.m
    double DNACorrFactor = 1.0; // not determined
    double RNACorrFactor = 1.0; // not determined
=======
    double n_DD = 0.0; // number of non-exchangeable (strongly bound) Deuterium (after perdeuteration)
    
    int DNA, RNA;
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548

    int NumberOfH = 0;
    int NumberOfAlternativeAtoms = 0;
    int NumberOfOtherLines = 0;
    int i = i_start;
    while(fgets(buffer,sizeof(buffer),PointerToFile)!=NULL){
<<<<<<< HEAD
        
        // read atom name:
        if (sscanf(buffer, "ATOM%*73c%s",&Atoms[i].Name) == 1 || sscanf(buffer, "HETATM%*71c%s",&Atoms[i].Name) == 1) {}
        // read alternative atom pos:
        if (sscanf(buffer, "ATOM%*12c%c",&AlternativeAtomPosition) == 1 || sscanf(buffer, "HETATM%*10c%c",&AlternativeAtomPosition) == 1) {}
        // ignore H and D:
        if (strcmp(&Atoms[i].Name,"H") == 0 || strcmp(&Atoms[i].Name,"D") == 0) {NumberOfH++;}
        // ignore alternative position atoms (only use position A):
        else if (AlternativeAtomPosition == 'B') {NumberOfAlternativeAtoms++;}
        // read all other atoms:
=======
        R = 0;
        
        // read atom name.
        if (sscanf(buffer, "ATOM%*73c%s",&Atoms[i].Name) == 1 || sscanf(buffer, "HETATM%*71c%s",&Atoms[i].Name) == 1) {}
        // read alternative atom pos
        if (sscanf(buffer, "ATOM%*12c%c",&AlternativeAtomPosition) == 1 || sscanf(buffer, "HETATM%*10c%c",&AlternativeAtomPosition) == 1) {}
        
        // if atom is H or D then ignore and if it is a repeated atom with alternative position then ignore (only use position A, else read atom.
        if (strcmp(&Atoms[i].Name,"H") == 0 || strcmp(&Atoms[i].Name,"D") == 0) {NumberOfH++;}
        else if (AlternativeAtomPosition == 'B') {NumberOfAlternativeAtoms++;}
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
        else if (sscanf(buffer,"ATOM%*18c%d%*4c%lf%lf%lf%*22c%s",&R,&Atoms[i].x,&Atoms[i].y,&Atoms[i].z,&AtomName)==5 ||
            sscanf(buffer,"HETATM%*16c%d%*4c%lf%lf%lf%*22c%s",&R,&Atoms[i].x,&Atoms[i].y,&Atoms[i].z,&AtomName)==5 ) {
            n_D = 0;
            n_H = 0;
            DNA = 0;
            RNA = 0;
<<<<<<< HEAD
            LIP = 0;
=======
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
            
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
<<<<<<< HEAD
                CopyPhysicalParameters(&Atoms[i], &element[UNKNOWN]); u++;
            }
            
            // Add H/D. According to http://www.rcsb.org/ligand
            if (sscanf(buffer,"ATOM%*13c%s", &AminoName) == 1 || sscanf(buffer,"HETATM%*11c%s", &AminoName) == 1)
            {
                // ADD H/D to AMINO ACIDS
=======
                CopyPhysicalParameters(&Atoms[i], &element[UNKNOWN]); u++;}
            
            // Add H/D to amino acids. Ordered alphabetically. Then for DNA and RNA
            if (sscanf(buffer,"ATOM%*13c%s", &AminoName) == 1)
            {
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
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
<<<<<<< HEAD
                // ADD H/D to DNA nucleotides
=======
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
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
<<<<<<< HEAD
                // ADD H/D to RNA nucleotides
=======
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
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
<<<<<<< HEAD
                // ADD H/D to RNA to other heteroatomic molecueles
                else if (strcmp(&AminoName,"HOH")==0) {
=======
            }
            // H/D added to HETATM according to http://www.rcsb.org/ligand
            else if (sscanf(buffer,"HETATM%*11c%s", &AminoName) == 1)
            {
                if (strcmp(&AminoName,"HOH")==0) {
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
                    n_D = 2; // HOH is a water atom in the crystal
                    Atoms[i].Volume = 30.0 - 2*5.15; // implicit solvent, therefore the volume of two H/D is subtracted
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
<<<<<<< HEAD
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
            // Update volume and scattering length with implicit H and D (adjust if DNA or RNA)
            Atoms[i].Volume += (n_D + n_H) * 5.15;
            
            // Volume corrections factors for DNA, RNA and LIP
            if (DNA == 1) {Atoms[i].Volume *= DNACorrFactor;}
            else if (RNA == 1) {Atoms[i].Volume *= RNACorrFactor;}
            else if (LIP == 1) {Atoms[i].Volume *= LipCorrFactor;}
            
            SumOfVolume += Atoms[i].Volume;
            
            if(SolventD2O < 0.0){
                Atoms[i].XRayScatteringLength += (n_D + n_H) * element[HD].XRayScatteringLength;
                Ba[i] = Atoms[i].XRayScatteringLength;
            }
            else {
                n_DD = Perdeuteration * n_H;
                n_H = n_H - n_DD;
                Atoms[i].NeutronScatteringLength += n_D * element[HD].NeutronScatteringLength - n_H * element[HYDROGEN].NeutronScatteringLength + n_DD * element[DEUTERIUM].NeutronScatteringLength;
=======
            }
            // Update volume and scattering length with implicit H and D (adjust if DNA or RNA)
            Atoms[i].Volume += (n_D + n_H) * 5.15;
            if (DNA == 1) {
                Atoms[i].Volume *= 1.0; // DNA volume correctioin
            }
            if (RNA == 1) {
                Atoms[i].Volume *= 1.0; /// RNA volume correctioin
            }
            
            SumOfVolume += Atoms[i].Volume;
            
            
            if(SolventD2O < 0.0){
                Atoms[i].XRayScatteringLength += (n_D + n_H) * element[DEUTERIUM].XRayScatteringLength;
                Ba[i] = Atoms[i].XRayScatteringLength;
            }
            else {
                n_D = n_D; // exchangeable hydrogens are the same
                n_DD = Perdeuteration * n_H;
                n_H = n_H - n_DD;
                Atoms[i].NeutronScatteringLength += n_D * element[DEUTERIUM].NeutronScatteringLength - n_H * element[HYDROGEN].NeutronScatteringLength + n_DD * 6.671e-13 ;
>>>>>>> 509c72e74235e58d805f89c2db77102f65256548
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
    if(i_start == 0){
        if(SolventD2O >= 0.0 && SolventD2O <= 1.0){
            sprintf(buf,"_N%1.0fP%1.0f_pr.dat",SolventD2O*100,Perdeuteration*100);
        } else {
            sprintf(buf,"_X%1.0f_pr.dat",PrcSucrose);
        }
    } else {sprintf(buf,"_pr.dat");}
    char *outputfilename = strcat(inputPDB, buf); // adds string to 1HEJ
    
    if (i_start == 0 && fopen(outputfilename,"r") && TOTAL_pr != 1) {
        printf("");
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
