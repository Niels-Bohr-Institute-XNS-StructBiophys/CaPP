void CalcCrossTerms(char *filename, int NumberOfAtomsPDB, int NumberOfAtomsWL, double *Ba, double *Bs, double *dB, struct Atom *Atoms, double WaterLayerContrast, double Delta_r, int PointsInPofr)
{
    // calculate p(r) for cross terms
    int i, j, index;
    double r, D, Dmax = 0.0;
    double *Pofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Gofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Hofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Jofr = (double *) calloc( PointsInPofr, sizeof(double));
    double *Kofr = (double *) calloc( PointsInPofr, sizeof(double));
    int TotalNumberOfAtoms = NumberOfAtomsPDB + NumberOfAtomsWL;
    
    for(i=NumberOfAtomsPDB;i<TotalNumberOfAtoms;i++){
        for(j=0;j<NumberOfAtomsPDB;j++){
            D = Distance(Atoms[i],Atoms[j]);
            index = floor(D/Delta_r);
            Pofr[index] +=       dB[i] * dB[j];
            Gofr[index] += 2.0 * Ba[i] * Ba[j]; // water layer - protein
            Hofr[index] += 2.0 * Ba[i] * Bs[j]; // water layer - solvent
            Jofr[index] += 2.0 * Bs[i] * Ba[j]; // solvent     - protein
            Kofr[index] += 2.0 * Bs[i] * Bs[j]; // solvent     - solvent
            if (D > Dmax){Dmax = D;}
        }
    }
    
    // write to outputfile
    FILE *outputFile;
    char *inputPDB = GetCStringBeforeLastDelimiter(filename, '.'); //extracts 1HEJ from 1HEJ.pdb
    char *outputfilename = strcat(inputPDB, "_cross_pr.dat"); // adds pr.dat to 1HEJ
    outputFile = fopen(outputfilename,"w");
    fprintf(outputFile, "#\n#\n#\n#\n#\n#\n");
    fprintf(outputFile, "# Dmax %8.2lf\n", Dmax);
    fprintf(outputFile, "#\n");
    fprintf(outputFile, "# Pair distance distribution function:\n");
    fprintf(outputFile, "#     r [AA]  P(r) [cm^2]   G(r) [cm^2]   H(r) [cm^2]   J(r) [cm^2]   K(r) [cm^2]\n");
    fprintf(outputFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for(i=0;i<PointsInPofr;i++)
    {
        r=(i+0.5)*Delta_r;
        fprintf(outputFile, "%12.4g  %-12.6g  %-12.6g  %-12.6g  %-12.6g  %-12.6g\n", r, Pofr[i], Gofr[i], Hofr[i], Jofr[i], Kofr[i]);
    }
    fclose(outputFile);
}
