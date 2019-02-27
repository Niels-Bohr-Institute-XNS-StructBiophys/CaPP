//subfunction
void WritePDB_Water(char filename[], double HalfBilayerThickness, double *Wx, double *Wy, double *Wz, int NW)
{
    FILE * fil;
    int i;
    fil=fopen(filename,"w");
    int j = 1;
    
    for(i=0;i<NW;i++)
    {
        if (fabs(Wz[i]) > HalfBilayerThickness){
            fprintf(fil,"ATOM  %5d  Q   WAT W%4d    %8.3lf%8.3lf%8.3lf  1.00 10.00           Q\n",j,j,Wx[i], Wy[i],Wz[i]);
            j++;}
    }
    fprintf(fil,"END\n");
    fclose(fil);
}

//mainfunction 1
char * GenerateWaterLayerPDB(char *inputPDB, double HalfBilayerThickness)
{
    // find number of residues in the PDB (CA atoms for protein + C4' atoms for DNA/RNA)
    int Nres = CheckNumberOfResidues(inputPDB);

    // alpha-carbon positions (C4' pos for DNA/RNA)
    double *CAx = (double *) calloc( Nres, sizeof(double));
    double *CAy = (double *) calloc( Nres, sizeof(double));
    double *CAz = (double *) calloc( Nres, sizeof(double));
    
    // water bead positions
    double *Wx = (double *) calloc( Nres, sizeof(double));
    double *Wy = (double *) calloc( Nres, sizeof(double));
    double *Wz = (double *) calloc( Nres, sizeof(double));
    
    // get CA positions
    FILE *fil;
    char *filename;
    char buffer[180];
    int i = 0;
    if( (fil = fopen(inputPDB,"r"))==0){
        filename = strcat(inputPDB, ".pdb"); // adds .pdb extension to inputPDB
    }
    if( (fil=fopen(inputPDB,"r"))==0){
        printf("Function GenerateWaterLayerPDB(): cannot open file:%s \n",filename);
        exit(1);
    }
    double x,y,z;
    while(fgets(buffer,sizeof(buffer),fil)!=NULL){
        if(sscanf(buffer,"ATOM%*9cCA%*2c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC5'%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC5%*2c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC19%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 ||sscanf(buffer,"ATOM%*9cC12%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC15%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC20%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 ||sscanf(buffer,"ATOM%*9cC24%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 ){
            CAx[i] = x;
            CAy[i] = y;
            CAz[i] = z;
            i++;
        }
    }
    fclose(fil);
    
    // find surface residues, place water beads (see Kynde's thesis)
    int NW = 0, N; // NW = number of water residues, N = Number of vacent residues (within 10 aa)
    double dx,dy,dz,D, dist, sumdx,sumdy,sumdz;
    for(int j=0; j<Nres; j++){
        N = 0;
        sumdx = 0; sumdy = 0; sumdz = 0;
        for(int i=0; i<Nres; i++){
            dx = CAx[j] - CAx[i];
            dy = CAy[j] - CAy[i];
            dz = CAz[j] - CAz[i];
            dist = sqrt(dx*dx + dy*dy + dz*dz);
            if(dist<10.0){
                sumdx += dx;
                sumdy += dy;
                sumdz += dz;
                N++;
            }
        }
        D = sqrt(sumdx*sumdx + sumdy*sumdy + sumdz*sumdz);
        //printf("D = %f \n",D);
        if(D>=0.1*N*N || D==0.0){
            Wx[NW] = CAx[j] + 5.0*sumdx/D;
            Wy[NW] = CAy[j] + 5.0*sumdy/D;
            Wz[NW] = CAz[j] + 5.0*sumdz/D;
            NW++;
        }
    }
    
    // make water pdb file
    char *inputPDB1 = GetCStringBeforeLastDelimiter(inputPDB, '.'); //extracts 1HEJ from 1HEJ.pdb
    char *waterfilename = strcat(inputPDB1, "_w_only.pdb");
    WritePDB_Water(waterfilename, HalfBilayerThickness, Wx, Wy, Wz, NW);
    return waterfilename;
}

//mainfunction 2
char * GenerateTotalPDB(char *filename)
{
    char *inputPDB1 = GetCStringBeforeLastDelimiter(filename, '.'); //extracts 1HEJ from 1HEJ.pdb
    char *waterfilename = strcat(inputPDB1, "_w_only.pdb");
    char *inputPDB2 = GetCStringBeforeLastDelimiter(filename, '.'); //extracts 1HEJ from 1HEJ.pdb
    char *totalfilename = strcat(inputPDB2, "_w.pdb");
    
    FILE *PointerToPDBFile;
    PointerToPDBFile = fopen(filename,"r");
    FILE *PointerToTotalFile;
    PointerToTotalFile = fopen(totalfilename,"w");
    
    char buffer[256];
    while(fgets(buffer,sizeof(buffer),PointerToPDBFile)!=NULL){
        if (strcmp(buffer,"END\n") != 0) {
            fprintf(PointerToTotalFile,"%s", buffer);
        }
    }
    fclose(PointerToPDBFile);
    
    FILE *PointerToWaterFile;
    PointerToWaterFile = fopen(waterfilename,"r");
    while(fgets(buffer,sizeof(buffer),PointerToWaterFile)!=NULL){
        fprintf(PointerToTotalFile,"%s",buffer);
    }
    fclose(PointerToWaterFile);
    fclose(PointerToTotalFile);
    return totalfilename;
}
