void CopyPhysicalParameters(struct Atom *Copy, struct Atom *Original){
    Copy->XRayScatteringLength    = Original->XRayScatteringLength;
    Copy->NeutronScatteringLength = Original->NeutronScatteringLength;
    Copy->Volume                  = Original->Volume;
    Copy->Weight                  = Original->Weight;
}

double max(double a, double b) {
    if (a>b)
        return a;
    else
        return b;
}

double SquareDistance(struct Atom pos1 , struct Atom pos2){
    return pow(pos1.x-pos2.x,2.0) + pow(pos1.y-pos2.y,2.0) + pow(pos1.z-pos2.z,2.0);
}

double Distance(struct Atom pos1 , struct Atom pos2){
    return sqrt(SquareDistance(pos1,pos2));
}

int CheckNumberOfAtomsInPDBFile(char *filename)
{
    FILE *PointerToFile;
    char buffer[256];
    int DummyID;
    int NumberOfAtoms = 0;

    if( (PointerToFile = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename
    }
    
    if( (PointerToFile = fopen(filename, "r")) == 0){
        printf("\n\nCannot open PDB file:%s \n",filename);
        exit(-1);
    }

    while(fgets(buffer, sizeof(buffer), PointerToFile) != NULL) {
        if(sscanf(buffer, "ATOM%*18c%d*4c", &DummyID) == 1 ||
           sscanf(buffer, "HETATM%*16c%d%*4c", &DummyID) == 1){
            ++NumberOfAtoms;
        }
    }

    fclose(PointerToFile);
    return NumberOfAtoms;
}


char *GetCStringBeforeLastDelimiter(char *str, char delim)
{
    size_t len;
    char *new_str;
    char *delim_pos = strrchr(str, delim);
    
    // new string is the length from the start of the old string to the delimiter
    len = delim_pos - str;
    new_str = malloc(len + 100);
    memcpy(new_str, str, len);
    new_str[len] = '\0'; // NULL terminates the new string
    
    return new_str;
}

static const char SwitchChar = '-';
static const char Unknown = '?';

int OptIndex = 1;       /* first option should be argv[1] */
char *OptArg = NULL;    /* global option argument pointer */

int options(int argc, char *argv[], const char *legal)
{
    static char *posn = "";  /* position in argv[OptIndex] */
    char *legal_index = NULL;
    int letter = 0;
    
    if(!*posn){
        /* no more args, no SwitchChar or no option letter ? */
        if((OptIndex >= argc) ||
           (*(posn = argv[OptIndex]) != SwitchChar) ||
           !*++posn)
            return -1;
        /* find double SwitchChar ? */
        if(*posn == SwitchChar){
            OptIndex++;
            return -1;
        }
    }
    letter = *posn++;
    if(!(legal_index = strchr(legal, letter))){
        if(!*posn)
            OptIndex++;
        return Unknown;
    }
    if(*++legal_index != ':'){
        /* no option argument */
        OptArg = NULL;
        if(!*posn)
            OptIndex++;
    } else {
        if(*posn)
        /* no space between opt and opt arg */
            OptArg = posn;
        else
            if(argc <= ++OptIndex){
                posn = "";
                return Unknown;
            } else
                OptArg = argv[OptIndex];
        posn = "";
        OptIndex++;
    }
    return letter;
}

int CheckNumberOfResidues(char filename[])
{
    FILE *fil;
    int ii = 0;
    char buffer[180];
    double x,y,z;
    
    if( (fil = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename
    }
    
    if( (fil=fopen(filename,"r"))==0){
        printf("Function CheckNumberOfResidues() cannot open file:%s \n",filename);
        exit(0);
    }
    // check for C-alpha (protein), C4' (RNA/DNA), or C5 (glycosylation)
    while(fgets(buffer,sizeof(buffer),fil)!=NULL){
        if(sscanf(buffer,"ATOM%*9cCA%*2c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC5'%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC5%*2c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC19%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 ||sscanf(buffer,"ATOM%*9cC12%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC15%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 || sscanf(buffer,"ATOM%*9cC20%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 ||sscanf(buffer,"ATOM%*9cC24%*1c%*s%*10c%lf%lf%lf",&x,&y,&z) == 3 ){
            ii++;
            //printf("Number of residues found: %d\n",ii);
        }
    }
    fclose(fil);
    return ii;
}

int CheckNumberOfAtoms(char filename[])
{
    FILE *fil;
    char buffer[180];
    int dum,ii=0;
    
    if( (fil = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename
    }
    if( (fil=fopen(filename,"r"))==0){
        printf("Function CheckNumberOfAtoms() cannot open file:%s \n",filename);
        exit(1);
    }
    while(fgets(buffer,sizeof(buffer),fil)!=NULL){
        if(sscanf(buffer,"ATOM %d",&dum)==1){
            ii++;
        }
        if(sscanf(buffer,"HETATM %d",&dum)==1){
            ii++;
        }
    }
    fclose(fil);
    return ii;
}

double CheckHalfBilayerThickness(char filename[])
{
    FILE *fil;
    char buffer[180];
    double T;
    char dum1[20], dum2[20], dum3[20], dum4[20], dum5[20];
    if( (fil=fopen(filename,"r"))==0)
    {
        printf("CheckHalfBilayerThickness cannot open file:%s \n",filename);
        exit(0);
    }
    // reads thickness from first line in file
    fgets(buffer,sizeof(buffer),fil);
    sscanf(buffer,"%s %s %s %s %s %lf", dum1, dum2, dum3, dum4, dum5, &T);
    
    fclose(fil);
    
    //check:
    if (strcmp(dum4, "bilayer")) {
        printf("\nERROR: Dear user,\nI don't think your PDB is from the OPM database (because word no. 4 in the first line is NOT \"bilayer\", but \"%s\") \nPlease give it another try with the manual option -m or with a PDB from the OPM database.\n\n", dum4);
        exit(-1);
    }
    return T;
}

void CopyFile(const char *copyFrom, const char *copyTo)
{
    /*
     This function copies content of copyFrom to copyTo (except the line END\n)
     */
    char content[80];
    char *new_str;
    
    //Step 1: Open text files and check that they open
    FILE *fp1, *fp2;
    fp1 = fopen(copyFrom,"r");
    fp2 = fopen(copyTo,"w");
    if(fp1 == NULL || fp2 == NULL){
        printf("\nError reading file, in function CopyFile()\n");
        exit(0);}
    
    //Step 2: Get text from original file. Stop before "END"
    while(fgets(content, sizeof(content), fp1) !=NULL)
    {
        new_str = malloc(80);
        memcpy(new_str, &content,3);
        new_str[3] = '\0'; // NULL terminates the new string
        if(strcmp(new_str,"END") != 0) {fprintf(fp2, "%s", content);}
    }
    
    //Step 3: Close both files and end program
    fclose(fp1);
    fclose(fp2);
}

double char2double(char *C)
{
    strtod(C,NULL);
    float D = atof(C);
    return D;
}
