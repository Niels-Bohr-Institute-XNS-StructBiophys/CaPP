///                           ///
///  Copyright 2015,          ///
///  University of Copenhagen ///
///                           ///
///  Andreas N Larsen         ///
///  based on                 ///
///  Søren Kynde's program    ///
///                           ///
///  anlarsen@nbi.ku.dk       ///
///                           ///

//Compile:
//gcc pr_implicit_HD.c -o pr


#define TRUE 1
#define FALSE 0
#define NC 2 /* Number of components */
#define PRO 0 /* Each components ID number*/
#define WAT 1
#define END (1.0/0.0)

struct Atom {
    // Coordinates
    double x;
    double y;
    double z;

    // Physical properties
    double XRayScatteringLength;
    double NeutronScatteringLength;
    double Volume;
    double Weight;

    // Book-keeping
    char Name;
};

void CopyPhysicalParameters(struct Atom *Copy, struct Atom *Original){
    Copy->XRayScatteringLength    = Original->XRayScatteringLength;
    Copy->NeutronScatteringLength = Original->NeutronScatteringLength;
    Copy->Volume                  = Original->Volume;
    Copy->Weight                  = Original->Weight;
};

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
        if(sscanf(buffer, "ATOM%*18c%d%*4c", &DummyID) == 1 ||
           sscanf(buffer, "HETATM%*16c%d%*4c", &DummyID) == 1){
            ++NumberOfAtoms;
        }
    }

    fclose(PointerToFile);
    return NumberOfAtoms;
}

char *ExtractString(char *str, char delim)
{
    size_t len;
    char *new_str;
    char *delim_pos = strchr(str, delim);
    
    // new string is the length from the start of the old string to the delimiter
    len = delim_pos - str;
    new_str = malloc(len + 1);
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

typedef double array;
array * Array(int N)
{
    array *px;
    px=calloc((size_t)(N+1),sizeof(double));
    px[N]= 1.0/0.0;
    return px;
}

struct bead
{
    double x,y,z;/*Space coordinates of center of gyration*/
    double xa,ya,za; /*Space coordinates of C-alpha atom.*/
    double delta;
    double r;   /*radius of gyration*/
    int N;  /*Number of neighbours*/
    int Neigh[10]; /*Array with names of neighbours*/
    array *F; /*Pointer to formfactor*/
    char *amin;
    int shadow;
    int oswich;
    int nagg;
    int visible;
    int w[20];
    double W;
    double Wa;
    double f;
};
typedef struct bead BEAD;

struct comp
{
    BEAD *B;
    /*double **rr;
     double **ra;*/
    double stat[12];
    double statexp[12];
    double staterr[12];
    double neigh1[20];
    double neigh1exp[20];
    double neigh1err[20];
    double neigh2[20];
    double neigh2exp[20];
    double neigh2err[20];
    double neigh3[20];
    double neigh3exp[20];
    double neigh3err[20];
    double lambda;
    double connectstrength;
    int N;
    int ID;
    int connect;
};
typedef struct comp COMP;

struct agg
{
    COMP C[NC];
    array  *Sx;
    array  *Sn;
    int sizeSx;
    int sizeSn;
};
typedef struct agg AGG;

int length(array *A)
{
    int s=0;
    do{
    }while(A[++s]!=END);
    return s;
}

void Component(COMP *comp,int id,int size)
{
    int i,j;
    comp->N=size;
    comp->ID=id;
    comp->B=calloc((size_t) size,sizeof(BEAD));
}

double adist(BEAD *a,BEAD *b)
{
    return sqrt(pow(a->xa-b->xa,2)+pow(a->ya-b->ya,2)+pow(a->za-b->za,2));
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

int CheckNumberOfResidues(char filename[])
{
    FILE *fil;
    int R,Rprev,ii=0;
    char buffer[180];
    char url[200];
    double dumx,dumy,dumz;
    
    if( (fil = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename
    }
    
    if( (fil=fopen(filename,"r"))==0){
        printf("CheckNumberOfResidues cannot open file:%s \n",filename);
        exit(0);
    }
    while(fgets(buffer,sizeof(buffer),fil)!=NULL){
        R=0;
        if(sscanf(buffer,"ATOM%*18c%d%*4c%lf%lf%lf",&R,&dumx,&dumy,&dumz)==4){
            if(R!=Rprev && R!=0)
                ii++;
            Rprev=R;;
        }
    }
    fclose(fil);
    return ii;
}

double We(char * amin, char aname, char anum, char abranch){
    double WH=-0.720,WC=0.509,WN=6.168,WO=4.951,WS=9.367,WP=13.03;
    double W=1e-5;
    switch(aname){
        case 'O':
            W=WO;
            switch(anum){
                case ' ':
                    return W;
                    break;
                case 'T':
                    return W;
                    break;
            };
            break;
        case 'N':
            W=WN;
            switch(anum){
                case ' ':
                    W+=WH;
                    break;
            };
            break;
        case 'C':
            W=WC;
            switch(anum){
                case 'A':
                    W+=WH;
                    break;
                case ' ':
                    return W;
            };
            break;
        case 'H':
            W=0;
            return W;
            break;
        default:
            ;
    }
    
    if(strcmp(amin,"ASN")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
            case 'E':
                switch(abranch){
                    case '2':
                        W+=2*WH;
                }
                break;
        }
        return W;
    }
    if(strcmp(amin,"GLN")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
            case 'E':
                switch(abranch){
                    case '2':
                        W+=2*WH;
                        break;
                }
                break;
        }
        return W;
    }
    if(strcmp(amin,"GLY")==0){
        switch(anum){
            case 'A':
                W+=WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"ALA")==0){
        switch(anum){
            case 'B':
                W+=3*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"VAL")==0){
        switch(anum){
            case 'B':
                W+=WH;
                break;
            case 'G': //Two branches
                W+=3*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"LEU")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=WH;
                break;
            case 'D': //Two branches
                W+=3*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"MET")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
            case 'D':
                W=WS;
                break;
            case 'E':
                W+=3*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"ILE")==0){
        switch(anum){
            case 'B':
                W+=WH;
                break;
            case 'G':
                switch(abranch){
                    case '1':
                        W+=2*WH;
                        break;
                    case '2':
                        W+=3*WH;
                        break;
                }
                break;
            case 'D':
                W+=3*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"PHE")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'D': //Two branches
                W+=WH;
                break;
            case 'E': //Two branches
                W+=WH;
                break;
            case 'Z':
                W+=WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"TYR")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'D': //Two branches
                W+=WH;
                break;
            case 'E': //Two branches
                W+=WH;
                break;
            case 'H':
                W+=WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"TRP")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'D':
                switch(abranch){
                    case '1':
                        W+=WH;
                }
                break;
            case 'E':
                switch(abranch){
                    case '1':
                        W+=WH;
                        break;
                    case '3':
                        W+=WH;
                        break;
                }
                break;
            case 'Z':   //2 branches
                W+=WH;
                break;
            case 'H':
                W+=WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"LYS")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
            case 'D':
                W+=2*WH;
                break;
            case 'E':
                W+=2*WH;
                break;
            case 'Z':
                W+=3*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"ARG")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
            case 'D':
                W+=2*WH;
                break;
            case 'E':
                W+=WH;
                break;
            case 'H': //Two branches
                W+=2*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"HIS")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'D':
                W+=WH; //Two branches
                break;
            case 'E':
                switch(abranch){
                    case '2':
                        W+=WH;
                        break;
                }
                break;
        }
        return W;
    }
    if(strcmp(amin,"SER")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"THR")==0){
        switch(anum){
            case 'B':
                W+=WH;
                break;
            case 'G':
                switch(abranch){
                    case '1':
                        W+=WH;
                        break;
                    case '2':
                        W+=3*WH;
                        break;
                }
                break;
        }
        return W;
    }
    if(strcmp(amin,"CYS")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W=WS;
                break;
        }
        return W;
    }
    if(strcmp(amin,"PRO")==0){
        switch(anum){
            case ' ':
                if(aname==*"N")
                    W-=WH;
                break;
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
            case 'D':
                W+=2*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"ASP")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
        }
        return W;
    }
    if(strcmp(amin,"GLU")==0){
        switch(anum){
            case 'B':
                W+=2*WH;
                break;
            case 'G':
                W+=2*WH;
                break;
        }
        return W;
    }
    return W;
    
}

void ReadAminoPDB(char filename[],COMP *co)
{
    int nagg=co->N;
    BEAD *agg=co->B;
    FILE *fil;
    char buffer[180];
    int i=0,Rprev=0,R=0;
    double weight=0,delta,W=0,Aweight,A=0,angle;
    double dx,dy,dz;
    double CA[3]={0,0,0};
    double COS[3]={0,0,0};
    double X=0,Y=0,Z=0,x,y,z,XA=0,YA=0,ZA=0;
    double WH=-0.720,WC=0.509,WN=6.168,WO=4.951,WS=9.367,WP=13.03;
    char anum,aname,abranch;
    char atom;
    char* scanline="hej";
    if( (fil=fopen(filename,"r"))==0){
        printf("ReadAminoPDB cannot open file:%s \n",filename);
        exit(1);
    }
    agg[i].amin=(char *) calloc(3,sizeof(char));
    while(scanline!=NULL){
        scanline=fgets(buffer,sizeof(buffer),fil);
        atom=0; R=0;
        if(sscanf(buffer,"\nATOM%*9c%c%*8c%d%*4c%lf%lf%lf",&atom,&R,&x,&y,&z)==5){
            if(R!=Rprev){
                if(Rprev!=0){
                    COS[0]=X/weight  ,COS[1]=Y/weight  ,COS[2]=Z/weight  ;
                    CA[0] =agg[i].xa ,CA[1] =agg[i].ya ,CA[2] =agg[i].za ;
                    
                    if(strcmp("pdb","CA")==0){ agg[i].x=CA[0]; agg[i].y=CA[1]; agg[i].z=CA[2];}
                    else{ agg[i].x=COS[0]; agg[i].y=COS[1]; agg[i].z=COS[2];}
                    
                    agg[i].delta=sqrt(pow(COS[0]-CA[0],2)+pow(COS[1]-CA[1],2)+pow(COS[2]-CA[2],2));
                    X=0;Y=0;Z=0;weight=0;XA=0;YA=0;ZA=0;Aweight=0;
                    i++;
                    if(i<nagg) {agg[i].amin=(char *) calloc(3,sizeof(char));}
                }
                Rprev=R;
            }
        }
        sscanf(buffer,"ATOM%*9cCA%*2c%*s%*10c%lf%lf%lf",&agg[i].xa,&agg[i].ya,&agg[i].za);
        sscanf(buffer,"ATOM%*13c%s",agg[i].amin);
        if(sscanf(buffer,"ATOM%*9c%c%c%c",&aname,&anum,&abranch))
            W=We(agg[i].amin,aname,anum,abranch);
        weight+=W;
        X+=W*x;
        Y+=W*y;
        Z+=W*z;
    }
    fclose(fil);
}

int CheckNumberOfAtoms(char filename[])
{
    FILE *fil;
    char buffer[180];
    double q,S,err;
    int dum,ii=0;
    
    if( (fil = fopen(filename,"r"))==0){
        filename = strcat(filename, ".pdb"); // adds .pdb extension to filename
    }
    
    if( (fil=fopen(filename,"r"))==0){
        printf("CheckNumberOfAtoms cannot open file:%s \n",filename);
        exit(1);
    }
    while(fgets(buffer,sizeof(buffer),fil)!=NULL){
        if(sscanf(buffer,"ATOM %d",&dum)==1){
            ii++;
        }
        if(sscanf(buffer,"HETATM %d",&dum)==1){
            ii++;
        }
        if(sscanf(buffer,"TER %d",&dum)==1){
            ii++;
        }
        
    }
    fclose(fil);
    return ii;
}

void WritePDB_Water(char filename[], double HalfBilayerThickness, COMP *comp)
{
    FILE * fil;
    int i, N=comp->N;
    BEAD *B=comp->B;
    fil=fopen(filename,"w");
    int j=1;
    
    for(i=0;i<N;i++)
    {
        if (fabs(B[i].z) > HalfBilayerThickness){
            fprintf(fil,"ATOM  %5d  Q   WAT W%4d    %8.3lf%8.3lf%8.3lf  1.00 10.00           Q\n",j,j,B[i].x, B[i].y,B[i].z);
            j++;}
    }
    fprintf(fil,"END");
    fclose(fil);
}

int CopyFile(const char *copyFrom, const char *copyTo)
{
    /*
     This function copies content of copyFrom to copyTo (except the line END\n)
     anlarsen@nbi.ku.dk
     Dec 2015
     */
    
    char content[80];
    char stringtestforEND;
    
    //Step 1: Open text files and check that they open
    FILE *fp1, *fp2;
    fp1 = fopen(copyFrom,"r");
    fp2 = fopen(copyTo,"w");
    
    if(fp1 == NULL || fp2 == NULL){
        printf("\nError reading file\n");
        exit(0);}
    
    //Step 2: Get text from original file. Stop before "END"
    while(fgets(content, sizeof(content), fp1) !=NULL)
    {
        sscanf(content,"%s",&stringtestforEND);
        if(strcmp(&stringtestforEND,"END") != 0) {fprintf(fp2, "%s", content);}
    }
    
    //Step 3: Close both files and end program
    fclose(fp1);
    fclose(fp2);
    return 0;
}

void WritePDB_ProteinAndWater(char inputfilename[], int Nato, double HalfBilayerThickness, char outputfilename[],COMP *comp)
{
    CopyFile(inputfilename, outputfilename);
    
    //open files
    FILE *outfil;
    outfil = fopen(outputfilename, "a");
    
    int j = 1;
    int i;
    int N = comp->N;
    fprint("Here 1");
    BEAD * B = comp->B;
    fprint("Here 2");
    for(i=0;i<N;i++)
    {
        fprint("Here 3");
        if (fabs(B[i].z) > HalfBilayerThickness){
            fprintf(outfil,"ATOM  %5d  Q   WAT W%4d    %8.3lf%8.3lf%8.3lf  1.00 10.00           Q\n",j+Nato,j,B[i].x, B[i].y,B[i].z);
            j++;}
    }
    fprintf(outfil,"END");
    fclose(outfil);
}

int PlaceWater(AGG *akk, int SIZE, array *S)
{
    int i,j,n;
    BEAD* mol=akk->C[PRO].B;
    BEAD* wat=akk->C[WAT].B;
    int nwat=akk->C[WAT].N;
    int nmol=akk->C[PRO].N;
    int N,s=0;
    double dx=0,dy=0,dz=0,D;
    double sumW;
    double shellvolume;
    double A,AA;
    if(nwat<2)
        return 1;
    for(n=0;n<nmol;n++){
        N=0;
        dx=0;dy=0;dz=0;
        for(i=0;i<nmol;i++){
            if(adist(&mol[n],&mol[i])<10.0){
                dx+=mol[n].xa-mol[i].xa;
                dy+=mol[n].ya-mol[i].ya;
                dz+=mol[n].za-mol[i].za;
                N++;
            }
        }
        D=sqrt(dx*dx+dy*dy+dz*dz);
        if(D>=0.10*N*N){
            wat[n].xa=mol[n].xa+5.*dx/D;
            wat[n].ya=mol[n].ya+5.*dy/D;
            wat[n].za=mol[n].za+5.*dz/D;
            wat[n].delta=1;
            wat[n].Wa=(double) 1;
            s++;
        }
        else {
            wat[n].delta=0;
            wat[n].Wa=0;
        }
    }
    shellvolume=3*37.2*s;
    printf("\n            Number of surface residues: %d\n", s);
    sumW=0;
    for(i=0;i<nmol;i++){
        sumW+=wat[i].Wa;
    }
    s=0;
    for(n=0;n<nmol;n++){
        wat[n].x=wat[n].xa;
        wat[n].y=wat[n].ya;
        wat[n].z=wat[n].za;
    }
    sumW=0;
    for(i=0;i<nmol;i++){
        if(wat[i].visible){
            sumW+=wat[i].W;
        }
    }
    return 1;
}

double char2double(char *C)
{
    double temp = strtod(C,NULL);
    float D = atof(C);
    return D;
}
