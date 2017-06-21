char * AddWaterLayerToPDB(char *inputPDB, double HalfBilayerThickness)
{
    int SIZE = 100;
  
    int Nres, Nato;
    Nres = CheckNumberOfResidues(inputPDB);
    Nato = CheckNumberOfAtoms(inputPDB);

    /*********Allocate memory***********/
    AGG  agg;
    agg.sizeSx=SIZE;
    agg.Sn=Array(SIZE);
    agg.Sx=Array(SIZE);
    Component(&agg.C[0],0,Nres);
    Component(&agg.C[1],1,Nres);
    /***********************************/
    
    ReadAminoPDB(inputPDB,&agg.C[0]);

    printf("\n            Total number of residues: %d\n",Nres);

    PlaceWater(&agg,SIZE,agg.Sx);
    
    char *inputPDB1 = ExtractString(inputPDB, '.'); //extracts 1HEJ from 1HEJ.pdb

    char *outputfilename = strcat(inputPDB1, "_w.pdb"); // adds w.pdb to 1HEJ
    WritePDB_ProteinAndWater(inputPDB, Nato, HalfBilayerThickness, outputfilename, &agg.C[WAT]);

    char *inputPDB2 = ExtractString(inputPDB, '.'); //extracts 1HEJ from 1HEJ.pdb again
    char *waterfilename = strcat(inputPDB2, "_w_only.pdb");
    WritePDB_Water(waterfilename, HalfBilayerThickness, &agg.C[WAT]);
    
    return outputfilename;
}
