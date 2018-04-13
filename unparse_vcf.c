#include "header.c"

void printHelp (FILE *fp)
{
    fprintf(fp,"\n\n\n");
    
    fprintf(fp," VCF_unparser_demo\n");
    fprintf(fp,"\t -input inputFolder\n");
    fprintf(fp,"\t -output outputFile\n");
    fprintf(fp,"\t -Wmin snipID\n");
    fprintf(fp,"\t -Wmax snipID\n");
    fprintf(fp,"\n\n");
    fprintf(fp,"\t-input <STRING>\t\tSpecifies the path of the input alignment parsed files\n");
    fprintf(fp,"\t-output <STRING>\t\tSpecifies the name of the output alignment file.\n");
    fprintf(fp,"\t      \t\t\tSupported file formats: VCF.\n\n");
    fprintf(fp,"\t-Wmin ID of the minimum snip to be included, minimum 1 (default)\n");
    fprintf(fp,"\t-Wmax ID of the maximum snip to be included, maximum total Snips (default)\n");
    fprintf(fp,"\n\n");
}

void commandLineParser(int argc, char** argv, 
                       char * inPath, 
                       char * outfile,
                       int * Wmin,
                       int * Wmax)
{
    int i, pathSet = 0, fileSet=0;
    
    for(i=1; i<argc; ++i)
    {
        if(!strcmp(argv[i], "-input")) 
        { 
            if (i!=argc-1)
            {           
                strcpy(inPath,argv[++i]);

                DIR* directory = opendir(inPath);

                if (directory == NULL) {
                	fprintf(stderr, "\n ERROR: Directory %s does not exist.\n\n",inPath);
                	exit(0);
                }
                else {
                    fprintf(stdout, "\n Directory %s opened.\n\n",inPath);
                    closedir(directory);
                }

                int pos = strlen(inPath)-1;
                if(inPath[pos] != '/') {
                    inPath[pos+1] = '/';
                    inPath[pos+2] = '\0';
                }
                pathSet=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-output")) 
        {
            if (i!=argc-1)
            {           
                strcpy(outfile,argv[++i]);
                fileSet=1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-Wmin")) 
        {
            if (i!=argc-1)
            {           
                *Wmin = atoi(argv[++i]);
            }
            continue;
        }

        if(!strcmp(argv[i], "-Wmax")) 
        {
            if (i!=argc-1)
            {           
                *Wmax = atoi(argv[++i]);
            }
            continue;
        }
        
        if(!strcmp(argv[i], "-help")||!strcmp(argv[i], "-h"))
        {             
            printHelp (stdout);
            
            exit(0);
        }
        
        fprintf(stderr, "\n ERROR: %s is not a valid command line parameter\n\n",argv[i]);
        exit(0);
    }
    
    if (pathSet==0)
    {
        fprintf(stderr, "\n ERROR: Please specify a path for the output with -output\n\n");
        exit(0);
    }
    
    if (fileSet==0)
    {
        fprintf(stderr, "\n ERROR: Please specify an alignment folder with -input\n\n");
        exit(0);
    }
}


void readHeaderFile(gzFile fpOut, char* inputPathName, char ** headerLine1, char ** headerLine2, char* allignmentId, int* snipsPerFile, int* snipSize, int* totalSnips,char **line, int* lineLength)
{
    int eol=0, eof=0, status, headerFound = 0,allignmentIdInt;
    char headerFile[INFILENAMESIZE+30];
    DIR *dp;
    struct dirent *entry;
    struct stat statbuf;
    if((dp = opendir(inputPathName)) == NULL) {
        fprintf(stderr,"cannot open directory: %s\n", inputPathName);
        assert(0);
    }
    status = chdir(inputPathName);
    while((entry = readdir(dp)) != NULL) {
        lstat(entry->d_name,&statbuf);
        if(!S_ISDIR(statbuf.st_mode)) {
            if(strstr(entry->d_name, "header") != NULL) {
                if(headerFound == 0) {
                    sscanf(entry->d_name,"%[^'_']",allignmentId);
                    sprintf(headerFile,"%s%s",inputPathName,entry->d_name);
                    printf("Found Header file %s\n",headerFile);
                    headerFound = 1;
                }
                else {
                    printf("NOTICE - Multiple Header files found!! : %s\nNOTICE - Processing first header file\n",entry->d_name);
                }
            }
        } 
    }
    status = chdir("..");
    closedir(dp);
    if(headerFound == 0) {
    	printf("ERROR: Could not find Header file\n");
    	exit(0);
    }
    printf("allignment %s found\n",allignmentId);
    gzFile fpHeader = gzopen(headerFile,"r");
    status =  getNextLine(fpHeader, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1) {
        *headerLine1 = (char *) malloc(sizeof(char)*(strlen(*line)));
        strcpy(*headerLine1,*line);
        gzprintf(fpOut,"%s\n",*line);
        while(eof == 0) {
            status =  getNextLine(fpHeader, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1) {
                writeLine(fpOut,line);
                if(strlen(*line) >=6) {
                    if((*line)[0] == '#' && (*line)[1] == 'C' && (*line)[2] == 'H' && (*line)[3] == 'R' && (*line)[4] == 'O' && (*line)[5] == 'M') {
                        *headerLine2 = (char *) malloc(sizeof(char)*(strlen(*line)));
                        strcpy(*headerLine2,*line);
                        status =  getNextLine(fpHeader, line, &eol, &eof, lineLength);
                        if(status == 1 && eol==1) {
                            sscanf(*line,"%d %d %d %d",&allignmentIdInt,snipsPerFile,snipSize, totalSnips);
                            if(allignmentIdInt != atoi(allignmentId)) {
                                printf("Different internal and external ID: %d != %s\n",allignmentIdInt, allignmentId);
                                assert(0);
                            }
                        }
                        else {
                            printf("wrong header format\n");
                            assert(0);
                        }
                    }

                }
            }
        }
    }
    else if(eof==1) {
        printf("Empty header file\n");
        assert(0);
    }
    return;
}

void readPart(gzFile fpOut, char* inputPathName, char *allignmentId, int * snipSize, int *counter, char ** headerLine1, char ** headerLine2, int *snipsPerFile, int *Wmin, int *Wmax, char **line, int* lineLength, char **word, int* wordLength)
{
    int tmpId=0,tmpMin=0,tmpMax=0, fileFound = 0;
    int eol=0, eof=0, status;
    char inputFileName[INFILENAMESIZE+30];
    DIR *dp;
    gzFile fpIn;
    struct dirent *entry;
    struct stat statbuf;
    if((dp = opendir(inputPathName)) == NULL) {
        fprintf(stderr,"cannot open directory: %s\n", inputPathName);
        assert(0);
    }
    status = chdir(inputPathName);
    while((entry = readdir(dp)) != NULL) {
        lstat(entry->d_name,&statbuf);
        if(!S_ISDIR(statbuf.st_mode)) {
            if(strstr(entry->d_name, allignmentId) != NULL && strstr(entry->d_name, "header") == NULL) {
                sscanf(entry->d_name,"%d_%d_%d",&tmpId,&tmpMin,&tmpMax);
                if((tmpMin <= (*Wmin + (*counter)))  && ((*Wmin + (*counter)) <= tmpMax) && (((tmpMax - tmpMin + 1) == *snipsPerFile) || tmpMax == *Wmax)) {
                    sprintf(inputFileName,"%s%s",inputPathName,entry->d_name);
                    fileFound = 1;
                    break;
                }
            }
        } 
    }
    status = chdir("..");
    closedir(dp);
    if(fileFound == 1) {
        printf("\033[A\33[2K");
        printf("Opening file %s, %d %%\n",inputFileName,((*counter)*100)/(*Wmax-*Wmin));
        fpIn = gzopen(inputFileName,"r");
        if(fpIn == NULL) {
            printf("failed to open file %s", inputFileName);
            assert(0);
        }

        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1) {
           if(strcmp(*line,*headerLine1)) {
                printf("Header 1 different\n");
                assert(0);
           }
           else {
                status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
                if(status == 1 && eol==1) {
                   if(strcmp(*line,*headerLine2)) {
                        printf("Header 2 different\n");
                        assert(0);
                   }
                   else {
                        while(eof == 0 && (*counter)< (*Wmax - *Wmin + 1)) {
                            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
                            if(status == 1 && eol==1 && tmpMin >= *Wmin && tmpMin <= *Wmax) {
                                writeLine(fpOut,line);
                                (*counter)++;
                            }
                            tmpMin++;
                        }
                   }
               }
           }
        }
        else {
            printf("file %s empty",inputFileName);
            assert(0);
        }
        gzclose(fpIn);
    }
    else {
        if(*counter > (*Wmax - *snipsPerFile))
            printf("File missing: probably %s_%d_%d.vcf.gz or %s_%d_%d.vcf.gz\n",allignmentId,(*Wmin + (*counter)),(*Wmin + (*counter) + *snipsPerFile-1),allignmentId,(*Wmin + (*counter)),*Wmax);
        else
            printf("File missing: probably %s_%d_%d.vcf.gz\n",allignmentId,(*Wmin + (*counter)),(*Wmin + (*counter) + *snipsPerFile-1));
        assert(0);
    }
    //printf("--closing file %s, %d %% \n",inputFileName,((*counter)*100)/(*Wmax-*Wmin));

}

int main(int argc, char** argv)
{
    int snipSize = 0,totalSnips = 0,snipsPerFile=0, counter = 0, Wmin = 0, Wmax = 0;
    
    double time1,time2;
    char inputPathName[INFILENAMESIZE];
    char outputFileName[INFILENAMESIZE];
    char allignmentId[INFILENAMESIZE];

    gzFile fpOut=NULL;

    char ** headerLine1, ** headerLine2;
    int  lineLength = STRINGLENGTH,wordLength = STRINGLENGTH;
    headerLine1 = (char **) malloc(sizeof(char*));
    headerLine2 = (char **) malloc(sizeof(char*));

    commandLineParser(argc, argv, inputPathName, outputFileName, &Wmin, &Wmax);
    
    char **line = (char**)malloc(sizeof(char*));
    *line = (char*)malloc(sizeof(char)*STRINGLENGTH);
    char **word = (char**)malloc(sizeof(char*));
    *word = (char*)malloc(sizeof(char)*STRINGLENGTH);
    fpOut = gzopen(outputFileName,"w");
    time1 = gettime();
    readHeaderFile(fpOut, inputPathName, headerLine1, headerLine2, allignmentId, &snipsPerFile, &snipSize, &totalSnips, line, &lineLength);
    if(Wmin == 0)
    	Wmin = 1;
    if(Wmax == 0)
    	Wmax = totalSnips;

    if(Wmin < 1 || Wmin > totalSnips) {
    	printf("ERROR: Wmin must be between 1 and totalSnips, setting Wmin to 1 (default)\n");
    	Wmin = 1;
    }

    if(Wmax < 1 || Wmax > totalSnips) {
    	printf("ERROR: Wmax must be between 1 and totalSnips, setting Wmax to total Snips (default)\n");
    	Wmax = totalSnips;
    }

    if(Wmin > Wmax) {
    	printf("ERROR: Wmin must be <= Wmax, setting Wmin = Wmax\n");
    	Wmin = Wmax;
    }
    if(snipsPerFile < 0 || snipsPerFile > totalSnips)
    	snipsPerFile = totalSnips;

    printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, Wmin: %d, Wmax: %d\n\n\n",totalSnips,snipSize,snipsPerFile, Wmin, Wmax);
    while(counter < (Wmax - Wmin + 1)) {
        readPart(fpOut, inputPathName, allignmentId, &snipSize, &counter, headerLine1, headerLine2, &snipsPerFile, &Wmin, &Wmax, line, &lineLength, word, &wordLength);
    }    
    gzclose(fpOut);
    time2 = gettime();
    printf("%d snips processed\n",counter);
    printf("Total Time: %fs \n", time2-time1);
    return 1;
}