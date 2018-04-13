#include "header.c"

void printHelp (FILE *fp)
{
    fprintf(fp,"\n\n\n");
    
    fprintf(fp," VCF_unparser_demo\n");
    fprintf(fp,"\t -input inputFolder\n");
    fprintf(fp,"\t -output outputFile\n");
    fprintf(fp,"\t -Wmin snip index\n");
    fprintf(fp,"\t -Wmax snip index\n");
    fprintf(fp,"\t -posWmin snip pos\n");
    fprintf(fp,"\t -posWmax snip pos\n");
    fprintf(fp,"\t -inputList inputFile\n");
    fprintf(fp,"\n\n");
    fprintf(fp,"\t-input <STRING>\t\tSpecifies the path of the input alignment parsed files\n");
    fprintf(fp,"\t-output <STRING>\t\tSpecifies the name of the output alignment file.\n");
    fprintf(fp,"\t      \t\t\tSupported file formats: VCF.\n\n");
    fprintf(fp,"\t-Wmin \t\t\tindex of the minimum snip to be included, minimum 1 (default)\n");
    fprintf(fp,"\t-Wmax \t\t\tindex of the maximum snip to be included, maximum total Snips (default)\n");
    fprintf(fp,"\t-posWmin \t\t\tpos of the minimum snip to be included, must be valid\n");
    fprintf(fp,"\t-posWmax \t\t\tpos of the maximum snip to be included, must be valid\n");
    fprintf(fp,"\t-inputList\t\t\tinput text file with the pos to keep\n");
    fprintf(fp,"\n\n");
}

void commandLineParser(int argc, char** argv, 
                       char * inPath, 
                       char * outfile,
                       int * Wmin,
                       int * Wmax,
                       int *Wset,
                       int * posWmin,
                       int * posWmax,
                       int *posWset,
                       char * inList, 
                       int * inListSet)
{
    int i, pathSet = 0, fileSet=0;
    gzFile fp;
    
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
        if(!strcmp(argv[i], "-inputList")) 
        { 
            if (i!=argc-1)
            {			
                strcpy(inList,argv[++i]);
                
                fp=gzopen(inList,"r");
                
                if (fp==NULL)
                {
                    fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",inList);
                    exit(0);
                }
                else
                {
                    gzclose(fp);
                    *inListSet=1;
                }
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
                *Wset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-Wmax")) 
        {
            if (i!=argc-1)
            {           
                *Wmax = atoi(argv[++i]);
                *Wset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-posWmin")) 
        {
            if (i!=argc-1)
            {           
                *posWmin = atoi(argv[++i]);
                *posWset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-posWmax")) 
        {
            if (i!=argc-1)
            {           
                *posWmax = atoi(argv[++i]);
                *posWset = 1;
            }
            continue;
        }
        
        if(!strcmp(argv[i], "-help")||!strcmp(argv[i], "-h")||!strcmp(argv[i], "--help"))
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

    if(((*Wset)==1 && (*posWset)==1) || ((*Wset)==1 && (*inListSet)==1) || ((*inListSet)==1 && (*posWset)==1))
    {
        fprintf(stderr, "\n ERROR: Please specify only an index window, or a pos window, or an input List\n\n");
        exit(0);
    }
}


void readHeaderFile(gzFile fpOut, char* inputPathName, char ** headerLine1, char ** headerLine2, char* allignmentId, int* snipsPerFile, int* snipSize, int* totalSnips, int *posMin,int *posMax,char **line, int* lineLength)
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
                            sscanf(*line,"%d %d %d %d %d %d",&allignmentIdInt,snipsPerFile,snipSize, totalSnips,posMin, posMax );
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

void readPart(gzFile fpOut, char* inputPathName, char *allignmentId, int * eop, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int Wmin, int Wmax, int Wset, int posWmin, int *tmpPosMinIndex, int posWmax, int posWset, char **line, int* lineLength, char **word, int* wordLength, int ** inList,int listSize, int inListSet,int *inListPos)
{
    int tmpId=0,tmpMin=0,tmpMax=0,tmpPosMin=0,tmpPosMax=0, fileFound = 0, pos;
    int eol=0, eof = 0, status;
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
                sscanf(entry->d_name,"%d_%d_%d_%d_%d",&tmpId,&tmpMin,&tmpMax,&tmpPosMin,&tmpPosMax);
                if(
                  	(
                  	 (tmpMin <= (Wmin + (*counter))  && (Wmin + (*counter)) <= tmpMax && Wset == 1) ||
                  	 (((tmpPosMin <= posWmin   && posWmin <= tmpPosMax && *tmpPosMinIndex == 0) || (tmpMin <= *tmpPosMinIndex  && *tmpPosMinIndex <= tmpMax)) && posWset == 1) || 
                  	 (tmpPosMin <= (*inList)[(*inListPos)] && (*inList)[(*inListPos)] <= tmpPosMax && inListSet == 1) ||
                  	 (tmpMin <= (1 + (*counter))  && (1 + (*counter)) <= tmpMax && Wset == 0 && posWset == 0 && inListSet == 0)
                  	) && 
                  	(((tmpMax - tmpMin + 1) == snipsPerFile) || tmpMax == totalSnips)
                  ) {
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
        if(Wset == 1)
        	printf("Opening file %s, %d %%\n",inputFileName,((*counter)*100)/(Wmax-Wmin+1));
        else if(posWset == 1)
        	printf("Opening file %s\n",inputFileName);
        else if(inListSet == 1)
        	printf("Opening file %s, %d %%\n",inputFileName,((*counter)*100)/listSize);
        else
        	printf("Opening file %s, %d %%\n",inputFileName,((*counter)*100)/totalSnips);
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
                   		do{
                            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
                            if(status == 1 && eol==1){
                            	int index = 0;
					            status = getWordFromString(*line, word, &eol, wordLength, &index);//chrom
					            assert(status);
					            status = getWordFromString(*line, word, &eol, wordLength, &index);//pos
					            assert(status);
						        pos = atoi(*word);
						        if((tmpMin >= Wmin && tmpMin <= Wmax && Wset ==1) || (pos >= posWmin && pos <= posWmax && posWset ==1) || (pos == (*inList)[(*inListPos)] && inListSet == 1) || (Wset ==0 && posWset == 0 && inListSet == 0)){
		                            writeLine(fpOut,line);
		                            (*counter)++;
        							(*inList)[(*inListPos)] = 0;
        							if(inListSet==1)
										(*inListPos)++;
		                        }

		                        if(inListSet==1 && pos > (*inList)[(*inListPos)])
									(*inListPos)++;

	                            if(((tmpPosMin > (*inList)[(*inListPos)] || (*inList)[(*inListPos)] > tmpPosMax) && inListSet == 1))
	                            	eof = 1;
		                        tmpMin++;
                            }
                            else {
                            	pos = 0;
                            }
                        }while((eof) == 0 && ( ((*counter)< (Wmax - Wmin + 1) && Wset == 1) || (pos < posWmax && posWset == 1) || ((*inListPos) < listSize &&  inListSet == 1)  || ((*counter)< totalSnips && Wset == 0 && posWset == 0 &&  inListSet == 0) ));
                        *tmpPosMinIndex = tmpMax+1;
                        if( ((*counter)>= (Wmax - Wmin + 1) && Wset == 1) || (pos >= posWmax && posWset == 1) || ((*inListPos) >= listSize &&  inListSet == 1) || ((*counter)>= totalSnips && Wset == 0 && posWset == 0  && inListSet == 0) )
                        	*eop = 1;
                   }
               }
           }
        }
        else {
            printf("file %s empty\n",inputFileName);
            assert(0);
        }
        gzclose(fpIn);
    }
    else {
    	if(inListSet==1 && (*inListPos) < listSize)
			(*inListPos)++;
		else if(inListSet==1 &&  (*inListPos) >= listSize)
			*eop = 1;
		else {
	        printf("File missing or wrong window\n");
	        assert(0);
	    }
    }
}

int main(int argc, char** argv)
{
    int snipSize = 0,totalSnips = 0,snipsPerFile=0, counter = 0, Wmin = 0, Wmax = 0, Wset = 0, posWmin = 0, posWmax = 0, posMin = 0, posMax = 0, posWset = 0, eop = 0,tmpPosMinIndex=0,inListSet = 0, listSize = 0,inListPos=0;
    
    double time1,time2;
    char inputPathName[INFILENAMESIZE];
    char inputListName[INFILENAMESIZE];
    char outputFileName[INFILENAMESIZE];
    char allignmentId[INFILENAMESIZE];

    gzFile fpOut=NULL;

    char ** headerLine1, ** headerLine2;
    int  lineLength = STRINGLENGTH,wordLength = STRINGLENGTH;
    headerLine1 = (char **) malloc(sizeof(char*));
    headerLine2 = (char **) malloc(sizeof(char*));

    commandLineParser(argc, argv, inputPathName, outputFileName, &Wmin, &Wmax, &Wset, &posWmin, &posWmax, &posWset,inputListName,&inListSet);
    char **line = (char**)malloc(sizeof(char*));
    *line = (char*)malloc(sizeof(char)*STRINGLENGTH);
    char **word = (char**)malloc(sizeof(char*));
    *word = (char*)malloc(sizeof(char)*STRINGLENGTH);
    int **inList = (int**)malloc(sizeof(int*));
    *inList = (int*)malloc(sizeof(int));
    fpOut = gzopen(outputFileName,"w");
    time1 = gettime();
    readHeaderFile(fpOut, inputPathName, headerLine1, headerLine2, allignmentId, &snipsPerFile, &snipSize, &totalSnips, &posMin, &posMax, line, &lineLength);
    
    if(Wset == 1) {
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
	}
	else if(posWset == 1) {
	    if(posWmin < posMin || posWmin > posMax)
	    	posWmin = posMin;
		
	    if(posWmax < posMax || posWmax > posMax)
	    	posWmax = posMax;
	}
	else if(inListSet == 1) {
		listSize = sortList(inputListName,inList, line, &lineLength);
	}

    if(snipsPerFile < 0 || snipsPerFile > totalSnips)
    	snipsPerFile = totalSnips;
    
	if(Wset == 1)
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, Wmin: %d, Wmax: %d\n\n\n",totalSnips,snipSize,snipsPerFile, Wmin, Wmax);
	else if(posWset == 1)
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, posWmin: %d, posWmax: %d\n\n\n",totalSnips,snipSize,snipsPerFile, posWmin, posWmax);
	else if(inListSet == 1)
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, using input List\n\n\n",totalSnips,snipSize,snipsPerFile);
	else
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d\n\n\n",totalSnips,snipSize,snipsPerFile);

    while(eop == 0) {
        readPart(fpOut, inputPathName, allignmentId, &eop, totalSnips, snipSize, &counter, headerLine1, headerLine2, snipsPerFile, Wmin, Wmax, Wset, posWmin, &tmpPosMinIndex, posWmax, posWset, line, &lineLength, word, &wordLength, inList, listSize, inListSet,&inListPos);
    }    
	if(inListSet == 1) {
		printf("\n");
		for(inListPos=0;inListPos<listSize;inListPos++)
			if((*inList)[inListPos] != 0)
        		printf("\tWARNING: list pos %d not found\n",(*inList)[inListPos]);
		printf("\n\n");
	}
    gzclose(fpOut);
	free(*headerLine1);
	free(*headerLine2);
	free(headerLine1);
	free(headerLine2);
    free(*line);
    free(line);
    free(*word);
    free(word);
    free(*inList);
    free(inList);
    time2 = gettime();
    printf("%d snips processed\n",counter);
    printf("Total Time: %fs \n", time2-time1);
    return 1;
}