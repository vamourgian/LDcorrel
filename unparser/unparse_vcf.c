#include "header.h"

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

int main(int argc, char** argv)
{
    int snipSize = 0,totalSnips = 0,snipsPerFile=0, counter = 0, Wmin = 0, Wmax = 0, Wset = 0, posWmin = 0, posWmax = 0, posMin = 0, posMax = 0, posWset = 0,inListSet = 0, listSize = 0;
    
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
		
	    if(posWmax < posMin || posWmax > posMax)
	    	posWmax = posMax;
	}
	else if(inListSet == 1) {
		listSize = sortList(inputListName,inList, line, &lineLength);
	}

    if(snipsPerFile < 0 || snipsPerFile > totalSnips)
    	snipsPerFile = totalSnips;
    
	if(Wset == 1)
	{
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, Wmin: %d, Wmax: %d\n\n\n",totalSnips,snipSize,snipsPerFile, Wmin, Wmax);
		readPartIndexW(fpOut, inputPathName, allignmentId, totalSnips, snipSize, &counter, headerLine1, headerLine2, snipsPerFile, Wmin, Wmax, line, &lineLength, word, &wordLength);
	}
	else if(posWset == 1)
	{
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, posWmin: %d, posWmax: %d\n\n\n",totalSnips,snipSize,snipsPerFile, posWmin, posWmax);
		readPartPosW(fpOut, inputPathName, allignmentId, totalSnips, snipSize, &counter, headerLine1, headerLine2, snipsPerFile, posWmin, posWmax, line, &lineLength, word, &wordLength);
	}
	else if(inListSet == 1)
	{
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, using input List\n\n\n",totalSnips,snipSize,snipsPerFile);
		readPartPosL(fpOut, inputPathName, allignmentId, totalSnips, snipSize, &counter, headerLine1, headerLine2, snipsPerFile, inList, listSize, line, &lineLength, word, &wordLength);
	}
	else
	{
		printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d\n\n\n",totalSnips,snipSize,snipsPerFile);
		readPart(fpOut, inputPathName, allignmentId, totalSnips, snipSize, &counter, headerLine1, headerLine2, snipsPerFile, line, &lineLength, word, &wordLength);
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