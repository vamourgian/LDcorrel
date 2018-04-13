#include "header.c"

int getFileFormat (gzFile fp)
{
    
    char tmp;
    tmp=gzgetc(fp);
    char macsFLAG[] = "COMMAND";
    int flaglength = 7,j;
    char vcfFLAG[] = "##fileformat=VCF";
    int flaglengthVCF = 16;
    int format = OTHER_FORMAT;
    while(tmp!=EOF)
    {
        if(tmp=='/')
        {
            tmp = gzgetc(fp);
            
            if(tmp=='/')
            {
                format = MS_FORMAT;
                break;
            }
            else
                tmp = gzgetc(fp);				
        }
        else
        {
            if(tmp=='>')
            {
                format = FASTA_FORMAT;
                break;
            }
            else
            {
                int counter = 0;
                
                while(counter < flaglength)
                {
                    if(tmp != macsFLAG[counter])
                        break;
                    
                    tmp = gzgetc(fp);
                    
                    ++counter;
                }
                
                j = counter;
                
                
                if(j == flaglength)
                {
                    format = MACS_FORMAT;
                    break;
                }
                
                else
                {
                    
                    gzseek(fp, -j - 1, SEEK_CUR);
                    
                    tmp = gzgetc( fp);
                    
                    int counter = 0;
                    
                    while(counter < flaglengthVCF)
                    {
                        if(tmp != vcfFLAG[counter])
                            break;
                        
                        tmp = gzgetc(fp);
                        
                        ++counter;
                    }
                    
                    j = counter;
                    
                    if(j == flaglengthVCF)
                    {
                        format = VCF_FORMAT;
                        break;
                    }
                    
                    else
                        tmp = gzgetc(fp);
                }
            }
        }		
    }
    
    return format;
}

void printHelp (FILE *fp)
{
    fprintf(fp,"\n\n\n");
    
    fprintf(fp," VCF_parser\n");
    fprintf(fp,"\t -input inputFile\n");
    fprintf(fp,"\t -output outputFolder\n");
    fprintf(fp,"\t -size parts_size\n");
    fprintf(fp,"\t -Wmin snip index\n");
    fprintf(fp,"\t -Wmax snip index\n");
    fprintf(fp,"\t -posWmin snip pos\n");
    fprintf(fp,"\t -posWmax snip pos\n");
    fprintf(fp,"\t -toSingleOutput\n");
    fprintf(fp,"\n\n");
    fprintf(fp,"\t-input <STRING>\t\tSpecifies the name of the input alignment file.\n");
    fprintf(fp,"\t-output <STRING>\tSpecifies the path of the output alignment files.\n");
    fprintf(fp,"\t-size <STRING>\t\tSpecifies the size of the memory footprint of the output alignment files in MB, if toSingleOutput is set, size is not needed.\n");
    fprintf(fp,"\t      \t\t\tSupported file formats: VCF.\n\n");
    fprintf(fp,"\t-Wmin \t\t\tindex of the minimum snip to be included, minimum 1 (default)\n");
    fprintf(fp,"\t-Wmax \t\t\tindex of the maximum snip to be included, maximum total Snips (default)\n");
    fprintf(fp,"\t-posWmin \t\t\tpos of the minimum snip to be included, must be valid\n");
    fprintf(fp,"\t-posWmax \t\t\tpos of the maximum snip to be included, must be valid\n");
    fprintf(fp,"\t-toSingleOutput\t\tUsed to generate a new VCF that is part of the input file, -Wmin and -Wmax mandatory with this command\n");
    fprintf(fp,"\n\n");
}

void commandLineParser(int argc, char** argv, 
                       char * infile, 
                       char * outPath, 
                       int * fileFormat,
                       int * size,
                       int * Wmin,
                       int * Wmax,
                       int *Wset,
                       int * posWmin,
                       int * posWmax,
                       int *posWset,
                       int * toSingleOutput)
{
    int i, pathSet = 0, fileSet=0, sizeSet=0;
    gzFile fp;
    
    for(i=1; i<argc; ++i)
    {
        if(!strcmp(argv[i], "-input")) 
        { 
            if (i!=argc-1)
            {			
                strcpy(infile,argv[++i]);
                
                fp=gzopen(infile,"r");
                
                if (fp==NULL)
                {
                    fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",infile);
                    exit(0);
                }
                else
                {
                    *fileFormat = getFileFormat (fp);
                    if (*fileFormat!=VCF_FORMAT)
                    {
                        fprintf(stderr, "\n ERROR: File %s is not VCF format.\n\n",infile);
                        exit(0);
                    }
                    
                    gzclose(fp);
                    
                    fileSet=1;
                }
            }
            continue;
        }
        if(!strcmp(argv[i], "-output")) 
        {
            if (i!=argc-1)
            {			
    			char outputPathName[INFILENAMESIZE];
                strcpy(outputPathName,argv[++i]);
                struct stat st = {0};
                
                if (stat(outputPathName, &st) == -1) {
                    int status = mkdir(outputPathName, 0700);
                    if (status == -1) {
                        fprintf(stderr, "\n ERROR: Directory %s could not be created.\n\n",outputPathName);
                        exit(0);
                    }
                    else
                        fprintf(stdout, "\n Directory %s created.\n\n",outputPathName);
                }
                else {
                    fprintf(stdout, "\n Directory %s exists, renaming.\n\n",outputPathName);
                    int i=1;
    				char outputPathNameNew[INFILENAMESIZE];
                    do{
                    	sprintf(outputPathNameNew,"%s_%d",outputPathName,i);
                    	i++;
                    }while(stat(outputPathNameNew, &st) != -1);

                	strcpy(outputPathName,outputPathNameNew);
                    int status = mkdir(outputPathName, 0700);
                    if (status == -1) {
                        fprintf(stderr, "\n ERROR: Directory %s could not be created.\n\n",outputPathName);
                        exit(0);
                    }
                    else
                        fprintf(stdout, "\n Directory %s created.\n\n",outputPathName);
                }
                int pos = strlen(outputPathName)-1;
                if(outputPathName[pos] != '/') {
                    outputPathName[pos+1] = '/';
                    outputPathName[pos+2] = '\0';
                }
                strcpy(outPath,outputPathName);
                pathSet=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-size"))
        {
            if (i!=argc-1)
            {
                *size = atoi(argv[++i]);
                fprintf(stdout, "\n Part size: %d MB.\n\n",*size);
                
                if(*size!=0)
                    sizeSet=1;
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

        if(!strcmp(argv[i], "-toSingleOutput")) 
        {
        	*toSingleOutput = 1;
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
        fprintf(stderr, "\n ERROR: Please specify an alignment with -input\n\n");
        exit(0);
    }
    
    if (( sizeSet==0) && *toSingleOutput == 0)
    {
        fprintf(stderr, "\n ERROR: Please specify a size in MB for the output with -size\n\n");
        exit(0);
    }

    if(*toSingleOutput == 1 && (*Wset)==0 && (*posWset)==0)
    {
        fprintf(stderr, "\n ERROR: Please specify an index window with -Wmin and/or Wmax, or a pos window with -posWmin and/or -posWmax\n\n");
        exit(0);
    }

    if((*Wset)==1 && (*posWset)==1)
    {
        fprintf(stderr, "\n ERROR: Please specify an index window or a pos window, not both\n\n");
        exit(0);
    }
}

gzFile createHeaderFile(gzFile fpIn, char* outputPathName, char ** headerLine1, char ** headerLine2, char **line, int* lineLength, char **word, int* wordLength)
{
    int eol=0, eof=0, status;
    char headerFile[INFILENAMESIZE+30];
    sprintf(headerFile,"%sheader.vcf.gz",outputPathName);
    
    gzFile fpOut=NULL;
    printf("Creating Header: %s\n\n",headerFile);
    fpOut = gzopen(headerFile,"w");
    
    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1) {
        *headerLine1 = (char *) malloc(sizeof(char)*(strlen(*line)));
        strcpy(*headerLine1,*line);
        gzprintf(fpOut,"%s\n",*line);
        while(1)
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1) {
                writeLine(fpOut,line);
                if(strlen(*line) >=6) {
                    if((*line)[0] == '#' && (*line)[1] == 'C' && (*line)[2] == 'H' && (*line)[3] == 'R' && (*line)[4] == 'O' && (*line)[5] == 'M') {
                        *headerLine2 = (char *) malloc(sizeof(char)*(strlen(*line)));
                        strcpy(*headerLine2,*line);
                        return fpOut;
                    }
                }
            }
            else if(eof==1)
                assert(0);
        }
    }
    else if(eof==1)
        assert(0);
    return NULL;
}

void createPart(gzFile fpIn, gzFile *fpOut, char* outputPathName, int *first, char *allignmentId, int * snipSize,int *eof, int *counter, char ** headerLine1, char ** headerLine2, int *size, int *maxcount, int *startingValue, int *endingValue, char **line, int* lineLength, char **word, int* wordLength, int lines, int Wmin, int Wmax, int Wset, int posWmin, int posWmax, int posWset, int *lineCounter, int *printDelay, int * firstPos)
{    
    char outputFile[INFILENAMESIZE+30];
    int status=-1, eol=0,index = 0,pos;
    if(*first == 1) {
        status =  getNextLine(fpIn, line, &eol, eof, lineLength);
        if(status == 1 && eol==1) {
            status = getWordFromString(*line, allignmentId, &eol, wordLength, &index);//chrom
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//pos
            assert(status);
	        pos = atoi(*word);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//id
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//ref
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//alt
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//qual
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//filter
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//info
            assert(status);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//format
            assert(status);
            while(index < strlen(*line)) {
                if ((*line)[index] == '1' || (*line)[index] == '0')
                    (*snipSize)++;
                index++;
            }
        	(*lineCounter)++;
            *maxcount = ((*size)*1024*1024)/((*snipSize)/8);

    		if(*maxcount < 0 || *maxcount > lines)
    			*maxcount = lines;

        	if(Wset == 1)
    			printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, Wmin: %d, Wmax: %d\n\n\n",lines,*snipSize,*maxcount, Wmin, Wmax);
    		else if(posWset == 1)
    			printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, posWmin: %d, posWmax: %d\n\n\n",lines,*snipSize,*maxcount, posWmin, posWmax);
    		else
    			printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d\n\n\n",lines,*snipSize,*maxcount);

    		if((Wset == 1 && *lineCounter >= Wmin) || (posWset == 1 && pos >= posWmin) || (Wset == 0 && posWset == 0)) {
            	sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,allignmentId,*lineCounter,(*lineCounter)+(*maxcount)-1,pos);
            	*startingValue = *lineCounter;
            	*endingValue = (*lineCounter)+(*maxcount)-1;
            	*firstPos = pos;

	        	if(Wset == 1)
                	printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/(Wmax-Wmin+1));
            	else if (posWset == 1) 
                	printf("Creating File: %s\n",outputFile);
	            else
                	printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/lines);
    		
            	*fpOut = gzopen(outputFile,"w");
	            if(*fpOut == NULL) {
	                printf("failed to open file %s", outputFile);
	                assert(0);
	            }
            	writeLine(*fpOut,headerLine1);
            	writeLine(*fpOut,headerLine2);
            	writeLine(*fpOut,line);
        		(*counter)++;
        		*first = 0;  
        	}
        	else
        		*first = 2;  

            if((Wset == 1 && *lineCounter >= Wmax) || (posWset == 1 && pos >= posWmax) || (Wset == 0 && posWset == 0 && *lineCounter >= lines)) {
            	*eof =1;
        		gzclose(*fpOut);
        		if((*counter) == 0) {
		        	if(Wset == 1)
		        		printf("ERROR: Less snips than expected found\n");
		        	else if(posWset == 1)
		        		printf("ERROR: pos window is wrong\n");
		        	exit(0);
        		}

	    		char outputFileNew[INFILENAMESIZE+30];
				printf("\033[A\33[2K");
            	printf("Creating File: %s, 100%%\n",outputFile);  
	            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,allignmentId,*startingValue,(*lineCounter),*firstPos,pos);
				printf("Renaming %s to %s\n",outputFile,outputFileNew);
				int ret = rename(outputFile, outputFileNew);
				if(ret != 0) {
					printf("Error: unable to rename the file\n");
	                assert(0);
				}
            }
        }
        else if(*eof ==1) {
        	printf("ERROR: No snips found\n");
			assert(0);
        }   
    }
    else if(*first == 2) {
        status =  getNextLine(fpIn, line, &eol, eof, lineLength);
        if(status == 1 && eol==1) {
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//chrom
            assert(status);
            char * allID = malloc(strlen(*word)*sizeof(char));
            strcpy(allID,*word);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//pos
            assert(status);
            pos = atoi(*word);
    		(*lineCounter)++;
            if(!strcmp(allID,allignmentId)) {
    			if((Wset == 1 && *lineCounter >= Wmin) || (posWset == 1 && pos >= posWmin) || (Wset == 0 && posWset == 0)) {
            	sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,allignmentId,*lineCounter,(*lineCounter)+(*maxcount)-1,pos);
            		*startingValue = *lineCounter;
            		*endingValue = (*lineCounter)+(*maxcount)-1;
            		*firstPos = pos;

	    			printf("\033[A\33[2K");
					if(Wset == 1) 
	                	printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/(Wmax-Wmin+1));
	            	else if (posWset == 1) 
	                	printf("Creating File: %s\n",outputFile);
		            else
	                	printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/lines);

	                *fpOut = gzopen(outputFile,"w");
	                if(*fpOut == NULL) {
	                    printf("failed to open file %s", outputFile);
	                    assert(0);
	                }
		            writeLine(*fpOut,headerLine1);
		            writeLine(*fpOut,headerLine2);
		            writeLine(*fpOut,line);
        			(*counter)++;
        			*first = 0;
		        }
		    	else
		    		*first = 2;  

                if((Wset == 1 && *lineCounter >= Wmax) || (posWset == 1 && pos >= posWmax) || (Wset == 0 && posWset == 0 && *lineCounter >= lines)) {
            		*eof =1;
            		gzclose(*fpOut);
	        		if((*counter) == 0) {
			        	if(Wset == 1)
			        		printf("ERROR: Less snips than expected found\n");
			        	else if(posWset == 1)
			        		printf("ERROR: pos window is wrong\n");
			        	exit(0);
	        		}

		    		char outputFileNew[INFILENAMESIZE+30];
					printf("\033[A\33[2K");
			        printf("Creating File: %s, 100%%\n",outputFile);   
	            	sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,allignmentId,*startingValue,(*lineCounter),*firstPos,pos);
					printf("Renaming %s to %s\n",outputFile,outputFileNew);
					int ret = rename(outputFile, outputFileNew);
					if(ret != 0) {
						printf("Error: unable to rename the file\n");
		                assert(0);
					}
                }
            }
            else
            	printf("WARNING: snip %d is allignment %s and will be skipped\n",*lineCounter,allID);
            free(allID);
        }
        else if(*eof ==1) {
        	return;
        }
    }
    else {
        status =  getNextLine(fpIn, line, &eol, eof, lineLength);
        if(status == 1 && eol==1) {
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//chrom
            assert(status);
            char * allID = malloc(strlen(*word)*sizeof(char));
            strcpy(allID,*word);
            status = getWordFromString(*line, *word, &eol, wordLength, &index);//pos
            assert(status);
            pos = atoi(*word);
    		(*lineCounter)++;

    		if(Wset == 1 && (((*counter)*100)/(Wmax-Wmin+1))%5 == 0 && (*printDelay) != (((*counter)*100)/(Wmax-Wmin+1))){
    			printf("\033[A\33[2K");
                printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/(Wmax-Wmin+1));   
                (*printDelay) = (((*counter)*100)/(Wmax-Wmin+1));    			
    		}
    		else if(Wset == 0 && posWset == 0 && (((*counter)*100)/lines)%5 == 0 && (*printDelay) != (((*counter)*100)/lines)) {
    			printf("\033[A\33[2K");
                printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/lines);  
                (*printDelay) = (((*counter)*100)/lines); 
    		}

            if(!strcmp(allID,allignmentId)) {
                (*counter)++;
	            writeLine(*fpOut,line);
                if((Wset == 1 && *lineCounter >= Wmax) || (posWset == 1 && pos >= posWmax) || (Wset == 0 && posWset == 0 && *lineCounter >= lines)) {
            		*eof =1;
            		gzclose(*fpOut);

		    		char outputFileNew[INFILENAMESIZE+30];
            		printf("Creating File: %s, 100%%\n",outputFile);   
    				sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,allignmentId,*startingValue,(*lineCounter),*firstPos,pos);
					printf("Renaming %s to %s\n",outputFile,outputFileNew);
					int ret = rename(outputFile, outputFileNew);
					if(ret != 0) {
						printf("Error: unable to rename the file\n");
		                assert(0);
					}
                }
                else if((*lineCounter) == (*endingValue)) {
                    *first = 2;
            		gzclose(*fpOut);

		    		char outputFileNew[INFILENAMESIZE+30];
		    		printf("\033[A\33[2K");                
            		sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,allignmentId,*startingValue,(*lineCounter),*firstPos,pos);
					printf("Renaming %s to %s\n",outputFile,outputFileNew);
					int ret = rename(outputFile, outputFileNew);
					if(ret != 0) {
						printf("Error: unable to rename the file\n");
		                assert(0);
					}                    
                }
            }
            else
            	printf("WARNING: snip %d is allignment %s and will be skipped\n",*lineCounter,allID);
            free(allID);
        }
    }
}

void createSingleFile(gzFile fpIn, gzFile *fpOut, char* outputPathName, int *counter, char **line, int* lineLength, char **word, int* wordLength, int lines, int Wmin, int Wmax, int Wset, int posWmin, int posWmax, int posWset)
{

	int eol=0,eof=0, status,index = 0,first = 1, snipSize = 0, lineCounter = 0, printDelay = 0,pos,printCount=0,i,firstPos=0;
	char allignmentId[INFILENAMESIZE+30];

    while(1)
    {
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1) {
            writeLine(*fpOut,line);
            if(strlen(*line) >=6) {
                if((*line)[0] == '#' && (*line)[1] == 'C' && (*line)[2] == 'H' && (*line)[3] == 'R' && (*line)[4] == 'O' && (*line)[5] == 'M')
                	break;
            }
        }
        else if(eof==1){
        	printf("ERROR: Header malformed\n");
            exit(0);
        }
    }
    printf("Header copied\n");
 	while(1)
    {
    	index = 0;
	    if(first == 1) {
	        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
	        if(status == 1 && eol==1) {
	            status = getWordFromString(*line, allignmentId, &eol, wordLength, &index);//chrom
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//pos
	            assert(status);
	            pos = atoi(*word);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//id
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//ref
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//alt
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//qual
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//filter
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//info
	            assert(status);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//format
	            assert(status);
	            while(index < strlen(*line)) {
	                if ((*line)[index] == '1' || (*line)[index] == '0')
	                    snipSize++;
	                index++;
	            }
	        	lineCounter++;
	        	if(Wset == 1)
	    			printf("Total snips: %d, Snip Size: %d bits, Wmin: %d, Wmax: %d\n\n\n",lines,snipSize, Wmin, Wmax);
	    		else if(posWset == 1)
	    			printf("Total snips: %d, Snip Size: %d bits, posWmin: %d, posWmax: %d\n\n\n",lines,snipSize, posWmin, posWmax);

	    		if((Wset == 1 && lineCounter >= Wmin) || (posWset == 1 && pos >= posWmin)) {
	    			if(posWset == 1 && firstPos == 0)
	    				firstPos = pos;
	            	writeLine(*fpOut,line);
	        		(*counter)++;
	        	}
	        	else if (Wset == 0 && posWset == 0) {
	            	writeLine(*fpOut,line);
	        		(*counter)++;
	        	}

	        	first = 2;
	        	if(Wset == 1)
	            	printf("Processing: %d%%\n",((*counter)*100)/(Wmax-Wmin+1));
            	else if (posWset == 1) {  
	            	printf("Processing");
	            	for(i=0;i<=printCount;i++)
	            		printf(".");
	            	printf("\n");
	            	printCount++;
	            	if(printCount == 5)
	            		printCount = 0;
	            }
	            else
	            	printf("Processing: %d%%\n",((*counter)*100)/lines);

	            if((Wset == 1 && lineCounter >= Wmax) || (posWset == 1 && pos >= posWmax) || lineCounter >= lines) {
	        		gzclose(*fpOut);
	        		if((*counter) == 0) {
			        	if(Wset == 1)
			        		printf("ERROR: Less snips than expected found\n");
			        	else if(posWset == 1)
			        		printf("ERROR: pos window is wrong\n");
			        	exit(0);
	        		}

		    		char outputFileNew[INFILENAMESIZE+30];
		    		char outputFile[INFILENAMESIZE+30];
    				printf("\033[A\33[2K");
                	printf("Processing: 100%%\n");
	        		if(Wset == 1)
		            	sprintf(outputFileNew,"%s%s_%d_%d.vcf.gz",outputPathName,allignmentId,Wmin,Wmax);
            		else if (posWset == 1) 
		            	sprintf(outputFileNew,"%s%s_%d_%d.vcf.gz",outputPathName,allignmentId,firstPos,pos);

		            sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
					printf("Renaming %s to %s\n",outputFile,outputFileNew);
					int ret = rename(outputFile, outputFileNew);
					if(ret != 0) {
						printf("Error: unable to rename the file\n");
		                assert(0);
					}
					break;
	            }
	        }
	        else if(eof ==1) {
	        	printf("ERROR: No snips found\n");
				exit(0);
	        }   
	    }
	    else if(first == 2) {
	        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
	        if(status == 1 && eol==1) {
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//chrom
	            assert(status);
	            char * allID = malloc(strlen(*word)*sizeof(char));
	            strcpy(allID,*word);
	            status = getWordFromString(*line, *word, &eol, wordLength, &index);//pos
	            assert(status);
	            pos = atoi(*word);
	    		lineCounter++;

	    		if(Wset == 1 && (((*counter)*100)/(Wmax-Wmin+1))%5 == 0 && printDelay != (((*counter)*100)/(Wmax-Wmin+1))){
	    			printf("\033[A\33[2K");
	                printf("Processing: %d%%\n",((*counter)*100)/(Wmax-Wmin+1));
	                printDelay = (((*counter)*100)/(Wmax-Wmin+1));    			
	    		}
	    		if(posWset == 1 && (((*counter)*100)/lines)%5 == 0 && printDelay != (((*counter)*100)/lines)) {
	    			printf("\033[A\33[2K");
	    			printf("Processing");
	            	for(i=0;i<=printCount;i++)
	    				printf(".");
	            	printf("\n");
	            	printCount++;
	            	if(printCount == 5)
	            		printCount = 0;
	                printDelay = (((*counter)*100)/lines);    			
	    		}
	            if(!strcmp(allID,allignmentId)) {
		    		if((Wset == 1 && lineCounter >= Wmin) || (posWset == 1 && pos >= posWmin)) {
		    			if(posWset == 1 && firstPos == 0)
		    				firstPos = pos;
		            	writeLine(*fpOut,line);
		        		(*counter)++;
		        	}

	            	if((Wset == 1 && lineCounter >= Wmax) || (posWset == 1 && pos >= posWmax) || lineCounter >= lines) {
		        		gzclose(*fpOut);
		        		if((*counter) == 0) {
				        	if(Wset == 1)
				        		printf("ERROR: Less snips than expected found\n");
				        	else if(posWset == 1)
				        		printf("ERROR: pos window is wrong\n");
				        	exit(0);
		        		}
			    		char outputFileNew[INFILENAMESIZE+30];
			    		char outputFile[INFILENAMESIZE+30];
	    				printf("\033[A\33[2K");
	                	printf("Processing: 100%%\n");
		        		if(Wset == 1)
			            	sprintf(outputFileNew,"%s%s_%d_%d.vcf.gz",outputPathName,allignmentId,Wmin,Wmax);
	            		else if (posWset == 1) 
			            	sprintf(outputFileNew,"%s%s_%d_%d.vcf.gz",outputPathName,allignmentId,firstPos,pos);
			            sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
						printf("Renaming %s to %s\n",outputFile,outputFileNew);
						int ret = rename(outputFile, outputFileNew);
						if(ret != 0) {
							printf("Error: unable to rename the file\n");
			                assert(0);
						}
						break;
	                }
	            }
	            else
	            	printf("WARNING: snip %d is allignment %s and will be skipped\n",lineCounter,allID);
	            free(allID);
	        }
	        else if(eof ==1) {
	        	if(Wset == 1)
	        		printf("ERROR: Less snips than expected found\n");
	        	else if(posWset == 1)
	        		printf("ERROR: pos windows is wrong\n");
	        	exit(0);
	        }
	    }
	}
}


int main(int argc, char** argv)
{
    int fileFormat=OTHER_FORMAT,size=0,first = 1,snipSize = 0,eof=0,counter = 0,maxcount=0, startingValue=0, endingValue=0, Wmin=0,Wmax=0, Wset = 0, posWmin=0, posWmax=0, posWset = 0,toSingleOutput=0, lineCounter = 0, printDelay = 0, firstPos;
    double time1,time2,time3,time4;
    char inputFileName[INFILENAMESIZE];
    char outputPathName[INFILENAMESIZE];
    char allignmentId[INFILENAMESIZE];
    char headerFileOld[INFILENAMESIZE+30];
    char headerFileNew[INFILENAMESIZE+30];
    gzFile fpIn=NULL, *fpOut=NULL, fpHeader=NULL;
    fpOut = (gzFile*)malloc(sizeof(gzFile));
    char ** headerLine1, ** headerLine2;
    int  lineLength = STRINGLENGTH,wordLength = STRINGLENGTH, lines = 0;
    commandLineParser(argc, argv, inputFileName, outputPathName, &fileFormat, &size, &Wmin, &Wmax, &Wset, &posWmin, &posWmax, &posWset, &toSingleOutput);

    char **line = (char**)malloc(sizeof(char*));
    *line = (char*)malloc(sizeof(char)*STRINGLENGTH);
    char **word = (char**)malloc(sizeof(char*));
    *word = (char*)malloc(sizeof(char)*STRINGLENGTH);

    fpIn = gzopen(inputFileName,"r");
    time1 = gettime();
    lines = countLines(fpIn, line, &lineLength);
    gzclose(fpIn);
    if(Wset == 1) {
	    if(Wmin < 1 || Wmin > lines) {
	    	printf("ERROR: Wmin must be between 1 and totalSnips, setting Wmin to 1 (default)\n");
	    	Wmin = 1;
	    }

	    if(Wmax < 1 || Wmax > lines) {
	    	printf("ERROR: Wmax must be between 1 and totalSnips, setting Wmax to total Snips (default)\n");
	    	Wmax = lines;
	    }

	    if(Wmin > Wmax) {
	    	printf("ERROR: Wmin must be <= Wmax, setting Wmin = Wmax\n");
	    	Wmin = Wmax;
	    }
	}
	else if(posWset == 1) {
		if(posWmax == 0)
			posWmax = 2147483646;
	}

    fpIn = gzopen(inputFileName,"r");
    time2 = gettime();
    if(toSingleOutput == 0) {
    	printf("Creating Multiple Files\n");
    	headerLine1 = (char **) malloc(sizeof(char*));
    	headerLine2 = (char **) malloc(sizeof(char*));
	    fpHeader = createHeaderFile(fpIn, outputPathName, headerLine1, headerLine2, line, &lineLength, word, &wordLength);
	    while(eof == 0)
	        createPart(fpIn, fpOut, outputPathName, &first, allignmentId, &snipSize,&eof, &counter, headerLine1, headerLine2, &size, &maxcount, &startingValue, &endingValue, line, &lineLength, word, &wordLength, lines, Wmin, Wmax, Wset, posWmin, posWmax, posWset, &lineCounter, &printDelay, &firstPos);
	    
	    time3 = gettime();
	    if(maxcount <0 || maxcount>counter)
	    	maxcount = counter;
	    gzprintf(fpHeader,"%s %d %d %d\n",allignmentId,maxcount,snipSize, counter);
	    gzclose(fpHeader);
	    sprintf(headerFileNew,"%s%s_header.vcf.gz",outputPathName,allignmentId);

	    sprintf(headerFileOld,"%sheader.vcf.gz",outputPathName);
		printf("Renaming %s to %s\n",headerFileOld, headerFileNew);
	    //printf("--closed file %s\n",headerFileNew);
	    int ret = rename(headerFileOld, headerFileNew);
		if(ret != 0) {
			printf("Error: unable to rename the file\n");
	        assert(0);
		}
    	free(*headerLine1);
    	free(*headerLine2);
    	free(headerLine1);
    	free(headerLine2);
	}
	else {
		char outputFile[INFILENAMESIZE+30];
        sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
	    printf("Creating Single File: %s\n",outputFile);
	    *fpOut = gzopen(outputFile,"w");
        if(*fpOut == NULL) {
            printf("failed to open file %s", outputFile);
            assert(0);
        }
		createSingleFile(fpIn, fpOut, outputPathName, &counter, line, &lineLength, word, &wordLength, lines, Wmin, Wmax, Wset, posWmin, posWmax, posWset);
	    time3 = gettime();
	}

    gzclose(fpIn);
    free(fpOut);
    free(*line);
    free(line);
    free(*word);
    free(word);
    time4 = gettime();
    printf("%d snips processed\n",counter);
    printf("get lines: %fs, create Parts: %fs, finalize: %fs, Total Time: %fs \n",time2 - time1, time3-time2, time4-time3, time4-time1);
    return 1;
}