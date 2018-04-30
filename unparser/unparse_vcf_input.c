#include "header.h"

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

void readPart(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, char **line, int* lineLength, char **word, int* wordLength)
{
    int tmpId,tmpMin,tmpMax,tmpPosMin,tmpPosMax, fileFound, eol, eof, status,eop=0,index;

    char inputFileName[INFILENAMESIZE+30];
    DIR *dp;
    gzFile fpIn;
    struct dirent *entry;
    struct stat statbuf;

    while(eop==0) {
        tmpId=0;
        tmpMin=0;
        tmpMax=0;
        tmpPosMin=0;
        tmpPosMax=0;
        fileFound=0;
        eol=0;
        eof=0;

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
                    if( (tmpMin <= (1 + (*counter))  && (1 + (*counter)) <= tmpMax) && 
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
                                    index = 0;
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//chrom
                                    assert(status);
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//pos
                                    assert(status);
                                    writeLine(fpOut,line);
                                    counter++;
                                }
                            }while(eof == 0 && (*counter)< totalSnips);
                            if((*counter)>= totalSnips)
                                eop = 1;
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
            printf("File missing or wrong window\n");
            assert(0);
        }
    }
}

void readPartIndexW(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int Wmin, int Wmax, char **line, int* lineLength, char **word, int* wordLength)
{
    int tmpId,tmpMin,tmpMax,tmpPosMin,tmpPosMax, fileFound, eol, eof , status,eop=0,index;

    char inputFileName[INFILENAMESIZE+30];
    DIR *dp;
    gzFile fpIn;
    struct dirent *entry;
    struct stat statbuf;

    while(eop==0) {
        tmpId=0;
        tmpMin=0;
        tmpMax=0;
        tmpPosMin=0;
        tmpPosMax=0;
        fileFound=0;
        eol=0;
        eof=0;

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
                    if( (tmpMin <= (Wmin + (*counter))  && (Wmin + (*counter)) <= tmpMax) && 
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
            printf("Opening file %s, %d %%\n",inputFileName,((*counter)*100)/(Wmax-Wmin+1));
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
                                    index = 0;
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//chrom
                                    assert(status);
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//pos
                                    assert(status);
                                    if(tmpMin >= Wmin && tmpMin <= Wmax){
                                        writeLine(fpOut,line);
                                        counter++;
                                    }
                                    tmpMin++;
                                }
                            }while(eof == 0 && (*counter)< (Wmax - Wmin + 1));
                            if((*counter)>= (Wmax - Wmin + 1))
                                eop = 1;
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
            printf("File missing or wrong window\n");
            assert(0);
        }
    }
}

void readPartPosW(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int posWmin, int posWmax, char **line, int* lineLength, char **word, int* wordLength)
{
    int tmpId,tmpMin,tmpMax,tmpPosMin,tmpPosMax, fileFound, pos, eol, eof , status,eop=0,tmpPosMinIndex=0,index;

    char inputFileName[INFILENAMESIZE+30];
    DIR *dp;
    gzFile fpIn;
    struct dirent *entry;
    struct stat statbuf;

    while(eop==0) {
        tmpId=0;
        tmpMin=0;
        tmpMax=0;
        tmpPosMin=0;
        tmpPosMax=0;
        fileFound=0;
        eol=0;
        eof=0;

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
                        ((tmpPosMin <= posWmin   && posWmin <= tmpPosMax && tmpPosMinIndex == 0) || (tmpMin <= tmpPosMinIndex  && tmpPosMinIndex <= tmpMax)) && 
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
            printf("Opening file %s\n",inputFileName);

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
                                    index = 0;
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//chrom
                                    assert(status);
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//pos
                                    assert(status);
                                    pos = atoi(*word);
                                    if(pos >= posWmin && pos <= posWmax){
                                        writeLine(fpOut,line);
                                        counter++;
                                    }
                                }
                                else {
                                    pos = 0;
                                }
                            }while(eof == 0 && pos < posWmax );
                            tmpPosMinIndex = tmpMax+1;
                            if(pos >= posWmax)
                                eop = 1;
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
            printf("File missing or wrong window\n");
            assert(0);
        }
    }
}

void readPartPosL(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int ** inList,int listSize, char **line, int* lineLength, char **word, int* wordLength)
{
    int tmpId,tmpMin,tmpMax,tmpPosMin,tmpPosMax, fileFound, pos, eol, eof , status,eop=0,inListPos=0,index;

    char inputFileName[INFILENAMESIZE+30];
    DIR *dp;
    gzFile fpIn;
    struct dirent *entry;
    struct stat statbuf;

    while(eop==0) {
        tmpId=0;
        tmpMin=0;
        tmpMax=0;
        tmpPosMin=0;
        tmpPosMax=0;
        fileFound=0;
        eol=0;
        eof=0;

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
                    if((tmpPosMin <= (*inList)[inListPos] && (*inList)[inListPos] <= tmpPosMax) && 
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
            printf("Opening file %s, %d %%\n",inputFileName,((*counter)*100)/listSize);

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
                                    index = 0;
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//chrom
                                    assert(status);
                                    status = getWordFromString(*line, word, &eol, wordLength, &index);//pos
                                    assert(status);
                                    pos = atoi(*word);
                                    if(pos == (*inList)[inListPos]){
                                        writeLine(fpOut,line);
                                        counter++;
                                        (*inList)[inListPos] = 0;
                                        inListPos++;
                                    }
                                    else if( pos > (*inList)[inListPos])
                                        inListPos++;

                                    if(tmpPosMin > (*inList)[inListPos] || (*inList)[inListPos] > tmpPosMax)
                                        eof = 1;
                                }
                                else {
                                    pos = 0;
                                }
                            }while(eof == 0 && inListPos < listSize);
                            if( inListPos >= listSize )
                                eop = 1;
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
            if(inListPos < listSize)
                inListPos++;
            else if(inListPos >= listSize)
                eop = 1;
        }
    }
    printf("\n");
    for(inListPos=0;inListPos<listSize;inListPos++)
        if((*inList)[inListPos] != 0)
            printf("\tWARNING: list pos %d not found\n",(*inList)[inListPos]);
    printf("\n\n");
}