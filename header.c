#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "zlib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <dirent.h>
#include <errno.h>

#define INFILENAMESIZE 1000
#define STRINGLENGTH 1000
#define VCF_HLENGTH 9 // number of fields in the VCF header line
#define MAX_CHROM_NAME_VCF 100

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3

double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime , NULL);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}


int getNextWord(gzFile fp, char * word, int *readEOL, int *readEOF, int *wordLength)
{
    *readEOL = *readEOF = 0;
    
    int i=0;
    word[i] = gzgetc(fp);
    
    if(word[i] == EOF)
    {
        *readEOF = 1;
        word[i] = '\0';
        return 0;
    }
    
    if(word[i] == 10 || word[i] == 13)
    {
        *readEOL = 1;
        word[i] = '\0';
        return 0;
    }
    while(!(word[i] == 9 || word[i] == 32) && !(word[i] == 10 || word[i] == 13)  && word[i] != EOF)
    {
        if(i+1 >= (*wordLength) )
        {
            (*wordLength) = (*wordLength) << 1; 
            
            word = realloc( word, (*wordLength) * sizeof(char) );       
        }
        
        i++;
        word[i] = gzgetc(fp);
    }
    
    if(word[i] == 9 || word[i] == 32) {
        while(word[i] == 9 || word[i] == 32)
        {
            if(i+1 >= (*wordLength) )
            {
                (*wordLength) = (*wordLength) << 1; 
                
                word = realloc( word, (*wordLength) * sizeof(char) );       
            }
            
            i++;
            word[i] = gzgetc(fp);                
        }
    }
    
    if(word[i] == 10 || word[i] == 13)
        *readEOL = 1;
    else
        gzungetc(word[i],fp);  
    
    word[i] = '\0';   
    return 1;    
}
//first try, a little slower than getNextLine
int getNextLine2(gzFile fp, char ** line, int *readEOL, int *readEOF, int *lineLength)
{
    
    *readEOL = *readEOF = 0;
    
    char ent = gzgetc(fp);
    
    int i=0;
    
    if(ent == EOF)
    {
        *readEOF = 1;
        (*line)[0] = '\0';
        return 0;
    }
    
    if(ent == 10 || ent == 13)
    {
        *readEOL = 1;
        (*line)[0] = '\0';
        return 0;
    }
    
    while(!(ent == 10 || ent == 13)  && ent != EOF)
    {
        
        if(i+1 >= (*lineLength) )
        {
            (*lineLength) = (*lineLength) << 1; 
            
            (*line) = realloc( (*line), (*lineLength) * sizeof(char) );     
        }
        
        (*line)[i++] = ent;
        
        ent = gzgetc(fp);
    }
    
    (*line)[i] = '\0';
    
    if(ent == 10 || ent == 13)
        *readEOL = 1;
    else
        *readEOF = 1;        
    
    return 1;
}
//second try, a bit better
int getNextLine(gzFile fpIn, char ** line, int *readEOL, int *readEOF, int *lineLength)
{
    *readEOL = *readEOF = 0;

    if(gzgets(fpIn, *line, *lineLength) != NULL) {
        if((*line)[strlen(*line)-1] != '\n') {
            while((*line)[strlen(*line)-1] != '\n'){
                char * tmpLine = malloc( (*lineLength) * sizeof(char));
                (*lineLength) = (*lineLength) << 1; 
                (*line) = realloc( (*line), (*lineLength) * sizeof(char) );
                if(gzgets(fpIn, tmpLine, (*lineLength)>>1) != NULL) {
                    strcat((*line), tmpLine);
                }
                else {
                    (*line)[strlen(*line)] = '\0';
                    *readEOF = 1;
                    if(strlen(*line) >= 3) {
                        *readEOL = 1;
                        return 1;
                    }
                    free(tmpLine);
                    return 0;
                }
                free(tmpLine);
            }
        }      
        *readEOL = 1;
        (*line)[strlen(*line)-1] = '\0';
        if(strlen(*line) >= 3) 
            return 1;
        else
            return 0;
    }
    else {
        *readEOF = 1; 
        return 0; 
    }
    return 1;
}

//first try, EXTREMELY slower than writeLine
void writeLine2(gzFile fpOut,char **line) {
    int i;
    for (i=0;i<strlen(*line);i++)
        gzputc(fpOut,(*line)[i]);
    gzputc(fpOut,'\n');
}
//second try, the best
void writeLine(gzFile fpOut,char **line) {
    if(gzputs(fpOut,*line) < 0) {
        printf("error writing to file");
        assert(0);
    }
    gzputc(fpOut,'\n');
}
//not really needed
int countLines(gzFile fpIn, char **line, int* lineLength)
{
  	int lines=0, notHeader = 0;
  	while(gzgets(fpIn, *line, *lineLength) != NULL)
	{
		if((*line)[strlen(*line)-1] != '\n') {
			while((*line)[strlen(*line)-1] != '\n'){
	            char * tmpLine = malloc( (*lineLength) * sizeof(char));
                (*lineLength) = (*lineLength) << 1; 
                (*line) = realloc( (*line), (*lineLength) * sizeof(char) );
                if(gzgets(fpIn, tmpLine, (*lineLength)>>1) != NULL) {
                    strcat((*line), tmpLine);
                }
                else {
                    int tmp = strlen(*line);
                    (*line)[tmp] = '\n';
                    (*line)[tmp+1] = '\0';
                }
	            free(tmpLine);
	        }
		}
        if(strlen(*line) >=6 ) {
            if(notHeader ==1)
                lines++;

            if((*line)[0] == '#' && (*line)[1] == 'C' && (*line)[2] == 'H' && (*line)[3] == 'R' && (*line)[4] == 'O' && (*line)[5] == 'M')
                notHeader = 1;
        }
	}
	return lines;
}

int getWordFromString(char* line, char * word, int *readEOL, int *wordLength, int *index)
{
    *readEOL = 0;
    
    int i=0;
    char ent = line[(*index)++];
    while(ent==' '|| ent == 9) // horizontal tab
        ent = line[(*index)++]; 
    if(ent == '\0')
    {
        *readEOL = 1;
        word[0] = '\0';
        return 0;
    }
    
    while(!(ent == 9 || ent == 32) && ent != '\0')
    {
        
        if(i+1 >= (*wordLength) )
        {
            (*wordLength) = (*wordLength) << 1; 
            
            word = realloc( word, (*wordLength) * sizeof(char) );       
        }
        
        word[i++] = ent;
        
        
        ent = line[(*index)++];
    }
    
    word[i] = '\0';
    
    if(ent == 10 || ent == 13 || ent =='\0')
        *readEOL = 1;
    
    return 1;
}