#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "zlib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <dirent.h>
#include <errno.h>
#include <pthread.h>
#include <math.h>

#define INFILENAMESIZE 1000
#define STRINGLENGTH 1000
#define VCF_HLENGTH 9 // number of fields in the VCF header line
#define MAX_CHROM_NAME_VCF 100
#define MAX_STATES_VCF 5
#define MAXAFLENGTH 9

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3

#define STATESALL 8
#define ZERO '0'
#define ONE  '1'
#define GAP '-'
#define AD 'A'
#define CY 'C'
#define GU 'G'
#define TH 'T'
#define UN 'N'
#define ad 'a'
#define cy 'c'
#define gu 'g'
#define th 't'

#define min(a,b)     ( (a) > (b) ? (b) : (a) )
#define max(a,b)     ( (a) <= (b) ? (b) : (a) )

double gettime(void);
int getNextWord(gzFile fp, char * word, int *readEOL, int *readEOF, int *wordLength);
int getNextLine(gzFile fpIn, char ** line, int *readEOL, int *readEOF, int *lineLength);
void writeLine(gzFile fpOut,char **line);
int countLines(gzFile fpIn, char **line, int* lineLength);
int getWordFromString(char* line, char ** word, int *readEOL, int *wordLength, int *index);
int cmpfunc (const void * a, const void * b);
int sortList(char* inputListName,int** list, char** line, int *lineLength);

void readHeaderFile(gzFile fpOut, char* inputPathName, char ** headerLine1, char ** headerLine2, char* allignmentId, int* snipsPerFile, int* snipSize, int* totalSnips, int *posMin,int *posMax,char **line, int* lineLength);
void readPart(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, char **line, int* lineLength, char **word, int* wordLength);
void readPartIndexW(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int Wmin, int Wmax, char **line, int* lineLength, char **word, int* wordLength);
void readPartPosW(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int posWmin, int posWmax, char **line, int* lineLength, char **word, int* wordLength);
void readPartPosL(gzFile fpOut, char* inputPathName, char *allignmentId, int totalSnips, int snipSize, int *counter, char ** headerLine1, char ** headerLine2, int snipsPerFile, int ** inList,int listSize, char **line, int* lineLength, char **word, int* wordLength);
