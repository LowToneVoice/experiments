#include <assert.h>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <map>

#include<TGaxis.h>

#include<TROOT.h>
#include<TCut.h>
#include<TFile.h>
#include<TMath.h>
#include<TLegend.h>
#include<THStack.h>
#include<TTreeFormula.h>
#include<TRandom.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TApplication.h>

#include <TRandom.h>/* header file for gRandom */
#include <TCanvas.h>/* header file for TCanvas */
#include <TMath.h>  /* header file for TMath */

#include <TH1.h>    /* header file for 1-d histogram */
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>

#include <TGraph.h>
#include <TGraph2D.h>
#include "TFile.h"
#include "TLegend.h"

#include "TTree.h"
#include "TVector3.h"
#include "TNtuple.h"

#include <time.h>

#define barn 1e-24
#define avogadro 6.02e23
#define cmTOfm   1e13
#define MevTOneV 1e15
#define mg 1e-3
#define g 1.
#define cm3 1.

int main(int argc, char *argv[]){

double OpenGVTime         = 3000000;
double StorageTime        = 50000000;
double ImportTime1       = 100000000;
double ImportTime2       = 200000000;
double ImportTime3       = 300000000;
double ImportTime4       = 400000000;
double ImportTime5       = 500000000;
double DAQTime            = 30000000;
double CleaningTime       = 70000000;
double DAQAndCleaningTime = DAQTime + CleaningTime;
double FirstCleaningTime  = 100000000;

double ImportTime        = ImportTime1;

FILE *fp;
if((fp=fopen("log.dat","w"))==NULL){
printf("file cannot be opened!\n");
exit(1);
}

clock_t start;
time_t timer;

int i=0,j=0,run=0;

time(&timer);
printf("Cleaning Start at %s",ctime(&timer));
fprintf(fp,"Cleaning Start @ %s",ctime(&timer));

system("GV10");
time(&timer);
printf("GV10 is performed @ %s",ctime(&timer));
fprintf(fp,"GV10 is performed @ %s",ctime(&timer));
start=clock();

while(1){

if(i==0){

//EQUIBRILIUM SECTION
//Now is at time limit to release the UCNs confined in storage chamber
//Next step is a time to realize equilibrium state
///////////////////////////////////////////////////
if( (run==0 && (clock()-start) >= FirstCleaningTime) || (run!=0&& (clock()-start)>=DAQAndCleaningTime ) ){
run++;
time(&timer);
printf("Run %d Start at %s",run,ctime(&timer));
fprintf(fp,"Run %d Start @ %s",run,ctime(&timer));
system("GV11");
time(&timer);
printf("GV11 is performed @ %s",ctime(&timer));
fprintf(fp,"GV11 is performed @ %s",ctime(&timer));
i++;
start=clock();
}

}else if(i==1){

if( (clock()-start) >= OpenGVTime ){
//ACCUMULATE SECTION
//Now is at time limit to make storage chamber equilibrium
//Next step is to accumulate UCN Flux in storage chamber
system("GV01");
time(&timer);
printf("GV01 is performed @ %s",ctime(&timer));
fprintf(fp,"GV01 is performed @ %s",ctime(&timer));
i++;
start=clock();

ImportTime = (j==0)?ImportTime1:(j==1)?ImportTime2:(j==2)?ImportTime3:(j==3)?ImportTime4:ImportTime5;

if(j==4){
j=0;
}else{
j++;
}

}//

}else if(i==2){

if( (clock()-start) >= ImportTime ){
//STORAGE SECTION
//Now is at time limit to accumulate UCNs in stogare chamber
//Next step is to confine UCNs in storage chamber
system("GV00");
time(&timer);
printf("GV00 is performed @ %s",ctime(&timer));
fprintf(fp,"GV00 is performed @ %s",ctime(&timer));
i++;
start=clock();

}//

}else if(i==3){

if( (clock()-start) >= StorageTime ){
//RELEASE SECTION
//Now is at time limit to confine UCNs in storage chamber
//Next step is to release the UCNs confined in storage chamber
system("GV10");
time(&timer);
printf("GV10 is performed @ %s",ctime(&timer));
fprintf(fp,"GV10 is performed @ %s",ctime(&timer));
i=0;
start=clock();
}//

}

if( run > 21)
break;

}//end while

time(&timer);
printf("Run End at %s", ctime(&timer) );

system("GV11");

fclose(fp);

return 0;
}
