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

int runmax = 599;

double FillingTime         = 40. * 1000000.;
double PreparateRVTime         = 2. * 1000000.;
double StorageTime1       = 25. * 1000000.;
double StorageTime2       = 45. * 1000000.;
double StorageTime3       = 75. * 1000000.;
double StorageTime4       = 105. * 1000000.;
double StorageTime5       = 75. * 1000000.;
double StorageTime6       = 125. * 1000000.;
double DAQTime            = 40. * 1000000.;
double CleaningTime       = 20. * 1000000.;
double DAQAndCleaningTime = DAQTime + CleaningTime;
double FirstCleaningTime  = 60. * 1000000.;

/*
//Valve operation test
double FillingTime         = 10. * 1000000.;
double PreparateRVTime         = 2. * 1000000.;
double StorageTime1       = 15. * 1000000.;
double StorageTime2       = 16. * 1000000.;
double StorageTime3       = 17. * 1000000.;
double StorageTime4       = 18. * 1000000.;
double StorageTime5       = 19. * 1000000.;
double StorageTime6       = 20. * 1000000.;
double DAQTime            = 10. * 1000000.;
double CleaningTime       = 10. * 1000000.;
double DAQAndCleaningTime = DAQTime + CleaningTime;
double FirstCleaningTime  = 10. * 1000000.;
*/
double StorageTime=StorageTime1;

if(argc==10){
printf("argc=%d,argv[0]=%s,argv[1]=%s,argv[2]=%s,argv[3]=%s,argv[4]=%s,argv[5]=%s,argv[6]=%s,argv[7]=%s,argv[8]=%s,argv[9]=%s\n",argc,argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
FillingTime         = (int)(atof(argv[1])*1e6);
PreparateRVTime         = (int)(atof(argv[2])*1e6);
StorageTime1       = (int)(atof(argv[3])*1e6);
StorageTime2       = (int)(atof(argv[4])*1e6);
StorageTime3       = (int)(atof(argv[5])*1e6);
StorageTime4       = (int)(atof(argv[6])*1e6);
DAQTime            = (int)(atof(argv[7])*1e6);
CleaningTime       = (int)(atof(argv[8])*1e6);
FirstCleaningTime  = (int)(atof(argv[9])*1e6);
DAQAndCleaningTime = DAQTime + CleaningTime;
}


FILE *fp;
if((fp=fopen("log2021.dat","a"))==NULL){
printf("file cannot be opened!\n");
exit(1);
}

clock_t start;
time_t timer;

int i=0,j=0,run=0;

time(&timer);
printf("Cleaning Start at %s",ctime(&timer));
fprintf(fp,"Cleaning Start @ %s",ctime(&timer));

system("../GV11");
time(&timer);
printf("GV11 (RV:Detection, GV:Open) is performed @ %s",ctime(&timer));
fprintf(fp,"GV11 (RV:Detection, GV:Open) is performed @ %s",ctime(&timer));
start=clock();

printf("\n");
fprintf(fp,"\n");

while(1){

if(i==0){

//ACCUMULATE SECTION
//Now is at time limit to release the UCNs confined in storage chamber
//Next step is to accumulate UCN Flux in storage chamber
///////////////////////////////////////////////////
if( (run==0 && (clock()-start) >= FirstCleaningTime) || (run!=0&& (clock()-start)>=DAQAndCleaningTime ) ){
run++;
time(&timer);
printf("Run %d Start at %s",run,ctime(&timer));
fprintf(fp,"Run %d Start @ %s",run,ctime(&timer));
system("../GV01");
time(&timer);
printf("GV01 (RV:Filling, GV:Open) is performed @ %s",ctime(&timer));
fprintf(fp,"GV01 (RV:Filling, GV:Open) is performed @ %s",ctime(&timer));
i++;
start=clock();
}

}else if(i==1){

if( (clock()-start) >= FillingTime ){
//STORAGE SECTION
//Now is at time limit to accumulate UCNs in stogare chamber
//Next step is to confine UCNs in storage chamber
system("../GV00");
time(&timer);
printf("GV00 (RV:Filling, GV:Close) is performed @ %s",ctime(&timer));
fprintf(fp,"GV00 (RV:Filling, GV:Close) is performed @ %s",ctime(&timer));
i++;
start=clock();
}//

}else if(i==2){

if( (clock()-start) >= PreparateRVTime ){
//STORAGE SECTION
//Now is at time limit to close stogare chamber
//Next step is to preparate Rotary valve to detection mode
system("../GV10");
time(&timer);
printf("GV10 (RV:Detection, GV:Close) is performed @ %s",ctime(&timer));
fprintf(fp,"GV10 (RV:Detection, GV:Close) is performed @ %s",ctime(&timer));
i++;
start=clock();
StorageTime = (j==0)?StorageTime1:(j==1)?StorageTime2:(j==2)?StorageTime3:(j==3)?StorageTime4:(j==4)?StorageTime5:StorageTime6;
if(StorageTime < PreparateRVTime) StorageTime = PreparateRVTime;

printf("Storage time is %d sec.\n", int(StorageTime/1000000.));
fprintf(fp,"Storage time is %d sec.\n", int(StorageTime/1000000.));

if(j==3){//loop control
j=0;
}else{
j++;
}

}//

}else if(i==3){

if( (clock()-start) >= (StorageTime - PreparateRVTime)){
//RELEASE SECTION
//Now is at time limit to confine UCNs in storage chamber
//Next step is to release the UCNs confined in storage chamber
system("../GV11");
time(&timer);
printf("GV11 (RV:Detection, GV:Open) is performed @ %s",ctime(&timer));
fprintf(fp,"GV11 (RV:Detection, GV:Open) is performed @ %s",ctime(&timer));
i=0;
start=clock();

printf("\n");
fprintf(fp,"\n");

}//

}

if( run > runmax)
break;

}//end while


time(&timer);
printf("Run End at %s", ctime(&timer) );

system("../GVoff");
printf("GVoff (RV:Filling, GV:Close) is performed @ %s",ctime(&timer));
fprintf(fp,"GVoff (RV:Filling, GV:Close) is performed @ %s",ctime(&timer));
printf("Measurement ends. @ %s",ctime(&timer));
fprintf(fp,"Measurement ends. @ %s",ctime(&timer));

fclose(fp);

return 0;
}
