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

//using namespace std;

int main(){

//printf("%d\n",CLOCS_PER_SEC);
FILE *fp;
if((fp=fopen("log.dat","w"))==NULL){
printf("file cannot be open!\n");
exit(1);
}

system("GV11");

clock_t start0,start;
time_t timer;

time(&timer);
printf("Run Start at %s", ctime(&timer) );
fprintf(fp,"Run Start @ %s",ctime(&timer));

int i=0;

system("GV10");
time(&timer);
printf("GV10 is performed @ %s",ctime(&timer));
fprintf(fp,"GV01 is performed @ %s",ctime(&timer));

start=clock();
start0=start;

while(1){

if( (clock()-start) == 100000000){
printf("now spend time(us):%lf\n",(double)(clock()-start)/1000000);
i++;

if(i==1){
system("GV01");
time(&timer);
printf("GV01 is performed @ %s",ctime(&timer));
fprintf(fp,"GV10 is performed @ %s",ctime(&timer));
}else{
system("GV10");
time(&timer);
printf("GV10 is performed @ %s",ctime(&timer));
fprintf(fp,"GV01 is performed @ %s",ctime(&timer));
i=0;
}

start=clock();
if( (clock()-start0) > 1000000000000)
break;

}//end if clock-start

}//end while

time(&timer);
printf("Run End at %s", ctime(&timer) );

fclose(fp);

return 0;
}
