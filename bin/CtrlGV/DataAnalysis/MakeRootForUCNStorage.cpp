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

#define barn 1e-24
#define avogadro 6.02e23
#define cmTOfm   1e13
#define MevTOneV 1e15
#define mg 1e-3
#define g 1.
#define cm3 1.

//using namespace std;

int main(int argc,char *argv[]){//macro start

	char filename[256]={0};
	strncpy(filename,argv[1],256);

	printf("filename=%s\n",filename);

	gStyle->SetOptStat(1111111);
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);

	FILE *fp;

	if(filename==NULL)exit(1);

	if((fp=fopen(filename,"r"))==NULL){
	printf("file cannt be open\n");
	exit(1);
	}

	TNtuple *ntuple = new TNtuple("ntuple","ntuple","e:ch:t:tof:gvt:gvid:run");
	TNtuple *GVT = new TNtuple("GVTiming","GVTiming","timing:GVID:run");

	char str[256];

	int i=0,run=0,k=0;

	bool GV1,GV2,key[4]={false,false,false,true};
	double xxx,e,ch,t,dummy;

	double t0=0;

	double t1[5]={0,0,0,0,0};

	bool cleaning=true;

	TH1I *KP[4];
	KP[0] = new TH1I("KP","KP",200,1,201);
	KP[1] = new TH1I("KP2","KP2",200,1,201);
	KP[2] = new TH1I("KP3","KP3",200,1,201);
	KP[3] = new TH1I("KP4","KP4",200,1,201);

	TH1I *TP[4];
	TP[0] = new TH1I("TP","TP",200,1,201);
	TP[1] = new TH1I("TP2","TP2",200,1,201);
	TP[2] = new TH1I("TP3","TP3",200,1,201);
	TP[3] = new TH1I("TP4","TP4",200,1,201);

	TH1I *BM[4];
	BM[0] = new TH1I("BM","BM",200,1,201);
	BM[1] = new TH1I("BM2","BM2",200,1,201);
	BM[2] = new TH1I("BM3","BM3",200,1,201);
	BM[3] = new TH1I("BM4","BM4",200,1,201);

	TH1F *T[4][5];

	char strstr[200],strstr2[200];

	for(int l=0;l<4;l++){

	memset(strstr,0,sizeof(strstr));

	     if(l==0) strcpy(strstr,"T");
	else if(l==1) strcpy(strstr,"T2");
	else if(l==2) strcpy(strstr,"T3");
	else if(l==3) strcpy(strstr,"T4");

	for(int m=0;m<5;m++){

	memset(strstr2,0,sizeof(strstr2));

	     if(m==0) sprintf(strstr2,"%s_0s",strstr);
	else if(m==1) sprintf(strstr2,"%s_1000s",strstr);
	else if(m==2) sprintf(strstr2,"%s_2000s",strstr);
	else if(m==3) sprintf(strstr2,"%s_3000s",strstr);
	else if(m==4) sprintf(strstr2,"%s_400s",strstr);

	T[l][m] = new TH1F(strstr2,strstr2,100,0,100e6);

	}//end m for

	}//end l for


	TH1I *Correction = new TH1I("Correction","Correction",1000,0,1000);
	TH1I *BG = new TH1I("BG","BG",200,1,201);
	TH1I *BG2 = new TH1I("BG2","BG2",200,1,201);
	TH1I *Run = new TH1I("Run","Run",200,1,201);

	TH1F *KPNormalized[5];
	KPNormalized[0] = new TH1F("KPNormalized_0s","KPNormalized_0s",200,0,200);
	KPNormalized[1] = new TH1F("KPNormalized_100s","KPNormalized_100s",200,0,200);
	KPNormalized[2] = new TH1F("KPNormalized_200s","KPNormalized_200s",200,0,200);
	KPNormalized[3] = new TH1F("KPNormalized_300s","KPNormalized_300s",200,0,200);
	KPNormalized[4] = new TH1F("KPNormalized_400s","KPNormalized_400s",200,0,200);

	TH1F *BMNormalized[5];
	BMNormalized[0] = new TH1F("BMNormalized_0s","BMNormalized_0s",200,0,200);
	BMNormalized[1] = new TH1F("BMNormalized_100s","BMNormalized_100s",200,0,200);
	BMNormalized[2] = new TH1F("BMNormalized_200s","BMNormalized_200s",200,0,200);
	BMNormalized[3] = new TH1F("BMNormalized_300s","BMNormalized_300s",200,0,200);
	BMNormalized[4] = new TH1F("BMNormalized_400s","BMNormalized_400s",200,0,200);



	//while loop start
	while((fgets(str,256,fp))!=NULL){
	sscanf(str,"%lf,%lf,%lf,%lf,%lf\n",&xxx,&ch,&e,&t,&dummy);

	i++;

	if(i%100000==0){
	printf("i=%d\n",i);
	printf("e=%lf:ch=%lf\n",e,ch);
	}
	
	if(ch==2) t0=t;

	if(ch==4) GV1=true;
	else if(ch==6) GV2=true;
	else if(ch==5) GV1=false;
	else if(ch==7) GV2=false;

	if(GV1&&GV2){


	if(key[0]){
	run++;
	t1[0]=t;
	printf("GV1=Open, GV2= Open, daq time=%e ,dt=%e\n",t/1e6,t-t1[3]);

	key[0]=false;
	key[1]=true;
	GVT->Fill(t1[0]-t1[3],0,run);
	}

	if(ch==15) Correction->Fill(run);
	if(ch==0) KP[0]->Fill(run);
	if(ch==1) TP[0]->Fill(run);
	if(ch==10) BM[0]->Fill(run);

	ntuple->Fill(e,ch,t,t-t0,t-t1[0],0,(float)run);

	}//end if GV1&&GV2

	if(!GV1&&GV2){

	if(key[1]){
	printf("GV1= Close, GV2= Open, daq time = %e ,dt=%e\n",t/1e6, t-t1[0]);
	t1[1]=t;
	key[1]=false;
	key[2]=true;
	GVT->Fill(t1[1]-t1[0],1,run);
	}

	if(ch==15) BG->Fill(run);
	if(ch==0) KP[1]->Fill(run);
	if(ch==1) TP[1]->Fill(run);
	if(ch==10) BM[1]->Fill(run);

	ntuple->Fill(e,ch,t,t-t0,t-t1[1],1,(float)run);

	}//end if !GV1&&GV2

	if(!GV1&&!GV2){


	if(key[2]){
	printf("GV1=Close, GV2= Close, daq time = %e, dt=%e\n",t/1e6, t-t1[1]);
	t1[2]=t;
	key[2]=false;
	key[3]=true;
	GVT->Fill(t1[2]-t1[1],2,run);
	}

	if(ch==15) BG2->Fill(run);
	if(ch==0) KP[2]->Fill(run);
	if(ch==1) TP[2]->Fill(run);
	if(ch==10) BM[2]->Fill(run);

	ntuple->Fill(e,ch,t,t-t0,t-t1[2],2,(float)run);

	}//end if !GV1&&!GV2

	if(GV1&&!GV2){


	if(key[3]){
	t1[3]=t;
	printf("GV1=Open, GV2= Close,daq time = %e, dt=%e\n",t/1e6, t-t1[2]);
	printf("--------------------------\n");
	if(run==0)t1[2]=0;
	key[3]=false;
	key[0]=true;
	GVT->Fill(t1[3]-t1[2],3,run);
	}

	if(ch==15) Run->Fill(run);
	if(ch==0) KP[3]->Fill(run);
	if(ch==1) TP[3]->Fill(run);
	if(ch==10) BM[3]->Fill(run);

	ntuple->Fill(e,ch,t,t-t0,t-t1[3],3,(float)run);

	}//end if GV1&&!GV2

	}//end while fgets


	char ROOTFILENAME[200];

	sscanf(filename,"%s",ROOTFILENAME);
	strcat(ROOTFILENAME,"_dataanalysis.root");

	bool ok[200];
	double ratio[200];
	bool first_id=true;
	double first;

	for(int i=1;i<201;i++){
	printf("The %d Bin Content of KP2=%lf\n",i, KP[1]->GetBinContent(i));
	printf("The %d Bin Content of BM2=%lf\n",i, BM[1]->GetBinContent(i));

	if(KP[1]->GetBinContent(i)>TP[1]->GetBinContent(i)*0.9){//1
	ok[i] = true;

	if(first_id){
	first = BM[1]->GetBinContent(i);
	printf("aaaaaaaaaaaaaaaaaaaaaaa\n");
	printf("aaaaaaaaaaaaaaaaaaaaaaa\n");
	printf("aaaaaaaaaaaaaaaaaaaaaaa\n");
	printf("first->GetBinContent(i=%lf)\n",BM[1]->GetBinContent(i));
	printf("aaaaaaaaaaaaaaaaaaaaaaa\n");
	printf("aaaaaaaaaaaaaaaaaaaaaaa\n");
	printf("aaaaaaaaaaaaaaaaaaaaaaa\n");
	first_id = false;
	}

	}//1
	else{//2
	//printf("i=%d\n",i);
	ok[i]=false;
	}//2

	}//end for i

	int j[5]={0,0,0,0,0};

	ratio[0]=0; ok[0]=false;

	for(int i=1;i<201;i++){

	if(ok[i]){
	ratio[i]=BM[1]->GetBinContent(i)/first;
	BMNormalized[(i+4)%5]->Fill(j[(i+4)%5],ratio[i]*Run->GetBinContent(i));
	printf("i=%d,(i+4)%5=%d:BMNormalized->GetBinContent(%d)=%lf\n",i,(i+4)%5,j[(i+4)%5]+1,BMNormalized[(i+4)%5]->GetBinContent(j[(i+4)%5]+1));
	j[(i+4)%5]+=1;
	}else{
	ratio[i]=0;
	}

	}//end for i

	float nt_t,nt_ch,nt_run,nt_gvt,nt_gvid;

	ntuple->SetBranchAddress("t",&nt_t);
	ntuple->SetBranchAddress("ch",&nt_ch);
	ntuple->SetBranchAddress("run",&nt_run);
	ntuple->SetBranchAddress("gvt",&nt_gvt);
	ntuple->SetBranchAddress("gvid",&nt_gvid);

	for(int i=0;i<ntuple->GetEntries();i++){

	ntuple->GetEntry(i);

	if(nt_ch==15){

	if(nt_gvid==3)
	T[3][((int)nt_run+4)%5]->Fill(nt_gvt,ratio[(int)nt_run]);
	else
	T[(int)nt_gvid][((int)nt_run+4)%5]->Fill(nt_gvt);

	}//end if nt_ch==10

	}//for int i

	T[3][0]->Scale(1/(double)j[0]);
	T[3][1]->Scale(1/(double)j[1]);
	T[3][2]->Scale(1/(double)j[2]);
	T[3][3]->Scale(1/(double)j[3]);
	T[3][4]->Scale(1/(double)j[4]);

	//double x[4]={40,80,120,160};
	//double x[5]={0,100,200,300,400};
	double x[5]={100,200,300,400,0};
	double y[5];

	for(int i=0;i<5;i++){
	y[i]=T[3][i]->Integral(0,30);
	}

	TGraph* tg = new TGraph(5,x,y);
	tg->SetMarkerStyle(22);
	tg->SetMarkerSize(2);
	tg->Draw("AP");
	tg->Fit("expo");
	
	TFile *File = TFile::Open(ROOTFILENAME,"RECREATE");
	ntuple->Write();
	Correction->Write();
	Run->Write();
	BG->Write();
	BG2->Write();
	KP[0]->Write();
	KP[1]->Write();
	KP[2]->Write();
	KP[3]->Write();
	GVT->Write();
	for(int l=0;l<4;l++){
	for(int m=0;m<5;m++){
	T[l][m]->Write();
	}}//end for l,m


}//End MakePhase
