#include "MakeNiki.h"
#include "DrawNikiPMTs.cxx"

extern void InitGui();
Bool_t drawflag = 1;
Bool_t force_correct_kp = 1; //correct wrong arrival time of kp

VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT api("nikiglass","nikiglass display",initfuncs);

int main(int argc, char *argv[])
{

  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1001111);//to delete number pallete
  gStyle->SetOptFit(1110);
  TGaxis::SetMaxDigits(3);

  // ******************************************************************************
  // timer for information about running time of code
  TStopwatch timer;
  timer.Start();

  std::ostringstream os, osopt;
  char DataFileName[128];

  if(argc==1 || argc>3){
    std::cerr<<"Usage: ./MakeNiki [DataFileName] [Option Make Both:0(default), Only RPMT:1, Only single:2 (including UCN case), FRP:3]" <<std::endl;
    return -1;
  } else {
    os <<argv[1];
    sprintf(DataFileName,os.str().c_str());
    std::cout <<"Making ROOT file from "<<DataFileName<<std::endl;
    if(argc==2){
      nOption = 0;
    }else if(argc==3){
      osopt <<argv[2];
      nOption = atoi(osopt.str().c_str());
      std::cout <<"Option " <<nOption <<" was chosen."<<std::endl;
    }
  }

  if(nOption==0)   std::cout <<"Create RPMT and Single channel trees."<<std::endl;
  if(nOption==1)   std::cout <<"Create Only a RPMT tree."<<std::endl;
  if(nOption==2)   std::cout <<"Create Only a Single channel tree."<<std::endl;
  if(nOption==3)   std::cout <<"Create FRP tree."<<std::endl;
  if(nOption>3) {
    std::cerr <<"No such a option."<<std::endl;
    exit(1);
  }

  MakeNiki(DataFileName); //for RPMT
  if((nOption==1 || nOption==3) && drawflag)  DrawNikiPMTs(DataFileName); //for RPMT
  // ******************************************************************************
  // timer for information about running time of code
  cout<<"*********************************"<<endl;
  cout<<"ROOT - Time at the end of job:"<<endl;
  cout<<"CpuTime = "<<timer.CpuTime()<<" seconds"<<endl;
  cout<<"RealTime = "<<timer.RealTime()<<" seconds"<<endl;
  cout<<"*********************************"<<endl<<flush;
  timer.Stop();

  return 0;
}


Int_t MakeNiki(char* DataFileName)
{
  Int_t i=0;
  Int_t board,ch,e,dummy;
  Int_t kp=0,tp=0,up=0, mp=0;
  Int_t MRflag=-1;
  Double_t tof,ucntof;
  Double_t tof_rpmt;
  long int nlines=0,nlines_tp=0,nlines_kp=0;
  long long int t=0;            // use long long t t!!
  long long int t_rpmt=0;       // use long long t t!!
  long long int tzero=0;        // use long long t t!!
  long long int tzeroucn=0;     // use long long t t!!
  long long int t_tp=0;         // use long long t t!!
  long long int t_kp=0;         // use long long t t!!
  long long int difft=0;        // use long long t t!!
  Int_t a=0, b=0, c=0, d=0, f=0;
  Double_t x=0, y=0;

  TString filestr = gSystem->BaseName(DataFileName);
  Int_t numstr =  filestr.Last('_') ;
  TString dummystr = filestr(numstr+1,3);
  Int_t filenum = dummystr.Atof();

  struct sCntl sRpmt;
  sRpmt=ClearCntl(sRpmt);

  TString title_str;
  if(nOption==0 || nOption==1) title_str = "treeRPMT";
  else if (nOption==3) title_str = "treeFRP";

  //RPMT root file
  TTree *tup = new TTree("T",title_str);
  tup->Branch("a",     &a,     "a/I");
  tup->Branch("b",     &b,     "b/I");
  tup->Branch("c",     &c,     "c/I");
  tup->Branch("d",     &d,     "d/I");
  tup->Branch("f",     &f,     "f/I");
  tup->Branch("t",     &t_rpmt,"t/L");
  tup->Branch("tof",   &tof_rpmt,"tof/D");
  tup->Branch("kp",    &kp,    "kp/I");
  tup->Branch("tp",    &tp,    "tp/I");
  tup->Branch("up",    &up,    "up/I");
  tup->Branch("mp",    &mp,    "mp/I");
  tup->Branch("MRflag",&MRflag,"MRflag/I");
  tup->Branch("x",     &x,      "x/D");
  tup->Branch("y",     &y,      "y/D");

  //other root file
  TTree *six = new  TTree("Six","");
  six->Branch("e",     &e,     "e/I");
  six->Branch("t",     &t,     "t/L");
  six->Branch("tof",   &tof,   "tof/D");
  six->Branch("ucntof",&ucntof,"ucntof/D");
  six->Branch("kp",    &kp,    "kp/I");
  six->Branch("tp",    &tp,    "tp/I");
  six->Branch("up",    &up,    "up/I");
  six->Branch("mp",    &mp,    "mp/I");
  six->Branch("MRflag",&MRflag,"MRflag/I");
  six->Branch("ch",    &ch,    "ch/I");
  six->Branch("board", &board, "board/I");

  FILE *fp;

  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }

  // data read start
  while(1){

    // read a line and set parameters
    Int_t retvalue = fscanf(fp, "%d,%d,%d,%lld,%d",&board,&ch,&e,&t,&dummy) ;

    //open new _xxx.dat file if the file is at the end.
    if(retvalue==EOF){
      filenum++;
      fclose(fp);
      TString newfilestr  = filestr(0,numstr+1) + Form("%03d.dat",filenum);
      fp = fopen(newfilestr.Data(), "r");
      if(fp){
	cout  << newfilestr.Data() <<" was opened continuously"<<endl;
      } else{
	cout  << newfilestr.Data() <<" was not opened. No more files."<<endl;
	break;
      }
    }

    //show line number
    if (nlines%1000000==0){
      cout <<" nlines =  "<<nlines<<", kp = "<<kp<<endl;
    }

    ///// process for RPMT or FRP /////
    //if time is later than sett+dt, then Fill data in tree
    if(t-sRpmt.sett>dt || !(ch==rpmt_ch1 || ch==rpmt_ch2 || ch==rpmt_ch3 || ch==rpmt_ch4) ){

      Double_t RPMT_X = -1., RPMT_Y = -1.;
      for(i=0;i<4;i++) if(sRpmt.iFlag[i]<2)sRpmt.iSumFlag+=sRpmt.iFlag[i];
      f=sRpmt.iSumFlag;
      //a=x1, b=x2, c=y1, d=y2
      a=sRpmt.iPulse[0]; b=sRpmt.iPulse[1]; c=sRpmt.iPulse[2]; d=sRpmt.iPulse[3];

      if(nOption==0 || nOption==1){       //for RPMT
	Int_t sum_ab = a + b;
	Int_t sum_cd = c + d;
	if(sum_ab!=0)  RPMT_X = Double_t(a)/Double_t(sum_ab);
	if(sum_cd!=0)  RPMT_Y = Double_t(c)/Double_t(sum_cd);
      } else if (nOption==3){             //for FRP
	Int_t sum_abcd = a + b + c + d;
	if(sum_abcd!=0)  RPMT_X = Double_t(b+d)/Double_t(sum_abcd);
	if(sum_abcd!=0)  RPMT_Y = Double_t(c+d)/Double_t(sum_abcd);
      }

      t_rpmt = sRpmt.sett;
      if(kp==0) tof_rpmt = -9999.; //ignore tof before first kp
      else tof_rpmt=sRpmt.sett-tzero;
      x=RPMT_X;      y=RPMT_Y;

      //fill tup tree for PMT fils
      if(nOption==0 || nOption==1 || nOption==3) tup->Fill();
      // initialize structure
      sRpmt = ClearCntl(sRpmt);
      //reseting sett
      sRpmt.sett=t;
    }
    ///// end process for RPMT or FRP /////

    //ch flag control
    if(multi_kp==1){
      for(Int_t i=0; i<num_kps;i++){
	if(ch==kp_chs[i] && (t - t_kp > next_pulse) ){
	  ch_flag = 0;
	  break;
	} // kp
      }
      for(Int_t i=0; i<num_tps;i++){
	if(ch==tp_chs[i] && (t - t_tp > next_pulse) ){
	  ch_flag = 1;
	  break;
	}
      }// tp
    } else {
      if(ch==kp_ch) ch_flag = 0; // kp
      else if(ch==tp_ch) ch_flag = 1; // tp
      else if(ch==up_ch) ch_flag = 2; // tp
    }

    //This counts the number of timing pulse.
    if(ch_flag==1 && e > kp_LLD && e<kp_HLD){
      tp++;
      t_tp=t;
      nlines_tp = nlines;
      ch_flag = -1;
    }
    //This counts the number of kicker pulse.
    else if(ch_flag==0 && e > kp_LLD && e<kp_HLD) {
      if(MRflag == 2) {
	MRflag = 0; //kp with frame overlap
      } else if(MRflag == 1) {
	MRflag = 2; //first kp after MR (no frame overlap)
	mp++;
      }
      kp++;
      t_kp=t;
      difft = t_kp-t_tp;
      nlines_kp = nlines;
      ch_flag = -1;
      //error rejection for delayed kp time
      if( (nlines_kp-nlines_tp==1)&&difft>2 ){
	cout << "Wrong time stamp in kp: "<<kp;
	if(force_correct_kp) {
	  cout << " Corrected to tp time: "<<t_tp<<endl;
	  t = t_tp;
	}
      }
      tzero=t;
    }
    //This counts the number of ucn phase pulse.
    else if(ch==up_ch && e > kp_LLD && e<kp_HLD){
      up++;
      tzeroucn=t;
      ch_flag = -1;
    }


    //for other ch events
    if(kp==0) tof = -9999.; //ignore tof before first kp
    else tof = (Double_t)(t-tzero);

    if(up==0) ucntof = -1.; //ignore tof before first kp
    else ucntof = (Double_t)(t-tzeroucn);

    if(tof*timebin>repitition) MRflag = 1; //for tof tail flag, in case no kp comes in the repitition.

    if(nOption==0 || nOption==2)six->Fill();

    //fill as RPMT data
    Int_t rpmtch = -1;
    if(ch==rpmt_ch1) rpmtch=0;
    else if(ch==rpmt_ch2) rpmtch=1;
    else if(ch==rpmt_ch3) rpmtch=2;
    else if(ch==rpmt_ch4) rpmtch=3;
    if(rpmtch >= 0){
      sRpmt.iFlag[rpmtch]++;
      sRpmt.iPulse[rpmtch]=e;
    }
    //next line
    nlines++;
  }
  //end of while

  //for last is not taken into account

#define CREATE_ROOT 1
#if CREATE_ROOT
  TString ROOTstr = filestr(0,19);
  ROOTstr += ".root";
  TFile file(ROOTstr.Data(),"RECREATE");
  if(nOption==0 || nOption==1)tup->Write();
  if(nOption==0 || nOption==2)six->Write();
  if(nOption==3) tup->Write();
  cout << ROOTstr.Data() << " has been created" << endl;
  file.Close();
#endif

  std::cout << "kp = " << kp << std::endl;


  return 0;

}
