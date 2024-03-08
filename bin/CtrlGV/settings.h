#ifndef _SETTINGS_H_
#define _SETTINGS_H_

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

#include <time.h>

double T_storage[] = {  //Duration of Cell valve close -> Cell valve open (seconds) for cell experiment
		      /* Used between June 2 and 3
                           15.,
      			   20.,
                           30.,
                           50.,
			   65.,
                           80.,
			   100.,
                           120.,
			   150.,
                           180.,
                           250.,
                           320.
		      */

		      // Planned to be used for the long runs
		      /*
                           15.,
                           30.,
                           50.,
                           80.,
			   110.,
                           140.,
                           180.,
                           220.,
			   260.,
			   300.,
                           340.
		      */
		      // For the noise test:
                          /* 15.,
			   400.,
			   30.,
			   340.,
			   50.,
			   300.,
			   80.,
			   260.,
			   110.,
			   220.,
			   140.,
			   180.,*/
                            15.,
                            400.,
                            50.,
                            300.,
                            110.,
                            220.,
                            140.
                            };

double T_filling[] = {  //Duration of UCN filling (seconds) for filling time search
                             20.,
                             40.,
                             60.,
			     80.,
                            100.,
                            120.,
                            200.
                            };

double T_FillingFixed        = 60.;  //Duration of UCN filling (seconds) for fixed span
double T_PreparateRV    = 1.;  //Duration of Rotary valve actuating for detection (seconds)
double T_DAQ            = 60.;  //Duration of UCN detection (seconds)
double T_Cleaning       = 20.;  //Duration of cell cleaning (seconds)
double T_FirstCleaning  = 80.;  //Duration of cell cleaning before first filling (seconds)
double T_StorageFixed = 20.;  //Duration of Cell valve close -> Cell valve open (seconds) for fixed span

char logfilename_cell[] = "./log2023-II__cell.txt";
char logfilename_fill[] = "./log2023-II_fill.txt";

#endif
