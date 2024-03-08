#include "settings.h"

int main(int argc, char *argv[]){

  int runmax = 599;

  double PreparateRVTime    = T_PreparateRV * 1000000.;
  double DAQTime            = T_DAQ * 1000000.;
  double CleaningTime       = T_Cleaning * 1000000.;
  double DAQAndCleaningTime = DAQTime + CleaningTime;
  double FirstCleaningTime  = T_FirstCleaning * 1000000.;
  int n_filling = (int) sizeof(T_filling) / sizeof(T_filling[0]);
  int ii;
  double StorageTime = T_StorageFixed * 1000000.;
  double FillingTime = T_filling[0] * 1000000.;

  clock_t start;
  time_t timer;
  
  FILE *fp;
  if((fp=fopen(logfilename_fill,"a"))==NULL){
    printf("file cannot be opened!\n");
    exit(1);
  }

  time_t now = time(NULL);
  struct tm *lctime = localtime(&now);

  printf("\n########################################################################\n");
  printf("####            ###  ######  #####      #####        ####    #####  ####\n");
  printf("#########  ########  ######  ###  #####  ###  ######  ###  #  ####  ####\n");
  printf("#########  ########  ######  ###  ##########  ######  ###  ##  ###  ####\n");
  printf("#########  ########  ######  ###  ##########          ###  ###  ##  ####\n");
  printf("#########  ########  ######  ###  #####  ###  ######  ###  ####  #  ####\n");
  printf("#########  #########        #####       ####  ######  ###  #####    ####\n");
  printf("########################################################################\n");
  printf("\nStorage measurement starts at %d/%02d/%02d %d:%02d:%02d (JST).\n",
      lctime->tm_year+1900,lctime->tm_mon+1,lctime->tm_mday,lctime->tm_hour,lctime->tm_min,lctime->tm_sec);
  printf("========================== Start controll log ==========================\n");
  printf("#\n");
  printf("# RV = 'Rotary valve',  CV = 'Cell valve'\n");
  printf("# Preset parameters: \n");
  printf("# Duration of the first cleaning: %f s.\n",FirstCleaningTime / 1000000.);
  printf("# Duration from (RV:Filling, CV:Open) to (RV:Filling, CV:Close) is ...\n");
  for(ii=0;ii<n_filling;ii++){
    printf("#    Case %d:  %f s.\n", ii+1, T_filling[ii]);
  }
  printf("#    Total case number of filling time: %d\n", n_filling);
  printf("# Duration from (RV:Filling, CV:Close) to (RV:Detection, CV:Close): %f s.\n",PreparateRVTime / 1000000.);
  printf("# Duration from (RV:Filling, CV:Close) to (RV:Detection, CV:Open): %f s.\n",StorageTime / 1000000.);
  printf("# Duration from (RV:Detection, CV:Open) to (RV:Filling, CV:Open): %f s.\n",DAQAndCleaningTime / 1000000.);
  printf("#   >> Duration of DAQ: %f s.\n",DAQTime / 1000000.);
  printf("#   >> Duration of cleaning after DAQ: %f s.\n",CleaningTime / 1000000.);
  printf("#\n");
  printf("\n========================== Start valve controll ==========================\n\n");


  fprintf(fp,"\n\n\n");
  fprintf(fp,"\n########################################################################\n");
  fprintf(fp,"####            ###  ######  #####      #####        ####    #####  ####\n");
  fprintf(fp,"#########  ########  ######  ###  #####  ###  ######  ###  #  ####  ####\n");
  fprintf(fp,"#########  ########  ######  ###  ##########  ######  ###  ##  ###  ####\n");
  fprintf(fp,"#########  ########  ######  ###  ##########          ###  ###  ##  ####\n");
  fprintf(fp,"#########  ########  ######  ###  #####  ###  ######  ###  ####  #  ####\n");
  fprintf(fp,"#########  #########        #####       ####  ######  ###  #####    ####\n");
  fprintf(fp,"########################################################################\n");
  fprintf(fp,"\nStorage measurement starts at %d/%02d/%02d %d:%02d:%02d (JST).\n",
      lctime->tm_year+1900,lctime->tm_mon+1,lctime->tm_mday,lctime->tm_hour,lctime->tm_min,lctime->tm_sec);
  fprintf(fp,"========================== Start controll log ==========================\n");
  fprintf(fp,"#\n");
  fprintf(fp,"# RV = 'Rotary valve',  CV = 'Cell valve'\n");
  fprintf(fp,"# Preset parameters: \n");
  fprintf(fp,"# Duration of the first cleaning: %f s.\n",FirstCleaningTime / 1000000.);
  fprintf(fp,"# Duration from (RV:Filling, CV:Open) to (RV:Filling, CV:Close) is ...\n");
  for(ii=0;ii<n_filling;ii++){
    fprintf(fp,"#    Case %d:  %f s.\n", ii+1, T_filling[ii]);
  }
  fprintf(fp,"#    Total case number of filling time: %d\n", n_filling);
  fprintf(fp,"# Duration from (RV:Filling, CV:Close) to (RV:Detection, CV:Close): %f s.\n",PreparateRVTime / 1000000.);
  fprintf(fp,"# Duration from (RV:Filling, CV:Close) to (RV:Detection, CV:Open): %f s.\n",StorageTime / 1000000.);
  fprintf(fp,"# Duration from (RV:Detection, CV:Open) to (RV:Filling, CV:Open): %f s.\n",DAQAndCleaningTime / 1000000.);
  fprintf(fp,"#   >> Duration of DAQ: %f s.\n",DAQTime / 1000000.);
  fprintf(fp,"#   >> Duration of cleaning after DAQ: %f s.\n",CleaningTime / 1000000.);
  fprintf(fp,"#\n");
  fprintf(fp,"\n========================== Start valve controll ==========================\n\n");



  int i=0,j=0,run=0;

  time(&timer);
  printf("Cleaning Start at %s",ctime(&timer));
  fprintf(fp,"Cleaning Start @ %s",ctime(&timer));

  system("../GV11");
  time(&timer);
  printf("GV11 (RV:Detection, CV:Open) is performed @ %s",ctime(&timer));
  fprintf(fp,"GV11 (RV:Detection, CV:Open) is performed @ %s",ctime(&timer));
  start=clock();

  printf("\n");
  fprintf(fp,"\n");
  fflush(fp);

  while(1){
    if(i==0){

      //ACCUMULATE SECTION
      //Now is at time limit to release the UCNs confined in storage chamber
      //Next step is to accumulate UCN Flux in storage chamber
      /////////////////////////////////////////////////////(font adjustment)
      if( (run==0 && (clock()-start) >= FirstCleaningTime) || (run!=0&& (clock()-start)>=DAQAndCleaningTime ) ){
        run++;
        time(&timer);
        printf("Run %d Start at %s",run,ctime(&timer));
        fprintf(fp,"Run %d Start @ %s",run,ctime(&timer));
        system("../GV01");
        time(&timer);
        printf("GV01 (RV:Filling, CV:Open) is performed @ %s",ctime(&timer));
        fprintf(fp,"GV01 (RV:Filling, CV:Open) is performed @ %s",ctime(&timer));
        fflush(fp);
        i++;
        start=clock();
      }//(font adjustment)

    }else if(i==1){

      if( (clock()-start) >= FillingTime ){
        //STORAGE SECTION
        //Now is at time limit to accumulate UCNs in stogare chamber
        //Next step is to confine UCNs in storage chamber
        system("../GV00");
        time(&timer);
        printf("GV00 (RV:Filling, CV:Close) is performed @ %s",ctime(&timer));
        fprintf(fp,"GV00 (RV:Filling, CV:Close) is performed @ %s",ctime(&timer));
        fflush(fp);
        i++;
        start=clock();

        printf("Filling time is %d sec.\n", int(FillingTime/1000000.));
        fprintf(fp,"Filling time is %d sec.\n", int(FillingTime/1000000.));
        fflush(fp);

      }//(font adjustment)

    }else if(i==2){

      if( (clock()-start) >= PreparateRVTime ){
        //STORAGE SECTION
        //Now is at time limit to close stogare chamber
        //Next step is to preparate Rotary valve to detection mode
        system("../GV10");
        time(&timer);
        printf("GV10 (RV:Detection, CV:Close) is performed @ %s",ctime(&timer));
        fprintf(fp,"GV10 (RV:Detection, CV:Close) is performed @ %s",ctime(&timer));
        i++;
        start=clock();
        if(StorageTime < PreparateRVTime) StorageTime = PreparateRVTime;

      }//(font adjustment)

    }else if(i==3){

      if( (clock()-start) >= (StorageTime - PreparateRVTime)){
        //RELEASE SECTION
        //Now is at time limit to confine UCNs in storage chamber
        //Next step is to release the UCNs confined in storage chamber
        system("../GV11");
        time(&timer);
        printf("GV11 (RV:Detection, CV:Open) is performed @ %s",ctime(&timer));
        fprintf(fp,"GV11 (RV:Detection, CV:Open) is performed @ %s",ctime(&timer));
        fflush(fp);
        i=0;
        start=clock();

        printf("\n");
        fprintf(fp,"\n");

        if(j==n_filling - 1){//loop control
          j=0;
        }else{
          j++;
        }//(font adjustment)
        FillingTime = T_filling[j] * 1000000.;

      }//(font adjustment)

    }//(font adjustment)

    if( run > runmax)
      break;

  }//end while

  time(&timer);
  printf("Run End at %s", ctime(&timer) );

  system("../GVoff");
  printf("GVoff (RV:Filling, CV:Close) is performed @ %s",ctime(&timer));
  fprintf(fp,"GVoff (RV:Filling, CV:Close) is performed @ %s",ctime(&timer));
  printf("Measurement ends. @ %s",ctime(&timer));
  fprintf(fp,"Measurement ends. @ %s",ctime(&timer));

  fclose(fp);

  return 0;
}
