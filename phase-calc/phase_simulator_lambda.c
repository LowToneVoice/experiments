#include <stdio.h>
#include <math.h>/*M_PIは円周率*/
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/*装置定数定義*/
#define D (50e-6+89e-6) /*ミラーのギャップ*/
#define Lambda 0.30e-9/*入射中性子波長*/
#define Lambda_max 0.80e-9
#define L 150e-3/*ミラー同士の距離*/
#define N 9.1e+28/*原子数密度*/
#define b_c 9.1e-15/*原子半径*/
#define m 1.67e-27/*中性子質量*/
/*物理定数定義*/
#define h 6.626e-34/*プランク定数*/
#define g 9.8/*重力定数*/

int main(void)
{
    /*変数宣言*/
        double theta=0, delta_theta=0,  dtheta=0;/*サンプリング角度、ミラーの角度誤差、サンプリング誤差*/
        double Delta_Psi[21];
        
        double Delta_Psi_g=0,Delta_Psi_a=0, Delta_Psi_a_sub=0, Delta_Psi_g_sub=0;/*位相ずれの変数*/
        double lambda=0;/*波長*/
        double U=0;/*ダミーの変数*/
        int i=0,j=0,rd=0,pm=0,gnuswitch=0;/*処理用の変数*/
        /*ファイル名指定の変数*/
        char filename_gnu[20];
        char filename_lambda[20];
        /*初期化、前処理*/
        for(j=0;j<21;j++){
            Delta_Psi[j]=0;
        }
        U=(h*h)*N*b_c/(2*M_PI*m);
    /*filename setting*/
        printf("ファイル名を指定せよ。\n");
        /*printf("全データのファイル名は\t");scanf("%s",filename_lambda);*/
        printf("波長を変えた全位相差の出力ファイルは\t");scanf("%s",filename_lambda);  
        printf("gnuplot入力用ファイルは\t");scanf("%s",filename_gnu);  
        printf("各ファイル名は,\t%s\t %s\n",filename_lambda,filename_gnu);
        
    /*乱数生成*/
            srand((unsigned) time(NULL));
    /*ファイル呼び出し/作成*/    
        FILE *fout_lambda;
        fout_lambda=fopen(filename_lambda,"w");
        FILE *fout_gnu;
        fout_gnu=fopen(filename_gnu, " w ");
        /*FILE *gp;
        gp=popen("gnuplot", "w");*/

    /*処理*/

        for(i=0;i<=1600;i++){
            /*theta, delta_thetaの処理*/
                j=0;             
                theta = M_PI/180*1e-3+theta;/*角度の離散化*/ 
                Delta_Psi[j]=theta*180/M_PI;
            for(j=0;j<21;j++){
                if(j>0){
                /*randomの処理*/
                    pm=rand()%2;
                    if(pm==0){
                        rd=-rand()%101;
                     }
                    else{
                        rd=rand()%101;
                    }
                    delta_theta = M_PI/180*1e-5*rd;/*誤差は離散間隔の1/100~1/1*/           
                /*lambdaの変更*/
                    lambda=Lambda + (Lambda_max-Lambda)/20*j;
                /*Delta_Psi_gの計算*/
                    Delta_Psi_g = -(2*M_PI*lambda*(m*m)*g)/(h*h)*((2*D*L)/tan(2*theta));
                    Delta_Psi_g_sub = 2*M_PI*lambda*(m*m)*g/(h*h)*(D/sin(theta+delta_theta))*(D/sin(theta+delta_theta))*(sin(delta_theta)/2);  
                /*Delta_Psi_aの計算*/
                    Delta_Psi_a=(4*M_PI*D*delta_theta)/lambda;
                    Delta_Psi_a_sub=-(Delta_Psi_a*m*U*lambda*lambda)/(h*h*theta*theta);
                /*Delta_Psi_allの計算*/
                    Delta_Psi[j]=(Delta_Psi_g+Delta_Psi_g_sub)-(Delta_Psi_a+Delta_Psi_a_sub);
                }
                /*ファイルへの書き込み*/    
                    if(j<20){
                        fprintf(fout_lambda, "%9.7f \t",Delta_Psi[j]);
                    }
                    else{
                        fprintf(fout_lambda,"%9.7f\n",Delta_Psi[j]);
                    }
            }    
        }        
    /*gnuplot設定ファイル*/
                fprintf(fout_gnu, "load \'%s\'\n",filename_lambda);
                fprintf(fout_gnu,"set tics font \'arial,20\'\n");
                fprintf(fout_gnu,"set xlabel \'degree\'\n");
                fprintf(fout_gnu,"set ylabel \'phase\'\n");
                fprintf(fout_gnu,"plot [0:][-1000:100]");
                for(j=2;j<21;j++){
                    fprintf(fout_gnu," \'%s\' u 1:%d,",filename_lambda,j);
                }
                fprintf(fout_gnu," \'%s\' u 1:%d\n",filename_lambda,j);
    /*ファイルの閉封/データ書き込み*/
        fclose(fout_lambda); 
        fclose(fout_gnu);
    /*gnuplotへの書き込み*/
        /*fprintf(gp, "load \'%s\'\n",filename_all);
        fprintf(gp,"set tics font \'arial,20\'\n");
        fprintf(gp,"plot [0:][-1000:100] \'%s\' u 1:3, \'%s\' u 1:4, \'%s\' u 1:5,   \'%s\' u 1:6\n",filename_all,filename_all,filename_all,filename_all);
        fprintf(gp,"set terminal png size 400,300 enhanced font \"Helvetica,20\"\n");
        fprintf(gp, "set output \'output.png\'\n");
        pclose(gp);*/
    
    
    return 0;
}

