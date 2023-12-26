#include <stdio.h>
#include <math.h>/*M_PI�͉~����*/
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/*���u�萔��`*/
#define D (50e-6+89e-6) /*�~���[�̃M���b�v*/
#define Lambda 0.30e-9/*���˒����q�g��*/
#define Lambda_max 0.80e-9
#define L 150e-3/*�~���[���m�̋���*/
#define N 9.1e+28/*���q�����x*/
#define b_c 9.1e-15/*���q���a*/
#define m 1.67e-27/*�����q����*/
/*�����萔��`*/
#define h 6.626e-34/*�v�����N�萔*/
#define g 9.8/*�d�͒萔*/

int main(void)
{
    /*�ϐ��錾*/
        double theta=0, delta_theta=0,  dtheta=0;/*�T���v�����O�p�x�A�~���[�̊p�x�덷�A�T���v�����O�덷*/
        double Delta_Psi[21];
        
        double Delta_Psi_g=0,Delta_Psi_a=0, Delta_Psi_a_sub=0, Delta_Psi_g_sub=0;/*�ʑ�����̕ϐ�*/
        double lambda=0;/*�g��*/
        double U=0;/*�_�~�[�̕ϐ�*/
        int i=0,j=0,rd=0,pm=0,gnuswitch=0;/*�����p�̕ϐ�*/
        /*�t�@�C�����w��̕ϐ�*/
        char filename_gnu[20];
        char filename_lambda[20];
        /*�������A�O����*/
        for(j=0;j<21;j++){
            Delta_Psi[j]=0;
        }
        U=(h*h)*N*b_c/(2*M_PI*m);
    /*filename setting*/
        printf("�t�@�C�������w�肹��B\n");
        /*printf("�S�f�[�^�̃t�@�C������\t");scanf("%s",filename_lambda);*/
        printf("�g����ς����S�ʑ����̏o�̓t�@�C����\t");scanf("%s",filename_lambda);  
        printf("gnuplot���͗p�t�@�C����\t");scanf("%s",filename_gnu);  
        printf("�e�t�@�C������,\t%s\t %s\n",filename_lambda,filename_gnu);
        
    /*��������*/
            srand((unsigned) time(NULL));
    /*�t�@�C���Ăяo��/�쐬*/    
        FILE *fout_lambda;
        fout_lambda=fopen(filename_lambda,"w");
        FILE *fout_gnu;
        fout_gnu=fopen(filename_gnu, " w ");
        /*FILE *gp;
        gp=popen("gnuplot", "w");*/

    /*����*/

        for(i=0;i<=1600;i++){
            /*theta, delta_theta�̏���*/
                j=0;             
                theta = M_PI/180*1e-3+theta;/*�p�x�̗��U��*/ 
                Delta_Psi[j]=theta*180/M_PI;
            for(j=0;j<21;j++){
                if(j>0){
                /*random�̏���*/
                    pm=rand()%2;
                    if(pm==0){
                        rd=-rand()%101;
                     }
                    else{
                        rd=rand()%101;
                    }
                    delta_theta = M_PI/180*1e-5*rd;/*�덷�͗��U�Ԋu��1/100~1/1*/           
                /*lambda�̕ύX*/
                    lambda=Lambda + (Lambda_max-Lambda)/20*j;
                /*Delta_Psi_g�̌v�Z*/
                    Delta_Psi_g = -(2*M_PI*lambda*(m*m)*g)/(h*h)*((2*D*L)/tan(2*theta));
                    Delta_Psi_g_sub = 2*M_PI*lambda*(m*m)*g/(h*h)*(D/sin(theta+delta_theta))*(D/sin(theta+delta_theta))*(sin(delta_theta)/2);  
                /*Delta_Psi_a�̌v�Z*/
                    Delta_Psi_a=(4*M_PI*D*delta_theta)/lambda;
                    Delta_Psi_a_sub=-(Delta_Psi_a*m*U*lambda*lambda)/(h*h*theta*theta);
                /*Delta_Psi_all�̌v�Z*/
                    Delta_Psi[j]=(Delta_Psi_g+Delta_Psi_g_sub)-(Delta_Psi_a+Delta_Psi_a_sub);
                }
                /*�t�@�C���ւ̏�������*/    
                    if(j<20){
                        fprintf(fout_lambda, "%9.7f \t",Delta_Psi[j]);
                    }
                    else{
                        fprintf(fout_lambda,"%9.7f\n",Delta_Psi[j]);
                    }
            }    
        }        
    /*gnuplot�ݒ�t�@�C��*/
                fprintf(fout_gnu, "load \'%s\'\n",filename_lambda);
                fprintf(fout_gnu,"set tics font \'arial,20\'\n");
                fprintf(fout_gnu,"set xlabel \'degree\'\n");
                fprintf(fout_gnu,"set ylabel \'phase\'\n");
                fprintf(fout_gnu,"plot [0:][-1000:100]");
                for(j=2;j<21;j++){
                    fprintf(fout_gnu," \'%s\' u 1:%d,",filename_lambda,j);
                }
                fprintf(fout_gnu," \'%s\' u 1:%d\n",filename_lambda,j);
    /*�t�@�C���̕�/�f�[�^��������*/
        fclose(fout_lambda); 
        fclose(fout_gnu);
    /*gnuplot�ւ̏�������*/
        /*fprintf(gp, "load \'%s\'\n",filename_all);
        fprintf(gp,"set tics font \'arial,20\'\n");
        fprintf(gp,"plot [0:][-1000:100] \'%s\' u 1:3, \'%s\' u 1:4, \'%s\' u 1:5,   \'%s\' u 1:6\n",filename_all,filename_all,filename_all,filename_all);
        fprintf(gp,"set terminal png size 400,300 enhanced font \"Helvetica,20\"\n");
        fprintf(gp, "set output \'output.png\'\n");
        pclose(gp);*/
    
    
    return 0;
}

