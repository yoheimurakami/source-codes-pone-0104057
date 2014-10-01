#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <bits/nan.h>
#include "function.h"
#include "mt19937ar.h"

#define EQ_N 6      /* タンパク質の種類 = 微分方程式数 */
#define PARA_N 16     /* パラメーター (速度定数) の数 */
#define PERTURB_PARA_N 8   /* 摂動を加えるパラメータ数 */
#define PARTICLE_N 10000 /* 粒子の数 */
#define RK_N 360001   /* Runge-Kutta法の最大反復回数 */
#define TIME_POINT_N 7


void initial_condition_set(double y[][2]);
void parameter_set(double p[]);
void function_set(double (*function[])());

void rk_search(double Ins, double dt, double (*func[EQ_N])(), double y[][2],  double p[], double time_point[]);

void parameter_generation(double parameter[], double ref_parameter[]);

double indicator_function(double time_point[], double ref_time_point[], double epsilon);




/* MAIN 関数 */

int main(void)
{
  double y[EQ_N][2];
  double (*func[EQ_N])();
  double parameter[PARA_N];
  double ref_parameter[PARA_N];
  double time_point[TIME_POINT_N];
  double ref_time_point[TIME_POINT_N];
  double epsilon = 0.01;
  int accept_particle_number = 0;

  int i, j, k, l, m;
  int mcmc_num, pa_num, integral_num, particle_num, rk_num, timepoint_num, parameter_num, equation_num;
  
  int x_flag = 0;
  double dt;

  int sampling_num;

  double unit_r, chosen_num;

  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fin;

  clock_t start, end;

  dt = 0.001;


  function_set(func);


/*** 計算開始 **********************************************************************************************************************************/

  /* 時間を計る */
  start = clock();
  
  /* 乱数の初期化 */
  init_genrand(8);

  /* ファイル設定 */
  fp1 = fopen("information8.data","w");

  /* 実験データの読み込み */
  fin = fopen("ref_time_point.data","r");

  for(timepoint_num=0;timepoint_num<TIME_POINT_N;timepoint_num++){
    fscanf(fin, "%lf", &ref_time_point[timepoint_num]);
  }

  fclose(fin);
  /* 実験データの読み込み終了 */


  /* 参照パラメータ値設定 */
  parameter_set(ref_parameter);

  /* ランダムサンプリングによる 粒子数個のパラメータセットの発生 **************************************************************************************/
  for(sampling_num=0;sampling_num<PARTICLE_N;sampling_num++){
    
    /* パラメータを発生させる */
    parameter_set(parameter);
    parameter_generation(parameter, ref_parameter);

    /* time point データの発生 */
    initial_condition_set(y);
    rk_search(1, dt, func, y, parameter, time_point);

    if(indicator_function(time_point, ref_time_point, epsilon)==1.0){
      accept_particle_number++;
    }
  }
  /* ランダムサンプリングによる初期パラメータ発生 終了 ************************************************************************************************/


  fprintf(fp1, "accept particle number = %d\n", accept_particle_number);

  end = clock();
  /* 計算時間の出力 */
  fprintf(fp1, "%f min\n", (double)(end - start)/CLOCKS_PER_SEC/60.0);

  fclose(fp1);

  return 0;
}














/* パラメーターを発生させる */

void parameter_generation(double parameter[], double ref_parameter[])
{
  double unit_r, exponent;
  int i;

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[0] = 1.0e-6*pow(10.0, exponent);

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[1] = 1.0e-6*pow(10.0, exponent);

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[2] = 1.0e-6*pow(10.0, exponent);

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[3] = 1.0e-6*pow(10.0, exponent);


  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[12] = 1.0e-6*pow(10.0, exponent);

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[13] = 1.0e-6*pow(10.0, exponent);

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[14] = 1.0e-6*pow(10.0, exponent);

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,10] の乱数に変換 */
  exponent = 10.0*unit_r;
  /* ランダムな値をログスケールで発生させる */
  parameter[15] = 1.0e-6*pow(10.0, exponent);
}





void rk_search(double Ins, double dt, double (*func[EQ_N])(), double y[][2],  double p[], double time_point[])
{
  int i, j;
  double k1[EQ_N], k2[EQ_N], k3[EQ_N], k4[EQ_N];
  double pro_IRcom, iAKT, imTOR;

  pro_IRcom = 8.8877e+2;
  iAKT = 1.1182;
  imTOR = 8.8025e+2;

  /* ルンゲクッタ法 */
  for(i=0;i<=RK_N;i++){

    /* 時間発展データの記録 */
    if(i==0) time_point[0] = y[4][0];
    if(i==10000) time_point[1] = y[4][0];
    if(i==30000) time_point[2] = y[4][0];
    if(i==60000) time_point[3] = y[4][0];
    if(i==120000) time_point[4] = y[4][0];
    if(i==240000) time_point[5] = y[4][0];
    if(i==360000) time_point[6] = y[4][0];

    for(j=0;j<EQ_N;j++){
      k1[j] = dt*(*func[j])(Ins, pro_IRcom, y[0][0], y[1][0], y[2][0], y[3][0], iAKT, y[4][0], imTOR, y[5][0], \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }
    for(j=0;j<EQ_N;j++){
      k2[j] = dt*(*func[j])(Ins, pro_IRcom, y[0][0] + k1[0]/2.0, y[1][0] + k1[1]/2.0, y[2][0] + k1[2]/2.0, \
				y[3][0] + k1[3]/2.0, iAKT, y[4][0] + k1[4]/2.0, imTOR, y[5][0] + k1[5]/2.0, \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }
    for(j=0;j<EQ_N;j++){
      k3[j] = dt*(*func[j])(Ins, pro_IRcom, y[0][0] + k2[0]/2.0, y[1][0] + k2[1]/2.0, y[2][0] + k2[2]/2.0, \
				y[3][0] + k2[3]/2.0, iAKT, y[4][0] + k2[4]/2.0, imTOR, y[5][0] + k2[5]/2.0, \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }
    for(j=0;j<EQ_N;j++){
      k4[j] = dt*(*func[j])(Ins, pro_IRcom, y[0][0] + k3[0], y[1][0] + k3[1], y[2][0] + k3[2], \
				y[3][0] + k3[3], iAKT, y[4][0] + k3[4], imTOR, y[5][0] + k3[5], \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }

    for(j=0;j<EQ_N;j++){
      y[j][1] = y[j][0] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0;
    }

    for(j=0;j<EQ_N;j++){
      y[j][0] = y[j][1];
    }
  }

}





double indicator_function(double time_point[], double ref_time_point[], double epsilon)
{
  int i;
  double sse = 0.0;
 
  for(i=0;i<TIME_POINT_N;i++){
    sse = sse + pow(time_point[i]-ref_time_point[i], 2);
  }

  if(sse <= epsilon){
    return 1.0;
  }
  else{
    return 0.0;
  }
  
}





/* パラメーターの設定 */

void parameter_set(double p[])
{
  p[0] = 5.5702; /* k1_synthesis */
  p[1] = 2.3969; /* k1_InsIRcom */
  p[2] = 2.0934e-2; /* k2_InsIRcom */
  p[3] = 1.2140e-5; /* k1_p1IRcomDeg */
  p[4] = 0; /* k1_p1IRcomPhos */
  p[5] = 0; /* k1_p1p2IRcomdePhos */
  p[6] = 0; /* k1_IRcomPhos */
  p[7] = 0; /* k1_p2IRcomdePhos */
  p[8] = 0; /* k1_p2IRcomDeg */
  p[9] = 0; /* k1_Insp2IRcom */
  p[10] = 0; /* k2_Insp2IRcom */
  p[11] = 0; /* k1_p1p2IRcomDeg */
  p[12] = 4.4307e-5; /* k1_AKTPhos */
  p[13] = 3.5211e-1; /* k1_pAKTdePhos */
  p[14] = 3.1385e-5; /* k1_mTORPhos */
  p[15] = 9.5934e-3; /* k1_pmTORdePhos */
}




/* 関数ポインタの設定 */

void function_set(double (*function[])())
{
  function[0] = fIRcom;
  function[1] = fp1IRcom;
  function[2] = fp2IRcom;
  function[3] = fp1p2IRcom;
  function[4] = fpAKT;
  function[5] = fpmTOR;
}





/* 初期量の設定 */

void initial_condition_set(double y[][2])
{
  y[0][0] = 8.8877e+2; /* IRcom */
  y[1][0] = 0; /* p1IRcom */
  y[2][0] = 0; /* p2IRcom */
  y[3][0] = 0; /* p1p2IRcom */
  y[4][0] = 0; /* pAKT */
  y[5][0] = 0; /* pmTOR */
}
  
