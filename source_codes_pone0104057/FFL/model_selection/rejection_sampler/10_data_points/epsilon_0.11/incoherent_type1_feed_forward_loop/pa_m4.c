#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <bits/nan.h>
#include "mt19937ar.h"

#define EQ_N 2      /* タンパク質の種類 = 微分方程式数 */
#define PARA_N 8     /* パラメーター (速度定数) の数 */
#define INTEGRAL_N 100000 /* 粒子の数 */
#define RK_N 1000   /* Runge-Kutta法の最大反復回数 */
#define TIME_POINT_N 10


void function_set(double (*func[EQ_N])());
void parameter_set(double *p);

double fY(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H);
double fZ(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H);

void rk_search(double X, double dt, double (*func[EQ_N])(), double y[][EQ_N],  double p[], double time_point[]);


void parameter_generation(double parameter[]);

void box_muller(double mean, double variance, double *nr);

double indicator_function(double time_point[], double ref_time_point[], double epsilon);
double likelihood_function(double time_point[], double ref_time_point[], double variance);






/* MAIN 関数 */

int main(void)
{
  double y[EQ_N][2];
  double (*func[EQ_N])();
  double parameter[PARA_N];
  double time_point[TIME_POINT_N];
  double ref_time_point[TIME_POINT_N];

  int i, j, k, l, m;
  int mcmc_num, pa_num, integral_num, particle_num, rk_num, timepoint_num, parameter_num, equation_num;
  
  int x_flag = 0;
  double X, dt;

  double total = 0.0;
  double epsilon = 0.11;
  int sampling_num;

  double unit_r, chosen_num;

  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fin;

  clock_t start, end;

  dt = 0.01;
  X = 1.0;


  function_set(func);


/*** 計算開始 **********************************************************************************************************************************/

  /* 時間を計る */
  start = clock();
  
  /* 乱数の初期化 */
  init_genrand(4);

  /* ファイル設定 */
  fp1 = fopen("information4.data","w");

  /* 実験データの読み込み */
  fin = fopen("ref_time_point.data","r");

  for(timepoint_num=0;timepoint_num<TIME_POINT_N;timepoint_num++){
    fscanf(fin, "%lf", &ref_time_point[timepoint_num]);
  }

  fclose(fin);
  /* 実験データの読み込み終了 */




  /* ランダムサンプリングによるパラメータセットの発生 **************************************************************************************/
  for(sampling_num=0;sampling_num<INTEGRAL_N;sampling_num++){

    /* パラメータ、合計濃度、変動係数を発生させる */
    parameter_set(parameter);
    parameter_generation(parameter);

    /* time point データの発生 */
    y[0][0] = 0.0;
    y[1][0] = 0.0;
    rk_search(X, dt, func, y, parameter, time_point);

    /* 尤度の計算 */
    total = total + indicator_function(time_point, ref_time_point, epsilon);

  }


  /* ランダムサンプリングによるパラメータ発生 終了 ************************************************************************************************/

  fprintf(fp1, "%12.10f\n", total/INTEGRAL_N);
  fprintf(fp1, "%d particles\n", INTEGRAL_N);

  end = clock();
  /* 計算時間の出力 */
  fprintf(fp1, "%f min\n", (double)(end - start)/CLOCKS_PER_SEC/60.0);



  fclose(fp1);


  return 0;
}














/* Functions */

/* パラメーターを発生させる */

void parameter_generation(double parameter[])
{
  double unit_r, exponent;

  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r; 
  /* 0.1 ~ 10.0 範囲のランダムな値をログスケールで発生させる */
  parameter[0] = 0.1*pow(10.0, exponent);


  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r; 
  /* 0.1 ~ 10.0 範囲のランダムな値をログスケールで発生させる */
  parameter[1] = 0.1*pow(10.0, exponent);


  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r; 
  /* 0.1 ~ 10.0 範囲のランダムな値をログスケールで発生させる */
  parameter[2] = 0.1*pow(10.0, exponent);


  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r;  
  /* 0.1 ~ 10.0 範囲のランダムな値をログスケールで発生させる */
  parameter[3] = 0.1*pow(10.0, exponent);

  
  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r;  
  /* 0.01 ~ 1.0 範囲のランダムな値をログスケールで発生させる */
  parameter[4] = 0.01*pow(10.0, exponent);


  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r;  
  /* 0.01 ~ 1.0 範囲のランダムな値をログスケールで発生させる */
  parameter[5] = 0.01*pow(10.0, exponent);


  /* [0,1] 単位乱数を用意する */
  unit_r = genrand_real1();
  /* [0,2] の乱数に変換 */
  exponent = 2.0*unit_r;  
  /* 0.05 ~ 5.0 範囲のランダムな値をログスケールで発生させる */
  parameter[6] = 0.05*pow(10.0, exponent);
}




/* box-muller 法で　正規乱数を得る */

void box_muller(double mean, double variance, double *nr)
{
  double r1, r2;

  /* [0,1] 単位乱数を2つ用意する */
  r1 = genrand_real1();
  r2 = genrand_real1();

  *nr = sqrt(variance) * sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2) + mean;
}








void rk_search(double X, double dt, double (*func[EQ_N])(), double y[][EQ_N],  double p[], double time_point[])
{
  int i, j;
  double k1[EQ_N], k2[EQ_N], k3[EQ_N], k4[EQ_N];

  /* ルンゲクッタ法 */
  for(i=0;i<=RK_N;i++){

    /* 時間発展データの記録 */
    if(i==100) time_point[0] = y[1][0];
    if(i==200) time_point[1] = y[1][0];
    if(i==300) time_point[2] = y[1][0];
    if(i==400) time_point[3] = y[1][0];
    if(i==500) time_point[4] = y[1][0];
    if(i==600) time_point[5] = y[1][0];
    if(i==700) time_point[6] = y[1][0];
    if(i==800) time_point[7] = y[1][0];
    if(i==900) time_point[8] = y[1][0];
    if(i==1000) time_point[9] = y[1][0];

    for(j=0;j<EQ_N;j++){
      k1[j] = dt*(*func[j])(X, y[0][0], y[1][0], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
    }
    for(j=0;j<EQ_N;j++){
      k2[j] = dt*(*func[j])(X, y[0][0] + k1[0]/2.0, y[1][0] + k1[1]/2.0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
    }
    for(j=0;j<EQ_N;j++){
      k3[j] = dt*(*func[j])(X, y[0][0] + k2[0]/2.0, y[1][0] + k2[1]/2.0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
    }
    for(j=0;j<EQ_N;j++){
      k4[j] = dt*(*func[j])(X, y[0][0] + k3[0], y[1][0] + k3[1], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
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



double likelihood_function(double time_point[], double ref_time_point[], double variance)
{
  int i;
  double product = 1.0;
  
  for(i=0;i<TIME_POINT_N;i++){
    product = product * 1/sqrt(2.0*M_PI*variance) * exp(-pow((time_point[i] - ref_time_point[i]), 2)/2/variance);
    //printf("%f\n", 1/sqrt(2.0*M_PI*variance) * exp(-pow((time_point[i] - ref_time_point[i]), 2)/2/variance));
  }

  //printf("product = %f\n", product);
  
  return product;
}
 





/* パラメーターの設定 */

void parameter_set(double *p)
{
  *(p+0) = 1.0; /* alpha_y */
  *(p+1) = 1.0; /* alpha_z */
  *(p+2) = 1.0; /* beta_y */
  *(p+3) = 1.0; /* beta_z */
  *(p+4) = 0.1; /* K_xy */
  *(p+5) = 0.1; /* K_xz */
  *(p+6) = 0.5; /* K_yz */
  *(p+7) = 2; /* H */
}







/* 関数ポインタの設定 */

void function_set(double (*func[EQ_N])())
{
  func[0] = fY;
  func[1] = fZ;
}




double fY(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H)
{
  return beta_y*pow(X/K_xy, H)/(1+pow(X/K_xy, H)) - alpha_y*Y;
}

double fZ(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H)
{
  return beta_z*pow(X/K_xz, H)/(1+pow(X/K_xz, H))*1/(1+pow(Y/K_yz, H)) - alpha_z*Z;
}





