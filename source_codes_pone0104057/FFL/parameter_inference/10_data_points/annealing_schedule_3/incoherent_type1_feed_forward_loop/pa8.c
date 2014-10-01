#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <bits/nan.h>
#include "mt19937ar.h"

#define EQ_N 2      /* タンパク質の種類 = 微分方程式数 */
#define PARA_N 8     /* パラメーター (速度定数) の数 */
#define PERTURB_PARA_N 7   /* 摂動を加えるパラメータ数 */
#define PARTICLE_N 10000 /* 粒子の数 */
#define POPULATION_N 5   /* ポピュレーションの数 (世代数) */
#define MCMC_MAX_STEP PERTURB_PARA_N*1 /* PA　アルゴリズム内の MCMC の最大ステップ数 */
#define STEP_PARA 0.5
#define RK_N 1000   /* Runge-Kutta法の最大反復回数 */
#define TIME_POINT_N 10


void function_set(double (*func[EQ_N])());
void parameter_set(double *p);

double fY(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H);
double fZ(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H);

void rk_search(double X, double dt, double (*func[EQ_N])(), double y[][EQ_N],  double p[], double time_point[]);

void parameter_sampler(double parameter[]);

double pdf_uniform_parameter(double parameter[]);

void parameter_generation(double parameter[]);

void box_muller(double mean, double variance, double *nr);

double indicator_function(double time_point[], double ref_time_point[], double epsilon);
double likelihood_function(double time_point[], double ref_time_point[], double variance);




/* MAIN 関数 */

int main(void)
{
  double y[EQ_N][2];
  double (*func[EQ_N])();
  double parameter[PARTICLE_N][PARA_N], parameter_memory[PARA_N], resampled_parameter[PARTICLE_N][PARA_N];
  double time_point[TIME_POINT_N];
  double ref_time_point[TIME_POINT_N];

  int i, j, k, l, m;
  int mcmc_num, pa_num, integral_num, particle_num, rk_num, timepoint_num, parameter_num, equation_num;
  
  int x_flag = 0;
  double X, dt;

  double probability_before, probability_after, probability_ratio;
  int move_count = 0;

  double weight[PARTICLE_N], weight_tmp[PARTICLE_N];
  double total_weight = 0, power_total_weight, partial_total_weight, upper_weight, lower_weight;
  double log_weight[PARTICLE_N], log_weight_tmp[PARTICLE_N];
  double total_likelihood[PARTICLE_N], total_likelihood_previous[PARTICLE_N];
  double log_likelihood[PARTICLE_N], total_log_likelihood;
  double epsilon[POPULATION_N] = {1, 0.5, 0.25, 0.15, 0.11};
  double ess;
  int sampling_num;
  int resampling_count = 0;

  int mcmc_step;

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
  init_genrand(8);

  /* ファイル設定 */
  fp1 = fopen("information8.data","w");
  fp2 = fopen("pa_parameter_set8.data","w");
  fp3 = fopen("weight8.data", "w");

  /* 実験データの読み込み */
  fin = fopen("ref_time_point.data","r");

  for(timepoint_num=0;timepoint_num<TIME_POINT_N;timepoint_num++){
    fscanf(fin, "%lf", &ref_time_point[timepoint_num]);
  }

  fclose(fin);
  /* 実験データの読み込み終了 */



  /* サンプル番号初期化 */
  sampling_num = 0;

  printf("0 th population\n");

  /* ランダムサンプリングによる 粒子数個のパラメータセットの発生 **************************************************************************************/
  for( ; ; ){

    /* パラメータ、合計濃度、変動係数を発生させる */
    parameter_set(parameter[sampling_num]);
    parameter_generation(parameter[sampling_num]);

    sampling_num++;
    
    if(sampling_num == PARTICLE_N) break; 
  }


  /* 重みを初期化 正規化 */
  for(particle_num=0;particle_num<PARTICLE_N;particle_num++){
    weight[particle_num] = 1.0/PARTICLE_N;
  }


  /* ランダムサンプリングによる初期パラメータ発生 終了 ************************************************************************************************/









  
  /* Population annealing 開始 *****************************************************************************************************************/
  for(pa_num=0;pa_num<POPULATION_N;pa_num++){

    printf("%d th population\n", pa_num+1);
    
    /* 各粒子ごとに尤度の計算 */
    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){
      total_likelihood[particle_num] = 0; /* 初期化 */
      total_likelihood_previous[particle_num] = 0; /* 初期化 */
    }


    /* 中間分布の計算 */
    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){

      /* time point データの発生 */
      y[0][0] = 0.0;
      y[1][0] = 0.0;
      rk_search(X, dt, func, y, parameter[particle_num], time_point);

      /* 尤度の計算 */
      total_likelihood[particle_num] = indicator_function(time_point, ref_time_point, epsilon[pa_num]);

      if(pa_num==0){
        total_likelihood_previous[particle_num] = 1;
      }
      else{
        total_likelihood_previous[particle_num] = indicator_function(time_point, ref_time_point, epsilon[pa_num-1]);
      }

    }
    /* 中間分布の計算終了 */


    /* 重みの計算 */
    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){
      weight[particle_num] = weight[particle_num] * total_likelihood[particle_num]/total_likelihood_previous[particle_num];
      if(total_likelihood[particle_num]==0 || total_likelihood_previous[particle_num]==0) weight[particle_num] = 0.0;
    }

    total_weight = 0.0;

    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){
      total_weight = total_weight + weight[particle_num];
    }
    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){
      weight[particle_num] = weight[particle_num]/total_weight;
    }






    /* 各粒子ごとに */
    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){

      /* 重みが0なら飛ばす */
      if(weight[particle_num]==0) continue;

      /* フラグの初期化 */
      x_flag = 0;


      /* MCMC 計算 */
      for(mcmc_num=0;mcmc_num<MCMC_MAX_STEP;mcmc_num++){
 
        /* 許容されなければ、元のままのパラメータ */
        if(x_flag == 1){
          for(parameter_num=0;parameter_num<PARA_N;parameter_num++){
            parameter[particle_num][parameter_num] = parameter_memory[parameter_num];
          }
        }


        /* 確率の初期化 */
        probability_before = 1.0;
        probability_after = 1.0;
    

        /* フラグの初期化 */
        x_flag = 0;


        /* 新候補の計算前の値を記憶 */
        for(parameter_num=0;parameter_num<PARA_N;parameter_num++){
          parameter_memory[parameter_num] = parameter[particle_num][parameter_num];
        }


        /* 移動前の確率を計算 */
        probability_before = probability_before*pdf_uniform_parameter(parameter[particle_num]); /* 事前分布を乗じる */

        
        total_likelihood[particle_num] = 0.0; /* 初期化 */

        /* time point データの発生 */
        y[0][0] = 0.0;
        y[1][0] = 0.0;
        rk_search(X, dt, func, y, parameter[particle_num], time_point);

        /* 尤度の計算 */
        total_likelihood[particle_num] = indicator_function(time_point, ref_time_point, epsilon[pa_num]);

        probability_before = probability_before * total_likelihood[particle_num]; /* 移動前の確率 */  


        /* 新候補を計算 パラメータの中から1つだけ動かす */
        unit_r = genrand_real2();
        chosen_num = PERTURB_PARA_N*unit_r;
        parameter_sampler(parameter[particle_num]);


        /* 移動後の確率を計算 : 一様事前分布を乗じる */
        probability_after = probability_after*pdf_uniform_parameter(parameter[particle_num]); /* 事前分布を乗じる */
        

        if(probability_after == 0.0){
          x_flag = 1;
          for(parameter_num=0;parameter_num<PARA_N;parameter_num++){
            parameter[particle_num][parameter_num] = parameter_memory[parameter_num];
          }
          continue;
        }



        /* 移動後の確率を計算 */
        total_likelihood[particle_num] = 0.0; /* 初期化 */

        /* time point データの発生 */
        y[0][0] = 0.0;
        y[1][0] = 0.0;
        rk_search(X, dt, func, y, parameter[particle_num], time_point);

        /* 尤度の計算 */
        total_likelihood[particle_num] = indicator_function(time_point, ref_time_point, epsilon[pa_num]);

        probability_after = probability_after * total_likelihood[particle_num]; /* 移動後の確率 */


        if(probability_after == 0.0){
          x_flag = 1;
          for(parameter_num=0;parameter_num<PARA_N;parameter_num++){
            parameter[particle_num][parameter_num] = parameter_memory[parameter_num];
          }
          continue;
        }



        /* 移動前と移動後の確率の比の計算 */
        probability_ratio = probability_after/probability_before;

        /* [0,1) 単位乱数の発生 */
        unit_r = genrand_real2();   


        /* 移動の判定 */
        if(probability_ratio > unit_r){ 
          move_count++;
        }
        else{
          x_flag = 1;
          for(parameter_num=0;parameter_num<PARA_N;parameter_num++){
            parameter[particle_num][parameter_num] = parameter_memory[parameter_num];
          }
          continue;
        }
    
      }
      /* particle_num 番目の粒子の MCMC 終了 */
    }
    /* 全粒子の MCMC 終了 */
 




 
    /* リサンプリング ********************************************************/
    power_total_weight = 0.0;
    ess = 0.0;

    for(particle_num=0;particle_num<PARTICLE_N;particle_num++){
      power_total_weight = power_total_weight + pow(weight[particle_num], 2);
    }


    /* Effective Sample Size の計算 */
    ess = 1.0/power_total_weight;
    //if(pa_num==4) ess = 0.0;

    fprintf(fp1, "%d th population", pa_num);
    fprintf(fp1, "ESS = %f\t particle_num/2 = %f\n\n", ess, PARTICLE_N/2.0);
    

    /* リサンプリングの判定と実行 */
    if(ess < PARTICLE_N/2){

      fprintf(fp1, "resampling\n");
      resampling_count++;
      
      /* 粒子数個リサンプリングする 重複可 */
      for(m=0;m<PARTICLE_N;m++){

        /* [0,1) 乱数発生 */
        unit_r = genrand_real2();

        partial_total_weight = 0.0;

        /* l 番目の粒子が選ばれる */
        for(l=0;l<PARTICLE_N;l++){
          partial_total_weight = partial_total_weight + weight[l];

          upper_weight = partial_total_weight;
          lower_weight = partial_total_weight - weight[l];

          if((unit_r >= lower_weight) && (unit_r < upper_weight)){
            for(i=0;i<PARA_N;i++){
              resampled_parameter[m][i] = parameter[l][i];
            }
            break;
          }
        }

      }
      /* 粒子数個リサンプリング終了 */



      /* リサンプル後のパラメータと活性化時間の再設定 */
      for(l=0;l<PARTICLE_N;l++){
        for(i=0;i<PARA_N;i++){
          parameter[l][i] = resampled_parameter[l][i];
        }
      }



      /* リサンプル後の重みを初期化 正規化 */
      for(l=0;l<PARTICLE_N;l++){
        weight[l] = 1.0/PARTICLE_N;
      }


    }
    /* リサンプリングの判定と実行 終了 */
    /* リサンプリング終了 ****************************************************/

  }
  /* Population annealing 終了 *****************************************************************************************************************/
  

  /* ファイルへ出力 */
  for(particle_num=0;particle_num<PARTICLE_N;particle_num++){

    /* パラメータ */
    for(parameter_num=0;parameter_num<PARA_N;parameter_num++){
      fprintf(fp2, "%f\t", parameter[particle_num][parameter_num]);
    }
    fprintf(fp2, "\n");

    /* 重み */
    fprintf(fp3, "%12.10f\n", weight[particle_num]);

  }


  /* 移動回数、 リサンプリング回数、 モデル数を出力 */
  fprintf(fp1, "%d times move\n", move_count);
  fprintf(fp1, "%d times resampling\n", resampling_count);
  end = clock();
  /* 計算時間の出力 */
  fprintf(fp1, "%f min\n", (double)(end - start)/CLOCKS_PER_SEC/60.0);



  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

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




/* 事前分布 */
double pdf_uniform_parameter(double parameter[])
{
  int probability_0_flag = 0;

  if(parameter[0] < 0.1 || parameter[0] > 10.0) probability_0_flag = 1;
  if(parameter[1] < 0.1 || parameter[1] > 10.0) probability_0_flag = 1;
  if(parameter[2] < 0.1 || parameter[2] > 10.0) probability_0_flag = 1;
  if(parameter[3] < 0.1 || parameter[3] > 10.0) probability_0_flag = 1;
  if(parameter[4] < 0.01 || parameter[4] > 1.0) probability_0_flag = 1;
  if(parameter[5] < 0.01 || parameter[5] > 1.0) probability_0_flag = 1;
  if(parameter[6] < 0.05 || parameter[6] > 5.0) probability_0_flag = 1;

  if(probability_0_flag == 1){
    return 0.0;
  }
  else{
    return 1.0;
  }
 
}





/* パラメーターに摂動を加える */

void parameter_sampler(double parameter[])
{
  double unit_r;
  int chosen_num, i;

  for( ; ; ){
    unit_r = genrand_real2();
    chosen_num = (int)(PARA_N*unit_r);
  
    if(chosen_num < 7){
      unit_r = genrand_real1();
      parameter[chosen_num] = parameter[chosen_num]*pow(10.0, (unit_r-0.5)*STEP_PARA);
      break;
    }
    else{
      continue;
    }
  }
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
  }
  
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





