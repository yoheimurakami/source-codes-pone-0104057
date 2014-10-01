#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "function.h"

#define VARIABLE_N 6 /* 変数の数 */
#define PARAMETER_N 16 /* パラメータの数 */
#define STEP_N 500000 /* ODEステップ数 */

void initial_condition_set(double y[][2]);
void parameter_set(double p[]);
void function_set(double (*function[])());
void runge_kutta(double y[][2], double (*function[])(), double p[], double dt);



/* MAIN 関数 */

int main(void)
{
  double y[VARIABLE_N][2];
  double parameter[PARAMETER_N];
  double (*function[VARIABLE_N])();
  double dt;
  
  dt = 0.001;

  parameter_set(parameter);
  function_set(function);
  initial_condition_set(y);
  
  runge_kutta(y, function, parameter, dt);
  
  return 0;
}





/* 4th order Runge-Kutta method */

void runge_kutta(double y[][2], double (*function[])(), double p[], double dt)
{
  int i, j;
  double k1[VARIABLE_N], k2[VARIABLE_N], k3[VARIABLE_N], k4[VARIABLE_N];
  FILE *fp;
  double Ins;
  double pro_IRcom, iAKT, imTOR;

  pro_IRcom = 8.8877e+2;
  iAKT = 1.1182;
  imTOR = 8.8025e+2;
  
  fp = fopen("runge_kutta.data","w");

  Ins = 1;

  for(i=0;i<=STEP_N;i++){

    if(i%10 == 0) fprintf(fp,"%f\t %f\n", dt*i, y[4][0]);

    for(j=0;j<VARIABLE_N;j++){
      k1[j] = dt*(*function[j])(Ins, pro_IRcom, y[0][0], y[1][0], y[2][0], y[3][0], iAKT, y[4][0], imTOR, y[5][0], \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }
    
    
    for(j=0;j<VARIABLE_N;j++){
      k2[j] = dt*(*function[j])(Ins, pro_IRcom, y[0][0] + k1[0]/2.0, y[1][0] + k1[1]/2.0, y[2][0] + k1[2]/2.0, \
				y[3][0] + k1[3]/2.0, iAKT, y[4][0] + k1[4]/2.0, imTOR, y[5][0] + k1[5]/2.0, \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
      
    }
    
    
    for(j=0;j<VARIABLE_N;j++){
      k3[j] = dt*(*function[j])(Ins, pro_IRcom, y[0][0] + k2[0]/2.0, y[1][0] + k2[1]/2.0, y[2][0] + k2[2]/2.0, \
				y[3][0] + k2[3]/2.0, iAKT, y[4][0] + k2[4]/2.0, imTOR, y[5][0] + k2[5]/2.0, \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }
    
    
    for(j=0;j<VARIABLE_N;j++){
      k4[j] = dt*(*function[j])(Ins, pro_IRcom, y[0][0] + k3[0], y[1][0] + k3[1], y[2][0] + k3[2], \
				y[3][0] + k3[3], iAKT, y[4][0] + k3[4], imTOR, y[5][0] + k3[5], \
				p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], \
				p[10], p[11], p[12], p[13], p[14], p[15]);
    }
    
    
    for(j=0;j<VARIABLE_N;j++){
      y[j][1] = y[j][0] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0;
    }
    
    for(j=0;j<VARIABLE_N;j++){
      y[j][0] = y[j][1];
    }
  }
  
  fclose(fp);
}




/* パラメーターの設定 */

void parameter_set(double p[])
{
  p[0] = 5.5702; /* k1_synthesis */
  p[1] = 2.3969; /* k1_InsIRcom */
  p[2] = 2.0934e-2; /* k2_InsIRcom */
  p[3] = 1.2140e-5; /* k1_p1IRcomDeg */
  p[4] = 2.7510e-1; /* k1_p1IRcomPhos */
  p[5] = 7.2509e-3; /* k1_p1p2IRcomdePhos */
  p[6] = 7.5812e+2; /* k1_IRcomPhos */
  p[7] = 9.1758e-1; /* k1_p2IRcomdePhos */
  p[8] = 4.1292e-2; /* k1_p2IRcomDeg */
  p[9] = 1.3032e-4; /* k1_Insp2IRcom */
  p[10] = 1.9200e-4; /* k2_Insp2IRcom */
  p[11] = 2.9311e-2; /* k1_p1p2IRcomDeg */
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
  
