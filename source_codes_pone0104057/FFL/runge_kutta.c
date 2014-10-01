#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define EQ_N 2      /* タンパク質の種類 = 微分方程式数 */
#define PARA_N 8     /* パラメーター (速度定数) の数 */
#define RK_N 1050

void rk(double dt, double y[][2], double (*func[EQ_N])(), double p[], double X);
void function_set(double (*func[EQ_N])());
void parameter_set(double *p);

double fY(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H);
double fZ(double X, double Y, double Z, double alpha_y, double alpha_z, double beta_y, double beta_z, double K_xy, double K_xz, double K_yz, double H);


/* MAIN 関数 */

int main(void)
{
  double y[EQ_N][2];
  double dt;
  double X;
  double (*func[EQ_N])();
  double p[PARA_N];
  int i;

  X = 1.0;

  dt = 0.01;
  function_set(func);

  /* Y と Z */
  for(i=0;i<EQ_N;i++){
    y[i][0] = 0.0;
  }

  parameter_set(p);
 
  rk(dt, y, func, p, X);

  return 0;
}



/* Functions */

/* 4th order Runge-Kutta method */

void rk(double dt, double y[][2], double (*func[EQ_N])(), double p[], double X)
{
  int i,j,k;
  FILE *fp;
  double k1[EQ_N],k2[EQ_N],k3[EQ_N],k4[EQ_N];

  fp = fopen("rk.data","w");

  for(i=0;i<=RK_N;i++){
    fprintf(fp,"%f\t %f\n", dt*i, y[1][0]);
    //if(i>=500) X = 0.0;


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

  
  fclose(fp);
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
  return beta_z*pow(X/K_xz, H)/(1+pow(X/K_xz, H))*pow(Y/K_yz, H)/(1+pow(Y/K_yz, H)) - alpha_z*Z;
}




