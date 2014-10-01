/* Insulin-AKT model ODE functions */

double fIRcom(double Ins, double pro_IRcom, double IRcom, double p1IRcom, double p2IRcom, double p1p2IRcom, \
	      double iAKT, double pAKT, double imTOR, double pmTOR,	\
	      double k1_synthesis, double k1_InsIRcom, double k2_InsIRcom, double k1_p1IRcomDeg, \
	      double k1_p1IRcomPhos, double k1_p1p2IRcomdePhos, double k1_IRcomPhos, double k1_p2IRcomdePhos, \
	      double k1_p2IRcomDeg, double k1_Insp2IRcom, double k2_Insp2IRcom, double k1_p1p2IRcomDeg, \
	      double k1_AKTPhos, double k1_pAKTdePhos, double k1_mTORPhos, double k1_pmTORdePhos)
{
  return k1_synthesis*(pro_IRcom - IRcom)		\
    - k1_InsIRcom*Ins*IRcom + k2_InsIRcom*p1IRcom	\
    - k1_IRcomPhos*pmTOR*IRcom			\
    + k1_p2IRcomdePhos*p2IRcom;
}


double fp1IRcom(double Ins, double pro_IRcom, double IRcom, double p1IRcom, double p2IRcom, double p1p2IRcom, \
		double iAKT, double pAKT, double imTOR, double pmTOR,	\
		double k1_synthesis, double k1_InsIRcom, double k2_InsIRcom, double k1_p1IRcomDeg, \
		double k1_p1IRcomPhos, double k1_p1p2IRcomdePhos, double k1_IRcomPhos, double k1_p2IRcomdePhos, \
		double k1_p2IRcomDeg, double k1_Insp2IRcom, double k2_Insp2IRcom, double k1_p1p2IRcomDeg, \
		double k1_AKTPhos, double k1_pAKTdePhos, double k1_mTORPhos, double k1_pmTORdePhos)
{
  return k1_InsIRcom*Ins*IRcom - k2_InsIRcom*p1IRcom	\
    - k1_p1IRcomDeg*p1IRcom				\
    - k1_p1IRcomPhos*pmTOR*p1IRcom			\
    + k1_p1p2IRcomdePhos*p1p2IRcom;
}


double fp2IRcom(double Ins, double pro_IRcom, double IRcom, double p1IRcom, double p2IRcom, double p1p2IRcom, \
		double iAKT, double pAKT, double imTOR, double pmTOR,	\
		double k1_synthesis, double k1_InsIRcom, double k2_InsIRcom, double k1_p1IRcomDeg, \
		double k1_p1IRcomPhos, double k1_p1p2IRcomdePhos, double k1_IRcomPhos, double k1_p2IRcomdePhos, \
		double k1_p2IRcomDeg, double k1_Insp2IRcom, double k2_Insp2IRcom, double k1_p1p2IRcomDeg, \
		double k1_AKTPhos, double k1_pAKTdePhos, double k1_mTORPhos, double k1_pmTORdePhos)
{
  return k1_IRcomPhos*pmTOR*IRcom		\
    - k1_p2IRcomdePhos*p2IRcom			\
    - k1_p2IRcomDeg*p2IRcom					\
    - k1_Insp2IRcom*Ins*p2IRcom + k2_Insp2IRcom*p1p2IRcom;
}


double fp1p2IRcom(double Ins, double pro_IRcom, double IRcom, double p1IRcom, double p2IRcom, double p1p2IRcom, \
		  double iAKT, double pAKT, double imTOR, double pmTOR,	\
		  double k1_synthesis, double k1_InsIRcom, double k2_InsIRcom, double k1_p1IRcomDeg, \
		  double k1_p1IRcomPhos, double k1_p1p2IRcomdePhos, double k1_IRcomPhos, double k1_p2IRcomdePhos, \
		  double k1_p2IRcomDeg, double k1_Insp2IRcom, double k2_Insp2IRcom, double k1_p1p2IRcomDeg, \
		  double k1_AKTPhos, double k1_pAKTdePhos, double k1_mTORPhos, double k1_pmTORdePhos)
{
  return k1_p1IRcomPhos*pmTOR*p1IRcom		\
    - k1_p1p2IRcomdePhos*p1p2IRcom			  \
    + k1_Insp2IRcom*Ins*p2IRcom - k2_Insp2IRcom*p1p2IRcom \
    - k1_p1p2IRcomDeg*p1p2IRcom;
}


double fpAKT(double Ins, double pro_IRcom, double IRcom, double p1IRcom, double p2IRcom, double p1p2IRcom, \
	     double iAKT, double pAKT, double imTOR, double pmTOR,	\
	     double k1_synthesis, double k1_InsIRcom, double k2_InsIRcom, double k1_p1IRcomDeg, \
	     double k1_p1IRcomPhos, double k1_p1p2IRcomdePhos, double k1_IRcomPhos, double k1_p2IRcomdePhos, \
	     double k1_p2IRcomDeg, double k1_Insp2IRcom, double k2_Insp2IRcom, double k1_p1p2IRcomDeg, \
	     double k1_AKTPhos, double k1_pAKTdePhos, double k1_mTORPhos, double k1_pmTORdePhos)
{
  return k1_AKTPhos*(iAKT - pAKT)*p1IRcom	\
    - k1_pAKTdePhos*pAKT;
}


double fpmTOR(double Ins, double pro_IRcom, double IRcom, double p1IRcom, double p2IRcom, double p1p2IRcom, \
	      double iAKT, double pAKT, double imTOR, double pmTOR,	\
	      double k1_synthesis, double k1_InsIRcom, double k2_InsIRcom, double k1_p1IRcomDeg, \
	      double k1_p1IRcomPhos, double k1_p1p2IRcomdePhos, double k1_IRcomPhos, double k1_p2IRcomdePhos, \
	      double k1_p2IRcomDeg, double k1_Insp2IRcom, double k2_Insp2IRcom, double k1_p1p2IRcomDeg, \
	      double k1_AKTPhos, double k1_pAKTdePhos, double k1_mTORPhos, double k1_pmTORdePhos)		       
{
  return k1_mTORPhos*(imTOR - pmTOR)*pAKT	\
    - k1_pmTORdePhos*pmTOR;
}



