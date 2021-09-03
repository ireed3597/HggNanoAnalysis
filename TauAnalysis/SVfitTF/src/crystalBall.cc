#include "TauAnalysis/SVfitTF/interface/crystalBall.h"
#include <math.h>
#include "TMath.h"
#include <iostream>

using namespace std;

double sum_gaus_fnc(double *x, double *par); // sum of 3gaus+pol7

double exp_fcn(double x, double norm, double slope); // exp

double gauss_fcn(double x, double norm, double mean, double sigma); // gauss

double pol_fcn(double x, double p0, double p1=0, double p2=0, double p3=0, double p4=0, double p5=0, double p6=0, double p7=0); // pol1


// function that return the analytic computation (from mathematica) of the CB
double normalizedCrystalBall(double *pp)
{
  double p1       = pp[6];
  //double q0       = pp[16]; // pol right
  double q1       = pp[17];
  double sexp     = pp[4];
  double ng1Left  = pp[0]; // norm gaus
  double ng2      = pp[10];
  double ng3      = pp[13];
  double mg1      = pp[1]; // mean gaus
  double mg2      = mg1+pp[11];
  double mg3      = mg2+pp[14];
  double sg1Left  = pp[2]; // sigma
  double sg1Right = pp[9];
  double sg2      = pp[12];
  double sg3      = pp[15];
  double r2       = mg1-0.5*sg1Left-pp[8]; // range
  double r1       = r2-pp[7];
  double r3       = pp[18];

  double nexp     = gauss_fcn(r2, ng1Left, mg1, sg1Left)/TMath::Exp(sexp*r2);
  double p0       = exp_fcn(r1, nexp, sexp) - p1*r1 ;
  double R3[]{r3};
  double q0       = sum_gaus_fnc(R3, pp)-q1*r3;

  // these are integrals from mathematica
  double polLeft  = 0.5*r1*(2*p0+p1*r1);

  double polRight = -0.5*(-2+r3)*(2*q0+q1*(2+r3));

  double expo     = ((-exp(r1*sexp) + exp(r2*sexp))*nexp)/sexp;

  double gaus1Left = ng1Left*sqrt(TMath::Pi()/2)*sg1Left*TMath::Erf((mg1 - r2)/(sqrt(2)*sg1Left));

  double ng1Right = ng1Left - ( ng2*exp(-0.5*pow((mg1 - mg2)/sg2, 2)) +
                    ng3*exp(-0.5*pow((mg1 - mg3)/sg3, 2)) ) ;

  double sumGaus =
    -sqrt(TMath::Pi()/2)*
    (ng1Right* sg1Right* TMath::Erf((mg1 - r3) /(sqrt(2)* sg1Right)) +
     ng2* sg2*           TMath::Erf((mg1 - mg2)/(sqrt(2)* sg2))      +
     ng2* sg2*           TMath::Erf((mg2 - r3) /(sqrt(2)* sg2))      +
     ng3* sg3*            TMath::Erf((mg1 - mg3)/(sqrt(2)* sg3))      +
     ng3* sg3*            TMath::Erf((mg3 - r3) /(sqrt(2)* sg3)))     ;

  return polLeft + expo + gaus1Left +  sumGaus + polRight;
} // normalizedCrystalBall


double crystalBall(double *xx, double *pp)
{
  double x = xx[0];

  double normGaus1     = pp[0];
  double meanGaus1     = pp[1];
  double sigGaus1_left = pp[2];
  double slopeExp_left = pp[4];
  double p1_left       = pp[6];
  double p1_right      = pp[17];
  double R2            = meanGaus1 - 0.5*sigGaus1_left - pp[8];
  double R1            = R2 - pp[7]; // min exp
  double R3            = pp[18];
  bool drawPolLeft     {true};
  bool drawExp         {true};
  bool drawGaus1       {true};
  bool drawGaus2       {true};
  bool drawGaus3       {true};
  bool drawPolRight    {true};
  bool drawAnyRight    = drawGaus1 || drawGaus2 || drawGaus3;

  // remove discontinuity at x=R2: solve exp=gauss for norm exp and replace R2 by (meanGaus1-R2*R2)
  double normExp_left = gauss_fcn(R2, normGaus1, meanGaus1, sigGaus1_left)/TMath::Exp(slopeExp_left*R2);

  // remove discontinuity at x=R1: solve pol1=exp for p0 left
  double p0_left = exp_fcn(R1, normExp_left, slopeExp_left) - p1_left*R1 ;

  // remove discontinuity at x=R3: solve pol1=sumGaus for p0 right
  double r3[]{R3};
  double p0_right = sum_gaus_fnc(r3, pp)-p1_right*R3;

  // compute probability
  double result = 0.;
  if (drawPolLeft  && x >= 0.        && x < R1       ) result += pol_fcn  (x, p0_left, p1_left);
  if (drawExp      && x >= R1        && x < R2       ) result += exp_fcn  (x, normExp_left, slopeExp_left);
  if (drawGaus1    && x >= R2        && x < meanGaus1) result += gauss_fcn(x, normGaus1, meanGaus1, sigGaus1_left);
  if (drawAnyRight && x >= meanGaus1 && x < R3       ) result += sum_gaus_fnc(xx, pp);
  if (drawPolRight && x >= R3        && x < 2.       ) result += pol_fcn  (x, p0_right, p1_right);

  return result;
} // crystalBall



//////////////////////
double exp_fcn(double x, double norm, double slope){
  return norm*TMath::Exp(slope*x);
}

double gauss_fcn(double x, double norm, double mean, double sigma){
  return norm*TMath::Exp(-0.5*TMath::Power((x - mean)/sigma, 2.));
}

double pol_fcn(double x, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7){
  return p0 + p1*x + p2*TMath::Power(x,2) +
    p3*TMath::Power(x,3) + p4*TMath::Power(x,4) +
    p5*TMath::Power(x,5) + p6*TMath::Power(x,6) +
    p7*TMath::Power(x,7);
}

double sum_gaus_fnc(double *xx, double *pp)
{
  double x = xx[0];

  double normGaus1Left  = pp[0];
  double meanGaus1      = pp[1];
  double valGaus1_at_meanGaus1 = normGaus1Left;

  double normGaus2      = pp[10];
  double meanGaus2      = meanGaus1 + pp[11];
  double sigGaus2       = pp[12];
  double valGaus2_at_meanGaus1 = gauss_fcn(meanGaus1, normGaus2, meanGaus2, sigGaus2);
  double valGaus2       = gauss_fcn(x, normGaus2, meanGaus2, sigGaus2);

  double normGaus3      = pp[13];
  double meanGaus3      = meanGaus2 + pp[14];
  double sigGaus3       = pp[15];
  double valGaus3_at_meanGaus1 = gauss_fcn(meanGaus1, normGaus3, meanGaus3, sigGaus3);
  double valGaus3       = gauss_fcn(x, normGaus3, meanGaus3, sigGaus3);

  double normGaus1Right = valGaus1_at_meanGaus1 - (valGaus2_at_meanGaus1 + valGaus3_at_meanGaus1);
  double sigGaus1_right = pp[9];
  double valGaus1_right = gauss_fcn(x, normGaus1Right, meanGaus1, sigGaus1_right);

  bool drawGaus1 {true};
  bool drawGaus2 {true};
  bool drawGaus3 {true};

  //if ( drawGaus1 && drawGaus2 && drawGaus3 && x > 0.99 && x < 1.01 ) {
  //  std::cout << "<sum_gaus_fnc>:" << std::endl;
  //  std::cout << " x = " << x << ": g1 = " << valGaus1_right << ", g2 = " << valGaus2 << ", g3 = " << valGaus3 << std::endl;
  //}

  double sum = 0.;
  if ( drawGaus1       ) sum += valGaus1_right;
  if ( drawGaus2       ) sum += valGaus2;
  if ( drawGaus3       ) sum += valGaus3;
  if ( normGaus1Right < 0. ) sum /= (1. + normGaus1Right*normGaus1Right);

  return sum;
}





