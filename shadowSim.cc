// Simulate the effect of shadowing of CRs by sun and moon
// Max Malacari - April 2018

// Classics
#include <iostream>
#include <sstream>

// ROOT headers
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TPaveStats.h"

using namespace std;

// Globals
const double sunMoonRad = 2.5; // degrees - angular radius of sun/moon
const double maxAngle = 25.; // degrees - max. angle to simulate from sun/moon (should be larger than fit range)
const double angRes = 2.; // degrees - detector angular resolution
const double maxFitTheta = 10.; // degrees - only fit out to this angle
const int nEvents = 500000; // how many events do you want within maxFitTheta?

// Function prototypes
double GetRandomPhi();
double GetRandomTheta();
void GaussFlucVector(double sigma, double theta1, double phi1, double &theta2, double &phi2); // fluctuate vector directly
void GaussFlucTangentPlane(double sigma, double theta1, double phi1, double &theta2, double &phi2); // fluctuate in tangent plane
static inline void loadBar(long long int x, long long int n, int r, int w);

int main(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetPalette(51);
  
  const int nThetaBins = 100;

  // Some diagnostic plots
  TH1F phiDist("phiDist","Polar angle around sun/moon",100,0.,360.);
  TH1F phiDist_smeared("phiDist_smeared","",100,0.,360.);
  TH1F thetaDist("thetaDist","Angular sep. from sun/moon",nThetaBins,0.,maxAngle+5.);
  TH1F thetaDist_smeared("thetaDist_smeared","",nThetaBins,0.,maxAngle+5.);
  TH1F smearPlot("smearPlot","Space angle",100,0.,5*angRes);
  TH2F tangentPlane("tangentPlane","Tangent plane w/o ang. res.",80,-maxAngle-5,maxAngle+5,80,-maxAngle-5,maxAngle+5);
  TH2F tangentPlane_smeared("tangentPlane_smeared","Tangent plane w/ ang. res.",80,-maxAngle-5,maxAngle+5,80,-maxAngle-5,maxAngle+5);
  
  cout << "Simulating " << nEvents << " events within " << maxFitTheta << " deg of the sun/moon..." << endl;
  

  double theta, phi, theta_smeared, phi_smeared;
  int i = 0;
  while (i < nEvents) { // generate events until we have nEvents within maxFitTheta deg

    theta = GetRandomTheta() * 180/TMath::Pi();
    phi = GetRandomPhi() * 180/TMath::Pi();

    //GaussFlucTangentPlane(angRes, theta, phi, theta_smeared, phi_smeared);

    // New smearing
    GaussFlucVector(angRes, theta, phi, theta_smeared, phi_smeared);

    if (theta_smeared <= maxFitTheta) i++;

    thetaDist.Fill(theta);
    phiDist.Fill(phi);
    tangentPlane.Fill(theta*cos(phi*TMath::Pi()/180.), theta*sin(phi*TMath::Pi()/180.));
    thetaDist_smeared.Fill(theta_smeared);
    phiDist_smeared.Fill(phi_smeared);
    tangentPlane_smeared.Fill(theta_smeared*cos(phi_smeared*TMath::Pi()/180.), theta_smeared*sin(phi_smeared*TMath::Pi()/180.));

    // Check the smearing angle
    TVector3 vec1(1.,1.,1.), vec2(1.,1.,1.);
    vec1.SetMagThetaPhi(1., theta*TMath::Pi()/180., phi*TMath::Pi()/180.);
    vec2.SetMagThetaPhi(1., theta_smeared*TMath::Pi()/180., phi_smeared*TMath::Pi()/180.);
    smearPlot.Fill( vec1.Angle(vec2)*180./TMath::Pi() );
    
    loadBar(i, nEvents, 100, 50);
  } // end for loop
  


  // Create the plot of interest and fit ------------

  TH1F eventDist("eventDist","",nThetaBins,0.,maxAngle+5.);
  
  // Fill histogram and normalize
  for (int i=1; i<=nThetaBins;i++) {
    double angleLo = eventDist.GetBinLowEdge(i) * TMath::Pi()/180.;
    double angleHi = angleLo + eventDist.GetBinWidth(i)*TMath::Pi()/180.;
    double dOmega = 2.*TMath::Pi()* (cos(angleLo)-cos(angleHi)); // sr
    double dDeg2 = TMath::Pi() * ( pow(angleHi*180./TMath::Pi(),2.) - pow(angleLo*180./TMath::Pi(),2.) ); // deg^2

    eventDist.SetBinContent(i, thetaDist_smeared.GetBinContent(i)/dOmega);
    eventDist.SetBinError(i, sqrt(thetaDist_smeared.GetBinContent(i))/dOmega);
  }

  // Flat function
  stringstream flatFuncStr;
  flatFuncStr << "[0]";
  TF1 flatFunc("flatFunc",flatFuncStr.str().c_str(),0.,180.); 
  flatFunc.SetParameter(0, eventDist.GetMaximum());
  eventDist.Fit("flatFunc","","", maxFitTheta, maxFitTheta+10.); // Fit where we expect a flat distribution
  
  // Ang. res. function
  stringstream fitFuncStr;
  fitFuncStr << "[0]*(1. - " << sunMoonRad*sunMoonRad << "/(2.*[1]^2)*exp(-x^2 / (2.*[1]^2)))";
  TF1 fitFunc("fitFunc",fitFuncStr.str().c_str(),0,180.); // fit only up to maxFitTheta
  //fitFunc.SetParameter(0, eventDist.GetMaximum());
  fitFunc.FixParameter(0, flatFunc.GetParameter(0)); // Fix this parameter to the flat func. fit
  fitFunc.SetParameter(1, 2);
  fitFunc.SetParName(0, "Norm.");
  fitFunc.SetParName(1, "Ang. res.");
  eventDist.Fit("fitFunc","","", 0., maxFitTheta); // fit this up to maxFitTheta

  // Draw and save plot
  TCanvas fitCanvas;
  eventDist.GetXaxis()->SetRangeUser(0,maxFitTheta);
  eventDist.SetXTitle("Angle to sun/moon [deg]");
  eventDist.SetYTitle("events / sr");
  eventDist.GetXaxis()->CenterTitle();
  eventDist.GetYaxis()->CenterTitle();
  eventDist.GetYaxis()->SetTitleOffset(1.5);
  eventDist.SetLineWidth(2);
  eventDist.SetLineColor(kBlack);
  eventDist.SetMarkerStyle(8);
  eventDist.SetMarkerSize(0.7);
  eventDist.Draw();
  flatFunc.SetLineWidth(2);
  flatFunc.SetLineColor(kBlue);
  flatFunc.Draw("same");
  fitFunc.SetLineWidth(2);
  fitFunc.SetLineColor(kRed);
  fitFunc.Draw("same");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)eventDist.FindObject("stats");
  st->SetX1NDC(0.534384);
  st->SetX2NDC(0.8954155);
  st->SetY1NDC(0.1033755);
  st->SetY2NDC(0.2236287);
  st->SetLineColor(kWhite);
  fitCanvas.SaveAs("./output/fit.pdf");
  
  // -----------------------------------------

  
  // Make some other interesting disagnostic plots
  TCanvas diagnosticCanvas("diagCanvas","",900,800);
  diagnosticCanvas.Divide(2,2);
  diagnosticCanvas.cd(1);
  thetaDist.SetMinimum(0);
  thetaDist.SetLineColor(kBlue);
  thetaDist.SetLineWidth(2);
  thetaDist_smeared.SetLineWidth(2);
  thetaDist_smeared.SetLineColor(kRed);
  thetaDist.SetXTitle("#theta [deg]");
  thetaDist.SetYTitle("events");
  thetaDist.GetXaxis()->CenterTitle();
  thetaDist.GetYaxis()->CenterTitle();
  thetaDist.GetYaxis()->SetTitleOffset(1.5);
  thetaDist.Draw();
  thetaDist_smeared.Draw("same");
  diagnosticCanvas.cd(2);
  // phiDist.SetMinimum(0);
  // phiDist.SetLineWidth(2);
  // phiDist.SetLineColor(kBlue);
  // phiDist_smeared.SetLineWidth(2);
  // phiDist_smeared.SetLineColor(kRed);
  // phiDist.SetXTitle("#phi [deg]");
  // phiDist.SetYTitle("Events");
  // phiDist.GetXaxis()->CenterTitle();
  // phiDist.GetYaxis()->CenterTitle();
  // phiDist.GetYaxis()->SetTitleOffset(1.5);
  // phiDist.Draw();
  // phiDist_smeared.Draw("same");
  smearPlot.SetMinimum(0);
  smearPlot.SetLineWidth(2);
  smearPlot.SetLineColor(kBlue);
  smearPlot.SetXTitle("Smearing angle [deg]");
  smearPlot.SetYTitle("events");
  smearPlot.GetXaxis()->CenterTitle();
  smearPlot.GetYaxis()->CenterTitle();
  smearPlot.GetYaxis()->SetTitleOffset(1.5);
  smearPlot.Draw();
  diagnosticCanvas.cd(3);
  tangentPlane.SetXTitle("X [deg]");
  tangentPlane.SetYTitle("Y [deg]");
  tangentPlane.GetXaxis()->CenterTitle();
  tangentPlane.GetYaxis()->CenterTitle();
  tangentPlane.Draw("COL4Z");
  diagnosticCanvas.cd(4);
  tangentPlane_smeared.SetXTitle("X [deg]");
  tangentPlane_smeared.SetYTitle("Y [deg]");
  tangentPlane_smeared.GetXaxis()->CenterTitle();
  tangentPlane_smeared.GetYaxis()->CenterTitle();
  tangentPlane_smeared.Draw("COL4Z");
  diagnosticCanvas.SaveAs("./output/diagCanvas.pdf");

  return 0;

}

// Generate a random polar angle around the sun/moon position
double GetRandomPhi() {
  static TRandom3 randGen(time(NULL));
  return randGen.Uniform(0., 2.*TMath::Pi());
}

// Generate a random angular separation between sunMoonRad and maxAngle
double GetRandomTheta() {
  static TRandom3 randGen(time(NULL)+1);
  static const double cosThetaMin = cos(sunMoonRad*TMath::Pi()/180.);
  static const double cosThetaMax = cos(maxAngle*TMath::Pi()/180.);
  const double randN = randGen.Uniform(0.,1.);
  return acos( (1.-randN)*cosThetaMin + randN*cosThetaMax );
}

// Smear the direction according to the detector angular resolution
void GaussFlucTangentPlane(double sigma, double theta1, double phi1, double &theta2, double &phi2){

  static TRandom3 randGen(time(NULL));
  
  double g1 = randGen.Uniform(0.,1.);
  double g2 = randGen.Uniform(0.,1.);

  double u1 = sigma*sqrt(-2.*log(g1))*sin(2.*TMath::Pi()*g2)*TMath::Pi()/180.;
  double u2 = sigma*sqrt(-2.*log(g1))*cos(2.*TMath::Pi()*g2)*TMath::Pi()/180.;
  
  double w1[3], w2[3], w3[3];
  // unit vector aligned the (theta_ref, phi_ref) direction
  w3[0] = sin(theta1*TMath::Pi()/180.)*cos(phi1*TMath::Pi()/180.);
  w3[1] = sin(theta1*TMath::Pi()/180.)*sin(phi1*TMath::Pi()/180.);
  w3[2] = cos(theta1*TMath::Pi()/180.);
  
  // unit vector perpendicular to (theta_ref, phi_ref) direction  
  w2[0] =  cos(theta1*TMath::Pi()/180.)*cos(phi1*TMath::Pi()/180.);
  w2[1] =  cos(theta1*TMath::Pi()/180.)*sin(phi1*TMath::Pi()/180.);
  w2[2] = -sin(theta1*TMath::Pi()/180.);

  // unit vector equal to w2 x w3
  w1[0] = w2[1]*w3[2] - w2[2]*w3[1];
  w1[1] = w2[2]*w3[0] - w2[0]*w3[2];
  w1[2] = w2[0]*w3[1] - w2[1]*w3[0];
  
  double x = w3[0] + u1*w1[0] + u2*w2[0];
  double y = w3[1] + u1*w1[1] + u2*w2[1];
  double z = w3[2] + u1*w1[2] + u2*w2[2];
  
  theta2 = acos(z/sqrt(x*x+y*y+z*z))*180./TMath::Pi();
  phi2   = atan2(y,x)*180./TMath::Pi();

  if (phi2 < 0)
    phi2 += 360.;

}

void GaussFlucVector(double sigma, double theta1, double phi1, double &theta2, double &phi2) {
  
  static TRandom3 randGen(time(NULL)+13);

  static TF1 angResFun("angResFunc","x*exp(-x^2/(2*[0]^2))",0,15.*TMath::Pi()/180.);
  static bool firstRun = true;
  if (firstRun) {
    angResFun.SetParameter(0,sigma*TMath::Pi()/180.);
    firstRun = false;
  }
  
  double randPhi = randGen.Uniform(0.,2*TMath::Pi());
  //double randTheta = randGen.Gaus(0.,sigma*TMath::Pi()/180.);
  double randTheta = angResFun.GetRandom();
  //cout << randTheta << endl;
  
  TVector3 flucVec(1.,1.,1.);
  flucVec.SetTheta(randTheta);
  flucVec.SetPhi(randPhi);
  flucVec.SetMag(1.);

  //cout << randTheta << " " << randPhi << " " << flucVec.Theta() << " " << flucVec.Phi() << endl;

  TVector3 eventVec(1.,1.,1.);
  eventVec.SetTheta(theta1*TMath::Pi()/180.);
  eventVec.SetPhi(phi1*TMath::Pi()/180.);
  eventVec.SetMag(1.);
  flucVec.RotateUz(eventVec);

  theta2 = flucVec.Theta()*180./TMath::Pi();
  phi2 = flucVec.Phi()*180./TMath::Pi();
}

// A nice loading bar to show simulation progress
static inline void loadBar(long long int x, long long int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r +1) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}
