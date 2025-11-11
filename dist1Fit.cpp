#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <TMinuit.h>
#include <TRandom2.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>

using namespace std;

// ---------- GLOBALS ----------
// Declare pointer to data as global
// The use of global variables is a disfavored coding practice, but
// in this case it will allow us to avoid writing more complex code
TH1F *hdata;
// also passing the fit function, to make the objective fcn more generic
TF1 *fparam;

//-------------------------------------------------------------------------
// gumbel distribution for fitting
double gumbel(double* xPtr, double par[]){        
  double x = *xPtr;
  double A = par[0];         // normalization
  double mu = par[1];       // mean of x
  double beta = par[2];
  double z = (x-mu)/beta;
  return A*TMath::Exp(-(z+TMath::Exp(-z)));
}

// double gaussian for fitting
double doubleGaus(double* xPtr, double par[]){        
  double x = *xPtr;
  double A1 = par[0];         // normalization
  double mu1 = par[1];       // mean of x
  double sig1 = par[2];
  double peak1 = A1*TMath::Exp(-0.5*(x-mu1)*(x-mu1)/sig1/sig1);

  double A2 = par[3];
  double mu2 = par[4];
  double sig2 = par[5];
  double peak2 = A2*TMath::Exp(-0.5*(x-mu2)*(x-mu2)/sig2/sig2);
  return peak1 + peak2;
}

//-------------------------------------------------------------------------
// This is our OBJECTIVE Function
// return NLL given a histogram and function

double calcNLL(TH1F* h, TF1* f){
  double nll=0;
  for (int i=1; i<=h->GetNbinsX(); i++){
    double x=h->GetBinCenter(i);
    int n=(int)(h->GetBinContent(i));
    double mu=f->Eval(x);
    if (mu<1e-10) mu=1e-10;    // avoid log(0) problems if we go outside a reasonable range!
    nll -= n * TMath::Log(mu) - mu  - TMath::LnGamma(n+1);
  }
  // cout << "nll "<< nll <<endl;
  return 2*nll;   // factor of 2 so the 1 sigma error contours follow the chi^2 convention
}


// return chi2 given a histogram and function

double calcchi2(TH1F* h, TF1* f){
  double chi2=0;
  for (int i=1; i<=h->GetNbinsX(); i++){
    double x=h->GetBinCenter(i); // bin center
    int n=(int)(h->GetBinContent(i)); // y bin
    double mu=f->Eval(x); // expected y bin
    double err=h->GetBinError(i); // bin error

    if (err==0) continue; // avoid division by zero

    double diff = mu-n;
    chi2 += diff*diff/(err*err);
  }

  return chi2;   
}


//-------------------------------------------------------------------------
// Minuit fcn: calculates value of the function to be minimized using
// the data and the model function
// This is the interface used to define our objective function.  We use
// only a subset of the input parameters below.
// npar: number of parameters
// par: array of parameter values
// f: the value of the objective function
// Minuit can also pass the gradient(deriv) of the objective function wrt the
// current parameters or flags that can trigger special operations (eg perform some initialization)

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  for (int i=0; i<npar; i++){
    fparam->SetParameter(i,par[i]);
  }

  // f = calcNLL(hdata,fparam);
  f = calcchi2(hdata,fparam);
 
}


//-------------------------------------------------------------------------

int main(int argc, char **argv) {


  TCanvas* canvas = new TCanvas("canvas","",1200,600);
  canvas->Divide(2,1);
  canvas->cd(1);

  auto *f=new TFile("distros.root");
  auto dist1 = (TH1F*) f->Get("dist1");
  


  // ***************************************
  // not important for this example
  // Set a bunch of parameters to make the plot look nice
  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);        
  canvas->SetFrameBorderMode(0);  
  //gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();
  gROOT->ForceStyle();

  //gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();
  // ***************************************

  
  //**************Double Gaussian**********************************//
  
  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 6;              // the number of parameters
 
  TMinuit minuit(npar);
  minuit.SetPrintLevel(0);         // set print level of minuit to normal
  minuit.SetFCN(fcn);              // the fcn to be minized (objective fcn) 
  // this calculates chi^2 or NLL, it is not the parameterization of the data!

  // Make a TF1, we can use this to fit the data and plot the results
  TF1* myfunc = new TF1("myfunc", doubleGaus, dist1->GetXaxis()->GetXmin(), dist1->GetXaxis()->GetXmax(), npar);  
  double par[npar];               // the free parameters
  double stepSize[npar];          // (limiting) step sizes for parameters
  double minVal[npar];            // minimum/maximum bounds on parameters 
  double maxVal[npar];            
  TString parName[npar];          // Optional names for nicer output!

  // Initial parameters MUST be set by some means to get things started
  par[0] = 0.9*dist1->GetMaximum(); // A1
  par[1] = 80; // mu1
  par[2] = 5; // sigma1
  par[3] = 0.2*par[0]; // A2
  par[4] = 90; //mu2
  par[5] = 10; //sigma2
  stepSize[0] = TMath::Abs(par[0]*0.1);   // usually 10% of initial guess is OK for starting step size, YMMV
  stepSize[1] = TMath::Abs(par[1]*0.1);   // step sizes MUST be positive!
  stepSize[2] = TMath::Abs(par[2]*0.1);
  stepSize[3] = TMath::Abs(par[3]*0.1);
  stepSize[4] = TMath::Abs(par[4]*0.1);
  stepSize[5] = TMath::Abs(par[5]*0.1);
  minVal[0] = 0;      // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 0;
  minVal[1] = 60; 
  maxVal[1] = 120;
  minVal[2] = 0; 
  maxVal[2] = 0;

  minVal[3] = 0; 
  maxVal[3] = 0;
  minVal[4] = 60; 
  maxVal[4] = 120;
  minVal[5] = 0; 
  maxVal[5] = 0;

  parName[0] = "A1";       // let's give our fit parameters useful names
  parName[1] = "mu1";
  parName[2] = "sig1";       
  parName[3] = "A2";
  parName[4] = "mu2";       
  parName[5] = "sig2";

  // initialize the parameters
  for (int i=0; i<npar; i++){
    minuit.DefineParameter(i, parName[i].Data(), 
			   par[i], stepSize[i], minVal[i], maxVal[i]);
  }

  // here we define the pointers to pass information to Minuit's fcn to calculate the
  // objective function
  // note: the use of global variables is discouraged in general, but for these
  // examples we'll let convenience and code simplicity outweigh the programming
  // style considerations.
  hdata=dist1;     // histogram to fit
  fparam=myfunc;  // our model

  // Do the minimization!
  minuit.Migrad();       // Minuit's best minimization algorithm

  // Get the result
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
    minuit.GetParameter(i,outpar[i],err[i]);
  }
  // // run a minos error analysis
  // // this perfoms a brute force scan around the minimum of the objective function
  // // to get better estimates of the errors on the fit parameters
  // // gMinuit->mnmnos();

  // store the fit parameters in our TF1, decorate our function
  myfunc->SetParameters(outpar);

  myfunc->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
  // myfunc->SetLineColor(1);             //  black (default)
  myfunc->SetLineWidth(2);

  myfunc->GetXaxis()->SetTitle("x");
  myfunc->GetYaxis()->SetTitle("Double Gaussian");

  // summarize the fitting results
  cout << "\n==========================\n"<<endl;
  cout << "Results of Double Gaussian"<< endl;
  double fmin, fedm, errdef;
  int npari, nparx, istat;  // see https://root.cern/doc/master/classTMinuit.html 
  minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  // cout << "minimum of NLL = " << fmin << endl;
  cout << "minimum of chi sqr = " << fmin << endl;
  cout << "fit status = " << istat << endl;
  cout << "best fit parameters\n" <<endl;
  for (int i=0; i<npar; ++i){
    cout << i << " : " << outpar[i] << " +- " << err[i] << endl;
  }

  int ndof = dist1->GetNbinsX()-npar;
  double redChiSqr = fmin/ndof;
  double p = TMath::Prob(fmin, ndof);
  cout << "degrees of freedom = " << ndof << endl;
  cout << "reduced chi sqr = " << redChiSqr << endl;
  cout << "p value = " << p << endl;


  //****************************Gumbel***********************//
  // Initialize minuit, set initial values etc. of parameters.
  const int nparG = 3;              // the number of parameters

  TMinuit minuitG(nparG);
  minuitG.SetPrintLevel(0);         // set print level of minuit to normal
  minuitG.SetFCN(fcn);              // the fcn to be minized (objective fcn) 
  // this calculates chi^2 or NLL, it is not the parameterization of the data!

  // Make a TF1, we can use this to fit the data and plot the results
  TF1* myfuncG = new TF1("myfuncG", gumbel, dist1->GetXaxis()->GetXmin(), dist1->GetXaxis()->GetXmax(), nparG);  
  double parG[nparG];               // the free parameters
  double stepSizeG[nparG];          // (limiting) step sizes for parameters
  double minValG[nparG];            // minimum/maximum bounds on parameters 
  double maxValG[nparG];            
  TString parNameG[nparG];          // Optional names for nicer output!

  // Initial parameters MUST be set by some means to get things started
  parG[0] = 2.8*dist1->GetMaximum(); // A1
  parG[1] = 80; // mu
  parG[2] = 6; // beta
  stepSizeG[0] = TMath::Abs(parG[0]*0.1);   // usually 10% of initial guess is OK for starting step size, YMMV
  stepSizeG[1] = TMath::Abs(parG[1]*0.1);   // step sizes MUST be positive!
  stepSizeG[2] = TMath::Abs(parG[2]*0.1);
  minValG[0] = 0;      // if min and max values = 0, parameter is unbounded.
  maxValG[0] = 4*dist1->GetMaximum();
  minValG[1] = 70; 
  maxValG[1] = 90;
  minValG[2] = 2; 
  maxValG[2] = 10;


  parNameG[0] = "A";       // let's give our fit parameters useful names
  parNameG[1] = "mu";
  parNameG[2] = "beta";

  // initialize the parameters
  for (int i=0; i<nparG; i++){
    minuitG.DefineParameter(i, parNameG[i].Data(), 
        parG[i], stepSizeG[i], minValG[i], maxValG[i]);
  }

  // here we define the pointers to pass information to Minuit's fcn to calculate the
  // objective function
  // note: the use of global variables is discouraged in general, but for these
  // examples we'll let convenience and code simplicity outweigh the programming
  // style considerations.
  hdata=dist1;     // histogram to fit
  fparam=myfuncG;  // our model

  // Do the minimization!
  minuitG.Migrad();       // Minuit's best minimization algorithm

  // Get the result
  double outparG[nparG], errG[nparG];
  for (int i=0; i<nparG; i++){
    minuitG.GetParameter(i,outparG[i],errG[i]);
  }

  // store the fit parameters in our TF1, decorate our function
  myfuncG->SetParameters(outparG);

  myfuncG->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
  // myfuncG->SetLineColor(1);             //  black (default)
  myfuncG->SetLineWidth(2);

  myfuncG->GetXaxis()->SetTitle("x");
  myfuncG->GetYaxis()->SetTitle("Gumbel");

  // summarize the fitting results
  cout << "\n==========================\n"<<endl;
  cout << "Results of Gumbel"<< endl;
  double fminG, fedmG, errdefG;
  int npariG, nparxG, istatG;  // see https://root.cern/doc/master/classTMinuit.html 
  minuit.mnstat(fminG, fedmG, errdefG, npariG, nparxG, istatG);
  // cout << "minimum of NLL = " << fminG << endl;
  cout << "minimum of chi sqr = " << fminG << endl;
  cout << "fit status = " << istatG << endl;
  cout << "best fit parameters\n" <<endl;
  for (int i=0; i<nparG; ++i){
    cout << i << " : " << outparG[i] << " +- " << errG[i] << endl;
  }

  int ndofG = dist1->GetNbinsX()-nparG;
  double redChiSqrG = fminG/ndofG;
  double pG = TMath::Prob(fminG, ndofG);
  cout << "degrees of freedom = " << ndofG << endl;
  cout << "reduced chi sqr = " << redChiSqrG << endl;
  cout << "p value = " << pG << endl;



  // plot

  myfunc->SetTitle("Double Gaussian;;");
  myfunc->Draw();
  dist1->Draw("e same");
  

  canvas->cd(2);
  myfuncG->SetTitle("Gumbel Distribution;;");
  myfuncG->Draw();
  dist1->Draw("e same");
  

  canvas->Update();
  canvas->SaveAs("ex1.pdf");
  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;

  return 0;

}
