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
TH1F *hdata1;
TH1F *hdata2;
// also passing the fit function, to make the objective fcn more generic
TF1 *fparam1;
TF1 *fparam2;

//-------------------------------------------------------------------------
// gaussian plus exponentially falling background
double gausExpBkg(double* xPtr, double par[]){        
  double x = *xPtr;
  double A1 = par[0];         // normalization
  double mu = par[1];       // mean of x
  double sigma = par[2];
  double Abkg = par[3];
  double lambda = par[4];
  
  double peak = A1*TMath::Exp(-0.5*(x-mu)*(x-mu)/sigma/sigma);
  return peak + Abkg * TMath::Exp(-x/lambda);
}

// gaussian plus power law background
double gausPowBkg(double* xPtr, double par[]){        
  double x = *xPtr;
  double A2 = par[5];         // normalization
  double mu = par[1];       // mean of x
  double sigma = par[2];
  double peak = A2*TMath::Exp(-0.5*(x-mu)*(x-mu)/sigma/sigma);

  double Abkg2 = par[6];
  double n = par[7];
  return peak + Abkg2 * TMath::Power(x, n);
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
  // calculate chi2 for both graphs
  for (int i=0; i<npar; i++){
    fparam1->SetParameter(i,par[i]);
    fparam2->SetParameter(i,par[i]);
  }

  // f = calcNLL(hdata,fparam);
  f = calcchi2(hdata1,fparam1)+calcchi2(hdata2,fparam2);
 
}


//-------------------------------------------------------------------------

int main(int argc, char **argv) {


  TCanvas* canvas = new TCanvas("canvas","",1200,600);
  canvas->Divide(2,1);
  canvas->cd(1);

  auto *f=new TFile("experiments.root");
  auto hexp1 = (TH1F*) f->Get("hexp1");
  auto hexp2 = (TH1F*) f->Get("hexp2");
  


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

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();
  // ***************************************

  
  //**************Double Gaussian**********************************//
  
  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 8;              // the number of parameters
 
  TMinuit minuit(npar);
  minuit.SetPrintLevel(0);         // set print level of minuit to normal
  minuit.SetFCN(fcn);              // the fcn to be minized (objective fcn) 
  // this calculates chi^2 or NLL, it is not the parameterization of the data!

  // Make a TF1, we can use this to fit the data and plot the results
  TF1* expbkg = new TF1("expbkg", gausExpBkg, hexp1->GetXaxis()->GetXmin(), hexp1->GetXaxis()->GetXmax(), npar);
  TF1* powbkg = new TF1("powbkg", gausPowBkg, hexp2->GetXaxis()->GetXmin(), hexp2->GetXaxis()->GetXmax(), npar);
  double par[npar];               // the free parameters
  double stepSize[npar];          // (limiting) step sizes for parameters
  double minVal[npar];            // minimum/maximum bounds on parameters 
  double maxVal[npar];            
  TString parName[npar];          // Optional names for nicer output!

  // Initial parameters MUST be set by some means to get things started
  par[0] = 40; // A1
  par[1] = 75; // mu
  par[2] = 5; // sigma
  par[3] = 1550;  // Abkg1
  par[4] = 15; // lambda
  par[5] = 20; // A2
  par[6] = 250000; // Abkg2
  par[7] = -2; // n

  stepSize[0] = TMath::Abs(par[0]*0.1);   // usually 10% of initial guess is OK for starting step size, YMMV
  stepSize[1] = TMath::Abs(par[1]*0.1);   // step sizes MUST be positive!
  stepSize[2] = TMath::Abs(par[2]*0.1);
  stepSize[3] = TMath::Abs(par[3]*0.1);
  stepSize[4] = TMath::Abs(par[4]*0.1);
  stepSize[5] = TMath::Abs(par[5]*0.1);
  stepSize[6] = TMath::Abs(par[6]*0.1);
  stepSize[7] = TMath::Abs(par[7]*0.1);

  minVal[0] = 10;      // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 70;
  minVal[1] = 70; 
  maxVal[1] = 80;
  minVal[2] = 0; 
  maxVal[2] = 10;
  minVal[3] = 0; 
  maxVal[3] = 0;
  minVal[4] = 5; 
  maxVal[4] = 30;
  minVal[5] = 10; 
  maxVal[5] = 50;
  minVal[6] = 0; 
  maxVal[6] = 0;
  minVal[7] = -3; 
  maxVal[7] = 0;

  parName[0] = "A1";       // let's give our fit parameters useful names
  parName[1] = "mu";
  parName[2] = "sigma";       
  parName[3] = "Abkg1";
  parName[4] = "lambda";       
  parName[5] = "A2";
  parName[6] = "Abkg2";
  parName[7] = "n";

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
  hdata1=hexp1;     // histogram to fit
  fparam1=expbkg;  // our model
  hdata2=hexp2;
  fparam2=powbkg;

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
  expbkg->SetParameters(outpar);
  powbkg->SetParameters(outpar);

  // myfunc->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
  // // myfunc->SetLineColor(1);             //  black (default)
  // myfunc->SetLineWidth(2);

  // myfunc->GetXaxis()->SetTitle("x");
  // myfunc->GetYaxis()->SetTitle("Double Gaussian");

  // summarize the fitting results
  cout << "\n==========================\n"<<endl;
  cout << "Results of double histogram chi sqr minimization"<< endl;
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

  
  int ndof = (hexp2->GetNbinsX()+hexp1->GetNbinsX())-npar; // combined dof
  double redChiSqr = fmin/ndof;
  double p = TMath::Prob(fmin, ndof);
  cout << "degrees of freedom = " << ndof << endl;
  cout << "reduced chi sqr = " << redChiSqr << endl;
  cout << "p value = " << p << endl;


  



  // plot

  // myfunc->SetTitle("Double Gaussian;;");
  hexp1->Draw("e");
  expbkg->Draw("same");
  

  canvas->cd(2);
  // myfuncG->SetTitle("Gumbel Distribution;;");
  hexp2->Draw("e");
  powbkg->Draw("same");
  

  canvas->Update();
  canvas->SaveAs("ex2.pdf");
  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;

  return 0;

}
