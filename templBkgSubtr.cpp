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
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>

using namespace std;

// ---------- GLOBALS ----------
// Declare pointer to data as global
// The use of global variables is a disfavored coding practice, but
// in this case it will allow us to avoid writing more complex code
TH2F *hdataGlobal;
TH2F *hbkgGlobal;
// also passing the fit function, to make the objective fcn more generic
TF2 *fparam;

//-------------------------------------------------------------------------
// 2D gaussian plus background
double gaus2DBkg(double *xPtr, double par[]){        
  double x = xPtr[0];
  double y = xPtr[1];
  double A= par[0];         // normalization
  double mux = par[1];       // mean of x
  double sigmax = par[2];
  double muy = par[3];       // mean of x
  double sigmay = par[4];
  double Abkg = par[5];
  
  double peak = A*TMath::Exp(-0.5*(x-mux)*(x-mux)/sigmax/sigmax)*TMath::Exp(-0.5*(y-muy)*(y-muy)/sigmay/sigmay);
  return peak + Abkg * hbkgGlobal->GetBinContent(hbkgGlobal->FindBin(x,y));
}



//-------------------------------------------------------------------------
// This is our OBJECTIVE Function
// return NLL given a histogram and function

// double calcNLL(TH1F* h, TF1* f){
//   double nll=0;
//   for (int i=1; i<=h->GetNbinsX(); i++){
//     double x=h->GetBinCenter(i);
//     int n=(int)(h->GetBinContent(i));
//     double mu=f->Eval(x);
//     if (mu<1e-10) mu=1e-10;    // avoid log(0) problems if we go outside a reasonable range!
//     nll -= n * TMath::Log(mu) - mu  - TMath::LnGamma(n+1);
//   }
//   // cout << "nll "<< nll <<endl;
//   return 2*nll;   // factor of 2 so the 1 sigma error contours follow the chi^2 convention
// }


// return chi2 given a histogram and function

double calcchi2(TH2F* h, TF2* f){
  double chi2=0;
  for (int i=1; i<=h->GetNbinsX(); i++){
    for (int j=1; j<h->GetNbinsY(); j++) {
      double x=h->GetXaxis()->GetBinCenter(i); // bin center in x
      double y=h->GetYaxis()->GetBinCenter(j); // bin center in y
      int n=(int)(h->GetBinContent(i,j)); // y bin
      double mu=f->Eval(x,y); // expected y bin
      double err=h->GetBinError(i,j); // bin error

      if (err==0) continue; // avoid division by zero

      double diff = mu-n;
      chi2 += diff*diff/(err*err);
    }
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
    fparam->SetParameter(i,par[i]);
  }

  // f = calcNLL(hdata,fparam);
  f = calcchi2(hdataGlobal,fparam);
 
}




//-------------------------------------------------------------------------

int main(int argc, char **argv) {


  TCanvas* canvas = new TCanvas("canvas","",1200,1200);
  canvas->Divide(2,2);
  canvas->cd(1);

  auto *f=new TFile("fitInputs.root");
  auto hdata = (TH2F*) f->Get("hdata");
  auto hbkg = (TH2F*) f->Get("hbkg");
  


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


  
  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 6;              // the number of parameters
 
  TMinuit minuit(npar);
  minuit.SetPrintLevel(0);         // set print level of minuit to normal
  minuit.SetFCN(fcn);              // the fcn to be minized (objective fcn) 
  // this calculates chi^2 or NLL, it is not the parameterization of the data!

  // Make a TF1, we can use this to fit the data and plot the results
  TF2* myfunc = new TF2("myfunc", gaus2DBkg, 0., 6., 0., 6., npar);
  myfunc->SetNpx(hdata->GetNbinsX());
  myfunc->SetNpy(hdata->GetNbinsY());
  double par[npar];               // the free parameters
  double stepSize[npar];          // (limiting) step sizes for parameters
  double minVal[npar];            // minimum/maximum bounds on parameters 
  double maxVal[npar];            
  TString parName[npar];          // Optional names for nicer output!

  // Initial parameters MUST be set by some means to get things started
  par[0] = 60; // A
  par[1] = 3.5; // mux
  par[2] = 0.5; // sigmax
  par[3] = 1.5;  // muy
  par[4] = 1; // sigmay
  par[5] = 0.25; // Abkg

  stepSize[0] = TMath::Abs(par[0]*0.1);   // usually 10% of initial guess is OK for starting step size, YMMV
  stepSize[1] = TMath::Abs(par[1]*0.1);   // step sizes MUST be positive!
  stepSize[2] = TMath::Abs(par[2]*0.1);
  stepSize[3] = TMath::Abs(par[3]*0.1);
  stepSize[4] = TMath::Abs(par[4]*0.1);
  stepSize[5] = TMath::Abs(par[5]*0.1);

  minVal[0] = 30;      // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 100;
  minVal[1] = 3; 
  maxVal[1] = 4;
  minVal[2] = 0; 
  maxVal[2] = 1;
  minVal[3] = 0; 
  maxVal[3] = 3;
  minVal[4] = 0; 
  maxVal[4] = 2;
  minVal[5] = 0; 
  maxVal[5] = 0.75;
 
  parName[0] = "A";       // let's give our fit parameters useful names
  parName[1] = "mux";
  parName[2] = "sigmax";       
  parName[3] = "muy";
  parName[4] = "sigmay";       
  parName[5] = "Abkg";

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
  hdataGlobal=hdata;     // 2D histogram to fit
  hbkgGlobal=hbkg; // 2D histogram model of background
  fparam=myfunc;

  // // Do the minimization!
  minuit.Migrad();       // Minuit's best minimization algorithm

  // Get the result
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
    minuit.GetParameter(i,outpar[i],err[i]);
  }
  // // // run a minos error analysis
  // // // this perfoms a brute force scan around the minimum of the objective function
  // // // to get better estimates of the errors on the fit parameters
  // // // gMinuit->mnmnos();

  // store the fit parameters in our TF1, decorate our function
  myfunc->SetParameters(outpar);

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

  
  int ndof = (hdata->GetNbinsX()*hdata->GetNbinsY())-npar; // combined dof
  double redChiSqr = fmin/ndof;
  double p = TMath::Prob(fmin, ndof);
  cout << "degrees of freedom = " << ndof << endl;
  cout << "reduced chi sqr = " << redChiSqr << endl;
  cout << "p value = " << p << endl;
  


  // Fill a TH2F with the fitted result
  const int nbinsx = hdata->GetNbinsX();
  const int nbinsy = hdata->GetNbinsY();
  TH2F *h2 = new TH2F("h2", "fit result;;;", nbinsx, 0, 6, nbinsy, 0, 6);

  // Fill histogram with TF2 values
  for (int i = 1; i <= h2->GetNbinsX(); i++) {
      for (int j = 1; j <= h2->GetNbinsY(); j++) {
          double x = h2->GetXaxis()->GetBinCenter(i);
          double y = h2->GetYaxis()->GetBinCenter(j);
          double z = myfunc->Eval(x, y);
          h2->SetBinContent(i, j, z);
      }
  }

  // TH2F for residuals
  TH2F *hres = new TH2F("hres", "residuals;;;", nbinsx, 0, 6, nbinsy, 0, 6);

  // Fill histogram with residuals (data-fit)
  for (int i = 1; i <= hres->GetNbinsX(); i++) {
      for (int j = 1; j <= hres->GetNbinsY(); j++) {
          double x = hres->GetXaxis()->GetBinCenter(i);
          double y = hres->GetYaxis()->GetBinCenter(j);
          double z = hdata->GetBinContent(hdata->FindBin(x,y)) - h2->GetBinContent(i,j);
          hres->SetBinContent(i, j, z);
      }
  }

  // TH2F for signal minus background
  TH2F *hsig = new TH2F("hsig", "data minus background;;;", nbinsx, 0, 6, nbinsy, 0, 6);

  // Fill histogram with signal minus background fit
  int nSig = 0;
  double nSigErr = 0;
  for (int i = 1; i <= hsig->GetNbinsX(); i++) {

      for (int j = 1; j <= hsig->GetNbinsY(); j++) {
        double x = hsig->GetXaxis()->GetBinCenter(i);
        double y = hsig->GetYaxis()->GetBinCenter(j);
        double bkgVal = hbkg->GetBinContent(hbkg->FindBin(x,y));
        double dataVal = hdata->GetBinContent(hdata->FindBin(x,y));
        double z = dataVal - outpar[5]* bkgVal;

        nSig += z;
        nSigErr += dataVal + TMath::Power(bkgVal*err[5], 2); // error of each bin will be sqrt(sum of error in data bin sqrd and error in bkg fit sqrd) -> error in total counts will be sqrt the sum of these sqrd

        if (z<0) z=0; // remove negative bins for plotting histogram
        hsig->SetBinContent(i, j, z);
      }
  }

  nSigErr = TMath::Power(nSigErr, 0.5);

  cout << "signal count = " << nSig << endl;
  cout << "signal count error = " << nSigErr << endl;

  // plot

  // myfunc->SetTitle("Double Gaussian;;");
  hdata->Draw("lego");
  

  canvas->cd(2);
  // myfuncG->SetTitle("Gumbel Distribution;;");
  // hbkg->Draw("lego");
  h2->Draw("lego");

  canvas->cd(3);
  // myfunc->Draw("lego");
  hres->Draw("colz");
  

  canvas->cd(4);
  // myfunc->Draw("lego");
  hsig->Draw("lego");
  

  canvas->Update();
  canvas->SaveAs("ex3.pdf");
  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;

  return 0;

}
