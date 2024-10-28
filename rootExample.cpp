// A simple C++ program to illustrate the use of some ROOT classes
// in your code

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TRandom2.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>

using namespace std;

void prettyPlot(TCanvas *canvas);
void rootExample();

int main(int argc, char **argv) {

  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this line can be removed
  TApplication theApp("App", &argc, argv);

  // optional:  here all the "root stuff" is put into its own function
  // this allows us to either compile the program in C++ as a stand alone
  // exe or to use the function in ROOT directly, eg.  root -l rootExample.cpp
  // If we build a shared library, we can call the function from Python as well
  rootExample();

  
  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program
  theApp.Run(true);
  return 0;
}

// here we create and display two histograms filled wit hrandom data
void rootExample(){
  TCanvas* canvas = new TCanvas();
  prettyPlot(canvas);

  // Generate some data in a histogram 
  TRandom2 r(0);

  // create a histogram and fill it w/ randomly distributed data
  // based on our function
  double xmin = 0.0;
  double xmax = 5.0;
  TH1F *hexp=new TH1F("hexp","exponential distribution;x;events",100,xmin,xmax);
  TH1F *hgaus=new TH1F("hgaus","normal distribution;x;events",100,xmin,xmax);

  for (int i=0; i<5000; i++){
    hexp->Fill(r.Exp(3.14159));
    hgaus->Fill(r.Gaus(2.5,0.5));
  }
  
  canvas->Divide(2,1);
  canvas->cd(1)->SetLogy();
  hexp->Draw();
  canvas->cd(2);
  hgaus->Draw();

}


void prettyPlot(TCanvas *canvas){
  // ***************************************
  // Set a bunch of parameters to make the plot look nice
  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);        
  canvas->SetFrameBorderMode(0);  
  gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();
  gROOT->ForceStyle();

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();
  // ***************************************
}
