
#include <TH1.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>

TH1* fitAndReplace(TH1* h, Double_t fitLow, Double_t fitHigh, Double_t stichPoint, Int_t func, Double_t highEnd, TF1** funcpp = NULL) {

  TH1* h_orig = (TH1*)h->Clone();
  
  for (Int_t bin = 0; bin < h->GetNbinsX()+1; bin++) {
    h->SetBinContent(bin,h->GetBinContent(bin)/h->GetBinWidth(bin));
    h->SetBinError(bin,h->GetBinError(bin)/h->GetBinWidth(bin));
    h->SetYTitle("Events/GeV");
  }

  
  TF1* fitFunc;

  if (func == 0) {

    fitFunc = new TF1("dijet","TMath::Exp(-1.0*[0])*TMath::Power(x,[1])*TMath::Power(x,[2]*TMath::Log(x))",fitLow,highEnd);
    fitFunc->SetParameter(0,7.6);
    fitFunc->SetParameter(1,2.1);
    fitFunc->SetParameter(2,-0.5);
  }
  else if (func == 1) {

    fitFunc = new TF1("modpower","[0]/pow(x+[1],[2])",fitLow,highEnd);
    fitFunc->SetParameter(0,6.8e14);
    fitFunc->SetParameter(1,1e2);
    fitFunc->SetParameter(2,5.2);
  }
  else if (func == 2) {

    fitFunc = new TF1("dilepton","([0]*(2.4952/(pow((91.1876-x),2)+pow(2.4952,2)))*pow((1-(x/13000)),[1])*pow((x/13000),([2]+[3]*log(x/13000)+[4]*pow(log(x/13000),2)+[5]*pow(log(x/13000),3))))",fitLow,highEnd);
    fitFunc->SetParameter(1,1.5);
    fitFunc->SetParameter(2,-12.38);
    fitFunc->SetParameter(3,-4.295);
    fitFunc->SetParameter(4,-.9191);
    fitFunc->SetParameter(5,-0.0845);
    //fitFunc->SetParameter(0,1780000./fitFunc->Integral(220,6000));
  }
  else {

    cout << "UNKNOWN FIT FUNCTION!\n";
    exit(1);
  }


  if (fitLow > stichPoint) {

    cout << "ERROR: fitLow > stichPoint\n";
    exit(1);
  }


  TVirtualFitter::Fitter(h)->SetMaxIterations(10000000);
    
  TFitResultPtr fitresult = h->Fit(fitFunc,"S","",fitLow,fitHigh);

  Int_t fitStatus = (Int_t)fitresult;
  
  if (fitStatus != 0) {

    cout << "FIT FAILED WITH STATUS " << fitStatus << endl;
    cout << fitLow << "  " << fitHigh << "  " << func << endl;
    exit(1);
  }


  Int_t Nbins = h_orig->GetNbinsX();
  
  Double_t* mTlows = new Double_t[Nbins];
  for (Int_t i = 1; i <= Nbins; i++)
    mTlows[i-1] = h_orig->GetBinLowEdge(i);
  

  Double_t* mTups = new Double_t[Nbins];
  for (Int_t i = 1; i <= Nbins; i++)
    mTups[i-1] = h_orig->GetBinLowEdge(i) + h_orig->GetBinWidth(i);


  for (Int_t bin = 1; bin <= Nbins; bin++) {

    Double_t fitInt = fitFunc->Integral(mTlows[bin-1],mTups[bin-1]);
    Double_t statErr = fitFunc->IntegralError(mTlows[bin-1],mTups[bin-1],fitresult->GetParams(), 
					      fitresult->GetCovarianceMatrix().GetMatrixArray());

    if (mTlows[bin-1] > stichPoint) {
      h_orig->SetBinContent(bin,fitInt);
      h_orig->SetBinError(bin,statErr);
    }
  }

  if (funcpp != NULL)
    *funcpp = fitFunc;

  delete[] mTlows;
  delete[] mTups;

  delete h;

  return h_orig;
}
