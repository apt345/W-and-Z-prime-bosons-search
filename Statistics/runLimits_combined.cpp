#include "FYS5555/Statistics/MC_Bayes.cpp+"

{

  Int_t Nmasses;
  Int_t Nchannels;
  Int_t mass;
  Double_t dummy;

  ifstream input("FYS5555/Statistics/inputs_bothChannels.txt");
  
  input >> Nmasses;
  input >> Nchannels;

  Double_t* masses = new Double_t[Nmasses];

  Double_t** bkg = new Double_t*[Nmasses];
  for (Int_t i = 0; i < Nmasses; i++)
    bkg[i] = new Double_t[Nchannels];
  
  Double_t** bkgUncertainty = new Double_t*[Nmasses];
  for (Int_t i = 0; i < Nmasses; i++)
    bkgUncertainty[i] = new Double_t[Nchannels];
  
  Double_t** efficiency = new Double_t*[Nmasses];
  for (Int_t i = 0; i < Nmasses; i++)
    efficiency[i] = new Double_t[Nchannels];
  
  Double_t** efficiencyUncertainty = new Double_t*[Nmasses];
  for (Int_t i = 0; i < Nmasses; i++)
    efficiencyUncertainty[i] = new Double_t[Nchannels];
  
  Double_t** Nsignal = new Double_t*[Nmasses];
  for (Int_t i = 0; i < Nmasses; i++)
    Nsignal[i] = new Double_t[Nchannels];
  
  Int_t** Nobs = new Int_t*[Nmasses];
  for (Int_t i = 0; i < Nmasses; i++)
    Nobs[i] = new Int_t[Nchannels];


  for (mass = 0; mass < Nmasses; mass++) {
    input >> masses[mass];
    input >> dummy;
    
    for (Int_t chan = 0; chan < Nchannels; chan++) {
      input >> efficiency[mass][chan];
      input >> efficiencyUncertainty[mass][chan];
      input >> Nsignal[mass][chan]; 
      input >> bkg[mass][chan];
      input >> bkgUncertainty[mass][chan];
      input >> Nobs[mass][chan];
    }
  }


  ofstream outputFile("FYS5555/Statistics/limitschannelcomb.txt");

  for (mass = 0; mass < Nmasses; mass++) {

    cout << "\n----------------------------\n";
    cout << "mass = " << masses[mass] << " GeV:" << endl;
    
    Double_t intLum = 10.6;
    Double_t intLumUncertainty = 0.4;

    cout << "\nInputs:\n";
    cout << "Int. luminosity = (" << intLum << " +/- " << intLumUncertainty << ")/fb\n";
    
      
    Double_t* background = new Double_t[Nchannels];
    for (Int_t chan = 0; chan < Nchannels; chan++)
      background[chan] = bkg[mass][chan]/intLum;
    
    Double_t* backgroundUncertainty = new Double_t[Nchannels];
    for (Int_t chan = 0; chan < Nchannels; chan++)
      backgroundUncertainty[chan] = bkgUncertainty[mass][chan]/intLum;
    
    Int_t* Nobserved = new Int_t[Nchannels];
    for (Int_t chan = 0; chan < Nchannels; chan++)
      Nobserved[chan] = Nobs[mass][chan];

    Double_t* signalEfficiency = new Double_t[Nchannels];
    for (Int_t chan = 0; chan < Nchannels; chan++)
      signalEfficiency[chan] = efficiency[mass][chan];
    
    Double_t* signalEfficiencyUncertainty = new Double_t[Nchannels];
    for (Int_t chan = 0; chan < Nchannels; chan++)
      signalEfficiencyUncertainty[chan] = efficiencyUncertainty[mass][chan];

    
    for (Int_t chan = 0; chan < Nchannels; chan++) {
      cout << "Channel " << chan << ":\n";
      cout << "   Background = " << bkg[mass][chan] << " +/- " << bkgUncertainty[mass][chan] << endl;
      cout << "   Observed events = " << Nobs[mass][chan] << endl;
      cout << "   Signal efficiency = " << efficiency[mass][chan] << " +/- " << efficiencyUncertainty[mass][chan] << endl;
    }

	  
    Int_t Nmc = 5000;


    MC_Bayes* mcb = new MC_Bayes(Nchannels, intLum, intLumUncertainty, background, backgroundUncertainty, Nobserved, Nmc, signalEfficiency, signalEfficiencyUncertainty, 
				   false,true,true,true);

    Double_t step = 1.e100;
    for (Int_t chan = 0; chan < Nchannels; chan++) {
      Double_t approxLimitEvts = 1.6 * TMath::Sqrt(bkgUncertainty[mass][chan]*bkgUncertainty[mass][chan] + bkg[mass][chan]);
      if (approxLimitEvts < 3.0)
	approxLimitEvts = 3.0;
      Double_t thisStep = approxLimitEvts/(efficiency[mass][chan]*intLum)/500.0;
      if (thisStep < step)
	step = thisStep;
    }
    
    Double_t limit = mcb->excludedSignal(step, false);
	
    cout << "\nObserved limit: " << limit << " fb" << endl;
    
    cout << "Observed limit / step = " << limit/step << endl;


    Int_t Npseudos = 1000; //number of background-only pseudo-experiments
	
    //The observed limit runs very fast, but when we repeat this calculation many times for the green and yellow
    //bands, running time can be long. Let's try some accuracy tuning to get things to run in reasonable time.
    //You can comment this out to run with "full" precision.
    step /= 2.0;
    mcb->Nmc = 1500;
    
    Bool_t highStat = false;
    Bool_t lowStat = true;
    for (Int_t chan = 0; chan < Nchannels; chan++) {
      if (bkg[mass][chan] > 2500.)
	highStat = true;
      if (bkg[mass][chan] > 5.)
	lowStat = false;
    }

    //For very high background levels, might be better to run fewer pseudos with higher precision
    if (highStat) {
      Npseudos /= 2;
      mcb->Nmc *= 2.5;
    }

    //For very low background levels, there are few different outcomes, and we can crank up the number of pseudos
    if (lowStat)
      Npseudos = 10000;
    

    Double_t* exp = mcb->expectedExclusion(step,Npseudos,false);

    cout << "\nExpected limit and bands [fb]:\n";
    cout << "-2sigma     -1sigma     median    +1sigma     +2sigma\n";
    cout << exp[0] << "     " << exp[1] << "     " << exp[2] << "     " << exp[3] << "     " << exp[4] << endl;

    
    outputFile << masses[mass] << " " << Nsignal[mass][0]/(efficiency[mass][0]*intLum) << " " << limit;
    for (Int_t k = 0; k < 5; k++)
      outputFile << " " << exp[k];
    outputFile << endl;

    delete mcb;
    delete[] exp;
  }
  cout << endl;

  outputFile.close();

  delete[] masses;
  delete[] bkg;
  delete[] bkgUncertainty;
  delete[] efficiency;
  delete[] efficiencyUncertainty;
  delete[] Nsignal;
  delete[] Nobs;
}
