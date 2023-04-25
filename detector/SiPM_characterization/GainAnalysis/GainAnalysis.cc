// Copyright 2023 Andrea De Vita - Nicolò Salimbeni

//root libraries
#include "TDirectory.h"
#include "TF1.h"
#include "TH1F.h"
#include "TSpectrum.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSpectrum.h"
#include <TFile.h>
#include <TTree.h>

//cpp libraries
#include <fstream>
#include <ostream>
#include <string>
#include <stdint.h>
#include <vector>

//custom libraries
#include "AnUtil.h"
#include "Event.h"
#include "InfoAcq.h"
#include "InfoAcq.cc"


/////////////// FUNCTION DECLARATIONS /////////////////
void FillOutFile(   const TString root_file, 
                    const char *dirname = "C:/root/folder/",
                    const char *ext     = ".root");

void Find_peaks(    const TString file_name,    const TString histogram_name,
                    int n_expected_peaks,       int fit_width, 
                    bool print_out = false,     const TString file_txt = "peaks_stats.txt");

/////////////// MAIN FUNCTION /////////////////////////
void GetGain(   const TString root_file,  const char *dirname, const char *ext,
                const bool out_val){

    //Prima di tutto devo creare un root file con all'interno tutti gli istogrammi
    FillOutFile(root_file, dirname, ext);

    //definisco il nome del txt_file che magari vorrò riempire (se out_val == true)
    TString txt_file;
    if(out_val ==true){
        
        //chiedo nome del txt file da creare
        std::cout << "What is the name of the txt file: ";
        std::cin >> txt_file;

        //apro txt file e ci metto intestazione
        std::ofstream file(txt_file, std::ios::out);
        file << "V[dV]\tGain[mV]\tSigmaV[dV]\tSigmaGain[mv]" << std::endl;

        //chiudo txt file
        file.close();
    }

    //
    // A questo punto ho prodotto 2 file: un root file con gli istogrammi e un txt file con intestazione
    //

    //apro il root_file e guardo i suoi contenuti (gli istogrammi)
    TFile *file = TFile::Open(root_file,"READ");
    TList* list = file->GetListOfKeys();

    // loop sugli istogrammi in modo tale da ricavare il gain per ognuno di essi 
    // (se out_val == true stampiamo anche i risultati nel txt_file)
    TIter next(list);
        TKey* key;
        while ((key = (TKey*)next())) {
            TString histo_name = key->GetName();
            Find_peaks(root_file, histo_name,3,7,out_val,txt_file);
        }
}


//////////// MINOR FUNCTION DEFINITIONS ///////////////////
float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue) {
  uint16_t inputRanges[12] = {10,   20,   50,   100,   200,   500,
                              1000, 2000, 5000, 10000, 20000, 50000};

  return (raw * inputRanges[rangeIndex]) * 1. / maxADCValue;
}

TH1F* ReadTree(const char *fileName, bool negative, int channel = chB,
                Long64_t nEvtMax = -1) {
  // dichiaro le struct
  InfoAcq::chSettings       chSet1;
  InfoAcq::chSettings       chSet2;
  InfoAcq::samplingSettings sampSet;

  // dichiaro le variabili dell'evento
  unsigned long long ID;
  int                samplesStored;
  long long          triggerInstant;
  short              timeUnit;
  short             *sample0;
  short             *sample1;

  // apro il file in sola lettura
  TFile *input_file = new TFile(fileName, "READ");

  // leggo i trees
  TTree *treeCh   = (TTree *)input_file->Get("Channels");
  TTree *treeSamp = (TTree *)input_file->Get("SampSets");
  TTree *treeEvt  = (TTree *)input_file->Get("Event");
  // TFile->Get() restituisce un oggetto generico che va
  // convertito esplicitamente anteponendo (TTree*)

  // prelevo i branch con le info e li associo alle struct
  treeCh->SetBranchAddress("Ch1", &chSet1.enabled);
  treeCh->SetBranchAddress("Ch2", &chSet2.enabled);
  treeSamp->SetBranchAddress("Settings", &sampSet.max_adc_value);

  // leggo le entries
  // dichiaro l'oggetto InfoAcq e lo riempio
  InfoAcq *info = new InfoAcq();
  cout << "Riempio l'oggetto INFO\n";
  treeCh->GetEntry(0);
  treeSamp->GetEntry(0);
  info->FillSettings(&chSet1, &chSet2, &sampSet);

  InfoAcq::chSettings chSet;
  chSet = channel == chA ? chSet1 : chSet2;

  if (chSet.enabled == false) {
    cout << " channel " << channel << " not enabled, stopping" << endl;
  }
  // imposto i branches per gli eventi
  sample0 = new short[sampSet.samplesStoredPerEvent];
  sample1 = new short[sampSet.samplesStoredPerEvent];
  treeEvt->SetBranchAddress("ID", &ID);
  treeEvt->SetBranchAddress("nSamp", &samplesStored);
  treeEvt->SetBranchAddress("Instant", &triggerInstant);
  treeEvt->SetBranchAddress("TimeUnit", &timeUnit);
  treeEvt->SetBranchAddress("WaveformsA", &sample0[0]);
  treeEvt->SetBranchAddress("WaveformsB", &sample1[0]);

  Long64_t nEvt = treeEvt->GetEntries();
  if (nEvtMax >= 0 && nEvtMax < nEvt) {
    nEvt = nEvtMax;
  }
  float maximum = 0.0;
  float minimum = 0.0;

  // spettro in energia

  // float xmin= negative?
  // adc_to_mv(sampSet.max_adc_value,chSet.range,-1*sampSet.max_adc_value) : 0 ;
  // float xmax = negative? 0 :
  // adc_to_mv(sampSet.max_adc_value,chSet.range,sampSet.max_adc_value) ;
  float xmin = -100;
  float xmax = 0;
  TH1F *spectrumMaximum =
      new TH1F("hMax", "Maxima Distribution Spectrum", 200, xmin, xmax);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //	TH1F* spectrumMaximum = new TH1F( "hMax", "Maximum Spectrum", 256,
  //xmin,xmax );
  cerr << " Range   : " << xmax - xmin << " mV " << endl;
  cerr << " #Events : " << nEvt << endl;
  cerr << " #Samples: " << sampSet.samplesStoredPerEvent << endl;
  cerr << " Timestamp: " << sampSet.timeIntervalNanoseconds << " ns" << endl;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (Long64_t index = 0; index < nEvt; index++) {
    treeEvt->GetEntry(index);
    for (int ii = 0; ii < sampSet.samplesStoredPerEvent; ii++) {
      short sample = channel == chA ? sample0[ii] : sample1[ii];
      float value  = adc_to_mv(sample, chSet.range, sampSet.max_adc_value);
      if (value > maximum) maximum = value;
      if (value < minimum) minimum = value;
    }
    spectrumMaximum->Fill(negative ? minimum : maximum);
    maximum = 0.0;
    minimum = 0.0;
  }
  return(spectrumMaximum);
}

void Find_peaks(const TString file_name, const TString histogram_name,
                int n_expected_peaks, int fit_width, 
                bool print_out = false, const TString file_txt = "peaks_stats.txt") {
 
    // Open the TFile
    TFile* file = TFile::Open(file_name);

    // Read the histogram from the file
    TH1F* histogram = dynamic_cast<TH1F*>(file->Get(histogram_name));

    // Create a TSpectrum object to find peaks
    TSpectrum* spectrum = new TSpectrum(n_expected_peaks);

    // Find the peaks in the histogram using TSpectrum
    int n_peaks = spectrum->Search(histogram, 2, "", 0.01);

    // Get the positions of the peaks
    double* peak_positions = spectrum->GetPositionX();

    // We don't now the order of the peaks in the array
    // I order them from the lower to the upper
    AnUtil::Reorder(peak_positions, n_peaks);

    // Fit gaussians
    TF1*     gaussians[n_peaks];  // 4 gaussians declaration
    Double_t bin_width = histogram->GetXaxis()->GetBinWidth(1);

    for (int i = 0; i < n_peaks; i++) {
        // take mean, and range values to fit
        Double_t mean      = peak_positions[i];
        Double_t range_inf = mean - fit_width * bin_width;  // 10 bin blow the mean
        Double_t range_max = mean + fit_width * bin_width;  // 10 bin above the mean

        gaussians[i] = new TF1(Form("gaussian%d", i), "gaus", range_inf, range_max);
        gaussians[i]->SetParameter(1, mean);
        histogram->Fit(gaussians[i], "R+");
    }

    // I save all the values of the gain in a vector
    std::vector<Double_t> gains;
    for (int i = 0; i < n_peaks - 1; i++) {
        Double_t gain = gaussians[i + 1]->GetParameter(1) - gaussians[i]->GetParameter(1);
        // remember: 1 is the parameter corresponding to the mean

        gains.push_back(gain);
    }

    // print the mean gain
    Double_t mean_gain = AnUtil::Mean(gains, gains.size());
    Double_t rms_gain  = AnUtil::Rms(gains, gains.size());
    std::cout << "Gain mean is:\t" << mean_gain << std::endl;
    std::cout << "Gain rms is:\t" << rms_gain << std::endl;

    //se print_out è true salviamo il risultato in un txt file già esistente e già inizializzato
    if (print_out == true){
        std::ofstream file(file_txt, std::ios::app);
        file <<histogram_name(10,13)<<"\t"<< mean_gain<<"\t"<< "3"<<"\t"<<rms_gain<< std::endl;
        file.close();
  }

}

void FillOutFile(   const TString root_file, 
                    const char *dirname ,
                    const char *ext){
    

    TSystemDirectory         dir(dirname, dirname);
    TList                   *files = dir.GetListOfFiles();

    std::vector<std::string> file_list;

    if (files) {
        TSystemFile *file;
        TString      fname;
        TIter        next(files);
        for (TObject *obj : *files) {
            file  = reinterpret_cast<TSystemFile *>(obj);
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
            file_list.push_back(fname.Data());
            }
        }
    }

    std::sort(file_list.begin(), file_list.end());
    for (int j = 0; j < file_list.size(); j++) {

        TString file_path = dirname + file_list[j];
        TH1F* Hist=(TH1F *)ReadTree(file_path,true,chA,-1)->Clone("Hist");

        TFile *f = TFile::Open(root_file,"UPDATE");
        TString histname((file_path)(5,13));
        
        Hist->Write(histname);
        f->Close();

    }
}