void computeEfficiencies(const std::string& regex) {
  std::cout << "Using '" << regex << "' to open files ..." << std::endl;
  std::string filenames = regex + ".root";

  std::cout << "Computing efficiencies and luminosity ..." << std::endl;

  double totalLumi = 0.;
  long long globalEvents = 0;
  long long globalPassedEvents = 0;

  TSystemDirectory dir("", ".");
  TList *files = dir.GetListOfFiles();
  TIter next(files);
  TRegexp compiledRegex(regex.c_str());
  TSystemFile* f = NULL;
  while (f = (TSystemFile *) next()) {
    TString name = f->GetName();
    if (name.Contains(compiledRegex)) {
      TFile * rootFile = TFile::Open(name);

      TParameter<double>* lumi = (TParameter<double>*) rootFile->Get("gammaJet/total_luminosity");
      TParameter<Long64_t>* totalEvents = (TParameter<Long64_t>*) rootFile->Get("gammaJet/total_events");
      TParameter<Long64_t>* passedEvents = (TParameter<Long64_t>*) rootFile->Get("gammaJet/passed_events");

      globalEvents += totalEvents->GetVal();
      globalPassedEvents += passedEvents->GetVal();

      if(lumi) {
        totalLumi += lumi->GetVal();
      } else
        std::cout << " WARNING! File '" << name << "' has no lumi information." << std::endl;

      std::cout << name << ": " << (double) passedEvents->GetVal() / (double) totalEvents->GetVal() * 100 << "%" << std::endl;

      rootFile->Close();
      delete rootFile;
    }
  }

  delete files;

  std::cout << "\tTotal lumi: " << totalLumi * 1e-9 << " fb-1" << std::endl;
  std::cout << "\Total efficiency: " << (double) globalPassedEvents / (double) globalEvents * 100 << "%" << std::endl;
}
