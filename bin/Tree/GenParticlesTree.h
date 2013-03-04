//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 28 16:00:10 2012 by ROOT version 5.32/00
// from TTree muons/muons tree
// found on file: output_mc.root
//////////////////////////////////////////////////////////

#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class GenParticlesTree {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    std::vector<int>     *pdg_id;
    std::vector<int>     *main_index;
    std::vector<int>     *status;
    std::vector<int>     *number_of_mothers;
    std::vector<int>     *index_of_first_mother;
    std::vector<int>     *index_of_second_mother;
    std::vector<float>   *energy;
    std::vector<float>   *pt;
    std::vector<float>   *px;
    std::vector<float>   *py;
    std::vector<float>   *pz;
    std::vector<float>   *eta;
    std::vector<float>   *phi;

    // List of branches
    TBranch        *b_pdg_id;   //!
    TBranch        *b_main_index;   //!
    TBranch        *b_status;   //!
    TBranch        *b_number_of_mothers;   //!
    TBranch        *b_index_of_first_mother;   //!
    TBranch        *b_index_of_second_mother;   //!
    TBranch        *b_energy;   //!
    TBranch        *b_pt;   //!
    TBranch        *b_px;   //!
    TBranch        *b_py;   //!
    TBranch        *b_pz;   //!
    TBranch        *b_eta;   //!
    TBranch        *b_phi;   //!

    GenParticlesTree();
    virtual ~GenParticlesTree();
    virtual Int_t    GetEntry(Long64_t entry);

    virtual void     Init(TTree *tree);
};

GenParticlesTree::GenParticlesTree() : fChain(0) 
{
}

GenParticlesTree::~GenParticlesTree()
{
  if (!fChain)
    return;

  delete fChain->GetCurrentFile();
}

Int_t GenParticlesTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain)
    return 0;

  return fChain->GetEntry(entry);
}

void GenParticlesTree::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree)
    return;

  // Set object pointer
  pdg_id = NULL;
  main_index = NULL;
  status = NULL;
  number_of_mothers = NULL;
  index_of_first_mother = NULL;
  index_of_second_mother = NULL;
  energy = NULL;
  pt = NULL;
  px = NULL;
  py = NULL;
  pz = NULL;
  eta = NULL;
  phi = NULL;

  fChain = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("pdg_id", &pdg_id, &b_pdg_id);
  fChain->SetBranchAddress("main_index", &main_index, &b_main_index);
  fChain->SetBranchAddress("status", &status, &b_status);
  fChain->SetBranchAddress("number_of_mothers", &number_of_mothers, &b_number_of_mothers);
  fChain->SetBranchAddress("index_of_first_mother", &index_of_first_mother, &b_index_of_first_mother);
  fChain->SetBranchAddress("index_of_second_mother", &index_of_second_mother, &b_index_of_second_mother);
  fChain->SetBranchAddress("energy", &energy, &b_energy);
  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("px", &px, &b_px);
  fChain->SetBranchAddress("py", &py, &b_py);
  fChain->SetBranchAddress("pz", &pz, &b_pz);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
}
