//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 28 15:43:50 2012 by ROOT version 5.32/00
// from TTree photon/photon tree
// found on file: output_mc.root
//////////////////////////////////////////////////////////

#ifndef BaseTree_h
#define BaseTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "BaseTree.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTree: public BaseTree {
  public :

    // Declaration of leaf types
    Float_t         jet_area;

    // List of branches
    TBranch        *b_jet_area;   //!

    JetTree();

    virtual void     Init(TTree *tree);
};

#endif

JetTree::JetTree():
  BaseTree::BaseTree()
{
}

void JetTree::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (! tree)
    return;

  BaseTree::Init(tree);

  fChain->SetBranchAddress("jet_area", &jet_area, &b_jet_area);
}
