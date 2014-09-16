// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <DataFormats/PatCandidates/interface/Photon.h>
#include "DataFormats/PatCandidates/interface/MET.h"

#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

using namespace std;

class PhotonIsolationProducer : public edm::EDProducer {
  public:
    explicit PhotonIsolationProducer(const edm::ParameterSet&);
    ~PhotonIsolationProducer();

  private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data --------------------------
    edm::InputTag src_;

    // Photon ID
    PFIsolationEstimator mPFIsolator;
};

PhotonIsolationProducer::PhotonIsolationProducer(const edm::ParameterSet& iConfig)
{
  src_= iConfig.getParameter<edm::InputTag>("src");

  mPFIsolator.initializePhotonIsolation(true);
  mPFIsolator.setConeSize(0.3);

  produces<edm::ValueMap<double>>("chargedHadronsIsolation");
  produces<edm::ValueMap<double>>("photonIsolation");
  produces<edm::ValueMap<double>>("neutralHadronsIsolation");
  produces<edm::ValueMap<bool>>("hasMatchedPromptElectron");
  //Footprint
  produces<edm::ValueMap<double>>("footchargediso");
  produces<edm::ValueMap<double>>("footchargedisoprimvtx");
  produces<edm::ValueMap<double>>("footneutraliso");
  produces<edm::ValueMap<double>>("footphotoniso");
  produces<edm::ValueMap<double>>("footchargedisorcone");
  produces<edm::ValueMap<double>>("footchargedisoprimvtxrcone");  
  produces<edm::ValueMap<double>>("footneutralisorcone");
  produces<edm::ValueMap<double>>("footphotonisorcone");
  produces<edm::ValueMap<double>>("footetarcone");
  produces<edm::ValueMap<double>>("footphircone");
  produces<edm::ValueMap<bool>>("footrconeisOK");
  //  produces<edm::ValueMap<std::vector<int>>>("footpfcandindexfootprint");
 produces<edm::ValueMap<double>>("footprintPx");
 produces<edm::ValueMap<double>>("footprintPy");
 produces<edm::ValueMap<double>>("footprintMExCorr");
 produces<edm::ValueMap<double>>("footprintMEyCorr");
//
 produces<double>("footprintMExraw");
 produces<double>("footprintMEyraw");
}


PhotonIsolationProducer::~PhotonIsolationProducer()
{

}

// ------------ method called to produce the data  ------------
void PhotonIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<edm::ValueMap<double>> chIsoMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler chFiller(*chIsoMap);

  std::auto_ptr<edm::ValueMap<double>> phIsoMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler phFiller(*phIsoMap);

  std::auto_ptr<edm::ValueMap<double>> nhIsoMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler nhFiller(*nhIsoMap);

  std::auto_ptr<edm::ValueMap<bool>> promptConvMap(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler promptConvFiller(*promptConvMap);

  //foot

  std::auto_ptr<edm::ValueMap<double>> footChIsoMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footchIsoFiller(*footChIsoMap);

  std::auto_ptr<edm::ValueMap<double>> footChIsoPVMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footchIsoPVFiller(*footChIsoPVMap);

  std::auto_ptr<edm::ValueMap<double>> footNIsoMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footnIsoFiller(*footNIsoMap);

  std::auto_ptr<edm::ValueMap<double>> footPhIsoMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footphIsoFiller(*footPhIsoMap);

  std::auto_ptr<edm::ValueMap<double>> footChIsoRconeMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footchIsoRconeFiller(*footChIsoRconeMap);

  std::auto_ptr<edm::ValueMap<double>> footChIsoPVRconeMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footchIsoPVRconeFiller(*footChIsoPVRconeMap);

  std::auto_ptr<edm::ValueMap<double>> footNIsoRconeMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footnIsoRconeFiller(*footNIsoRconeMap);

  std::auto_ptr<edm::ValueMap<double>> footPhIsoRconeMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footphIsoRconeFiller(*footPhIsoRconeMap);

  std::auto_ptr<edm::ValueMap<double>> footEtaRconeMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footEtaRconeFiller(*footEtaRconeMap);

  std::auto_ptr<edm::ValueMap<double>> footPhiRconeMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler footPhiRconeFiller(*footPhiRconeMap);

  std::auto_ptr<edm::ValueMap<bool>> footRconeOKMap(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler footRconeOKFiller(*footRconeOKMap);

  std::auto_ptr<edm::ValueMap<double>> FootPxMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler FootPxFiller(*FootPxMap);

  std::auto_ptr<edm::ValueMap<double>> FootPyMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler FootPyFiller(*FootPyMap);

  std::auto_ptr<edm::ValueMap<double>> FootprintMExMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler FootprintMExFiller(*FootprintMExMap);

  std::auto_ptr<edm::ValueMap<double>> FootprintMEyMap(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler FootprintMEyFiller(*FootprintMEyMap);

  std::auto_ptr<double> MExraw(new double(0.));
  std::auto_ptr<double> MEyraw(new double(0.));
/*
std::auto_ptr<float> px(new float(0.));
std::auto_ptr<float> py(new float(0.));
std::auto_ptr<float> FootprintMEx(new float(0.));
std::auto_ptr<float> FootprintMEy(new float(0.));
*/
  /*  std::auto_ptr<edm::ValueMap<std::vector<int>>> footPFcandIdxFPrintMap(new edm::ValueMap<std::vector<int>>());
      edm::ValueMap<std::vector<int>>::Filler footPFcIdxFPFiller(*footPFcandIdxFPrintMap);*/

  edm::Handle<pat::PhotonCollection> photonsHandle;
  iEvent.getByLabel(src_, photonsHandle);


  edm::Handle<pat::METCollection> rawMets;
//  iEvent.getByLabel(std::string("patPFMet" + ((*it == "AK5Calo") ? "" : *it)), rawMets);
  iEvent.getByLabel(std::string("patPFMetPFlowAK5chs"), rawMets);
  const pat::MET& rawMet = rawMets->at(0);

 

 
  // First, conversion safe electron veto
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle;

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);

  edm::Handle<reco::GsfElectronCollection> hElectrons;
  iEvent.getByLabel("gsfElectrons", hElectrons);

  edm::Handle<reco::PFCandidateCollection> hPFCandidates;
  iEvent.getByLabel("particleFlow", hPFCandidates);
  const reco::PFCandidateCollection& pfCandidates = *hPFCandidates;

  edm::Handle<reco::VertexCollection>  vertexCollection;
  iEvent.getByLabel("goodOfflinePrimaryVertices", vertexCollection);
  reco::VertexRef vertexRef(vertexCollection, 0);



  if (vertexCollection->empty())
    return;

  std::vector<double> chIsoValues;
  std::vector<double> phIsoValues;
  std::vector<double> nhIsoValues;
  std::vector<bool>   promptConvValues;
  chIsoValues.reserve(photonsHandle->size());
  phIsoValues.reserve(photonsHandle->size());
  nhIsoValues.reserve(photonsHandle->size());
  promptConvValues.reserve(photonsHandle->size());

  //foot
  std::vector<double> footChIsovalues;
  std::vector<double> footChIsoPVValues;
  std::vector<double> footNIsoValues;
  std::vector<double> footPhIsoValues;
  std::vector<double> footChIsoRconeValues;
  std::vector<double> footChIsoPVRconeValues;
  std::vector<double> footNIsoRconeValues;
  std::vector<double> footPhIsoRconeValues;
  std::vector<double> footEtaRconeValues;
  std::vector<double> footPhiRconeValues;
  std::vector<bool> footRconeOKValues;
  std::vector<double> footPxValues;
  std::vector<double> footPyValues;
  std::vector<double> footFootprintMExValues;
  std::vector<double> footFootprintMEyValues;

  //  std::vector<std::vector<int>>  footPFcandIdxFPrint;
  footChIsovalues.reserve(photonsHandle->size());
  footChIsoPVValues.reserve(photonsHandle->size());
  footNIsoValues.reserve(photonsHandle->size());
  footPhIsoValues.reserve(photonsHandle->size());
  footChIsoRconeValues.reserve(photonsHandle->size());
  footChIsoPVRconeValues.reserve(photonsHandle->size());
  footNIsoRconeValues.reserve(photonsHandle->size());
  footPhIsoRconeValues.reserve(photonsHandle->size());
  footEtaRconeValues.reserve(photonsHandle->size());
  footPhiRconeValues.reserve(photonsHandle->size());
  footRconeOKValues.reserve(photonsHandle->size());
  footPxValues.reserve(photonsHandle->size());
  footPyValues.reserve(photonsHandle->size());
  footFootprintMExValues.reserve(photonsHandle->size());
  footFootprintMEyValues.reserve(photonsHandle->size());

 //  footPFcandIdxFPrint.reserve(photonsHandle->size());

  SuperClusterFootprintRemoval remover(iEvent,iSetup);
  PFIsolation_struct PFIso_struct;

  double px;
  double py;
  double FootprintMEx;
  double FootprintMEy;

/*
  edm::EventID eventId = iEvent.id();
  int event = eventId.event();
  int run = eventId.run();
  int lumiBlock = eventId.luminosityBlock();
*/

  pat::PhotonCollection::const_iterator it = photonsHandle->begin();
  for (; it != photonsHandle->end(); ++it) {

  *MExraw=0.;
  *MEyraw=0.;
  px=0.;
  py=0.;
  FootprintMEx=0.;
  FootprintMEy=0.;
 //
   pat::Photon photon = *it;
    reco::SuperClusterRef scref = photon.superCluster();

    promptConvValues.push_back(ConversionTools::hasMatchedPromptElectron(photon.superCluster(), hElectrons, hConversions, beamspot.position()));

    mPFIsolator.fGetIsolation(&photon, &pfCandidates, vertexRef, vertexCollection);

    chIsoValues.push_back(mPFIsolator.getIsolationCharged());
    phIsoValues.push_back(mPFIsolator.getIsolationPhoton());
    nhIsoValues.push_back(mPFIsolator.getIsolationNeutral());

   PFIso_struct = remover.PFIsolation(scref,edm::Ptr<Vertex>(vertexCollection,0));

   footChIsovalues.push_back(PFIso_struct.chargediso);
   footChIsoPVValues.push_back(PFIso_struct.chargediso_primvtx);
   footNIsoValues.push_back(PFIso_struct.neutraliso);
   footPhIsoValues.push_back(PFIso_struct.photoniso);
   footChIsoRconeValues.push_back(PFIso_struct.chargediso_rcone);
   footChIsoPVRconeValues.push_back(PFIso_struct.chargediso_primvtx_rcone);
   footNIsoRconeValues.push_back(PFIso_struct.neutraliso_rcone);
   footPhIsoRconeValues.push_back(PFIso_struct.photoniso_rcone);
   footEtaRconeValues.push_back(PFIso_struct.eta_rcone);
   footPhiRconeValues.push_back(PFIso_struct.phi_rcone);
   footRconeOKValues.push_back(PFIso_struct.rcone_isOK);
// initialize alternative calculation of raw met

//calculate pfcandidates to remove from met calculation
//already stored in the footprint
std::set<int> PFCandsToRemove;
for (unsigned int i=0; i<PFIso_struct.pfcandindex_footprint.size(); i++) PFCandsToRemove.insert(PFIso_struct.pfcandindex_footprint.at(i)); 
//look if there is an overlapping photon in the pfCands list
reco::SuperCluster* phsc = (reco::SuperCluster*)(photon.superCluster()).get();
reco::SuperCluster* pfCsc=0;
reco::SuperCluster* elsc=0;
edm::Ptr<reco::CaloCluster> pfCcc;
edm::Ptr<reco::CaloCluster> phcc = phsc->seed();
for(unsigned int i=0; i<pfCandidates.size(); ++i ) {
   if ((pfCandidates)[i].particleId()==reco::PFCandidate::gamma){
   if ((pfCandidates)[i].mva_nothing_gamma()>0.){
   if(((pfCandidates)[i].superClusterRef()).isNonnull()) {
     pfCsc = (reco::SuperCluster*)((pfCandidates)[i].superClusterRef()).get();
     pfCcc = pfCsc->seed();
//       if(sqrt(pow((pfCcc->eta()-phcc->eta()),2)+pow((pfCcc->phi()-phcc->phi()),2))<0.001 && fabs(pfCsc->rawEnergy()-phsc->rawEnergy())<0.001) {
       if( fabs(pfCsc->rawEnergy())/cosh(pfCcc->eta())>20. && sqrt(pow((pfCcc->eta()-phcc->eta()),2)+pow((pfCcc->phi()-phcc->phi()),2))<0.2 && fabs((pfCsc->rawEnergy()-phsc->rawEnergy())/phsc->rawEnergy())<0.5) {
        PFCandsToRemove.insert(i);
       }
      }//sc nonnull
    }
   }
  }

//look if there is an overlapping electron
bool foundEgSC = false;
edm::Ptr<GsfElectron> elPtrSl;
reco::GsfElectron  selected_electron;
reco::GsfElectronCollection::const_iterator elit = hElectrons->begin();
 for (; elit != hElectrons->end(); ++elit) {
   reco::GsfElectron electron = *elit;
  if(electron.superCluster().isNonnull()) {
     elsc = (reco::SuperCluster*)(electron.superCluster()).get();
     edm::Ptr<reco::CaloCluster>  elcc = elsc->seed();
//    if(sqrt(pow((elcc->eta()-phcc->eta()),2)+pow((elcc->phi()-phcc->phi()),2))<0.001 && (elsc->rawEnergy()-phsc->rawEnergy())<0.001){
    if( fabs(elsc->rawEnergy())/cosh(elcc->eta())>20. && sqrt(pow((elcc->eta()-phcc->eta()),2)+pow((elcc->phi()-phcc->phi()),2))<0.2 && fabs((elsc->rawEnergy()-phsc->rawEnergy())/phsc->rawEnergy())<0.50){
     selected_electron=electron;
     foundEgSC=true;
     break;
     }
    }//sc nonnull
    }
//now look for it in the pfcands list
if (foundEgSC){
 double MVACut_ = -1.;
 for(unsigned int i=0; i< pfCandidates.size(); ++i ) {
  if ((pfCandidates)[i].particleId()==reco::PFCandidate::e && (pfCandidates)[i].gsfTrackRef().isNull()==false && (pfCandidates)[i].mva_e_pi()>MVACut_){ 
   if((pfCandidates)[i].gsfTrackRef()==selected_electron.gsfTrack()) {
  PFCandsToRemove.insert(i);
   }
  }
 }
}

//re-loop over all pfcandidates
for(unsigned int i=0; i<pfCandidates.size(); ++i ) {
//fill alternative calculation of raw met
*MExraw += -1.*(pfCandidates)[i].px();
*MEyraw += -1.*(pfCandidates)[i].py();
bool UseThisPfC=true;
for(unsigned int v: PFCandsToRemove) {
if(i==v) UseThisPfC=false;
}
if(UseThisPfC) {
//calculatio of FootPrintMET
FootprintMEx += -1.*(pfCandidates)[i].px();
FootprintMEy += -1.*(pfCandidates)[i].py();
}
}

//now loop over candidates to remove and sum up their px and py
for(int v: PFCandsToRemove) {
px+=(pfCandidates)[v].px();
py+=(pfCandidates)[v].py();
}
//fill alternative met calc variables
  footPxValues.push_back(px);
  footPyValues.push_back(py);
  footFootprintMExValues.push_back(FootprintMEx);
  footFootprintMEyValues.push_back(FootprintMEy);
}


  chFiller.insert(photonsHandle, chIsoValues.begin(), chIsoValues.end());
  chFiller.fill();
  
  phFiller.insert(photonsHandle, phIsoValues.begin(), phIsoValues.end());
  phFiller.fill();
  
  nhFiller.insert(photonsHandle, nhIsoValues.begin(), nhIsoValues.end());
  nhFiller.fill();
  
  promptConvFiller.insert(photonsHandle, promptConvValues.begin(), promptConvValues.end());
  promptConvFiller.fill();

  footchIsoFiller.insert(photonsHandle, footChIsovalues.begin(), footChIsovalues.end());
  footchIsoFiller.fill();

  footchIsoPVFiller.insert(photonsHandle, footChIsovalues.begin(), footChIsovalues.end());
  footchIsoPVFiller.fill();

  footnIsoFiller.insert(photonsHandle, footNIsoValues.begin(), footNIsoValues.end());
  footnIsoFiller.fill();

  footphIsoFiller.insert(photonsHandle, footPhIsoValues.begin(), footPhIsoValues.end());
  footphIsoFiller.fill();

  footchIsoRconeFiller.insert(photonsHandle, footChIsoRconeValues.begin(), footChIsoRconeValues.end());
  footchIsoRconeFiller.fill();

  footchIsoPVRconeFiller.insert(photonsHandle, footChIsoPVRconeValues.begin(), footChIsoPVRconeValues.end());
  footchIsoPVRconeFiller.fill();

  footnIsoRconeFiller.insert(photonsHandle, footNIsoRconeValues.begin(), footNIsoRconeValues.end());
  footnIsoRconeFiller.fill();

  footphIsoRconeFiller.insert(photonsHandle, footPhIsoRconeValues.begin(), footPhIsoRconeValues.end());
  footphIsoRconeFiller.fill();

  footEtaRconeFiller.insert(photonsHandle, footEtaRconeValues.begin(), footEtaRconeValues.end());
  footEtaRconeFiller.fill();

  footPhiRconeFiller.insert(photonsHandle, footPhiRconeValues.begin(), footPhiRconeValues.end());
  footPhiRconeFiller.fill();

  footRconeOKFiller.insert(photonsHandle, footRconeOKValues.begin(), footRconeOKValues.end());
  footRconeOKFiller.fill();


  FootPxFiller.insert(photonsHandle, footPxValues.begin(), footPxValues.end());
  FootPxFiller.fill();

  FootPyFiller.insert(photonsHandle, footPyValues.begin(), footPyValues.end());
  FootPyFiller.fill();

  FootprintMExFiller.insert(photonsHandle, footFootprintMExValues.begin(), footFootprintMExValues.end());
  FootprintMExFiller.fill();

  FootprintMEyFiller.insert(photonsHandle, footFootprintMEyValues.begin(), footFootprintMEyValues.end());
  FootprintMEyFiller.fill();


 
  iEvent.put(chIsoMap, "chargedHadronsIsolation");
  iEvent.put(phIsoMap, "photonIsolation");
  iEvent.put(nhIsoMap, "neutralHadronsIsolation");
  iEvent.put(promptConvMap, "hasMatchedPromptElectron");

  iEvent.put(footChIsoMap, "footchargediso");
  iEvent.put(footChIsoPVMap, "footchargedisoprimvtx");
  iEvent.put(footNIsoMap, "footneutraliso");
  iEvent.put(footPhIsoMap, "footphotoniso");
  iEvent.put(footChIsoRconeMap, "footchargedisorcone");
  iEvent.put(footChIsoPVRconeMap, "footchargedisoprimvtxrcone");
  iEvent.put(footNIsoRconeMap, "footneutralisorcone");
  iEvent.put(footPhIsoRconeMap, "footphotonisorcone");
  iEvent.put(footEtaRconeMap, "footetarcone");
  iEvent.put(footPhiRconeMap, "footphircone");
  iEvent.put(footRconeOKMap, "footrconeisOK");
  iEvent.put(FootPxMap, "footprintPx");
  iEvent.put(FootPyMap, "footprintPy");
  iEvent.put(MExraw, "footprintMExraw");
  iEvent.put(MEyraw, "footprintMEyraw");
  iEvent.put(FootprintMExMap, "footprintMExCorr");
  iEvent.put(FootprintMEyMap, "footprintMEyCorr");
}

// ------------ method called once each job just before starting event loop  ------------
void
PhotonIsolationProducer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhotonIsolationProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonIsolationProducer);
