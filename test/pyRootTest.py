import ROOT
from DataFormats.FWLite import Events, Handle

events = Events('patTuple_PF2PAT.root')
handle = Handle('std::vector<pat::Photon>')
label = ("selectedPatPhotons")

for event in events:
    event.getByLabel(label,handle)


    photons = handle.product()
    for photon in photons:
        print photon.userFloat('eleNewEnergiesProducer:energySCEleJoshPhoSemiParamV5ecorr')
        #print photon.userFloat('eCorrections_:regression1Energy')

        
