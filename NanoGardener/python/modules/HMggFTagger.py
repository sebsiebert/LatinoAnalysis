import ROOT
import os
import math

from ZZMatrixElement.MELA.mela import Mela, SimpleParticleCollection_t, TVar
SimpleParticle_t = ROOT.SimpleParticle_t

ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class HMggFTaggerClass(Module):
    '''
    put this file in LatinoAnalysis/NanoGardener/python/modules/
    Add extra variables to NANO tree
    '''
    def __init__(self):
        self.LHCsqrts = 13.
        self.mh = 125
        self.mela_con = 0.01
        verbosity = ROOT.TVar.SILENT

        cmssw_base = os.getenv('CMSSW_BASE')
        cmssw_arch = os.getenv('SCRAM_ARCH')
        # MELA
        # Stuff needed:
        # Directory ZZMatrixElement: https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git / https://twiki.cern.ch/twiki/bin/view/CMS/MELAProject
        # Directory MelaAnalytics: https://github.com/usarica/MelaAnalytics.git
        ROOT.gSystem.AddIncludePath("-I"+cmssw_base+"/interface/")
        ROOT.gSystem.AddIncludePath("-I"+cmssw_base+"/src/")
        ROOT.gSystem.Load("libZZMatrixElementMELA.so")
        ROOT.gSystem.Load("libMelaAnalyticsCandidateLOCaster.so")
        ROOT.gSystem.Load(cmssw_base+"/src/ZZMatrixElement/MELA/data/"+cmssw_arch+"/libmcfm_705.so")

        # self.mela = ROOT.Mela(LHCsqrts_, mh_)
        self.mela = Mela(self.LHCsqrts, self.mh, verbosity)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ggF_MELA_var", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        IsBoosted  = getattr(event, 'IsBoosted')
        IsResolved = getattr(event, 'IsBoosted')
        IsVbfFat   = getattr(event, 'IsVbfFat')
        IsVbfjj    = getattr(event, 'IsVbfjj')

        Lept_col    = Collection(event, 'Lepton')
        CJet_col    = Collection(event, 'CleanJet')
        CFatJet_col = Collection(event, 'CleanFatJet')
        WJet_col    = Collection(event, 'CleanFatJetPassMBoosted')

        Wlep_pt_Puppi   = getattr(event, "Wlep_pt_Puppi")
        Wlep_eta_Puppi  = getattr(event, "Wlep_eta_Puppi")
        Wlep_phi_Puppi  = getattr(event, "Wlep_phi_Puppi")
        Wlep_mass_Puppi = getattr(event, "Wlep_mass_Puppi")

        Wjj_ClJet0_idx  = getattr(event, "idx_j1")
        Wjj_ClJet1_idx  = getattr(event, "idx_j2")

        IsggF = False



        daughter_coll   = SimpleParticleCollection_t()
        associated_coll = SimpleParticleCollection_t()

        # daughter particles
        #######################
        ## leptonic part
        #######################
        chargedLeptonID = Lept_col[0]['pdgId']
        chargedLeptonLV = ROOT.TLorentzVector()
        chargedLeptonLV.SetPtEtaPhiM(
                Lept_col[0]['pt'],
                Lept_col[0]['eta'],
                Lept_col[0]['phi'],
                0.0005 if abs(chargedLeptonID)==11 else 0.1057
        )
        chargedLepton = SimpleParticle_t(chargedLeptonID, chargedLeptonLV)
        neutralLeptonID = chargedLeptonID
        neutralLeptonID += neutralLeptonID / abs(neutralLeptonID) # increase abs(ID)
        neutralLeptonLV = ROOT.TLorentzVector()
        neutralLeptonLV.SetPtEtaPhiM(
                Wlep_pt_Puppi,
                Wlep_eta_Puppi,
                Wlep_phi_Puppi,
                Wlep_mass_Puppi
        )
        neutralLeptonLV = neutralLeptonLV - chargedLeptonLV
        neutralLepton = SimpleParticle_t(neutralLeptonID, neutralLeptonLV)

        daughter_coll.push_back(chargedLepton)
        daughter_coll.push_back(neutralLepton)

        # fill remaining leptons in associated_coll
        for i in range(2, len(Lept_col)):
            lepID = Lept_col[i]['pdgId']
            lepLV = ROOT.TLorentzVector()
            lepLV.SetPtEtaPhiM(
                    Lept_col[i]['pt'],
                    Lept_col[i]['eta'],
                    Lept_col[i]['phi'],
                    0.0005 if abs(chargedLeptonID)==11 else 0.1057
            )
            associated_coll.push_back(SimpleParticle_t(chargedLeptonID, chargedLeptonLV))

        ######################
        ## hadronic part
        ######################
        if IsBoosted and IsVbfFat:
            jetID = 0
            jetLV = ROOT.TLorentzVector()
            jetLV.SetPtEtaPhiM(
                    WJet_col[0]['pt'],
                    WJet_col[0]['eta'],
                    WJet_col[0]['phi'],
                    WJet_col[0]['mass']
            )
            daughter_coll.push_back(SimpleParticle_t(jetID, jetLV))

            # fill all CleanJetNotFat into associated_coll
            for i in range(2, len(CFatJet_col)):
                if i == WJet_col[0]['CFatJetIdx']: continue # don't add W jet twice
                jetLV = ROOT.TLorentzVector()
                jetLV.SetPtEtaPhiM(
                        CFatJet_col[i]['pt'],
                        CFatJet_col[i]['eta'],
                        CFatJet_col[i]['phi'],
                        CFatJet_col[i]['mass']
                )
                associated_coll.push_back(SimpleParticle_t(0, jetLV))


        elif IsResolved and IsVbfjj:
            jetID = 0
            jet0LV = ROOT.TLorentzVector()
            jet0LV.SetPtEtaPhiM(
                    CJet_col[Wjj_ClJet0_idx]['pt'],
                    CJet_col[Wjj_ClJet0_idx]['eta'],
                    CJet_col[Wjj_ClJet0_idx]['phi'],
                    CJet_col[Wjj_ClJet0_idx]['mass']
            )
            jet1LV = ROOT.TLorentzVector()
            jet1LV.SetPtEtaPhiM(
                    CJet_col[Wjj_ClJet1_idx]['pt'],
                    CJet_col[Wjj_ClJet1_idx]['eta'],
                    CJet_col[Wjj_ClJet1_idx]['phi'],
                    CJet_col[Wjj_ClJet1_idx]['mass']
            )
            daughter_coll.push_back(SimpleParticle_t(jetID, jet0LV))
            daughter_coll.push_back(SimpleParticle_t(jetID, jet1LV))

            # fill all other CleanJet into associated_coll
            for i in range(2, len(CJet_col)):
                if i == WWjj_ClJet0_idx or i == WWjj_ClJet1_idx: continue # don't add twice
                jetLV = ROOT.TLorentzVector()
                jetLV.SetPtEtaPhiM(
                        CJet_col[i]['pt'],
                        CJet_col[i]['eta'],
                        CJet_col[i]['phi'],
                        CJet_col[i]['mass']
                )
                associated_coll.push_back(SimpleParticle_t(0, jetLV))



        # compute discriminator
        self.mela.setCandidateDecayMode(TVar.CandidateDecay_WW)
        self.mela.setInputEvent(daughter_coll, associated_coll, 0, 0)
        self.mela.setCurrentCandidateFromIndex(0)

        self.mela.setProcess(TVar.bkgWW, TVar.MCFM, TVar.ZZQQB)
        p_qqWW = self.mela.computeP(False)
        self.mela.setProcess(TVar.HSMHiggs, TVar.MCFM, TVar.ZZGG)
        p_ggH = self.mela.computeP(False)


        melaDiscriminant = 1/(1+self.mela_con*(p_qqWW/p_ggH)) if p_ggH > 0 else 0
        self.mela.resetInputEvent()


        self.out.fillBranch("ggF_MELA_var", melaDiscriminant)
        return True

# to avoid having them loaded when not needed, define modules using the syntax:
# 'name = lambda : constructor'
HMggFTagger = lambda : HMggFTaggerClass()

"""
float mala_con = .01;

// calculate ME probabiities in event loop

SimpleParticleCollection_t daughter_coll, associated_coll;

// jets/subjets (WJ1 and WJ2) ordered by Pt
daughter_coll.push_back(SimpleParticle_t(lepid, *Lep));
daughter_coll.push_back(SimpleParticle_t(nuid, *Nu));
daughter_coll.push_back(SimpleParticle_t(0, *WJ1));
daughter_coll.push_back(SimpleParticle_t(0, *WJ2));

mela->setCandidateDecayMode(TVar::CandidateDecay_WW);
mela->setInputEvent(&daughter_coll, &associated_coll, 0, 0);
mela->setCurrentCandidateFromIndex(0);

mela->setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZQQB);
mela->computeP(p_qqWW, false);
mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
mela->computeP(p_ggH, false);

mela_var = 1/(1+mela_con*(p_qqWW/p_ggH));

mela->resetInputEvent();
"""
