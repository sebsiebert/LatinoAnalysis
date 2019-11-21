import ROOT
import os

from ZZMatrixElement.MELA.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class HMggFTaggerClass(Module):
    '''
    put this file in LatinoAnalysis/NanoGardener/python/modules/
    Add extra variables to NANO tree
    '''
    def __init__(self):
        self.LHCsqrts_ = 13.
        self.mh_ = 125

        cmssw_base = os.getenv('CMSSW_BASE')
        cmssw_arch = os.getenv('SCRAM_ARCH')
        # MELA
        # Stuff needed:
        # Directory ZZMatrixElement: https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git / https://twiki.cern.ch/twiki/bin/view/CMS/MELAProject
        # Directory MelaAnalytics: https://github.com/usarica/MelaAnalytics.git
        # ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/interface/")
        # ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/src/")
        # ROOT.gSystem.Load("libZZMatrixElementMELA.so")
        # ROOT.gSystem.Load("libMelaAnalyticsCandidateLOCaster.so")
        # ROOT.gSystem.Load(self.cmssw_base+"/src/ZZMatrixElement/MELA/data/"+cmssw_arch+"/libmcfm_705.so")

        # self.mela = ROOT.Mela(LHCsqrts_, mh_)
        self.mela = Mela(self.LHCsqrts_, self.mh_)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("IsggF", "O")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        IsBoosted   = getattr(event, 'IsBoosted')
        IsResolved  = getattr(event, 'IsBoosted')
        IsTopTagged = getattr(event, 'IsTopTagged')
        IsVbfFat    = getattr(event, 'IsVbfFat')
        IsVbfjj     = getattr(event, 'IsVbfjj')

        Lept_col           = Collection(event, 'Lepton')
        CFatJet_col        = Collection(event, 'CleanFatJet')
        CJet_col           = Collection(event, 'CleanJet')
        CleanJetNotFat_col = Collection(event, 'CleanJetNotFat')
        Jet_col            = Collection(event, 'Jet')

        IsggF = False

        daughter_coll   = SimpleParticleCollection_t()
        associated_coll = SimpleParticleCollection_t()

        if IsBoosted and IsVbfFat:
            pass

        elif IsResolved and IsVbfjj:
            pass

        self.out.fillBranch("IsggF", IsggF)
        return True

# define modules using the syntax
# 'name = lambda : constructor'
# to avoid having them loaded when not needed
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
