# instead of doing a delta R matching of the last copies of dark particles
# to the jet, we want to match the constituents of the dark jets to the last
# dark descendants of these last copies of dark particles.
# Something to investigate: why are there no SM jets when we replace the
# Delta R matching with constituent matching.
import ROOT as rt
import inspect
import sys
import numpy as np
import pandas as pd
from array import array
from DataFormats.FWLite import Events, Handle
import tdrstyle
import os
import matplotlib.pyplot as plt

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()
rt.gROOT.ProcessLine(".L loader.C+")

DarkMediatorID_ = [4900001,4900002,4900003,4900004,4900005,4900006] # 4900023 is schannel Zprime
DarkQuarkID_ = [4900101,4900102]
DarkHadronID_ = [4900111,4900113,4900211,4900213]
DarkGluonID_ = [4900021]
DarkStableID_ = [51,52,53]
SMQuarkID_ = [1,2,3,4,5,6,7,8]
pTCut = 100
coneSize = 0.8
dpFractCut = 0.7
AllDarkID = DarkMediatorID_ + DarkQuarkID_ + DarkHadronID_ + DarkGluonID_

def deltaPhi(jetphiL,metPhiL):
    dphi = jetphiL - metPhiL
    if dphi < -np.pi:
        dphi += 2*np.pi
    if dphi > np.pi:
        dphi -= 2*np.pi
    return abs(dphi)

def deltaR(jet1,jet2):
    phi1 = jet1.phi()
    phi2 = jet2.phi()
    eta1 = jet1.eta()
    eta2 = jet2.eta()
    dp = abs(phi1 - phi2)
    if dp > np.pi:
        dp -= 2*np.pi
    deltaR2 = (eta1 - eta2) * (eta1 - eta2) + dp * dp
    return np.sqrt(deltaR2)

def darkTruth(darklist,pId):
    dt = False
    for dark in darklist:
        if abs(pId) == dark:
            dt = True
    return dt

def medDecay(fdau,fmom,fdQPartFM,fdGPartFM,fSMqPart):
    for mddi in range(fdau.numberOfDaughters()):
        mdd = fdau.daughter(mddi)
        mddId = mdd.pdgId()
        if darkTruth(DarkMediatorID_,mddId):
            medDecay(mdd,fmom,fdQPartFM,fdGPartFM,fSMqPart)
        elif darkTruth(DarkQuarkID_,mddId):
            fdQPartFM.append([mdd,fmom])
        elif darkTruth(DarkGluonID_,mddId):
            fdGPartFM.append([mdd,fmom])
        elif darkTruth(SMQuarkID_,mddId):
            fSMqPart.append([mdd,fmom])

def firstDark(partid,part_i,fdMPart,fdQPart,fdGPart,fdQPartFM,fdGPartFM,fSMqPart):
    DarkIDs = DarkMediatorID_ + DarkQuarkID_ + DarkGluonID_
    if darkTruth(DarkIDs,partid):
        fmom = part_i.mother(0)
        momId = abs(fmom.pdgId())
        if darkTruth(DarkIDs,momId):
            firstDark(momId,fmom,fdMPart,fdQPart,fdGPart,fdQPartFM,fdGPartFM,fSMqPart)
        else:
            for dau_i in range(part_i.numberOfDaughters()):
                fdau = part_i.daughter(dau_i)
                fdauId = fdau.pdgId()
                if darkTruth(DarkMediatorID_,fdauId):
                    fdMPart.append(fdau)
                    medDecay(fdau,fdau,fdQPartFM,fdGPartFM,fSMqPart)
                if darkTruth(DarkQuarkID_,fdauId):
                    fdQPart.append(fdau)
                if darkTruth(DarkGluonID_,fdauId):
                    fdGPart.append(fdau)

def lastDau(lastD,dHad): # find the last daughters
    for idH in range(dHad.numberOfDaughters()):
        dauH = dHad.daughter(idH)
        dauHID = abs(dauH.pdgId())
        if dauH.numberOfDaughters() == 0 and not darkTruth(DarkStableID_,dauHID): # last particle that is not stable dark hadron
            lastD.append(dauH)
        elif dauH.numberOfDaughters() > 0:
            lastDau(lastD,dauH)

def lastDarkGQ(DarkIDList,partid,part_i,lastD):
    if darkTruth(DarkIDList,partid) and part_i.isLastCopy():
        lastDau(lastD,part_i)

def checkDark(ijet,lastD,idValue):
    jCatID = 0
    for i in range(ijet.numberOfDaughters()):
        if ijet.daughter(i) in lastD:
            jCatID = idValue
            break
    return jCatID

def darkPt(darkMList): # calculate the total dark pT
    darkMomTLSum = rt.TLorentzVector(0,0,0,0) # sum of the TLorentzVectors of the dark mothers
    for imom in darkMList:
        darkMomTLSum += rt.TLorentzVector(imom.px(),imom.py(),imom.pz(),imom.energy())
    return darkMomTLSum.Pt()

def storTLV(H,jets):
    H.clear()
    for jet in np.array(jets):
        H.push_back(rt.TLorentzVector(jet.px(),jet.py(),jet.pz(),jet.energy()))

def isfj(fPart,uj,idValue):
    isJ = 0
    for fP in fPart:
        if deltaR(fP,uj) < coneSize:
            isJ = idValue
            break
    return isJ

def jetSub_Calc(jetsList,Jgirth,JaxisMj,JaxisMn,JptD):
    for ijt in range(len(jetsList)):
        i_jet = jetsList[ijt]
        if i_jet.pt() < 200:
            continue
        momentGirth = 0.0
        sumPt = 0.0; sumPt2 = 0.0
        sumDeta = 0.0; sumDphi = 0.0; sumDeta2 = 0.0; sumDphi2 = 0.0; sumDetaDphi = 0.0
        for k in range(i_jet.numberOfDaughters()):
            i_part = i_jet.daughterPtr(k)
            if(i_part.isNonnull() and i_part.isAvailable()):
                dphi = deltaPhi(i_jet.phi(),i_part.phi())
                deta = i_part.eta() - i_jet.eta()
                dR = deltaR(i_jet,i_part)
                pT = i_part.pt()
                pT2 = pT*pT

                momentGirth += pT*dR
                sumPt += pT
                sumPt2 += pT2
                sumDeta += deta*pT2
                sumDphi += dphi*pT2
                sumDeta2 += deta*deta*pT2
                sumDphi2 += dphi*dphi*pT2
                sumDetaDphi += deta*dphi*pT2

        # finish axis calculations (eigenvectors)
        sumDeta /= sumPt2
        sumDphi /= sumPt2
        sumDeta2 /= sumPt2
        sumDphi2 /= sumPt2
        sumDetaDphi /= sumPt2
        a = 0.0; b = 0.0; c = 0.0; d = 0.0
        a = sumDeta2 - sumDeta*sumDeta
        b = sumDphi2 - sumDphi*sumDphi
        c = sumDeta*sumDphi - sumDetaDphi
        d = np.sqrt(abs((a-b)*(a-b)+4*c*c))
        if (a+b+d)>0:
            axis1 = np.sqrt(0.5*(a+b+d))
        else:
            axis1 = 0

        if (a+b-d)>0:
            axis2 = np.sqrt(0.5*(a+b-d))
        else:
            axis2 = 0.0

        # store values
        girthvl = momentGirth/i_jet.pt()
        axisMjvl = axis1
        axisMnvl = axis2
        ptDvl = np.sqrt(sumPt2)/sumPt

        Jgirth[ijt] = girthvl
        JaxisMj[ijt] = axisMjvl
        JaxisMn[ijt] = axisMnvl
        JptD[ijt] = ptDvl

def darkPtFract(adj,lastD,stableD):
    tot4vec = rt.TLorentzVector(0.,0.,0.,0.)
    dark4vec = rt.TLorentzVector(0.,0.,0.,0.)
    anum = adj.numberOfDaughters()
    for ani in range(anum):
        adj_dau = adj.daughter(ani)
        prt = rt.TLorentzVector(adj_dau.px(),adj_dau.py(),adj_dau.pz(),adj_dau.energy())
        tot4vec += prt
        if adj_dau in lastD+stableD:
            dark4vec += prt
# accounting for stable dark hadrons
    jetDarkF = dark4vec.Pt()/tot4vec.Pt();
    return jetDarkF

def nConstPerJet(JetsList,nPBranch):
    for bj in range(len(JetsList)):
        bjet = JetsList[bj]
        npjet = bjet.numberOfDaughters()

        nPBranch[bj] = npjet

def uniqueList(alist):
    uList = []
    for a in alist:
        if a not in uList:
            uList.append(a)

    return np.array(uList)

# args, in order mMed, mDark, Rinv, Alpha
mMed, mDark, rInv, Alpha, MadGraph = sys.argv[1:]
MadGraph = bool(int(MadGraph))
# example, 3000 20 0.3 peak 0
if not MadGraph:
    name = "s1_mMed-{}_mDark-{}_rinv-{}_alpha-{}_Pythia_newMT2.root".format(mMed,mDark,rInv,Alpha)
else:
    name = "s1_mMed-{}_mDark-{}_rinv-{}_alpha-{}_newMT2.root".format(mMed,mDark,rInv,Alpha)
outf = rt.TFile( "Skims/" + name, 'RECREATE' )
cdGMA = outf.mkdir("GenMassAnalyzer")
cdGMA.cd()

tr = rt.TTree( 'tree', 'GenMassAnalyzer Tree' )

nfdMPart = array( 'i', [ 0 ] )
nAK8Jets = array( 'i', [ 0 ] )

AK8_Jets = rt.std.vector(rt.TLorentzVector)()
AK8_MedPair1 = rt.std.vector(rt.TLorentzVector)()
AK8_MedPair2 = rt.std.vector(rt.TLorentzVector)()
AK8_MT2dQM1 = rt.std.vector(rt.TLorentzVector)()
AK8_MT2SMM1 = rt.std.vector(rt.TLorentzVector)()
AK8_MT2dQM2 = rt.std.vector(rt.TLorentzVector)()
AK8_MT2SMM2 = rt.std.vector(rt.TLorentzVector)()
AK8_MTPair = rt.std.vector(rt.TLorentzVector)()
AK8_MTdQM = rt.std.vector(rt.TLorentzVector)()
AK8_MTSMM = rt.std.vector(rt.TLorentzVector)()
MET = rt.std.vector(rt.TLorentzVector)()
# jet substructure
Nmax = 50 # maximum number of jet in an event
AK8_girth = array( 'd', Nmax*[ 0. ] )
AK8_axisMj = array( 'd', Nmax*[ 0. ] )
AK8_axisMn = array( 'd', Nmax*[ 0. ] )
AK8_ptD = array( 'd', Nmax*[ 0. ] )
# dark pT fraction
svj_t_jetCat = array( 'i', Nmax*[ 0 ] )
AK8_dptFrac = array( 'd', Nmax*[ 0. ] )
# number of particles for each jet
npAK8Jets = array( 'i', Nmax*[ 0 ] )

tr.Branch( 'nAK8Jets', nAK8Jets, 'nAK8Jets/I' )
tr.Branch("AK8Jets","vector<TLorentzVector>",AK8_Jets)
tr.Branch("AK8MedPair1","vector<TLorentzVector>",AK8_MedPair1)
tr.Branch("AK8MedPair2","vector<TLorentzVector>",AK8_MedPair2)
tr.Branch("AK8MT2dQM1","vector<TLorentzVector>",AK8_MT2dQM1)
tr.Branch("AK8MT2SMM1","vector<TLorentzVector>",AK8_MT2SMM1)
tr.Branch("AK8MT2dQM2","vector<TLorentzVector>",AK8_MT2dQM2)
tr.Branch("AK8MT2SMM2","vector<TLorentzVector>",AK8_MT2SMM2)
tr.Branch("AK8MTPair","vector<TLorentzVector>",AK8_MTPair)
tr.Branch("AK8MTdQM","vector<TLorentzVector>",AK8_MTdQM)
tr.Branch("AK8MTSMM","vector<TLorentzVector>",AK8_MTSMM)
tr.Branch("MET","vector<TLorentzVector>",MET)
# jet substructure
tr.Branch( 'AK8_girth', AK8_girth, 'girth[nAK8Jets]/D' )
tr.Branch( 'AK8_axisMj', AK8_axisMj, 'axisMj[nAK8Jets]/D' )
tr.Branch( 'AK8_axisMn', AK8_axisMn, 'axisMn[nAK8Jets]/D' )
tr.Branch( 'AK8_ptD', AK8_ptD, 'ptD[nAK8Jets]/D' )
# dark pT fraction
tr.Branch( 'svj_t_jetCat', svj_t_jetCat, 'svj_t_jetCat[nAK8Jets]/I' )
tr.Branch( 'AK8_dptFrac', AK8_dptFrac, 'AK8_dptFrac[nAK8Jets]/D' )
# number of particles in each jet
tr.Branch( 'npAK8Jets', npAK8Jets, 'npAK8Jets[nAK8Jets]/I' )

pbegin = 1    # run all: 1
pend = 101    # run all: 101
# MadGraph
if MadGraph == True:
    fileStart = "root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020/GEN/" # output directory of step1_GEN
    listOfFiles = [fileStart+"step1_LHE-GEN_t-channel_mMed-{}_mDark-{}_rinv-{}_alpha-{}_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_part-{}.root".format(mMed,mDark,rInv,Alpha,iPart) for iPart in range(pbegin,pend)]#101
# Pythia Only
else:
    fileStart = "root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020/PythiaOnly/GEN/" # output directory of step1_GEN
    listOfFiles = [fileStart+"step1_GEN_t-channel_mMed-{}_mDark-{}_rinv-{}_alpha-{}_yukawa-1_13TeV-pythia8_n-1000_part-{}.root".format(mMed,mDark,rInv,Alpha,iPart) for iPart in range(pbegin,pend)]#101

# testing purposes (only run over first file)
# fileStart = "root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020/GEN/" # output directory of step1_GEN
# listOfFiles = [fileStart+"step1_LHE-GEN_t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_part-1.root"]

events = Events ( listOfFiles )

handleJets = Handle("std::vector<reco::GenJet>")
handleMET = Handle("std::vector<reco::GenMET>")
handleParts = Handle("std::vector<reco::GenParticle>")

labelJets = ("ak8GenJets")
labelMET = ("genMetTrue")
labelPart = ("genParticles")

count = 0
MedEvent = 0 # number of events with mediators
stabOLEvent = 0 # number of events where more than one jet shares the same stable dark hadron.

numdJ = 0

totJets = 0
totDJets = 0
totSMJets = 0
nMT2Event = 0
MT2Event_MedJet = 0

for event in events:
    count += 1
    # if count == 1000: # which event number to stop the code
    #     break
    if count % 1000 == 0:
        print "Event number " + str(count)
    event.getByLabel(labelJets,handleJets)
    event.getByLabel(labelMET,handleMET)
    event.getByLabel(labelPart,handleParts)
    jets = handleJets.product()
    hMET = handleMET.product()
    hPart = handleParts.product()

    lastD = [] # last daughters of all the last dark hadrons

    jCatList = []
    fdQPart = [] # first dark quarks
    fdGPart = [] # first dark gluons
    fdMPart = [] # first dark mediators
    fdQPartFM = [] # first dark quarks from the dark mediator
    fdGPartFM = [] # first dark gluons from the dark mediator
    fSMqPart = []  # first quarks from the dark mediator decay
    AK8Jets = []
    stableD = [] # stable dark hadrons

    if len(jets) < 0:
        continue

    part_ind = 0
    for part_i in handleParts.product():
        part_ind += 1
        partid = abs(part_i.pdgId())
        lastDarkGQ(DarkHadronID_,partid,part_i,lastD)
        firstDark(partid,part_i,fdMPart,fdQPart,fdGPart,fdQPartFM,fdGPartFM,fSMqPart)
        if darkTruth(DarkStableID_,partid):
            stableD.append(part_i)

    # getting rid of duplicates; tried using np.unique(), but it gives very weird results for fdQPart
    fdQPart = uniqueList(fdQPart)
    fdGPart = uniqueList(fdGPart)
    fdMPart = uniqueList(fdMPart)
    fdQPartFM = uniqueList(fdQPartFM)
    fdGPartFM = uniqueList(fdGPartFM)
    fSMqPart = uniqueList(fSMqPart)

# finding jets that contain descendants from the last dark hadrons, these should be dark jets
    i = 0
    for i in range(len(jets)):
        ijet = jets[i]
        if ijet.pt() > pTCut:
            jCatLabel = 0
            AK8Jets.append(ijet)
            totJets += 1 # test
            jCatLabel += checkDark(ijet,lastD,1)
            jCatLabel += isfj(fdQPart,ijet,2)
            jCatLabel += isfj(fdGPart,ijet,4)
            if len(fdQPartFM) > 0:
                jCatLabel += isfj(fdQPartFM[:,0],ijet,8)
            if len(fSMqPart) > 0:
                jCatLabel += isfj(fSMqPart[:,0],ijet,16)
            if jCatLabel %2 != 0:
                totDJets += 1
            elif jCatLabel %2 == 0:
                totSMJets += 1
            jCatList.append(jCatLabel)
            svj_t_jetCat[i] = jCatLabel
            AK8_dptFrac[i] = darkPtFract(ijet,lastD,stableD)
    nAK8Jets[0] = len(AK8Jets)
    AK8Jets = np.array(AK8Jets)
    jCatList = np.array(jCatList)
#
    storTLV(AK8_Jets,AK8Jets)
    storTLV(MET,[hMET[0]])
    jetSub_Calc(AK8Jets,AK8_girth,AK8_axisMj,AK8_axisMn,AK8_ptD)
    nConstPerJet(AK8Jets,npAK8Jets)

    # figuring out which pair of jets came from 1 mediator; useful for MT comparison
    dQMJs = []
    SMMJs = []
    if len(fdMPart) == 1:
        dQM = fdQPartFM[fdQPartFM[:,1] == fdMPart[0]].flatten()[0] # dark quark from mediator
        SMM = fSMqPart[fSMqPart[:,1] == fdMPart[0]].flatten()[0] # SM quark from mediator
        manyParticlesPerJet = False
        for dj in AK8Jets:
            nPartPerJet = 0
            if deltaR(dj, dQM)<0.8:
                dQMJs.append(dj)
                nPartPerJet += 1
            if deltaR(dj, SMM)<0.8:
                SMMJs.append(dj)
                nPartPerJet += 1
            if nPartPerJet > 1:
                manyParticlesPerJet = True
        oneJetPerParticle = len(dQMJs) == 1 and len(SMMJs) == 1

        if (not oneJetPerParticle) or manyParticlesPerJet:
            dQMJs = []
            SMMJs = []
    storTLV(AK8_MTdQM,dQMJs)
    storTLV(AK8_MTSMM,SMMJs)
    storTLV(AK8_MTPair,dQMJs+SMMJs)

    # group the first quarks and dark quarks from mediator according to the mediator mothers.
    # Each mediator gives one quark and one dark quark.
    # This will be useful for calculating MT2.
    dQM1Js = []
    dQM2Js = []
    SMM1Js = []
    SMM2Js = []
    if len(fdMPart) == 2:
        dQM1 = fdQPartFM[fdQPartFM[:,1] == fdMPart[0]].flatten()[0] # dark quark from mediator 1
        dQM2 = fdQPartFM[fdQPartFM[:,1] == fdMPart[1]].flatten()[0] # dark quark from mediator 2
        SMM1 = fSMqPart[fSMqPart[:,1] == fdMPart[0]].flatten()[0] # SM quark from mediator 1
        SMM2 = fSMqPart[fSMqPart[:,1] == fdMPart[1]].flatten()[0] # SM quark from mediator 2
        manyParticlesPerJet = False
        for djc in AK8Jets:
            nPartPerJet = 0
            if deltaR(djc, dQM1)<0.8:
                dQM1Js.append(djc)
                nPartPerJet += 1
            if deltaR(djc, SMM1)<0.8:
                SMM1Js.append(djc)
                nPartPerJet += 1
            if deltaR(djc, dQM2)<0.8:
                dQM2Js.append(djc)
                nPartPerJet += 1
            if deltaR(djc, SMM2)<0.8:
                SMM2Js.append(djc)
                nPartPerJet += 1
            if nPartPerJet > 1:
                manyParticlesPerJet = True
                break
        oneJetPerParticle = len(dQM1Js) == 1 and len(dQM2Js) == 1 and len(SMM1Js) == 1 and len(SMM2Js) == 1
        oneParticlePerJet = not manyParticlesPerJet
        if (not oneJetPerParticle) or manyParticlesPerJet:
            dQM1Js = []
            dQM2Js = []
            SMM1Js = []
            SMM2Js = []
    storTLV(AK8_MT2dQM1,dQM1Js)
    storTLV(AK8_MT2SMM1,SMM1Js)
    storTLV(AK8_MT2dQM2,dQM2Js)
    storTLV(AK8_MT2SMM2,SMM2Js)
    storTLV(AK8_MedPair1,dQM1Js+SMM1Js)
    storTLV(AK8_MedPair2,dQM2Js+SMM2Js)

    tr.Fill()

print("totJets: {}".format(totJets))
print("totDJets: {}".format(totDJets))
print("totSMJets: {}".format(totSMJets))
print("nMT2Event: {}".format(nMT2Event))
tr.AutoSave()
outf.Close()
# saving the file to eos
# inputFold = "root://cmseos.fnal.gov//store/user/keanet/CondorOutput/tchannel/SVJP_Master_08272020/Skims/"
# os.system("xrdcp -f Skims/" + name + " " + inputFold + name)
# os.system("rm Skims/" + name)
