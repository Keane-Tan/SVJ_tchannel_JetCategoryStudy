import ROOT as rt
import inspect
import sys
import numpy as np
import pandas as pd
from array import array
from DataFormats.FWLite import Events, Handle
import tdrstyle

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()
rt.gROOT.ProcessLine(".L loader.C+")

DarkMediatorID_ = [4900001,4900002,4900003,4900004,4900005,4900006] # 4900023 is schannel Zprime
DarkQuarkID_ = [4900101,4900102]
DarkHadronID_ = [4900111,4900113,4900211,4900213]
DarkGluonID_ = [4900021]
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
	dt = 0
	for dark in darklist:
		if abs(pId) == dark:
			dt = 1
	return dt

def medDecay(fdau,fmom):
	for mddi in range(fdau.numberOfDaughters()):
		mdd = fdau.daughter(mddi)
		mddId = mdd.pdgId()
		if darkTruth(DarkMediatorID_,mddId):
			medDecay(mdd,fmom)
		elif darkTruth(DarkQuarkID_,mddId):
			fdQPartFM.append([mdd,fmom])
		elif darkTruth(DarkGluonID_,mddId):
			fdGPartFM.append([mdd,fmom])
		else:
			fSMqPart.append([mdd,fmom])

def firstDark(partid,part_ind):
	DarkIDs = DarkMediatorID_ + DarkQuarkID_ + DarkGluonID_
	if part_ind < 20:
		if darkTruth(DarkIDs,partid) and part_i.numberOfMothers()==2:
			for dau_i in range(part_i.numberOfDaughters()):
				fdau = part_i.daughter(dau_i)
				fdauId = fdau.pdgId()
				if darkTruth(DarkMediatorID_,fdauId):
					fdMPart.append(fdau)
					medDecay(fdau,fdau)
				if darkTruth(DarkQuarkID_,fdauId):
					fdQPart.append(fdau)
				if darkTruth(DarkGluonID_,fdauId):
					fdGPart.append(fdau)

def lastDarkGQ(DarkIDList,partid,dPartList):
	if darkTruth(DarkIDList,partid) and part_i.isLastCopy():
		dPartList.append(part_i)

def checkDark(dPart,ijet,dJets):
	for dqCan in dPart:
		dR = deltaR(ijet, dqCan)
		if dR<0.8:
			global isDark
			isDark = 1
			dPartSet.append(dqCan) # collecting all the dark particles in the dark jet
	if isDark == 1:       # if more than one quarks/gluons fall into the cone of one jet we still have only one dark jet
		dPartTrack.append(np.array(dPartSet))
		dJets.append(ijet)

def darkPt(darkMList): # calculate the total dark pT
	darkMomTLSum = rt.TLorentzVector(0,0,0,0) # sum of the TLorentzVectors of the dark mothers
	for imom in darkMList:
		darkMomTLSum += rt.TLorentzVector(imom.px(),imom.py(),imom.pz(),imom.energy())
	return darkMomTLSum.Pt()

def storTLV(H,jets):
	H.clear()
	for jet in jets:
		H.push_back(rt.TLorentzVector(jet.px(),jet.py(),jet.pz(),jet.energy()))

def isfj(fPart,uj,label):
	isJ = 0
	for fP in fPart:
		if deltaR(fP,uj) < 0.8:
			isJ = 1
	if isJ == 1:
		jlabels.append(label)

def lastDau(lastD,dHad): # find the last daughters
	for idH in range(dHad.numberOfDaughters()):
		dauH = dHad.daughter(idH)
		dauHID = dauH.pdgId()
		if dauH.numberOfDaughters() == 0 and (dauHID != 51 or dauHID != 52 or dauHID != 53): # last particle that is not stable dark hadron
			lastD.append(dauH)
		elif dauH.numberOfDaughters() > 0:
			lastDau(lastD,dauH)

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

def darkPtFract(JetsList,wdpTFrac):
	iptf = 0
	for adj in JetsList:
		tot4vec = rt.TLorentzVector(0.,0.,0.,0.)
		dark4vec = rt.TLorentzVector(0.,0.,0.,0.)
		anum = adj.numberOfDaughters()
		adauL = [] # list of all the daughters for each jet
		for ani in range(anum):
			adauL.append(adj.daughter(ani))

		adarkp = [prt for prt in adauL if prt in lastD] # list of all the dark particles

		for mpa in adarkp:
			dark4vec += rt.TLorentzVector(mpa.px(),mpa.py(),mpa.pz(),mpa.energy())

		for apa in adauL:
			tot4vec += rt.TLorentzVector(apa.px(),apa.py(),apa.pz(),apa.energy())

# accounting for stable dark hadrons
		totpT = tot4vec.Pt()
		darpT = dark4vec.Pt()

		jetDarkF = darpT/totpT

		wdpTFrac[iptf] = jetDarkF
		iptf += 1

		if JetsList == dJets:
			dJPts.append(jetDarkF)

def darkPtFract_Indi(jet):
	tot4vec = rt.TLorentzVector(0.,0.,0.,0.)
	dark4vec = rt.TLorentzVector(0.,0.,0.,0.)
	anum = jet.numberOfDaughters()
	adauL = [] # list of all the daughters for each jet
	for ani in range(anum):
		adauL.append(jet.daughter(ani))

	adarkp = [prt for prt in adauL if prt in lastD] # list of all the dark particles

	for mpa in adarkp:
		dark4vec += rt.TLorentzVector(mpa.px(),mpa.py(),mpa.pz(),mpa.energy())

	for apa in adauL:
		tot4vec += rt.TLorentzVector(apa.px(),apa.py(),apa.pz(),apa.energy())

# accounting for stable dark hadrons
	totpT = tot4vec.Pt()
	darpT = dark4vec.Pt()

	jetDarkF = darpT/totpT
	return jetDarkF

def nPartPerJet(JetsList,nPBranch):
	for bj in range(len(JetsList)):
		bjet = JetsList[bj]
		npjet = bjet.numberOfDaughters()

		nPBranch[bj] = npjet


# args, in order mMed, mDark, Rinv, Alpha
mMed, mDark, rInv, Alpha = sys.argv[1:]
# print("mMed mDark rInv Alpha")
# print("{} {} {} {}".format(mMed, mDark, rInv, Alpha))

pbegin = 1			# run all: 1
pend = 101	# run all: 101

outf = rt.TFile( "Skims/s1_mMed-{}_mDark-{}_rinv-{}_alpha-{}.root".format(mMed,mDark,rInv,Alpha), 'RECREATE' )
cdGMA = outf.mkdir("GenMassAnalyzer")
cdGMA.cd()

tr = rt.TTree( 'tree', 'GenMassAnalyzer Tree' )

nfdMPart = array( 'i', [ 0 ] )
nAK8Jets = array( 'i', [ 0 ] )
nDJets = array( 'i', [ 0 ] )
nSM_Jets = array( 'i', [ 0 ] )
nSMM_Jets = array( 'i', [ 0 ] )
nG_Jets = array( 'i', [ 0 ] )
nQM_Jets = array( 'i', [ 0 ] )
nQ_Jets = array( 'i', [ 0 ] )
nQM_QJets = array( 'i', [ 0 ] )
nQM_GJets = array( 'i', [ 0 ] )
nQ_GJets = array( 'i', [ 0 ] )
nG_SMJets = array( 'i', [ 0 ] )
nQM_SMJets = array( 'i', [ 0 ] )
nQ_SMJets = array( 'i', [ 0 ] )
nLD_lowDFJets = array( 'i', [ 0 ] )
nLD_highDFJets = array( 'i', [ 0 ] )
nLD_SMJets = array( 'i', [ 0 ] )

AK8Jets = rt.std.vector(rt.TLorentzVector)()
AK8_DarkJets = rt.std.vector(rt.TLorentzVector)()
AK8_SM_Jets = rt.std.vector(rt.TLorentzVector)()
AK8_SMM_Jets = rt.std.vector(rt.TLorentzVector)()
AK8_G_Jets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only first dark gluons
AK8_QM_Jets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only first dark quarks from mediator
AK8_Q_Jets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only first dark quarks not from mediator
AK8_QM_QJets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only first dark quarks from and not from mediator
AK8_QM_GJets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only first dark quarks from mediator and first dark gluons
AK8_Q_GJets = rt.std.vector(rt.TLorentzVector)()		# jets that have only first dark quarks not from mediator and first dark gluons
AK8_G_SMJets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only first dark gluons and SM quark from mediator
AK8_QM_SMJets = rt.std.vector(rt.TLorentzVector)()	# jets that have only first dark quarks from mediator and SM quark from mediator
AK8_Q_SMJets = rt.std.vector(rt.TLorentzVector)()	# jets that have only first dark quarks not from mediator and SM quark from mediator
AK8_LD_lowDFJets = rt.std.vector(rt.TLorentzVector)()	# jets that have only last dark particles
AK8_LD_highDFJets = rt.std.vector(rt.TLorentzVector)()	# jets that have only last dark particles
AK8_LD_SMJets = rt.std.vector(rt.TLorentzVector)() 	# jets that have only last dark particles and SM quarks from mediator
AK8_MedPair1 = rt.std.vector(rt.TLorentzVector)()
AK8_MedPair2 = rt.std.vector(rt.TLorentzVector)()
AK8_MTPair = rt.std.vector(rt.TLorentzVector)()
MET = rt.std.vector(rt.TLorentzVector)()
# jet substructure
Nmax = 50 # maximum number of jet in an event
dgirth = array( 'd', Nmax*[ 0. ] )
daxisMj = array( 'd', Nmax*[ 0. ] )
daxisMn = array( 'd', Nmax*[ 0. ] )
dptD = array( 'd', Nmax*[ 0. ] )
AK8_girth = array( 'd', Nmax*[ 0. ] )
AK8_axisMj = array( 'd', Nmax*[ 0. ] )
AK8_axisMn = array( 'd', Nmax*[ 0. ] )
AK8_ptD = array( 'd', Nmax*[ 0. ] )
SM_girth = array( 'd', Nmax*[ 0. ] )
SM_axisMj = array( 'd', Nmax*[ 0. ] )
SM_axisMn = array( 'd', Nmax*[ 0. ] )
SM_ptD = array( 'd', Nmax*[ 0. ] )
SMM_girth = array( 'd', Nmax*[ 0. ] )
SMM_axisMj = array( 'd', Nmax*[ 0. ] )
SMM_axisMn = array( 'd', Nmax*[ 0. ] )
SMM_ptD = array( 'd', Nmax*[ 0. ] )
G_girth = array( 'd', Nmax*[ 0. ] )
G_axisMj = array( 'd', Nmax*[ 0. ] )
G_axisMn = array( 'd', Nmax*[ 0. ] )
G_ptD = array( 'd', Nmax*[ 0. ] )
QM_girth = array( 'd', Nmax*[ 0. ] )
QM_axisMj = array( 'd', Nmax*[ 0. ] )
QM_axisMn = array( 'd', Nmax*[ 0. ] )
QM_ptD = array( 'd', Nmax*[ 0. ] )
Q_girth = array( 'd', Nmax*[ 0. ] )
Q_axisMj = array( 'd', Nmax*[ 0. ] )
Q_axisMn = array( 'd', Nmax*[ 0. ] )
Q_ptD = array( 'd', Nmax*[ 0. ] )
QM_Qgirth = array( 'd', Nmax*[ 0. ] )
QM_QaxisMj = array( 'd', Nmax*[ 0. ] )
QM_QaxisMn = array( 'd', Nmax*[ 0. ] )
QM_QptD = array( 'd', Nmax*[ 0. ] )
QM_Ggirth = array( 'd', Nmax*[ 0. ] )
QM_GaxisMj = array( 'd', Nmax*[ 0. ] )
QM_GaxisMn = array( 'd', Nmax*[ 0. ] )
QM_GptD = array( 'd', Nmax*[ 0. ] )
Q_Ggirth = array( 'd', Nmax*[ 0. ] )
Q_GaxisMj = array( 'd', Nmax*[ 0. ] )
Q_GaxisMn = array( 'd', Nmax*[ 0. ] )
Q_GptD = array( 'd', Nmax*[ 0. ] )
G_SMgirth = array( 'd', Nmax*[ 0. ] )
G_SMaxisMj = array( 'd', Nmax*[ 0. ] )
G_SMaxisMn = array( 'd', Nmax*[ 0. ] )
G_SMptD = array( 'd', Nmax*[ 0. ] )
QM_SMgirth = array( 'd', Nmax*[ 0. ] )
QM_SMaxisMj = array( 'd', Nmax*[ 0. ] )
QM_SMaxisMn = array( 'd', Nmax*[ 0. ] )
QM_SMptD = array( 'd', Nmax*[ 0. ] )
Q_SMgirth = array( 'd', Nmax*[ 0. ] )
Q_SMaxisMj = array( 'd', Nmax*[ 0. ] )
Q_SMaxisMn = array( 'd', Nmax*[ 0. ] )
Q_SMptD = array( 'd', Nmax*[ 0. ] )
LD_lowDFgirth = array( 'd', Nmax*[ 0. ] )
LD_lowDFaxisMj = array( 'd', Nmax*[ 0. ] )
LD_lowDFaxisMn = array( 'd', Nmax*[ 0. ] )
LD_lowDFptD = array( 'd', Nmax*[ 0. ] )
LD_highDFgirth = array( 'd', Nmax*[ 0. ] )
LD_highDFaxisMj = array( 'd', Nmax*[ 0. ] )
LD_highDFaxisMn = array( 'd', Nmax*[ 0. ] )
LD_highDFptD = array( 'd', Nmax*[ 0. ] )
LD_SMgirth = array( 'd', Nmax*[ 0. ] )
LD_SMaxisMj = array( 'd', Nmax*[ 0. ] )
LD_SMaxisMn = array( 'd', Nmax*[ 0. ] )
LD_SMptD = array( 'd', Nmax*[ 0. ] )
# dark pT fraction
dpTFrac = array( 'd', Nmax*[ 0. ] )
G_dpTFrac = array( 'd', Nmax*[ 0. ] )
QM_dpTFrac = array( 'd', Nmax*[ 0. ] )
Q_dpTFrac = array( 'd', Nmax*[ 0. ] )
QM_QdpTFrac = array( 'd', Nmax*[ 0. ] )
QM_GdpTFrac = array( 'd', Nmax*[ 0. ] )
Q_GdpTFrac = array( 'd', Nmax*[ 0. ] )
G_SMdpTFrac = array( 'd', Nmax*[ 0. ] )
QM_SMdpTFrac = array( 'd', Nmax*[ 0. ] )
Q_SMdpTFrac = array( 'd', Nmax*[ 0. ] )
LD_lowDFdpTFrac = array( 'd', Nmax*[ 0. ] )
LD_highDFdpTFrac = array( 'd', Nmax*[ 0. ] )
LD_SMdpTFrac = array( 'd', Nmax*[ 0. ] )
# number of particles for each jet
npAK8Jets = array( 'i', Nmax*[ 0 ] )
npDJets = array( 'i', Nmax*[ 0 ] )
npSM_Jets = array( 'i', Nmax*[ 0 ] )
npSMM_Jets = array( 'i', Nmax*[ 0 ] )
npG_Jets = array( 'i', Nmax*[ 0 ] )
npQM_Jets = array( 'i', Nmax*[ 0 ] )
npQ_Jets = array( 'i', Nmax*[ 0 ] )
npQM_QJets = array( 'i', Nmax*[ 0 ] )
npQM_GJets = array( 'i', Nmax*[ 0 ] )
npQ_GJets = array( 'i', Nmax*[ 0 ] )
npG_SMJets = array( 'i', Nmax*[ 0 ] )
npQM_SMJets = array( 'i', Nmax*[ 0 ] )
npQ_SMJets = array( 'i', Nmax*[ 0 ] )
npLD_lowDFJets = array( 'i', Nmax*[ 0 ] )
npLD_highDFJets = array( 'i', Nmax*[ 0 ] )
npLD_SMJets = array( 'i', Nmax*[ 0 ] )

tr.Branch( 'nDJets', nDJets, 'nDJets/I' )
tr.Branch( 'nfdMPart', nfdMPart, 'nfdMPart/I' )
tr.Branch( 'nAK8Jets', nAK8Jets, 'nAK8Jets/I' )
tr.Branch( 'nSM_Jets', nSM_Jets, 'nSM_Jets/I' )
tr.Branch( 'nSMM_Jets', nSMM_Jets, 'nSMM_Jets/I' )
tr.Branch( 'nG_Jets', nG_Jets, 'nG_Jets/I' )
tr.Branch( 'nQM_Jets', nQM_Jets, 'nQM_Jets/I' )
tr.Branch( 'nQ_Jets', nQ_Jets, 'nQ_Jets/I' )
tr.Branch( 'nQM_QJets', nQM_QJets, 'nQM_QJets/I' )
tr.Branch( 'nQM_GJets', nQM_GJets, 'nQM_GJets/I' )
tr.Branch( 'nQ_GJets', nQ_GJets, 'nQ_GJets/I' )
tr.Branch( 'nG_SMJets', nG_SMJets, 'nG_SMJets/I' )
tr.Branch( 'nQM_SMJets', nQM_SMJets, 'nQM_SMJets/I' )
tr.Branch( 'nQ_SMJets', nQ_SMJets, 'nQ_SMJets/I' )
tr.Branch( 'nLD_lowDFJets', nLD_lowDFJets, 'nLD_lowDFJets/I' )
tr.Branch( 'nLD_highDFJets', nLD_highDFJets, 'nLD_highDFJets/I' )
tr.Branch( 'nLD_SMJets', nLD_SMJets, 'nLD_SMJets/I' )
tr.Branch("AK8Jets","vector<TLorentzVector>",AK8Jets)
tr.Branch("AK8_DarkJets","vector<TLorentzVector>",AK8_DarkJets)
tr.Branch("AK8_SM_Jets","vector<TLorentzVector>",AK8_SM_Jets)
tr.Branch("AK8_SMM_Jets","vector<TLorentzVector>",AK8_SMM_Jets)
tr.Branch("AK8_G_Jets","vector<TLorentzVector>",AK8_G_Jets)
tr.Branch("AK8_QM_Jets","vector<TLorentzVector>",AK8_QM_Jets)
tr.Branch("AK8_Q_Jets","vector<TLorentzVector>",AK8_Q_Jets)
tr.Branch("AK8_QM_QJets","vector<TLorentzVector>",AK8_QM_QJets)
tr.Branch("AK8_QM_GJets","vector<TLorentzVector>",AK8_QM_GJets)
tr.Branch("AK8_Q_GJets","vector<TLorentzVector>",AK8_Q_GJets)
tr.Branch("AK8_G_SMJets","vector<TLorentzVector>",AK8_G_SMJets)
tr.Branch("AK8_QM_SMJets","vector<TLorentzVector>",AK8_QM_SMJets)
tr.Branch("AK8_Q_SMJets","vector<TLorentzVector>",AK8_Q_SMJets)
tr.Branch("AK8_LD_lowDFJets","vector<TLorentzVector>",AK8_LD_lowDFJets)
tr.Branch("AK8_LD_highDFJets","vector<TLorentzVector>",AK8_LD_highDFJets)
tr.Branch("AK8_LD_SMJets","vector<TLorentzVector>",AK8_LD_SMJets)
tr.Branch("AK8_MedPair1","vector<TLorentzVector>",AK8_MedPair1)
tr.Branch("AK8_MedPair2","vector<TLorentzVector>",AK8_MedPair2)
tr.Branch("AK8_MTPair","vector<TLorentzVector>",AK8_MTPair)
tr.Branch("MET","vector<TLorentzVector>",MET)
# jet substructure
tr.Branch( 'AK8_girth', AK8_girth, 'girth[nAK8Jets]/D' )
tr.Branch( 'AK8_axisMj', AK8_axisMj, 'axisMj[nAK8Jets]/D' )
tr.Branch( 'AK8_axisMn', AK8_axisMn, 'axisMn[nAK8Jets]/D' )
tr.Branch( 'AK8_ptD', AK8_ptD, 'ptD[nAK8Jets]/D' )
tr.Branch( 'Dark_girth', dgirth, 'girth[nDJets]/D' )
tr.Branch( 'Dark_axisMj', daxisMj, 'axisMj[nDJets]/D' )
tr.Branch( 'Dark_axisMn', daxisMn, 'axisMn[nDJets]/D' )
tr.Branch( 'Dark_ptD', dptD, 'ptD[nDJets]/D' )
tr.Branch( 'SM_girth', SM_girth, 'girth[nSM_Jets]/D' )
tr.Branch( 'SM_axisMj', SM_axisMj, 'axisMj[nSM_Jets]/D' )
tr.Branch( 'SM_axisMn', SM_axisMn, 'axisMn[nSM_Jets]/D' )
tr.Branch( 'SM_ptD', SM_ptD, 'ptD[nSM_Jets]/D' )
tr.Branch( 'SMM_girth', SMM_girth, 'girth[nSMM_Jets]/D' )
tr.Branch( 'SMM_axisMj', SMM_axisMj, 'axisMj[nSMM_Jets]/D' )
tr.Branch( 'SMM_axisMn', SMM_axisMn, 'axisMn[nSMM_Jets]/D' )
tr.Branch( 'SMM_ptD', SMM_ptD, 'ptD[nSMM_Jets]/D' )
tr.Branch( 'G_girth', G_girth, 'girth[nG_Jets]/D' )
tr.Branch( 'G_axisMj', G_axisMj, 'axisMj[nG_Jets]/D' )
tr.Branch( 'G_axisMn', G_axisMn, 'axisMn[nG_Jets]/D' )
tr.Branch( 'G_ptD', G_ptD, 'ptD[nG_Jets]/D' )
tr.Branch( 'QM_girth', QM_girth, 'girth[nQM_Jets]/D' )
tr.Branch( 'QM_axisMj', QM_axisMj, 'axisMj[nQM_Jets]/D' )
tr.Branch( 'QM_axisMn', QM_axisMn, 'axisMn[nQM_Jets]/D' )
tr.Branch( 'QM_ptD', QM_ptD, 'ptD[nQM_Jets]/D' )
tr.Branch( 'Q_girth', Q_girth, 'girth[nQ_Jets]/D' )
tr.Branch( 'Q_axisMj', Q_axisMj, 'axisMj[nQ_Jets]/D' )
tr.Branch( 'Q_axisMn', Q_axisMn, 'axisMn[nQ_Jets]/D' )
tr.Branch( 'Q_ptD', Q_ptD, 'ptD[nQ_Jets]/D' )
tr.Branch( 'QM_Qgirth', QM_Qgirth, 'girth[nQM_QJets]/D' )
tr.Branch( 'QM_QaxisMj', QM_QaxisMj, 'axisMj[nQM_QJets]/D' )
tr.Branch( 'QM_QaxisMn', QM_QaxisMn, 'axisMn[nQM_QJets]/D' )
tr.Branch( 'QM_QptD', QM_QptD, 'ptD[nQM_QJets]/D' )
tr.Branch( 'QM_Ggirth', QM_Ggirth, 'girth[nQM_GJets]/D' )
tr.Branch( 'QM_GaxisMj', QM_GaxisMj, 'axisMj[nQM_GJets]/D' )
tr.Branch( 'QM_GaxisMn', QM_GaxisMn, 'axisMn[nQM_GJets]/D' )
tr.Branch( 'QM_GptD', QM_GptD, 'ptD[nQM_GJets]/D' )
tr.Branch( 'Q_Ggirth', Q_Ggirth, 'girth[nQ_GJets]/D' )
tr.Branch( 'Q_GaxisMj', Q_GaxisMj, 'axisMj[nQ_GJets]/D' )
tr.Branch( 'Q_GaxisMn', Q_GaxisMn, 'axisMn[nQ_GJets]/D' )
tr.Branch( 'Q_GptD', Q_GptD, 'ptD[nQ_GJets]/D' )
tr.Branch( 'G_SMgirth', G_SMgirth, 'girth[nG_SMJets]/D' )
tr.Branch( 'G_SMaxisMj', G_SMaxisMj, 'axisMj[nG_SMJets]/D' )
tr.Branch( 'G_SMaxisMn', G_SMaxisMn, 'axisMn[nG_SMJets]/D' )
tr.Branch( 'G_SMptD', G_SMptD, 'ptD[nG_SMJets]/D' )
tr.Branch( 'QM_SMgirth', QM_SMgirth, 'girth[nQM_SMJets]/D' )
tr.Branch( 'QM_SMaxisMj', QM_SMaxisMj, 'axisMj[nQM_SMJets]/D' )
tr.Branch( 'QM_SMaxisMn', QM_SMaxisMn, 'axisMn[nQM_SMJets]/D' )
tr.Branch( 'QM_SMptD', QM_SMptD, 'ptD[nQM_SMJets]/D' )
tr.Branch( 'Q_SMgirth', Q_SMgirth, 'girth[nQ_SMJets]/D' )
tr.Branch( 'Q_SMaxisMj', Q_SMaxisMj, 'axisMj[nQ_SMJets]/D' )
tr.Branch( 'Q_SMaxisMn', Q_SMaxisMn, 'axisMn[nQ_SMJets]/D' )
tr.Branch( 'Q_SMptD', Q_SMptD, 'ptD[nQ_SMJets]/D' )
tr.Branch( 'LD_lowDFgirth', LD_lowDFgirth, 'girth[nLD_lowDFJets]/D' )
tr.Branch( 'LD_lowDFaxisMj', LD_lowDFaxisMj, 'axisMj[nLD_lowDFJets]/D' )
tr.Branch( 'LD_lowDFaxisMn', LD_lowDFaxisMn, 'axisMn[nLD_lowDFJets]/D' )
tr.Branch( 'LD_lowDFptD', LD_lowDFptD, 'ptD[nLD_lowDFJets]/D' )
tr.Branch( 'LD_highDFgirth', LD_highDFgirth, 'girth[nLD_highDFJets]/D' )
tr.Branch( 'LD_highDFaxisMj', LD_highDFaxisMj, 'axisMj[nLD_highDFJets]/D' )
tr.Branch( 'LD_highDFaxisMn', LD_highDFaxisMn, 'axisMn[nLD_highDFJets]/D' )
tr.Branch( 'LD_highDFptD', LD_highDFptD, 'ptD[nLD_highDFJets]/D' )
tr.Branch( 'LD_SMgirth', LD_SMgirth, 'girth[nLD_SMJets]/D' )
tr.Branch( 'LD_SMaxisMj', LD_SMaxisMj, 'axisMj[nLD_SMJets]/D' )
tr.Branch( 'LD_SMaxisMn', LD_SMaxisMn, 'axisMn[nLD_SMJets]/D' )
tr.Branch( 'LD_SMptD', LD_SMptD, 'ptD[nLD_SMJets]/D' )
# dark pT fraction
tr.Branch( 'dpTFrac', dpTFrac, 'dpTFrac[nDJets]/D' )
tr.Branch( 'G_dpTFrac', G_dpTFrac, 'G_dpTFrac[nG_Jets]/D' )
tr.Branch( 'QM_dpTFrac', QM_dpTFrac, 'QM_dpTFrac[nQM_Jets]/D' )
tr.Branch( 'Q_dpTFrac', Q_dpTFrac, 'Q_dpTFrac[nQ_Jets]/D' )
tr.Branch( 'QM_QdpTFrac', QM_QdpTFrac, 'QM_QdpTFrac[nQM_QJets]/D' )
tr.Branch( 'QM_GdpTFrac', QM_GdpTFrac, 'QM_GdpTFrac[nQM_GJets]/D' )
tr.Branch( 'Q_GdpTFrac', Q_GdpTFrac, 'Q_GdpTFrac[nQ_GJets]/D' )
tr.Branch( 'G_SMdpTFrac', G_SMdpTFrac, 'G_SMdpTFrac[nG_SMJets]/D' )
tr.Branch( 'QM_SMdpTFrac', QM_SMdpTFrac, 'QM_SMdpTFrac[nQM_SMJets]/D' )
tr.Branch( 'Q_SMdpTFrac', Q_SMdpTFrac, 'Q_SMdpTFrac[nQ_SMJets]/D' )
tr.Branch( 'LD_lowDFdpTFrac', LD_lowDFdpTFrac, 'LD_lowDFdpTFrac[nLD_lowDFJets]/D' )
tr.Branch( 'LD_highDFdpTFrac', LD_highDFdpTFrac, 'LD_highDFdpTFrac[nLD_highDFJets]/D' )
tr.Branch( 'LD_SMdpTFrac', LD_SMdpTFrac, 'LD_SMdpTFrac[nLD_SMJets]/D' )
# number of particles in each jet
tr.Branch( 'npAK8Jets', npAK8Jets, 'npAK8Jets[nAK8Jets]/I' )
tr.Branch( 'npDJets', npDJets, 'npDJets[nDJets]/I' )
tr.Branch( 'npSM_Jets', npSM_Jets, 'npSM_Jets[nSM_Jets]/I' )
tr.Branch( 'npSMM_Jets', npSMM_Jets, 'npSMM_Jets[nSMM_Jets]/I' )
tr.Branch( 'npG_Jets', npG_Jets, 'npG_Jets[nG_Jets]/I' )
tr.Branch( 'npQM_Jets', npQM_Jets, 'npQM_Jets[nQM_Jets]/I' )
tr.Branch( 'npQ_Jets', npQ_Jets, 'npQ_Jets[nQ_Jets]/I' )
tr.Branch( 'npQM_QJets', npQM_QJets, 'npQM_QJets[nQM_QJets]/I' )
tr.Branch( 'npQM_GJets', npQM_GJets, 'npQM_GJets[nQM_GJets]/I' )
tr.Branch( 'npQ_GJets', npQ_GJets, 'npQ_GJets[nQ_GJets]/I' )
tr.Branch( 'npG_SMJets', npG_SMJets, 'npG_SMJets[nG_SMJets]/I' )
tr.Branch( 'npQM_SMJets', npQM_SMJets, 'npQM_SMJets[nQM_SMJets]/I' )
tr.Branch( 'npQ_SMJets', npQ_SMJets, 'npQ_SMJets[nQ_SMJets]/I' )
tr.Branch( 'npLD_lowDFJets', npLD_lowDFJets, 'npLD_lowDFJets[nLD_lowDFJets]/I' )
tr.Branch( 'npLD_highDFJets', npLD_highDFJets, 'npLD_highDFJets[nLD_highDFJets]/I' )
tr.Branch( 'npLD_SMJets', npLD_SMJets, 'npLD_SMJets[nLD_SMJets]/I' )

fileStart = "root://cmseos.fnal.gov//store/user/keanet/tchannel/GEN/" # output directory of step1_GEN
listOfFiles = [fileStart+"step1_LHE-GEN_t-channel_mMed-{}_mDark-{}_rinv-{}_alpha-{}_13TeV-madgraphMLM-pythia8_n-1000_part-{}.root".format(mMed,mDark,rInv,Alpha,iPart) for iPart in range(pbegin,pend)]#101

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

for event in events:
	count += 1
	print "Event number " + str(count)
	event.getByLabel(labelJets,handleJets)
	event.getByLabel(labelMET,handleMET)
	event.getByLabel(labelPart,handleParts)
	jets = handleJets.product()
	hMET = handleMET.product()
	hPart = handleParts.product()

	dJets = [] # dark quark jets
	SMJets = [] # Standard Model jets

	dQPart = [] # last dark quarks
	dGPart = [] # last dark gluons
	dHPart = [] # last dark hadrons

	fdQPart = [] # first dark quarks
	fdGPart = [] # first dark gluons
	fdMPart = [] # first dark mediators

	fdQPartFM = [] # first dark quarks from the dark mediator
	fdGPartFM = [] # first dark gluons from the dark mediator
	fSMqPart = []  # first quarks from the dark mediator decay

	stableD = [] # stable dark hadrons

	G_ = []
	QM_ = []
	Q_ = []
	QM_Q = []
	QM_G = []
	Q_G = []
	G_SM = []
	QM_SM = []
	Q_SM = []
	LD_lowDF = []
	LD_highDF = []
	LD_SM = []
	SMM_ = [] # SM jets that contain SM quarks from mediator
	SM_ = []

	if len(jets) < 1:
		continue

# collecting all the last dark quarks and dark gluons
	part_ind = 0
	for part_i in handleParts.product():
		part_ind += 1
		partid = part_i.pdgId()
		lastDarkGQ(DarkQuarkID_,partid,dQPart)
		lastDarkGQ(DarkGluonID_,partid,dGPart)
		lastDarkGQ(DarkHadronID_,partid,dHPart)

		firstDark(partid,part_ind)

		if partid == 51 or partid == 52 or partid == 53:
			stableD.append(part_i)

# finding jets that contain dark gluons/quarks/both
	dPart = dQPart + dGPart
	dPartTrack = [] # keeping track of which dark jet contains which dark particles

	fjets = []
	for ijet in jets:
		isDark = 0
		if ijet.pt() < 200:
			continue
		fjets.append(ijet)
		dPartSet = []
		checkDark(dPart,ijet,dJets)
		if isDark == 0:
			SMJets.append(ijet)

	dJPts = [] # dark jets' pT fraction
	djnum = 0
	# print "Number of last dark daughters = " + str(len(lastD))

# every dark quark and gluon eventually turns into dark hadrons that may or may not decay.
# Here, we want to create a list of last particles (numberOfDaughters=0) that came from these dark hadrons.
	lastD = [] # last daughters of all the last dark hadrons
	for dmom in dHPart:
		lastDau(lastD,dmom)

	lastD += stableD
	darkPtFract(dJets,dpTFrac)

# pick the dark jet with the largest fraction of dark pT to associate with the dark quarks/gluons within the jet
	col = ["Dark pT Frac","Dark Jets","Particles in Jet"]

	dTab = pd.DataFrame(np.transpose(np.array([dJPts, dJets, dPartTrack])),columns=col)
	dTab = dTab.sort_values("Dark pT Frac",ascending=False)

	dJlist = np.array(dTab["Dark Jets"])
	dPlist = np.array(dTab["Particles in Jet"])

	for i in range(len(dPlist)):
		for j in range(len(dPlist)):
			if i < j:
				comval = np.intersect1d(dPlist[i], dPlist[j])

				for val in comval:
					dPlist[j] = np.delete(dPlist[j], np.argwhere(dPlist[j] == val))


	indDel = []
	for ai in range(len(dPlist)):
		arri = dPlist[ai]
		if len(arri) == 0:
			indDel.append(ai)

	for idD in indDel:
		SMJets.append(dJlist[idD])
	dJlist = np.delete(dJlist, indDel)
	dJets = list(dJlist)

# The dark jets identified so far are still not "pure" enough.

# Here we are attempting to further classify and reclassify these dark jets.
	relab = ["dq","dg","dqm","SMqm","lD"]
	dJetswL = [] # dark jets with corresponding labels

	fdQPFM_flat = []
	fSMqPFM_flat = []

	for ifd in fdQPartFM:
		fdQPFM_flat.append(ifd[0])
	for ifsm in fSMqPart:
		fSMqPFM_flat.append(ifsm[0])

	# reclassifying the SM jets into mediator SM and non-mediator SM
	for smj in SMJets:

		SMM_tr = 0

		for smmq in fSMqPFM_flat:
			if deltaR(smmq,smj) < 0.8:
				SMM_tr = 1
				break

		if SMM_tr == 1:
			SMM_.append(smj)
		else:
			SM_.append(smj)

	for udj in dJets:
		jlabels =[]
		isfj(fdQPart,udj,"dq")
		isfj(fdGPart,udj,"dg")
		isfj(fdQPFM_flat,udj,"dqm")
		isfj(fSMqPFM_flat,udj,"SMqm")
		isfj(dQPart+dGPart,udj,"lD")
		dJetswL.append([udj,jlabels])

	for djl in dJetswL:
		djjet = djl[0]
		djlab = djl[1]

		if ("dg" in djlab) and ("dq" not in djlab and "dqm" not in djlab and "SMqm" not in djlab):
			G_.append(djjet)
		elif ("dqm" in djlab) and ("dq" not in djlab and "dg" not in djlab and "SMqm" not in djlab):
			QM_.append(djjet)
		elif ("dq" in djlab) and ("dqm" not in djlab and "dg" not in djlab and "SMqm" not in djlab):
			Q_.append(djjet)
		elif ("dqm" in djlab and "dq" in djlab) and ("dg" not in djlab and "SMqm" not in djlab):
			QM_Q.append(djjet)
		elif ("dqm" in djlab and "dg" in djlab) and ("dq" not in djlab and "SMqm" not in djlab):
			QM_G.append(djjet)
		elif ("dq" in djlab and "dg" in djlab) and ("dqm" not in djlab and "SMqm" not in djlab):
			Q_G.append(djjet)
		elif ("dg" in djlab and "SMqm" in djlab) and ("dqm" not in djlab and "dq" not in djlab):
			G_SM.append(djjet)
		elif ("dqm" in djlab and "SMqm" in djlab) and ("dg" not in djlab and "dq" not in djlab):
			QM_SM.append(djjet)
		elif ("dq" in djlab and "SMqm" in djlab) and ("dg" not in djlab and "dqm" not in djlab):
			Q_SM.append(djjet)
		elif "lD" in djlab and ("dg" not in djlab and "dq" not in djlab and "dqm" not in djlab and "SMqm" not in djlab):
			if darkPtFract_Indi(djjet) < 0.7:
				LD_lowDF.append(djjet)
			elif darkPtFract_Indi(djjet) >= 0.7:
				LD_highDF.append(djjet)
		elif "lD" in djlab and "SMqm" in djlab and ("dq" not in djlab and "dqm" not in djlab and "dg" not in djlab):
			LD_SM.append(djjet)

	numdJ += len(dJets)

	nAK8Jets[0] = len(fjets)
	nDJets[0] = len(dJets)
	nSM_Jets[0] = len(SM_)
	nSMM_Jets[0] = len(SMM_)
	nG_Jets[0] = len(G_)
	nQM_Jets[0] = len(QM_)
	nQ_Jets[0] = len(Q_)
	nQM_QJets[0] = len(QM_Q)
	nQM_GJets[0] = len(QM_G)
	nQ_GJets[0] = len(Q_G)
	nG_SMJets[0] = len(G_SM)
	nQM_SMJets[0] = len(QM_SM)
	nQ_SMJets[0] = len(Q_SM)
	nLD_lowDFJets[0] = len(LD_lowDF)
	nLD_highDFJets[0] = len(LD_highDF)
	nLD_SMJets[0] = len(LD_SM)

	## This is for calculating dark pT fraction with genjetnonu
	# double = 0
	# allJets = dJets + SM_ + SMM_
	# for sthad in stableD:
	# 	matched = 0
	# 	for aj in allJets:
	# 		if deltaR(aj,sthad) < 0.8:
	# 			matched += 1
	#
	# 	if matched > 1:
	# 		double = 1
	# 		break
	#
	# if double == 1:
	# 	continue # skip events that contain shared stable dark hadrons
	#
	# stabOLEvent += 1

	storTLV(AK8Jets,fjets)
	storTLV(AK8_DarkJets,dJets)
	storTLV(AK8_SM_Jets,SM_)
	storTLV(AK8_SMM_Jets,SMM_)
	storTLV(AK8_G_Jets,G_)
	storTLV(AK8_QM_Jets,QM_)
	storTLV(AK8_Q_Jets,Q_)
	storTLV(AK8_QM_QJets,QM_Q)
	storTLV(AK8_QM_GJets,QM_G)
	storTLV(AK8_Q_GJets,Q_G)
	storTLV(AK8_G_SMJets,G_SM)
	storTLV(AK8_QM_SMJets,QM_SM)
	storTLV(AK8_Q_SMJets,Q_SM)
	storTLV(AK8_LD_lowDFJets,LD_lowDF)
	storTLV(AK8_LD_highDFJets,LD_highDF)
	storTLV(AK8_LD_SMJets,LD_SM)
	storTLV(MET,[hMET[0]])

	jetSub_Calc(fjets,AK8_girth,AK8_axisMj,AK8_axisMn,AK8_ptD)
	jetSub_Calc(dJets,dgirth,daxisMj,daxisMn,dptD)
	jetSub_Calc(SM_,SM_girth,SM_axisMj,SM_axisMn,SM_ptD)
	jetSub_Calc(SMM_,SMM_girth,SMM_axisMj,SMM_axisMn,SMM_ptD)
	jetSub_Calc(G_,G_girth,G_axisMj,G_axisMn,G_ptD)
	jetSub_Calc(QM_,QM_girth,QM_axisMj,QM_axisMn,QM_ptD)
	jetSub_Calc(Q_,Q_girth,Q_axisMj,Q_axisMn,Q_ptD)
	jetSub_Calc(QM_Q,QM_Qgirth,QM_QaxisMj,QM_QaxisMn,QM_QptD)
	jetSub_Calc(QM_G,QM_Ggirth,QM_GaxisMj,QM_GaxisMn,QM_GptD)
	jetSub_Calc(Q_G,Q_Ggirth,Q_GaxisMj,Q_GaxisMn,Q_GptD)
	jetSub_Calc(G_SM,G_SMgirth,G_SMaxisMj,G_SMaxisMn,G_SMptD)
	jetSub_Calc(QM_SM,QM_SMgirth,QM_SMaxisMj,QM_SMaxisMn,QM_SMptD)
	jetSub_Calc(Q_SM,Q_SMgirth,Q_SMaxisMj,Q_SMaxisMn,Q_SMptD)
	jetSub_Calc(LD_lowDF,LD_lowDFgirth,LD_lowDFaxisMj,LD_lowDFaxisMn,LD_lowDFptD)
	jetSub_Calc(LD_highDF,LD_highDFgirth,LD_highDFaxisMj,LD_highDFaxisMn,LD_highDFptD)
	jetSub_Calc(LD_SM,LD_SMgirth,LD_SMaxisMj,LD_SMaxisMn,LD_SMptD)

	darkPtFract(dJets,dpTFrac)
	darkPtFract(G_,G_dpTFrac)
	darkPtFract(QM_,QM_dpTFrac)
	darkPtFract(Q_,Q_dpTFrac)
	darkPtFract(QM_Q,QM_QdpTFrac)
	darkPtFract(QM_G,QM_GdpTFrac)
	darkPtFract(Q_G,Q_GdpTFrac)
	darkPtFract(G_SM,G_SMdpTFrac)
	darkPtFract(QM_SM,QM_SMdpTFrac)
	darkPtFract(Q_SM,Q_SMdpTFrac)
	darkPtFract(LD_lowDF,LD_lowDFdpTFrac)
	darkPtFract(LD_highDF,LD_highDFdpTFrac)
	darkPtFract(LD_SM,LD_SMdpTFrac)

	nPartPerJet(fjets,npAK8Jets)
	nPartPerJet(dJets,npDJets)
	nPartPerJet(SM_,npSM_Jets)
	nPartPerJet(SMM_,npSMM_Jets)
	nPartPerJet(G_,npG_Jets)
	nPartPerJet(QM_,npQM_Jets)
	nPartPerJet(Q_,npQ_Jets)
	nPartPerJet(QM_Q,npQM_QJets)
	nPartPerJet(QM_G,npQM_GJets)
	nPartPerJet(Q_G,npQ_GJets)
	nPartPerJet(G_SM,npG_SMJets)
	nPartPerJet(QM_SM,npQM_SMJets)
	nPartPerJet(Q_SM,npQ_SMJets)
	nPartPerJet(LD_lowDF,npLD_lowDFJets)
	nPartPerJet(LD_highDF,npLD_highDFJets)
	nPartPerJet(LD_SM,npLD_SMJets)

	# figuring out which pair of jets came from 1 mediator; useful for MT comparison
	dSMpair0 = []

	nfdMPart[0] = len(fdMPart)

	if len(fdMPart) == 1:
		pair0 = []
		for fdq0 in fdQPartFM:
			for fsm0 in fSMqPart:
				if fdq0[1] == fsm0[1]:
					pair0.append(fdq0[0])
					pair0.append(fsm0[0])
				elif fdq0[1] != fsm0[1]:
					continue

		for dj0 in dJets:
			if deltaR(dj0, pair0[0])<0.8:
				dSMpair0.append(dj0)

		SMpair0 = []
		for sm0 in SMJets:
			if deltaR(sm0, pair0[1])<0.8:
				SMpair0.append(sm0)

		def j_pt(x):
			return x.pt()

		SMpair0.sort(key=j_pt,reverse=True) # choose the SM jet close to the first quark from mediator based on pT

		if len(SMpair0) != 0:
			dSMpair0.append(SMpair0[0])

	if len(dSMpair0) != 2:
		dSMpair0 = []
	storTLV(AK8_MTPair,dSMpair0) # remember dSMpair0[0] is the dark jet, dSMpair0[1] is the SM jet.

	# group the first quarks and dark quarks from mediator according to the mediator mothers.
	# Each mediator gives one quark and one dark quark.
	# This will be useful for calculating MT2.

	dSMpair1 = []
	dSMpair2 = []
	if len(fdMPart) == 2:
		PartFM = fdQPartFM + fSMqPart
		pair1 = []
		pair2 = []
		for fdq in fdQPartFM:
			for fsm in fSMqPart:
				if fdq[1] == fsm[1] and len(pair1) == 0:
					pair1.append(fdq[0])
					pair1.append(fsm[0])
				elif fdq[1] == fsm[1] and len(pair1) > 0:
					pair2.append(fdq[0])
					pair2.append(fsm[0])
				elif fdq[1] != fsm[1]:
					continue

		for djc in dJets:
			if deltaR(djc, pair1[0])<0.8:
				dSMpair1.append(djc)
			elif deltaR(djc, pair2[0])<0.8:
				dSMpair2.append(djc)

		SMpair1 = []
		SMpair2 = []
		for smc in SMJets:
			if deltaR(smc, pair1[1])<0.8:
				SMpair1.append(smc)
			elif deltaR(smc, pair2[1])<0.8:
				SMpair2.append(smc)

		def j_pt(x):
			return x.pt()

		SMpair1.sort(key=j_pt,reverse=True) # choose the SM jet close to the first quark from mediator based on pT
		SMpair2.sort(key=j_pt,reverse=True)

		if len(SMpair1) != 0:
			dSMpair1.append(SMpair1[0])
		if len(SMpair2) != 0:
			dSMpair2.append(SMpair2[0])

	if len(dSMpair1) != 2 or len(dSMpair2) != 2:
		dSMpair1 = []
		dSMpair2 = []
	storTLV(AK8_MedPair1,dSMpair1) # remember dSMpair1[0] is the dark jet, dSMpair1[1] is the SM jet, same for dSMpair2.
	storTLV(AK8_MedPair2,dSMpair2)

	tr.Fill()
	# if count == 1000: # which event number to stop the code
	# 	break

# some statistics related to the SusJets (suspicious jets - SM jets misidentified as dark)
# print ("Total number of events = " + str(count))
# print ("Total number of events without shared dark hadrons = " + str(stabOLEvent))
tr.AutoSave()
outf.Close()
