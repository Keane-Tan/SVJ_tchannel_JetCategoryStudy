# -*- coding: utf-8 -*-

# modified from v5
# Testing the hypothesis: Q in Dark, LD_lowDF in Mix, and a new general category
# Darker_Mix: QM_SM, Q_SM, QM_Q

import ROOT as rt
import inspect
import sys
import numpy as np
import itertools
import pandas as pd
import matplotlib.pyplot as plt

rt.gROOT.LoadMacro("Interface/lester_mt2_bisect.h")
rt.gROOT.ProcessLine(".L loader.C+")
rt.gROOT.SetBatch(True)
rt.gROOT.SetStyle("Plain")

def VecToList(jetsVec):
	jetsList = []
	for jet in jetsVec:
		jetsList.append(jet)
	return jetsList

def deltaPhi(jetphiL,metPhiL):
	dphi = jetphiL - metPhiL
	if dphi < -np.pi:
		dphi += 2*np.pi
	if dphi > np.pi:
		dphi -= 2*np.pi
	return abs(dphi)

def deltaRforTL(jet1,jet2):
	phi1 = jet1.Phi()
	phi2 = jet2.Phi()
	eta1 = jet1.Eta()
	eta2 = jet2.Eta()
	dp = abs(phi1 - phi2)
	if dp > np.pi:
		dp -= 2*np.pi
	deltaR2 = (eta1 - eta2) * (eta1 - eta2) + dp * dp
	return np.sqrt(deltaR2)

def trans_mass_Njet(jet0, jet1, met, metPhi):
	visible = jet0 + jet1
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	mT_first2 = rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))
	delR = deltaRforTL(jet0,jet1)
	delPhi = deltaPhi(jet0.Phi(),jet1.Phi())
	delEta = abs(jet0.Eta() - jet1.Eta())
	return [mT_first2,delR,delPhi,delEta]

def allJetCom(pTReq=None):
	List4jets = list(itertools.combinations(range(len(jets)),4))
	List4jets3Com = []
	for _4jets in List4jets:
		_4j1 = _4jets[0]
		_4j2 = _4jets[1]
		_4j3 = _4jets[2]
		_4j4 = _4jets[3]

		if pTReq != None:
			if jets[_4j1].Pt() < 500 or jets[_4j2].Pt() < 500 or jets[_4j3].Pt() < 500 or jets[_4j4].Pt() < 500:
				continue

		List4jets3Com.append([_4j1,_4j2,_4j3,_4j4])
		List4jets3Com.append([_4j1,_4j3,_4j2,_4j4])
		List4jets3Com.append([_4j1,_4j4,_4j3,_4j2])

	return List4jets3Com

def allWrongCom(allCom,rightCom):
	for _4j3c in allCom:
		comval = np.intersect1d(_4j3c, rightCom)
		if len(comval) == 4:
			f2 = _4j3c[:2]
			l2 = _4j3c[2:]
			tf2 = rightCom[:2]
			if len(np.intersect1d(f2, tf2)) == 2 or len(np.intersect1d(l2, tf2)) == 2:
				allCom.remove(_4j3c)
	return allCom

def M_2J(j1,j2):
	totJets = rt.TLorentzVector(j1.Px(),j1.Py(),j1.Pz(),j1.Energy())
	totJets += rt.TLorentzVector(j2.Px(),j2.Py(),j2.Pz(),j2.Energy())
	return totJets.M()

def SumJet(j1,j2):
	totJets = rt.TLorentzVector(j1.Px(),j1.Py(),j1.Pz(),j1.Energy())
	totJets += rt.TLorentzVector(j2.Px(),j2.Py(),j2.Pz(),j2.Energy())
	return totJets

def MPCom(comList):
	diffList = []
	comJetList = []
	for com in comList:
		jc1 = jets[com[0]]
		jc2 = jets[com[1]]
		jc3 = jets[com[2]]
		jc4 = jets[com[3]]
		comJetList.append([jc1,jc2,jc3,jc4])
		diffList.append(abs(M_2J(jc1,jc2) - M_2J(jc3,jc4)))
	diffList = np.array(diffList)
	mdI = np.argmin(diffList)
	return comJetList[mdI]

def MT2DRCal(FDjet0,FSMjet0,FDjet1,FSMjet1,met4p):
	Fjet0 = FDjet0 + FSMjet0
	Fjet1 = FDjet1 + FSMjet1

	METx = met4p.Px()
	METy = met4p.Py()
	MT2v = rt.asymm_mt2_lester_bisect.get_mT2(
	Fjet0.M(), Fjet0.Px(), Fjet0.Py(),
	Fjet1.M(), Fjet1.Px(), Fjet1.Py(),
	METx, METy, 0.0, 0.0, 0
	)
	delR = deltaRforTL(Fjet0,Fjet1)
	delPhi = deltaPhi(Fjet0.Phi(),Fjet1.Phi())
	delEta = abs(Fjet0.Eta() - Fjet1.Eta())

	return [MT2v,delR,delPhi,delEta]

def MT_mpCalc(jetComList):
	mpCom = MPCom(jetComList)
	MT2 = MT2DRCal(mpCom[0],mpCom[1],mpCom[2],mpCom[3],met4p)
	return [MT2,mpCom]

def getPtWeight(PtHist):
	pTwList = []
	for ibin in range(42):
		pTw = PtHist.GetBinContent(ibin)
		pTwList.append(pTw)
	print pTwList

def eW(jet,HistNorm):
	for ibin in range(42):
		iLow = HistNorm.GetXaxis().GetBinLowEdge(ibin)
		iW = HistNorm.GetXaxis().GetBinWidth(ibin)

		if jet.Pt() >= 3500:
			evtWeight = HistNorm.GetBinContent(41);
			isNormal = 1
			break
		elif iLow <= jet.Pt() < iLow + iW:
			evtWeight = HistNorm.GetBinContent(ibin)
			isNormal = 1
			break

	return evtWeight

def dSfillJetVar(dSJetsList,dSjetxMult,dSjetxPt,dSjetxEta,dSjetxPhi,dSdeltaPhix,
metPhi,dSjetxPt2D,dSjetxEta2D,dSjetxPhi2D,dSdeltaPhix2D,dSdarkFrac): # similar to above but for SM and Dark jets instead of jet[N]
	for dsi in range(len(dSJetsList)):
		daf = dSdarkFrac[dsi]
		dSJets = dSJetsList[dsi]
		dSjetxPt.Fill(dSJets.Pt())
		dSjetxEta.Fill(dSJets.Eta())
		dSjetxPhi.Fill(dSJets.Phi())
		dSdeltaPhix.Fill(deltaPhi(dSJets.Phi(),metPhi))
		dSjetxPt2D.Fill(dSJets.Pt(),daf)
		dSjetxEta2D.Fill(dSJets.Eta(),daf)
		dSjetxPhi2D.Fill(dSJets.Phi(),daf)
		dSdeltaPhix2D.Fill(deltaPhi(dSJets.Phi(),metPhi),daf)

	dSjetxMult.Fill(len(dSJetsList))

def dSfillJetVar_no2D(dSJetsList,dSjetxMult,dSjetxPt,dSjetxEta,dSjetxPhi,dSdeltaPhix,
metPhi,dSdarkFrac): # similar to above but for SM and Dark jets instead of jet[N]
	for dsi in range(len(dSJetsList)):
		daf = dSdarkFrac[dsi]
		dSJets = dSJetsList[dsi]
		dSjetxPt.Fill(dSJets.Pt())
		dSjetxEta.Fill(dSJets.Eta())
		dSjetxPhi.Fill(dSJets.Phi())
		dSdeltaPhix.Fill(deltaPhi(dSJets.Phi(),metPhi))

	dSjetxMult.Fill(len(dSJetsList))

def fillSDPt(JetsList,Hist):
	for jet in JetsList:
		Hist.Fill(jet.Pt())

def fillM(JetsList,dSMx,dSMx2D,dSdarkFrac,HistNorm=None):
	for dsi in range(len(JetsList)):
		daf = dSdarkFrac[dsi]
		jet = JetsList[dsi]

		if HistNorm == None:
			w = 1
		else:
			w = eW(jet,HistNorm)

		dSMx.Fill(jet.M(),w)
		dSMx2D.Fill(jet.M(),daf,w)

def makejsubplot(kJets,kaMj,kaMn,kgir,kpD,knp,
_axisMajor,_axisMinor,_momentGirth,_ptD,_np,
_axisMajor2D,_axisMinor2D,_momentGirth2D,_ptD2D,_np2D,dSdarkFrac,HistNorm=None):
	for ja in range(len(kJets)):

		if HistNorm == None:
			w = 1
		else:
			jet = kJets[ja]
			w = eW(jet,HistNorm)

		daF = dSdarkFrac[ja]
		am = kaMj[ja]
		an = kaMn[ja]
		gi = kgir[ja]
		ptd = kpD[ja]
		npv = knp[ja]
		_axisMajor.Fill(am,w)
		_axisMinor.Fill(an,w)
		_momentGirth.Fill(gi,w)
		_ptD.Fill(ptd,w)
		_np.Fill(npv,w)
		_axisMajor2D.Fill(am,daF,w)
		_axisMinor2D.Fill(an,daF,w)
		_momentGirth2D.Fill(gi,daF,w)
		_ptD2D.Fill(ptd,daF,w)
		_np2D.Fill(npv,daF,w)

def tauplot(kJets,katau1,katau2,katau3,
Hist_tau1,Hist_tau2,Hist_tau3,Hist_tau21,Hist_tau32,
Hist_tau1_2D,Hist_tau2_2D,Hist_tau3_2D,Hist_tau21_2D,Hist_tau32_2D,dSdarkFrac,HistNorm=None):
	for ja in range(len(kJets)):

		if HistNorm == None:
			w = 1
		else:
			jet = kJets[ja]
			w = eW(jet,HistNorm)

		daF = dSdarkFrac[ja]
		tau1 = katau1[ja]
		tau2 = katau2[ja]
		tau3 = katau3[ja]
		Hist_tau1.Fill(tau1,w)
		Hist_tau2.Fill(tau2,w)
		Hist_tau3.Fill(tau3,w)
		if tau1 != 0:
			Hist_tau21.Fill(tau2/tau1,w)
			Hist_tau21_2D.Fill(tau2/tau1,daF,w)
		else:
			Hist_tau21.Fill(0.0,w)
			Hist_tau21_2D.Fill(0.0,daF,w)
		if tau2 != 0:
			Hist_tau32.Fill(tau3/tau2,w)
			Hist_tau32_2D.Fill(tau3/tau2,daF,w)
		else:
			Hist_tau32.Fill(0.0,w)
			Hist_tau32_2D.Fill(0.0,daF,w)

		Hist_tau1_2D.Fill(tau1,daF,w)
		Hist_tau2_2D.Fill(tau2,daF,w)
		Hist_tau3_2D.Fill(tau3,daF,w)

def SDFill(SDJetsList,SDMx,SDMx2D,SDdarkFrac,HistNorm=None): # similar to above but for SM and Dark jets instead of jet[N]
	for dsi in range(len(SDJetsList)):
		daf = SDdarkFrac[dsi]
		SDJet = SDJetsList[dsi]

		if HistNorm == None:
			w = 1
		else:
			w = eW(SDJet,HistNorm)

		SDMx.Fill(SDJet.M(),w)
		SDMx2D.Fill(SDJet.M(),daf,w)

def fillDPT(dF,dFHist,djets,HistNorm=None):
	for fi in range(len(dF)):
		frac = dF[fi]
		djet = djets[fi]

		if HistNorm == None:
			w = 1
		else:
			w = eW(djet,HistNorm)

		dFHist.Fill(frac,w)

def ST_Val(METv,jets):
	evtST = METv
	for jet in jets:
		evtST += jet.Pt()
	return evtST

def fill_DPhi_PP(jetCats,DPhi_PP,metPhi,DPhi_JJ_PP=None,DPhi_JM_N=None,DPhi_JM_F=None,):
	if DPhi_JJ_PP != None and len(jetCats) == 2:
		DPhiList = []
		for jet in jetCats:
			DPhiList.append(deltaPhi(jet.Phi(),metPhi))
			DPhi_PP.Fill(deltaPhi(jet.Phi(),metPhi))
		DPhiList.sort()
		DPhi_JM_N.Fill(DPhiList[0])
		DPhi_JM_F.Fill(DPhiList[1])
		DPhi_JJ_PP.Fill(deltaPhi(jetCats[0].Phi(),jetCats[1].Phi()))
	elif DPhi_JJ_PP == None:
		for jet in jetCats:
			DPhi_PP.Fill(deltaPhi(jet.Phi(),metPhi))

def closestToMET(jets,SDJets,DPhiMin,Pt_Cl,Eta_Cl,Phi_Cl,DPhi_Cl,M_Cl,
axisMajor_Cl,axisMinor_Cl,momentGirth_Cl,ptD_Cl,
tau1_Cl,tau2_Cl,tau3_Cl,tau21_Cl,tau32_Cl,mult_Cl,pmult_Cl,SDM_Cl,
axisMj,axisMn,girth,ptD,tau1,tau2,tau3,pmult):
	if len(jets) > 0:
		DPhiList = []
		for ijet in jets:
			DPhi = deltaPhi(ijet.Phi(),METPhiv)
			DPhiList.append(DPhi)
		DPhiMin.Fill(min(DPhiList))
		Cli = DPhiList.index(min(DPhiList))
		Pt_Cl.Fill(jets[Cli].Pt())
		Eta_Cl.Fill(jets[Cli].Eta())
		Phi_Cl.Fill(jets[Cli].Phi())
		DPhi_Cl.Fill(deltaPhi(jets[Cli].Phi(),METPhiv))
		M_Cl.Fill(jets[Cli].M())
		axisMajor_Cl.Fill(axisMj[Cli])
		axisMinor_Cl.Fill(axisMn[Cli])
		momentGirth_Cl.Fill(girth[Cli])
		ptD_Cl.Fill(ptD[Cli])
		tau1_Cl.Fill(tau1[Cli])
		tau2_Cl.Fill(tau2[Cli])
		tau3_Cl.Fill(tau3[Cli])
		if tau1[Cli] > 0:
			tau21_Cl.Fill(tau2[Cli]/tau1[Cli])
		if tau2[Cli] > 0:
			tau32_Cl.Fill(tau3[Cli]/tau2[Cli])
		mult_Cl.Fill(Cli+1)
		pmult_Cl.Fill(pmult[Cli])

	if len(SDJets) > 0:
		SDDPhiL = []
		for sjet in SDJets:
			DPhi = deltaPhi(sjet.Phi(),METPhiv)
			SDDPhiL.append(DPhi)
		Cli = SDDPhiL.index(min(SDDPhiL))
		SDM_Cl.Fill(SDJets[Cli].M())

def jNFill(jets,Pt_jN,Eta_jN,Phi_jN,DPhi_jN,M_jN,
axisMajor_jN,axisMinor_jN,momentGirth_jN,ptD_jN,
tau1_jN,tau2_jN,tau3_jN,tau21_jN,tau32_jN,mult_jN,pmult_jN,
axisMj,axisMn,girth,ptD,tau1,tau2,tau3,pmult,jNumList):
	for jN in jNumList:
		Pt_jN.Fill(jets[jN].Pt())
		Eta_jN.Fill(jets[jN].Eta())
		Phi_jN.Fill(jets[jN].Phi())
		DPhi_jN.Fill(deltaPhi(jets[jN].Phi(),METPhiv))
		M_jN.Fill(jets[jN].M())
		axisMajor_jN.Fill(axisMj[jN])
		axisMinor_jN.Fill(axisMn[jN])
		momentGirth_jN.Fill(girth[jN])
		ptD_jN.Fill(ptD[jN])
		tau1_jN.Fill(tau1[jN])
		tau2_jN.Fill(tau2[jN])
		tau3_jN.Fill(tau3[jN])
		if tau1[jN] > 0:
			tau21_jN.Fill(tau2[jN]/tau1[jN])
		if tau2[jN] > 0:
			tau32_jN.Fill(tau3[jN]/tau2[jN])
		mult_jN.Fill(jN+1)
		pmult_jN.Fill(pmult[jN])

def jNSDFill(SDJets,SDM_jN,jNumList):
	for jN in jNumList:
		SDM_jN.Fill(SDJets[jN].M())

# args, in order mMed, mDark, Rinv, Alpha
mMed, mDark, rInv, Alpha = sys.argv[1:]

print("mMed mDark rInv Alpha")
print("{} {} {} {}".format(mMed, mDark, rInv, Alpha))

msd_max = 200
ptsd_max = 3500
metRx = 0.45

SM_Pt = rt.TH1F("SM_Pt","SM_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
SM_Eta = rt.TH1F("SM_Eta","SM_Eta;#eta(j);Events",15, -7, 7)
SM_Phi = rt.TH1F("SM_Phi","SM_Phi;#phi(j);Events",15, -3.2, 3.2)
SM_DPhi = rt.TH1F("SM_DPhi","SM_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
SM_M = rt.TH1F("SM_M","SM_M; M(j);Events",40,0,1200)

SMM_Pt = rt.TH1F("SMM_Pt","SMM_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
SMM_Eta = rt.TH1F("SMM_Eta","SMM_Eta;#eta(j);Events",15, -7, 7)
SMM_Phi = rt.TH1F("SMM_Phi","SMM_Phi;#phi(j);Events",15, -3.2, 3.2)
SMM_DPhi = rt.TH1F("SMM_DPhi","SMM_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
SMM_M = rt.TH1F("SMM_M","SMM_M; M(j);Events",40,0,1200)

G_Pt = rt.TH1F("G_Pt","G_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
G_Eta = rt.TH1F("G_Eta","G_Eta;#eta(j);Events",15, -7, 7)
G_Phi = rt.TH1F("G_Phi","G_Phi;#phi(j);Events",15, -3.2, 3.2)
G_DPhi = rt.TH1F("G_DPhi","G_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
G_M = rt.TH1F("G_M","G_M; M(j);Events",40,0,1200)

QM_Pt = rt.TH1F("QM_Pt","QM_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
QM_Eta = rt.TH1F("QM_Eta","QM_Eta;#eta(j);Events",15, -7, 7)
QM_Phi = rt.TH1F("QM_Phi","QM_Phi;#phi(j);Events",15, -3.2, 3.2)
QM_DPhi = rt.TH1F("QM_DPhi","QM_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
QM_M = rt.TH1F("QM_M","QM_M; M(j);Events",40,0,1200)

Q_Pt = rt.TH1F("Q_Pt","Q_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
Q_Eta = rt.TH1F("Q_Eta","Q_Eta;#eta(j);Events",15, -7, 7)
Q_Phi = rt.TH1F("Q_Phi","Q_Phi;#phi(j);Events",15, -3.2, 3.2)
Q_DPhi = rt.TH1F("Q_DPhi","Q_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Q_M = rt.TH1F("Q_M","Q_M; M(j);Events",40,0,1200)

QM_QPt = rt.TH1F("QM_QPt","QM_QPt;p_{T}(j) [GeV];Events",40, 0, 3500)
QM_QEta = rt.TH1F("QM_QEta","QM_QEta;#eta(j);Events",15, -7, 7)
QM_QPhi = rt.TH1F("QM_QPhi","QM_QPhi;#phi(j);Events",15, -3.2, 3.2)
QM_QDPhi = rt.TH1F("QM_QDPhi","QM_QDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
QM_QM = rt.TH1F("QM_QM","QM_QM; M(j);Events",40,0,1200)

QM_GPt = rt.TH1F("QM_GPt","QM_GPt;p_{T}(j) [GeV];Events",40, 0, 3500)
QM_GEta = rt.TH1F("QM_GEta","QM_GEta;#eta(j);Events",15, -7, 7)
QM_GPhi = rt.TH1F("QM_GPhi","QM_GPhi;#phi(j);Events",15, -3.2, 3.2)
QM_GDPhi = rt.TH1F("QM_GDPhi","QM_GDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
QM_GM = rt.TH1F("QM_GM","QM_GM; M(j);Events",40,0,1200)

Q_GPt = rt.TH1F("Q_GPt","Q_GPt;p_{T}(j) [GeV];Events",40, 0, 3500)
Q_GEta = rt.TH1F("Q_GEta","Q_GEta;#eta(j);Events",15, -7, 7)
Q_GPhi = rt.TH1F("Q_GPhi","Q_GPhi;#phi(j);Events",15, -3.2, 3.2)
Q_GDPhi = rt.TH1F("Q_GDPhi","Q_GDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Q_GM = rt.TH1F("Q_GM","Q_GM; M(j);Events",40,0,1200)

G_SMPt = rt.TH1F("G_SMPt","G_SMPt;p_{T}(j) [GeV];Events",40, 0, 3500)
G_SMEta = rt.TH1F("G_SMEta","G_SMEta;#eta(j);Events",15, -7, 7)
G_SMPhi = rt.TH1F("G_SMPhi","G_SMPhi;#phi(j);Events",15, -3.2, 3.2)
G_SMDPhi = rt.TH1F("G_SMDPhi","G_SMDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
G_SMM = rt.TH1F("G_SMM","G_SMM; M(j);Events",40,0,1200)

QM_SMPt = rt.TH1F("QM_SMPt","QM_SMPt;p_{T}(j) [GeV];Events",40, 0, 3500)
QM_SMEta = rt.TH1F("QM_SMEta","QM_SMEta;#eta(j);Events",15, -7, 7)
QM_SMPhi = rt.TH1F("QM_SMPhi","QM_SMPhi;#phi(j);Events",15, -3.2, 3.2)
QM_SMDPhi = rt.TH1F("QM_SMDPhi","QM_SMDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
QM_SMM = rt.TH1F("QM_SMM","QM_SMM; M(j);Events",40,0,1200)

Q_SMPt = rt.TH1F("Q_SMPt","Q_SMPt;p_{T}(j) [GeV];Events",40, 0, 3500)
Q_SMEta = rt.TH1F("Q_SMEta","Q_SMEta;#eta(j);Events",15, -7, 7)
Q_SMPhi = rt.TH1F("Q_SMPhi","Q_SMPhi;#phi(j);Events",15, -3.2, 3.2)
Q_SMDPhi = rt.TH1F("Q_SMDPhi","Q_SMDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Q_SMM = rt.TH1F("Q_SMM","Q_SMM; M(j);Events",40,0,1200)

LD_lowDFPt = rt.TH1F("LD_lowDFPt","LD_lowDFPt;p_{T}(j) [GeV];Events",40, 0, 3500)
LD_lowDFEta = rt.TH1F("LD_lowDFEta","LD_lowDFEta;#eta(j);Events",15, -7, 7)
LD_lowDFPhi = rt.TH1F("LD_lowDFPhi","LD_lowDFPhi;#phi(j);Events",15, -3.2, 3.2)
LD_lowDFDPhi = rt.TH1F("LD_lowDFDPhi","LD_lowDFDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
LD_lowDFM = rt.TH1F("LD_lowDFM","LD_lowDFM; M(j);Events",40,0,1200)

LD_highDFPt = rt.TH1F("LD_highDFPt","LD_highDFPt;p_{T}(j) [GeV];Events",40, 0, 3500)
LD_highDFEta = rt.TH1F("LD_highDFEta","LD_highDFEta;#eta(j);Events",15, -7, 7)
LD_highDFPhi = rt.TH1F("LD_highDFPhi","LD_highDFPhi;#phi(j);Events",15, -3.2, 3.2)
LD_highDFDPhi = rt.TH1F("LD_highDFDPhi","LD_highDFDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
LD_highDFM = rt.TH1F("LD_highDFM","LD_highDFM; M(j);Events",40,0,1200)

LD_SMPt = rt.TH1F("LD_SMPt","LD_SMPt;p_{T}(j) [GeV];Events",40, 0, 3500)
LD_SMEta = rt.TH1F("LD_SMEta","LD_SMEta;#eta(j);Events",15, -7, 7)
LD_SMPhi = rt.TH1F("LD_SMPhi","LD_SMPhi;#phi(j);Events",15, -3.2, 3.2)
LD_SMDPhi = rt.TH1F("LD_SMDPhi","LD_SMDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
LD_SMM = rt.TH1F("LD_SMM","LD_SMM; M(j);Events",40,0,1200)

jet_SM_dpTFrac = rt.TH1F("SM_dpTFrac","SM_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_SMM_dpTFrac = rt.TH1F("SMM_dpTFrac","SMM_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_G_dpTFrac = rt.TH1F("G_dpTFrac","G_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_QM_dpTFrac = rt.TH1F("QM_dpTFrac","QM_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_Q_dpTFrac = rt.TH1F("Q_dpTFrac","Q_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_QM_QdpTFrac = rt.TH1F("QM_QdpTFrac","QM_QdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_QM_GdpTFrac = rt.TH1F("QM_GdpTFrac","QM_GdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_Q_GdpTFrac = rt.TH1F("Q_GdpTFrac","Q_GdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_G_SMdpTFrac = rt.TH1F("G_SMdpTFrac","G_SMdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_QM_SMdpTFrac = rt.TH1F("QM_SMdpTFrac","QM_SMdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_Q_SMdpTFrac = rt.TH1F("Q_SMdpTFrac","Q_SMdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_LD_lowDFdpTFrac = rt.TH1F("LD_lowDFdpTFrac","LD_lowDFdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_LD_highDFdpTFrac = rt.TH1F("LD_highDFdpTFrac","LD_highDFdpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_LD_SMdpTFrac = rt.TH1F("LD_SMdpTFrac","LD_SMdpTFrac;Dark pT Fraction;Events",40, 0, 1)

SM_mult = rt.TH1F("SM_mult","SM_mult;Number of Jets;Events",7,0,7)
SMM_mult = rt.TH1F("SMM_mult","SMM_mult;Number of Jets;Events",7,0,7)
G_mult = rt.TH1F("G_mult","G_mult;Number of Jets;Events",7,0,7)
QM_mult = rt.TH1F("QM_mult","QM_mult;Number of Jets;Events",7,0,7)
Q_mult = rt.TH1F("Q_mult","Q_mult;Number of Jets;Events",7,0,7)
QM_Qmult = rt.TH1F("QM_Qmult","QM_Qmult;Number of Jets;Events",7,0,7)
QM_Gmult = rt.TH1F("QM_Gmult","QM_Gmult;Number of Jets;Events",7,0,7)
Q_Gmult = rt.TH1F("Q_Gmult","Q_Gmult;Number of Jets;Events",7,0,7)
G_SMmult = rt.TH1F("G_SMmult","G_SMmult;Number of Jets;Events",7,0,7)
QM_SMmult = rt.TH1F("QM_SMmult","QM_SMmult;Number of Jets;Events",7,0,7)
Q_SMmult = rt.TH1F("Q_SMmult","Q_SMmult;Number of Jets;Events",7,0,7)
LD_lowDFmult = rt.TH1F("LD_lowDFmult","LD_lowDFmult;Number of Jets;Events",7,0,7)
LD_highDFmult = rt.TH1F("LD_highDFmult","LD_highDFmult;Number of Jets;Events",7,0,7)
LD_SMmult = rt.TH1F("LD_SMmult","LD_SMmult;Number of Jets;Events",7,0,7)

SM_axisMajor = rt.TH1F("SM_axisMajor","SM_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
SM_axisMinor = rt.TH1F("SM_axisMinor","SM_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
SM_momentGirth = rt.TH1F("SM_momentGirth","SM_momentGirth;girth(j);Events",40, 0, 0.5)
SM_ptD = rt.TH1F("SM_ptD","SM_ptD;p_{T}D(j);Events",40, 0, 1.2)

SMM_axisMajor = rt.TH1F("SMM_axisMajor","SMM_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
SMM_axisMinor = rt.TH1F("SMM_axisMinor","SMM_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
SMM_momentGirth = rt.TH1F("SMM_momentGirth","SMM_momentGirth;girth(j);Events",40, 0, 0.5)
SMM_ptD = rt.TH1F("SMM_ptD","SMM_ptD;p_{T}D(j);Events",40, 0, 1.2)

G_axisMajor = rt.TH1F("G_axisMajor","G_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
G_axisMinor = rt.TH1F("G_axisMinor","G_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
G_momentGirth = rt.TH1F("G_momentGirth","G_momentGirth;girth(j);Events",40, 0, 0.5)
G_ptD = rt.TH1F("G_ptD","G_ptD;p_{T}D(j);Events",40, 0, 1.2)

QM_axisMajor = rt.TH1F("QM_axisMajor","QM_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
QM_axisMinor = rt.TH1F("QM_axisMinor","QM_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_momentGirth = rt.TH1F("QM_momentGirth","QM_momentGirth;girth(j);Events",40, 0, 0.5)
QM_ptD = rt.TH1F("QM_ptD","QM_ptD;p_{T}D(j);Events",40, 0, 1.2)

Q_axisMajor = rt.TH1F("Q_axisMajor","Q_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
Q_axisMinor = rt.TH1F("Q_axisMinor","Q_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
Q_momentGirth = rt.TH1F("Q_momentGirth","Q_momentGirth;girth(j);Events",40, 0, 0.5)
Q_ptD = rt.TH1F("Q_ptD","Q_ptD;p_{T}D(j);Events",40, 0, 1.2)

QM_QaxisMajor = rt.TH1F("QM_QaxisMajor","QM_QaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
QM_QaxisMinor = rt.TH1F("QM_QaxisMinor","QM_QaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_QmomentGirth = rt.TH1F("QM_QmomentGirth","QM_QmomentGirth;girth(j);Events",40, 0, 0.5)
QM_QptD = rt.TH1F("QM_QptD","QM_QptD;p_{T}D(j);Events",40, 0, 1.2)

QM_GaxisMajor = rt.TH1F("QM_GaxisMajor","QM_GaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
QM_GaxisMinor = rt.TH1F("QM_GaxisMinor","QM_GaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_GmomentGirth = rt.TH1F("QM_GmomentGirth","QM_GmomentGirth;girth(j);Events",40, 0, 0.5)
QM_GptD = rt.TH1F("QM_GptD","QM_GptD;p_{T}D(j);Events",40, 0, 1.2)

Q_GaxisMajor = rt.TH1F("Q_GaxisMajor","Q_GaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
Q_GaxisMinor = rt.TH1F("Q_GaxisMinor","Q_GaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
Q_GmomentGirth = rt.TH1F("Q_GmomentGirth","Q_GmomentGirth;girth(j);Events",40, 0, 0.5)
Q_GptD = rt.TH1F("Q_GptD","Q_GptD;p_{T}D(j);Events",40, 0, 1.2)

G_SMaxisMajor = rt.TH1F("G_SMaxisMajor","G_SMaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
G_SMaxisMinor = rt.TH1F("G_SMaxisMinor","G_SMaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
G_SMmomentGirth = rt.TH1F("G_SMmomentGirth","G_SMmomentGirth;girth(j);Events",40, 0, 0.5)
G_SMptD = rt.TH1F("G_SMptD","G_SMptD;p_{T}D(j);Events",40, 0, 1.2)

QM_SMaxisMajor = rt.TH1F("QM_SMaxisMajor","QM_SMaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
QM_SMaxisMinor = rt.TH1F("QM_SMaxisMinor","QM_SMaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_SMmomentGirth = rt.TH1F("QM_SMmomentGirth","QM_SMmomentGirth;girth(j);Events",40, 0, 0.5)
QM_SMptD = rt.TH1F("QM_SMptD","QM_SMptD;p_{T}D(j);Events",40, 0, 1.2)

Q_SMaxisMajor = rt.TH1F("Q_SMaxisMajor","Q_SMaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
Q_SMaxisMinor = rt.TH1F("Q_SMaxisMinor","Q_SMaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
Q_SMmomentGirth = rt.TH1F("Q_SMmomentGirth","Q_SMmomentGirth;girth(j);Events",40, 0, 0.5)
Q_SMptD = rt.TH1F("Q_SMptD","Q_SMptD;p_{T}D(j);Events",40, 0, 1.2)

LD_lowDFaxisMajor = rt.TH1F("LD_lowDFaxisMajor","LD_lowDFaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
LD_lowDFaxisMinor = rt.TH1F("LD_lowDFaxisMinor","LD_lowDFaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
LD_lowDFmomentGirth = rt.TH1F("LD_lowDFmomentGirth","LD_lowDFmomentGirth;girth(j);Events",40, 0, 0.5)
LD_lowDFptD = rt.TH1F("LD_lowDFptD","LD_lowDFptD;p_{T}D(j);Events",40, 0, 1.2)

LD_highDFaxisMajor = rt.TH1F("LD_highDFaxisMajor","LD_highDFaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
LD_highDFaxisMinor = rt.TH1F("LD_highDFaxisMinor","LD_highDFaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
LD_highDFmomentGirth = rt.TH1F("LD_highDFmomentGirth","LD_highDFmomentGirth;girth(j);Events",40, 0, 0.5)
LD_highDFptD = rt.TH1F("LD_highDFptD","LD_highDFptD;p_{T}D(j);Events",40, 0, 1.2)

LD_SMaxisMajor = rt.TH1F("LD_SMaxisMajor","LD_SMaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
LD_SMaxisMinor = rt.TH1F("LD_SMaxisMinor","LD_SMaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
LD_SMmomentGirth = rt.TH1F("LD_SMmomentGirth","LD_SMmomentGirth;girth(j);Events",40, 0, 0.5)
LD_SMptD = rt.TH1F("LD_SMptD","LD_SMptD;p_{T}D(j);Events",40, 0, 1.2)

## taus
SM_tau1 = rt.TH1F("SM_tau1","SM_tau1;#tau_{1}(j);Events",40, 0, 0.8)
SM_tau2 = rt.TH1F("SM_tau2","SM_tau2;#tau_{2}(j);Events",40, 0, 0.65)
SM_tau3 = rt.TH1F("SM_tau3","SM_tau3;#tau_{3}(j);Events",40, 0, 0.35)
SM_tau21 = rt.TH1F("SM_tau21","SM_tau21;#tau_{21}(j);Events",40, 0, 1.3)
SM_tau32 = rt.TH1F("SM_tau32","SM_tau32;#tau_{32}(j);Events",40, 0, 1.3)

SMM_tau1 = rt.TH1F("SMM_tau1","SMM_tau1;#tau_{1}(j);Events",40, 0, 0.8)
SMM_tau2 = rt.TH1F("SMM_tau2","SMM_tau2;#tau_{2}(j);Events",40, 0, 0.65)
SMM_tau3 = rt.TH1F("SMM_tau3","SMM_tau3;#tau_{3}(j);Events",40, 0, 0.35)
SMM_tau21 = rt.TH1F("SMM_tau21","SMM_tau21;#tau_{21}(j);Events",40, 0, 1.3)
SMM_tau32 = rt.TH1F("SMM_tau32","SMM_tau32;#tau_{32}(j);Events",40, 0, 1.3)

G_tau1 = rt.TH1F("G_tau1","G_tau1;#tau_{1}(j);Events",40, 0, 0.8)
G_tau2 = rt.TH1F("G_tau2","G_tau2;#tau_{2}(j);Events",40, 0, 0.65)
G_tau3 = rt.TH1F("G_tau3","G_tau3;#tau_{3}(j);Events",40, 0, 0.35)
G_tau21 = rt.TH1F("G_tau21","G_tau21;#tau_{21}(j);Events",40, 0, 1.3)
G_tau32 = rt.TH1F("G_tau32","G_tau32;#tau_{32}(j);Events",40, 0, 1.3)

QM_tau1 = rt.TH1F("QM_tau1","QM_tau1;#tau_{1}(j);Events",40, 0, 0.8)
QM_tau2 = rt.TH1F("QM_tau2","QM_tau2;#tau_{2}(j);Events",40, 0, 0.65)
QM_tau3 = rt.TH1F("QM_tau3","QM_tau3;#tau_{3}(j);Events",40, 0, 0.35)
QM_tau21 = rt.TH1F("QM_tau21","QM_tau21;#tau_{21}(j);Events",40, 0, 1.3)
QM_tau32 = rt.TH1F("QM_tau32","QM_tau32;#tau_{32}(j);Events",40, 0, 1.3)

Q_tau1 = rt.TH1F("Q_tau1","Q_tau1;#tau_{1}(j);Events",40, 0, 0.8)
Q_tau2 = rt.TH1F("Q_tau2","Q_tau2;#tau_{2}(j);Events",40, 0, 0.65)
Q_tau3 = rt.TH1F("Q_tau3","Q_tau3;#tau_{3}(j);Events",40, 0, 0.35)
Q_tau21 = rt.TH1F("Q_tau21","Q_tau21;#tau_{21}(j);Events",40, 0, 1.3)
Q_tau32 = rt.TH1F("Q_tau32","Q_tau32;#tau_{32}(j);Events",40, 0, 1.3)

QM_Qtau1 = rt.TH1F("QM_Qtau1","QM_Qtau1;#tau_{1}(j);Events",40, 0, 0.8)
QM_Qtau2 = rt.TH1F("QM_Qtau2","QM_Qtau2;#tau_{2}(j);Events",40, 0, 0.65)
QM_Qtau3 = rt.TH1F("QM_Qtau3","QM_Qtau3;#tau_{3}(j);Events",40, 0, 0.35)
QM_Qtau21 = rt.TH1F("QM_Qtau21","QM_Qtau21;#tau_{21}(j);Events",40, 0, 1.3)
QM_Qtau32 = rt.TH1F("QM_Qtau32","QM_Qtau32;#tau_{32}(j);Events",40, 0, 1.3)

QM_Gtau1 = rt.TH1F("QM_Gtau1","QM_Gtau1;#tau_{1}(j);Events",40, 0, 0.8)
QM_Gtau2 = rt.TH1F("QM_Gtau2","QM_Gtau2;#tau_{2}(j);Events",40, 0, 0.65)
QM_Gtau3 = rt.TH1F("QM_Gtau3","QM_Gtau3;#tau_{3}(j);Events",40, 0, 0.35)
QM_Gtau21 = rt.TH1F("QM_Gtau21","QM_Gtau21;#tau_{21}(j);Events",40, 0, 1.3)
QM_Gtau32 = rt.TH1F("QM_Gtau32","QM_Gtau32;#tau_{32}(j);Events",40, 0, 1.3)

Q_Gtau1 = rt.TH1F("Q_Gtau1","Q_Gtau1;#tau_{1}(j);Events",40, 0, 0.8)
Q_Gtau2 = rt.TH1F("Q_Gtau2","Q_Gtau2;#tau_{2}(j);Events",40, 0, 0.65)
Q_Gtau3 = rt.TH1F("Q_Gtau3","Q_Gtau3;#tau_{3}(j);Events",40, 0, 0.35)
Q_Gtau21 = rt.TH1F("Q_Gtau21","Q_Gtau21;#tau_{21}(j);Events",40, 0, 1.3)
Q_Gtau32 = rt.TH1F("Q_Gtau32","Q_Gtau32;#tau_{32}(j);Events",40, 0, 1.3)

G_SMtau1 = rt.TH1F("G_SMtau1","G_SMtau1;#tau_{1}(j);Events",40, 0, 0.8)
G_SMtau2 = rt.TH1F("G_SMtau2","G_SMtau2;#tau_{2}(j);Events",40, 0, 0.65)
G_SMtau3 = rt.TH1F("G_SMtau3","G_SMtau3;#tau_{3}(j);Events",40, 0, 0.35)
G_SMtau21 = rt.TH1F("G_SMtau21","G_SMtau21;#tau_{21}(j);Events",40, 0, 1.3)
G_SMtau32 = rt.TH1F("G_SMtau32","G_SMtau32;#tau_{32}(j);Events",40, 0, 1.3)

QM_SMtau1 = rt.TH1F("QM_SMtau1","QM_SMtau1;#tau_{1}(j);Events",40, 0, 0.8)
QM_SMtau2 = rt.TH1F("QM_SMtau2","QM_SMtau2;#tau_{2}(j);Events",40, 0, 0.65)
QM_SMtau3 = rt.TH1F("QM_SMtau3","QM_SMtau3;#tau_{3}(j);Events",40, 0, 0.35)
QM_SMtau21 = rt.TH1F("QM_SMtau21","QM_SMtau21;#tau_{21}(j);Events",40, 0, 1.3)
QM_SMtau32 = rt.TH1F("QM_SMtau32","QM_SMtau32;#tau_{32}(j);Events",40, 0, 1.3)

Q_SMtau1 = rt.TH1F("Q_SMtau1","Q_SMtau1;#tau_{1}(j);Events",40, 0, 0.8)
Q_SMtau2 = rt.TH1F("Q_SMtau2","Q_SMtau2;#tau_{2}(j);Events",40, 0, 0.65)
Q_SMtau3 = rt.TH1F("Q_SMtau3","Q_SMtau3;#tau_{3}(j);Events",40, 0, 0.35)
Q_SMtau21 = rt.TH1F("Q_SMtau21","Q_SMtau21;#tau_{21}(j);Events",40, 0, 1.3)
Q_SMtau32 = rt.TH1F("Q_SMtau32","Q_SMtau32;#tau_{32}(j);Events",40, 0, 1.3)

LD_lowDFtau1 = rt.TH1F("LD_lowDFtau1","LD_lowDFtau1;#tau_{1}(j);Events",40, 0, 0.8)
LD_lowDFtau2 = rt.TH1F("LD_lowDFtau2","LD_lowDFtau2;#tau_{2}(j);Events",40, 0, 0.65)
LD_lowDFtau3 = rt.TH1F("LD_lowDFtau3","LD_lowDFtau3;#tau_{3}(j);Events",40, 0, 0.35)
LD_lowDFtau21 = rt.TH1F("LD_lowDFtau21","LD_lowDFtau21;#tau_{21}(j);Events",40, 0, 1.3)
LD_lowDFtau32 = rt.TH1F("LD_lowDFtau32","LD_lowDFtau32;#tau_{32}(j);Events",40, 0, 1.3)

LD_highDFtau1 = rt.TH1F("LD_highDFtau1","LD_highDFtau1;#tau_{1}(j);Events",40, 0, 0.8)
LD_highDFtau2 = rt.TH1F("LD_highDFtau2","LD_highDFtau2;#tau_{2}(j);Events",40, 0, 0.65)
LD_highDFtau3 = rt.TH1F("LD_highDFtau3","LD_highDFtau3;#tau_{3}(j);Events",40, 0, 0.35)
LD_highDFtau21 = rt.TH1F("LD_highDFtau21","LD_highDFtau21;#tau_{21}(j);Events",40, 0, 1.3)
LD_highDFtau32 = rt.TH1F("LD_highDFtau32","LD_highDFtau32;#tau_{32}(j);Events",40, 0, 1.3)

LD_SMtau1 = rt.TH1F("LD_SMtau1","LD_SMtau1;#tau_{1}(j);Events",40, 0, 0.8)
LD_SMtau2 = rt.TH1F("LD_SMtau2","LD_SMtau2;#tau_{2}(j);Events",40, 0, 0.65)
LD_SMtau3 = rt.TH1F("LD_SMtau3","LD_SMtau3;#tau_{3}(j);Events",40, 0, 0.35)
LD_SMtau21 = rt.TH1F("LD_SMtau21","LD_SMtau21;#tau_{21}(j);Events",40, 0, 1.3)
LD_SMtau32 = rt.TH1F("LD_SMtau32","LD_SMtau32;#tau_{32}(j);Events",40, 0, 1.3)

## softdrop
SM_SDM = rt.TH1F("SM_SDM","SM_SDM; m_{SD}(j);Events",40,0,msd_max)
SMM_SDM = rt.TH1F("SMM_SDM","SMM_SDM; m_{SD}(j);Events",40,0,msd_max)
G_SDM = rt.TH1F("G_SDM","G_SDM; m_{SD}(j);Events",40,0,msd_max)
QM_SDM = rt.TH1F("QM_SDM","QM_SDM; m_{SD}(j);Events",40,0,msd_max)
Q_SDM = rt.TH1F("Q_SDM","Q_SDM; m_{SD}(j);Events",40,0,msd_max)
QM_QSDM = rt.TH1F("QM_QSDM","QM_QSDM; m_{SD}(j);Events",40,0,msd_max)
QM_GSDM = rt.TH1F("QM_GSDM","QM_GSDM; m_{SD}(j);Events",40,0,msd_max)
Q_GSDM = rt.TH1F("Q_GSDM","Q_GSDM; m_{SD}(j);Events",40,0,msd_max)
G_SMSDM = rt.TH1F("G_SMSDM","G_SMSDM; m_{SD}(j);Events",40,0,msd_max)
QM_SMSDM = rt.TH1F("QM_SMSDM","QM_SMSDM; m_{SD}(j);Events",40,0,msd_max)
Q_SMSDM = rt.TH1F("Q_SMSDM","Q_SMSDM; m_{SD}(j);Events",40,0,msd_max)
LD_lowDFSDM = rt.TH1F("LD_lowDFSDM","LD_lowDFSDM; m_{SD}(j);Events",40,0,msd_max)
LD_highDFSDM = rt.TH1F("LD_highDFSDM","LD_highDFSDM; m_{SD}(j);Events",40,0,msd_max)
LD_SMSDM = rt.TH1F("LD_SMSDM","LD_SMSDM; m_{SD}(j);Events",40,0,msd_max)

SM_SDPt = rt.TH1F("SM_SDPt","SM_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
SMM_SDPt = rt.TH1F("SMM_SDPt","SMM_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
G_SDPt = rt.TH1F("G_SDPt","G_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
QM_SDPt = rt.TH1F("QM_SDPt","QM_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
Q_SDPt = rt.TH1F("Q_SDPt","Q_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
QM_QSDPt = rt.TH1F("QM_QSDPt","QM_QSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
QM_GSDPt = rt.TH1F("QM_GSDPt","QM_GSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
Q_GSDPt = rt.TH1F("Q_GSDPt","Q_GSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
G_SMSDPt = rt.TH1F("G_SMSDPt","G_SMSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
QM_SMSDPt = rt.TH1F("QM_SMSDPt","QM_SMSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
Q_SMSDPt = rt.TH1F("Q_SMSDPt","Q_SMSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
LD_lowDFSDPt = rt.TH1F("LD_lowDFSDPt","LD_lowDFSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
LD_highDFSDPt = rt.TH1F("LD_highDFSDPt","LD_highDFSDPt; pT_{SD}(j);Events",40,0,ptsd_max)
LD_SMSDPt = rt.TH1F("LD_SMSDPt","LD_SMSDPt; pT_{SD}(j);Events",40,0,ptsd_max)

# number of particles in each jet
SM_pmult = rt.TH1F("SM_pmult","SM_pmult;Number of Particles in a Jet;Events",40,0,400)
SMM_pmult = rt.TH1F("SMM_pmult","SMM_pmult;Number of Particles in a Jet;Events",40,0,400)
G_pmult = rt.TH1F("G_pmult","G_pmult;Number of Particles in a Jet;Events",40,0,400)
QM_pmult = rt.TH1F("QM_pmult","QM_pmult;Number of Particles in a Jet;Events",40,0,400)
Q_pmult = rt.TH1F("Q_pmult","Q_pmult;Number of Particles in a Jet;Events",40,0,400)
QM_Qpmult = rt.TH1F("QM_Qpmult","QM_Qpmult;Number of Particles in a Jet;Events",40,0,400)
QM_Gpmult = rt.TH1F("QM_Gpmult","QM_Gpmult;Number of Particles in a Jet;Events",40,0,400)
Q_Gpmult = rt.TH1F("Q_Gpmult","Q_Gpmult;Number of Particles in a Jet;Events",40,0,400)
G_SMpmult = rt.TH1F("G_SMpmult","G_SMpmult;Number of Particles in a Jet;Events",40,0,400)
QM_SMpmult = rt.TH1F("QM_SMpmult","QM_SMpmult;Number of Particles in a Jet;Events",40,0,400)
Q_SMpmult = rt.TH1F("Q_SMpmult","Q_SMpmult;Number of Particles in a Jet;Events",40,0,400)
LD_lowDFpmult = rt.TH1F("LD_lowDFpmult","LD_lowDFpmult;Number of Particles in a Jet;Events",40,0,400)
LD_highDFpmult = rt.TH1F("LD_highDFpmult","LD_highDFpmult;Number of Particles in a Jet;Events",40,0,400)
LD_SMpmult = rt.TH1F("LD_SMpmult","LD_SMpmult;Number of Particles in a Jet;Events",40,0,400)

# 2D plots (kinematic vs dark pT fraction)
SM_Pt_2D = rt.TH2F("SM_Pt_2D","SM_Pt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
SM_Eta_2D = rt.TH2F("SM_Eta_2D","SM_Eta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
SM_Phi_2D = rt.TH2F("SM_Phi_2D","SM_Phi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
SM_DPhi_2D = rt.TH2F("SM_DPhi_2D","SM_DPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

SMM_Pt_2D = rt.TH2F("SMM_Pt_2D","SMM_Pt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
SMM_Eta_2D = rt.TH2F("SMM_Eta_2D","SMM_Eta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
SMM_Phi_2D = rt.TH2F("SMM_Phi_2D","SMM_Phi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
SMM_DPhi_2D = rt.TH2F("SMM_DPhi_2D","SMM_DPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

G_Pt_2D = rt.TH2F("G_Pt_2D","G_Pt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
G_Eta_2D = rt.TH2F("G_Eta_2D","G_Eta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
G_Phi_2D = rt.TH2F("G_Phi_2D","G_Phi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
G_DPhi_2D = rt.TH2F("G_DPhi_2D","G_DPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

QM_Pt_2D = rt.TH2F("QM_Pt_2D","QM_Pt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
QM_Eta_2D = rt.TH2F("QM_Eta_2D","QM_Eta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
QM_Phi_2D = rt.TH2F("QM_Phi_2D","QM_Phi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
QM_DPhi_2D = rt.TH2F("QM_DPhi_2D","QM_DPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

Q_Pt_2D = rt.TH2F("Q_Pt_2D","Q_Pt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
Q_Eta_2D = rt.TH2F("Q_Eta_2D","Q_Eta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
Q_Phi_2D = rt.TH2F("Q_Phi_2D","Q_Phi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
Q_DPhi_2D = rt.TH2F("Q_DPhi_2D","Q_DPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

QM_QPt_2D = rt.TH2F("QM_QPt_2D","QM_QPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
QM_QEta_2D = rt.TH2F("QM_QEta_2D","QM_QEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
QM_QPhi_2D = rt.TH2F("QM_QPhi_2D","QM_QPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
QM_QDPhi_2D = rt.TH2F("QM_QDPhi_2D","QM_QDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

QM_GPt_2D = rt.TH2F("QM_GPt_2D","QM_GPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
QM_GEta_2D = rt.TH2F("QM_GEta_2D","QM_GEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
QM_GPhi_2D = rt.TH2F("QM_GPhi_2D","QM_GPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
QM_GDPhi_2D = rt.TH2F("QM_GDPhi_2D","QM_GDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

Q_GPt_2D = rt.TH2F("Q_GPt_2D","Q_GPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
Q_GEta_2D = rt.TH2F("Q_GEta_2D","Q_GEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
Q_GPhi_2D = rt.TH2F("Q_GPhi_2D","Q_GPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
Q_GDPhi_2D = rt.TH2F("Q_GDPhi_2D","Q_GDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

G_SMPt_2D = rt.TH2F("G_SMPt_2D","G_SMPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
G_SMEta_2D = rt.TH2F("G_SMEta_2D","G_SMEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
G_SMPhi_2D = rt.TH2F("G_SMPhi_2D","G_SMPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
G_SMDPhi_2D = rt.TH2F("G_SMDPhi_2D","G_SMDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

QM_SMPt_2D = rt.TH2F("QM_SMPt_2D","QM_SMPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
QM_SMEta_2D = rt.TH2F("QM_SMEta_2D","QM_SMEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
QM_SMPhi_2D = rt.TH2F("QM_SMPhi_2D","QM_SMPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
QM_SMDPhi_2D = rt.TH2F("QM_SMDPhi_2D","QM_SMDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

Q_SMPt_2D = rt.TH2F("Q_SMPt_2D","Q_SMPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
Q_SMEta_2D = rt.TH2F("Q_SMEta_2D","Q_SMEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
Q_SMPhi_2D = rt.TH2F("Q_SMPhi_2D","Q_SMPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
Q_SMDPhi_2D = rt.TH2F("Q_SMDPhi_2D","Q_SMDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

LD_lowDFPt_2D = rt.TH2F("LD_lowDFPt_2D","LD_lowDFPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
LD_lowDFEta_2D = rt.TH2F("LD_lowDFEta_2D","LD_lowDFEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
LD_lowDFPhi_2D = rt.TH2F("LD_lowDFPhi_2D","LD_lowDFPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
LD_lowDFDPhi_2D = rt.TH2F("LD_lowDFDPhi_2D","LD_lowDFDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

LD_highDFPt_2D = rt.TH2F("LD_highDFPt_2D","LD_highDFPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
LD_highDFEta_2D = rt.TH2F("LD_highDFEta_2D","LD_highDFEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
LD_highDFPhi_2D = rt.TH2F("LD_highDFPhi_2D","LD_highDFPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
LD_highDFDPhi_2D = rt.TH2F("LD_highDFDPhi_2D","LD_highDFDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

LD_SMPt_2D = rt.TH2F("LD_SMPt_2D","LD_SMPt_2D;p_{T}(j) [GeV];Events",40, 0, 3500, 40, 0, 1)
LD_SMEta_2D = rt.TH2F("LD_SMEta_2D","LD_SMEta_2D;#eta(j);Events",15, -7, 7, 40, 0, 1)
LD_SMPhi_2D = rt.TH2F("LD_SMPhi_2D","LD_SMPhi_2D;#phi(j);Events",15, -3.2, 3.2, 40, 0, 1)
LD_SMDPhi_2D = rt.TH2F("LD_SMDPhi_2D","LD_SMDPhi_2D; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5, 40, 0, 1)

SM_axisMajor_2D = rt.TH2F("SM_axisMajor_2D","SM_axisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
SM_axisMinor_2D = rt.TH2F("SM_axisMinor_2D","SM_axisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
SM_momentGirth_2D = rt.TH2F("SM_momentGirth_2D","SM_momentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
SM_ptD_2D = rt.TH2F("SM_ptD_2D","SM_ptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

SMM_axisMajor_2D = rt.TH2F("SMM_axisMajor_2D","SMM_axisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
SMM_axisMinor_2D = rt.TH2F("SMM_axisMinor_2D","SMM_axisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
SMM_momentGirth_2D = rt.TH2F("SMM_momentGirth_2D","SMM_momentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
SMM_ptD_2D = rt.TH2F("SMM_ptD_2D","SMM_ptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

G_axisMajor_2D = rt.TH2F("G_axisMajor_2D","G_axisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
G_axisMinor_2D = rt.TH2F("G_axisMinor_2D","G_axisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
G_momentGirth_2D = rt.TH2F("G_momentGirth_2D","G_momentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
G_ptD_2D = rt.TH2F("G_ptD_2D","G_ptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

QM_axisMajor_2D = rt.TH2F("QM_axisMajor_2D","QM_axisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
QM_axisMinor_2D = rt.TH2F("QM_axisMinor_2D","QM_axisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
QM_momentGirth_2D = rt.TH2F("QM_momentGirth_2D","QM_momentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
QM_ptD_2D = rt.TH2F("QM_ptD_2D","QM_ptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

Q_axisMajor_2D = rt.TH2F("Q_axisMajor_2D","Q_axisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
Q_axisMinor_2D = rt.TH2F("Q_axisMinor_2D","Q_axisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
Q_momentGirth_2D = rt.TH2F("Q_momentGirth_2D","Q_momentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
Q_ptD_2D = rt.TH2F("Q_ptD_2D","Q_ptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

QM_QaxisMajor_2D = rt.TH2F("QM_QaxisMajor_2D","QM_QaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
QM_QaxisMinor_2D = rt.TH2F("QM_QaxisMinor_2D","QM_QaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
QM_QmomentGirth_2D = rt.TH2F("QM_QmomentGirth_2D","QM_QmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
QM_QptD_2D = rt.TH2F("QM_QptD_2D","QM_QptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

QM_GaxisMajor_2D = rt.TH2F("QM_GaxisMajor_2D","QM_GaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
QM_GaxisMinor_2D = rt.TH2F("QM_GaxisMinor_2D","QM_GaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
QM_GmomentGirth_2D = rt.TH2F("QM_GmomentGirth_2D","QM_GmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
QM_GptD_2D = rt.TH2F("QM_GptD_2D","QM_GptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

Q_GaxisMajor_2D = rt.TH2F("Q_GaxisMajor_2D","Q_GaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
Q_GaxisMinor_2D = rt.TH2F("Q_GaxisMinor_2D","Q_GaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
Q_GmomentGirth_2D = rt.TH2F("Q_GmomentGirth_2D","Q_GmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
Q_GptD_2D = rt.TH2F("Q_GptD_2D","Q_GptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

G_SMaxisMajor_2D = rt.TH2F("G_SMaxisMajor_2D","G_SMaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
G_SMaxisMinor_2D = rt.TH2F("G_SMaxisMinor_2D","G_SMaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
G_SMmomentGirth_2D = rt.TH2F("G_SMmomentGirth_2D","G_SMmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
G_SMptD_2D = rt.TH2F("G_SMptD_2D","G_SMptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

QM_SMaxisMajor_2D = rt.TH2F("QM_SMaxisMajor_2D","QM_SMaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
QM_SMaxisMinor_2D = rt.TH2F("QM_SMaxisMinor_2D","QM_SMaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
QM_SMmomentGirth_2D = rt.TH2F("QM_SMmomentGirth_2D","QM_SMmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
QM_SMptD_2D = rt.TH2F("QM_SMptD_2D","QM_SMptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

Q_SMaxisMajor_2D = rt.TH2F("Q_SMaxisMajor_2D","Q_SMaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
Q_SMaxisMinor_2D = rt.TH2F("Q_SMaxisMinor_2D","Q_SMaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
Q_SMmomentGirth_2D = rt.TH2F("Q_SMmomentGirth_2D","Q_SMmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
Q_SMptD_2D = rt.TH2F("Q_SMptD_2D","Q_SMptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

LD_lowDFaxisMajor_2D = rt.TH2F("LD_lowDFaxisMajor_2D","LD_lowDFaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
LD_lowDFaxisMinor_2D = rt.TH2F("LD_lowDFaxisMinor_2D","LD_lowDFaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
LD_lowDFmomentGirth_2D = rt.TH2F("LD_lowDFmomentGirth_2D","LD_lowDFmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
LD_lowDFptD_2D = rt.TH2F("LD_lowDFptD_2D","LD_lowDFptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

LD_highDFaxisMajor_2D = rt.TH2F("LD_highDFaxisMajor_2D","LD_highDFaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
LD_highDFaxisMinor_2D = rt.TH2F("LD_highDFaxisMinor_2D","LD_highDFaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
LD_highDFmomentGirth_2D = rt.TH2F("LD_highDFmomentGirth_2D","LD_highDFmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
LD_highDFptD_2D = rt.TH2F("LD_highDFptD_2D","LD_highDFptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

LD_SMaxisMajor_2D = rt.TH2F("LD_SMaxisMajor_2D","LD_SMaxisMajor_2D;#sigma_{major}(j);Events",40, 0, 0.5, 40, 0, 1)
LD_SMaxisMinor_2D = rt.TH2F("LD_SMaxisMinor_2D","LD_SMaxisMinor_2D;#sigma_{minor}(j);Events",40, 0, 0.3, 40, 0, 1)
LD_SMmomentGirth_2D = rt.TH2F("LD_SMmomentGirth_2D","LD_SMmomentGirth_2D;girth(j);Events",40, 0, 0.5, 40, 0, 1)
LD_SMptD_2D = rt.TH2F("LD_SMptD_2D","LD_SMptD_2D;p_{T}D(j);Events",40, 0, 1.2, 40, 0, 1)

SM_M_2D = rt.TH2F("SM_M_2D","SM_M_2D; M(j);Events",40,0,1200, 40, 0, 1)
SMM_M_2D = rt.TH2F("SMM_M_2D","SMM_M_2D; M(j);Events",40,0,1200, 40, 0, 1)
G_M_2D = rt.TH2F("G_M_2D","G_M_2D; M(j);Events",40,0,1200, 40, 0, 1)
QM_M_2D = rt.TH2F("QM_M_2D","QM_M_2D; M(j);Events",40,0,1200, 40, 0, 1)
Q_M_2D = rt.TH2F("Q_M_2D","Q_M_2D; M(j);Events",40,0,1200, 40, 0, 1)
QM_QM_2D = rt.TH2F("QM_QM_2D","QM_QM_2D; M(j);Events",40,0,1200, 40, 0, 1)
QM_GM_2D = rt.TH2F("QM_GM_2D","QM_GM_2D; M(j);Events",40,0,1200, 40, 0, 1)
Q_GM_2D = rt.TH2F("Q_GM_2D","Q_GM_2D; M(j);Events",40,0,1200, 40, 0, 1)
G_SMM_2D = rt.TH2F("G_SMM_2D","G_SMM_2D; M(j);Events",40,0,1200, 40, 0, 1)
QM_SMM_2D = rt.TH2F("QM_SMM_2D","QM_SMM_2D; M(j);Events",40,0,1200, 40, 0, 1)
Q_SMM_2D = rt.TH2F("Q_SMM_2D","Q_SMM_2D; M(j);Events",40,0,1200, 40, 0, 1)
LD_lowDFM_2D = rt.TH2F("LD_lowDFM_2D","LD_lowDFM_2D; M(j);Events",40,0,1200, 40, 0, 1)
LD_highDFM_2D = rt.TH2F("LD_highDFM_2D","LD_highDFM_2D; M(j);Events",40,0,1200, 40, 0, 1)
LD_SMM_2D = rt.TH2F("LD_SMM_2D","LD_SMM_2D; M(j);Events",40,0,1200, 40, 0, 1)

SM_tau1_2D = rt.TH2F("SM_tau1_2D","SM_tau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
SM_tau2_2D = rt.TH2F("SM_tau2_2D","SM_tau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
SM_tau3_2D = rt.TH2F("SM_tau3_2D","SM_tau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
SM_tau21_2D = rt.TH2F("SM_tau21_2D","SM_tau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
SM_tau32_2D = rt.TH2F("SM_tau32_2D","SM_tau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
SM_SDM_2D = rt.TH2F("SM_SDM_2D","SM_SDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

SMM_tau1_2D = rt.TH2F("SMM_tau1_2D","SMM_tau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
SMM_tau2_2D = rt.TH2F("SMM_tau2_2D","SMM_tau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
SMM_tau3_2D = rt.TH2F("SMM_tau3_2D","SMM_tau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
SMM_tau21_2D = rt.TH2F("SMM_tau21_2D","SMM_tau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
SMM_tau32_2D = rt.TH2F("SMM_tau32_2D","SMM_tau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
SMM_SDM_2D = rt.TH2F("SMM_SDM_2D","SMM_SDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

G_tau1_2D = rt.TH2F("G_tau1_2D","G_tau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
G_tau2_2D = rt.TH2F("G_tau2_2D","G_tau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
G_tau3_2D = rt.TH2F("G_tau3_2D","G_tau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
G_tau21_2D = rt.TH2F("G_tau21_2D","G_tau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
G_tau32_2D = rt.TH2F("G_tau32_2D","G_tau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
G_SDM_2D = rt.TH2F("G_SDM_2D","G_SDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

QM_tau1_2D = rt.TH2F("QM_tau1_2D","QM_tau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
QM_tau2_2D = rt.TH2F("QM_tau2_2D","QM_tau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
QM_tau3_2D = rt.TH2F("QM_tau3_2D","QM_tau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
QM_tau21_2D = rt.TH2F("QM_tau21_2D","QM_tau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_tau32_2D = rt.TH2F("QM_tau32_2D","QM_tau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_SDM_2D = rt.TH2F("QM_SDM_2D","QM_SDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

Q_tau1_2D = rt.TH2F("Q_tau1_2D","Q_tau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
Q_tau2_2D = rt.TH2F("Q_tau2_2D","Q_tau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
Q_tau3_2D = rt.TH2F("Q_tau3_2D","Q_tau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
Q_tau21_2D = rt.TH2F("Q_tau21_2D","Q_tau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
Q_tau32_2D = rt.TH2F("Q_tau32_2D","Q_tau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
Q_SDM_2D = rt.TH2F("Q_SDM_2D","Q_SDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

QM_Qtau1_2D = rt.TH2F("QM_Qtau1_2D","QM_Qtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
QM_Qtau2_2D = rt.TH2F("QM_Qtau2_2D","QM_Qtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
QM_Qtau3_2D = rt.TH2F("QM_Qtau3_2D","QM_Qtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
QM_Qtau21_2D = rt.TH2F("QM_Qtau21_2D","QM_Qtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_Qtau32_2D = rt.TH2F("QM_Qtau32_2D","QM_Qtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_QSDM_2D = rt.TH2F("QM_QSDM_2D","QM_QSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

QM_Gtau1_2D = rt.TH2F("QM_Gtau1_2D","QM_Gtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
QM_Gtau2_2D = rt.TH2F("QM_Gtau2_2D","QM_Gtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
QM_Gtau3_2D = rt.TH2F("QM_Gtau3_2D","QM_Gtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
QM_Gtau21_2D = rt.TH2F("QM_Gtau21_2D","QM_Gtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_Gtau32_2D = rt.TH2F("QM_Gtau32_2D","QM_Gtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_GSDM_2D = rt.TH2F("QM_GSDM_2D","QM_GSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

Q_Gtau1_2D = rt.TH2F("Q_Gtau1_2D","Q_Gtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
Q_Gtau2_2D = rt.TH2F("Q_Gtau2_2D","Q_Gtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
Q_Gtau3_2D = rt.TH2F("Q_Gtau3_2D","Q_Gtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
Q_Gtau21_2D = rt.TH2F("Q_Gtau21_2D","Q_Gtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
Q_Gtau32_2D = rt.TH2F("Q_Gtau32_2D","Q_Gtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
Q_GSDM_2D = rt.TH2F("Q_GSDM_2D","Q_GSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

G_SMtau1_2D = rt.TH2F("G_SMtau1_2D","G_SMtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
G_SMtau2_2D = rt.TH2F("G_SMtau2_2D","G_SMtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
G_SMtau3_2D = rt.TH2F("G_SMtau3_2D","G_SMtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
G_SMtau21_2D = rt.TH2F("G_SMtau21_2D","G_SMtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
G_SMtau32_2D = rt.TH2F("G_SMtau32_2D","G_SMtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
G_SMSDM_2D = rt.TH2F("G_SMSDM_2D","G_SMSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

QM_SMtau1_2D = rt.TH2F("QM_SMtau1_2D","QM_SMtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
QM_SMtau2_2D = rt.TH2F("QM_SMtau2_2D","QM_SMtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
QM_SMtau3_2D = rt.TH2F("QM_SMtau3_2D","QM_SMtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
QM_SMtau21_2D = rt.TH2F("QM_SMtau21_2D","QM_SMtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_SMtau32_2D = rt.TH2F("QM_SMtau32_2D","QM_SMtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
QM_SMSDM_2D = rt.TH2F("QM_SMSDM_2D","QM_SMSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

Q_SMtau1_2D = rt.TH2F("Q_SMtau1_2D","Q_SMtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
Q_SMtau2_2D = rt.TH2F("Q_SMtau2_2D","Q_SMtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
Q_SMtau3_2D = rt.TH2F("Q_SMtau3_2D","Q_SMtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
Q_SMtau21_2D = rt.TH2F("Q_SMtau21_2D","Q_SMtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
Q_SMtau32_2D = rt.TH2F("Q_SMtau32_2D","Q_SMtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
Q_SMSDM_2D = rt.TH2F("Q_SMSDM_2D","Q_SMSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

LD_lowDFtau1_2D = rt.TH2F("LD_lowDFtau1_2D","LD_lowDFtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
LD_lowDFtau2_2D = rt.TH2F("LD_lowDFtau2_2D","LD_lowDFtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
LD_lowDFtau3_2D = rt.TH2F("LD_lowDFtau3_2D","LD_lowDFtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
LD_lowDFtau21_2D = rt.TH2F("LD_lowDFtau21_2D","LD_lowDFtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
LD_lowDFtau32_2D = rt.TH2F("LD_lowDFtau32_2D","LD_lowDFtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
LD_lowDFSDM_2D = rt.TH2F("LD_lowDFSDM_2D","LD_lowDFSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

LD_highDFtau1_2D = rt.TH2F("LD_highDFtau1_2D","LD_highDFtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
LD_highDFtau2_2D = rt.TH2F("LD_highDFtau2_2D","LD_highDFtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
LD_highDFtau3_2D = rt.TH2F("LD_highDFtau3_2D","LD_highDFtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
LD_highDFtau21_2D = rt.TH2F("LD_highDFtau21_2D","LD_highDFtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
LD_highDFtau32_2D = rt.TH2F("LD_highDFtau32_2D","LD_highDFtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
LD_highDFSDM_2D = rt.TH2F("LD_highDFSDM_2D","LD_highDFSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

LD_SMtau1_2D = rt.TH2F("LD_SMtau1_2D","LD_SMtau1_2D;#tau_{1}(j);Events",40, 0, 0.8, 40, 0, 1)
LD_SMtau2_2D = rt.TH2F("LD_SMtau2_2D","LD_SMtau2_2D;#tau_{2}(j);Events",40, 0, 0.65, 40, 0, 1)
LD_SMtau3_2D = rt.TH2F("LD_SMtau3_2D","LD_SMtau3_2D;#tau_{3}(j);Events",40, 0, 0.35, 40, 0, 1)
LD_SMtau21_2D = rt.TH2F("LD_SMtau21_2D","LD_SMtau21_2D;#tau_{21}(j);Events",40, 0, 1.3, 40, 0, 1)
LD_SMtau32_2D = rt.TH2F("LD_SMtau32_2D","LD_SMtau32_2D;#tau_{32}(j);Events",40, 0, 1.3, 40, 0, 1)
LD_SMSDM_2D = rt.TH2F("LD_SMSDM_2D","LD_SMSDM_2D;m_{SD}(j);Events",40,0,msd_max, 40, 0, 1)

SM_pmult_2D = rt.TH2F("SM_pmult_2D","SM_pmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
SMM_pmult_2D = rt.TH2F("SMM_pmult_2D","SMM_pmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
G_pmult_2D = rt.TH2F("G_pmult_2D","G_pmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
QM_pmult_2D = rt.TH2F("QM_pmult_2D","QM_pmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
Q_pmult_2D = rt.TH2F("Q_pmult_2D","Q_pmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
QM_Qpmult_2D = rt.TH2F("QM_Qpmult_2D","QM_Qpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
QM_Gpmult_2D = rt.TH2F("QM_Gpmult_2D","QM_Gpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
Q_Gpmult_2D = rt.TH2F("Q_Gpmult_2D","Q_Gpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
G_SMpmult_2D = rt.TH2F("G_SMpmult_2D","G_SMpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
QM_SMpmult_2D = rt.TH2F("QM_SMpmult_2D","QM_SMpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
Q_SMpmult_2D = rt.TH2F("Q_SMpmult_2D","Q_SMpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
LD_lowDFpmult_2D = rt.TH2F("LD_lowDFpmult_2D","LD_lowDFpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
LD_highDFpmult_2D = rt.TH2F("LD_highDFpmult_2D","LD_highDFpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)
LD_SMpmult_2D = rt.TH2F("LD_SMpmult_2D","LD_SMpmult_2D;Number of Particles in a Jet;Events",40,0,400, 40, 0, 1)

## DPhi for pair production processes only
SM_DPhi_PP = rt.TH1F("SM_DPhi_PP","SM_DPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
SMM_DPhi_PP = rt.TH1F("SMM_DPhi_PP","SMM_DPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
G_DPhi_PP = rt.TH1F("G_DPhi_PP","G_DPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
QM_DPhi_PP = rt.TH1F("QM_DPhi_PP","QM_DPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
Q_DPhi_PP = rt.TH1F("Q_DPhi_PP","Q_DPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
QM_QDPhi_PP = rt.TH1F("QM_QDPhi_PP","QM_QDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
QM_GDPhi_PP = rt.TH1F("QM_GDPhi_PP","QM_GDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
Q_GDPhi_PP = rt.TH1F("Q_GDPhi_PP","Q_GDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
G_SMDPhi_PP = rt.TH1F("G_SMDPhi_PP","G_SMDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
QM_SMDPhi_PP = rt.TH1F("QM_SMDPhi_PP","QM_SMDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
Q_SMDPhi_PP = rt.TH1F("Q_SMDPhi_PP","Q_SMDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
LD_lowDFDPhi_PP = rt.TH1F("LD_lowDFDPhi_PP","LD_lowDFDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
LD_highDFDPhi_PP = rt.TH1F("LD_highDFDPhi_PP","LD_highDFDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
LD_SMDPhi_PP = rt.TH1F("LD_SMDPhi_PP","LD_SMDPhi_PP; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)

QM_Dhi_JJ_PP = rt.TH1F("QM_Dhi_JJ_PP","QM_Dhi_JJ_PP; #Delta#phi(j_{1},j_{2}) (PP);Events",30,0,3.5)
QM_Dhi_JMET_Near = rt.TH1F("QM_Dhi_JMET_Near","QM_Dhi_JMET_Near; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
QM_Dhi_JMET_Far = rt.TH1F("QM_Dhi_JMET_Far","QM_Dhi_JMET_Far; #Delta#phi(j,#slash{E}_{T}) (PP);Events",30,0,3.5)
# Looser, more general categories, easier to interpret and use

# SM_Only: SM_ + SMM_
# Mix_: G_SM, LD_SM
# Dark_: QM_, QM_G, Q_G
# SoftDark_: G_, LD_highDF

SM_OnlyPt = rt.TH1F("SM_OnlyPt","SM_OnlyPt;p_{T}(j) [GeV];Events",40, 0, 3500)
SM_OnlyEta = rt.TH1F("SM_OnlyEta","SM_OnlyEta;#eta(j);Events",15, -7, 7)
SM_OnlyPhi = rt.TH1F("SM_OnlyPhi","SM_OnlyPhi;#phi(j);Events",15, -3.2, 3.2)
SM_OnlyDPhi = rt.TH1F("SM_OnlyDPhi","SM_OnlyDPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
SM_OnlyDPhiMin = rt.TH1F("SM_OnlyDPhiMin","SM_OnlyDPhiMin; #Delta#phi_{min}(j,#slash{E}_{T});Events",30,0,3.5)
SM_OnlyM = rt.TH1F("SM_OnlyM","SM_OnlyM; M(j);Events",40,0,1200)
SM_Onlymult = rt.TH1F("SM_Onlymult","SM_Onlymult;Number of Jets;Events",7,0,7)
SM_OnlyaxisMajor = rt.TH1F("SM_OnlyaxisMajor","SM_OnlyaxisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
SM_OnlyaxisMinor = rt.TH1F("SM_OnlyaxisMinor","SM_OnlyaxisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
SM_OnlymomentGirth = rt.TH1F("SM_OnlymomentGirth","SM_OnlymomentGirth;girth(j);Events",40, 0, 0.5)
SM_OnlyptD = rt.TH1F("SM_OnlyptD","SM_OnlyptD;p_{T}D(j);Events",40, 0, 1.2)
SM_Onlytau1 = rt.TH1F("SM_Onlytau1","SM_Onlytau1;#tau_{1}(j);Events",40, 0, 0.8)
SM_Onlytau2 = rt.TH1F("SM_Onlytau2","SM_Onlytau2;#tau_{2}(j);Events",40, 0, 0.65)
SM_Onlytau3 = rt.TH1F("SM_Onlytau3","SM_Onlytau3;#tau_{3}(j);Events",40, 0, 0.35)
SM_Onlytau21 = rt.TH1F("SM_Onlytau21","SM_Onlytau21;#tau_{21}(j);Events",40, 0, 1.3)
SM_Onlytau32 = rt.TH1F("SM_Onlytau32","SM_Onlytau32;#tau_{32}(j);Events",40, 0, 1.3)
SM_OnlySDM = rt.TH1F("SM_OnlySDM","SM_OnlySDM; m_{SD}(j);Events",40,0,msd_max)
SM_OnlySDPt = rt.TH1F("SM_OnlySDPt","SM_OnlySDPt; pT_{SD}(j);Events",40,0,ptsd_max)
SM_Onlypmult = rt.TH1F("SM_Onlypmult","SM_Onlypmult;Number of Particles in a Jet;Events",40,0,400)

Mix_Pt = rt.TH1F("Mix_Pt","Mix_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
Mix_Eta = rt.TH1F("Mix_Eta","Mix_Eta;#eta(j);Events",15, -7, 7)
Mix_Phi = rt.TH1F("Mix_Phi","Mix_Phi;#phi(j);Events",15, -3.2, 3.2)
Mix_DPhi = rt.TH1F("Mix_DPhi","Mix_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Mix_DPhiMin = rt.TH1F("Mix_DPhiMin","Mix_DPhiMin; #Delta#phi_{min}(j,#slash{E}_{T});Events",30,0,3.5)
Mix_M = rt.TH1F("Mix_M","Mix_M; M(j);Events",40,0,1200)
Mix_mult = rt.TH1F("Mix_mult","Mix_mult;Number of Jets;Events",7,0,7)
Mix_axisMajor = rt.TH1F("Mix_axisMajor","Mix_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
Mix_axisMinor = rt.TH1F("Mix_axisMinor","Mix_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
Mix_momentGirth = rt.TH1F("Mix_momentGirth","Mix_momentGirth;girth(j);Events",40, 0, 0.5)
Mix_ptD = rt.TH1F("Mix_ptD","Mix_ptD;p_{T}D(j);Events",40, 0, 1.2)
Mix_tau1 = rt.TH1F("Mix_tau1","Mix_tau1;#tau_{1}(j);Events",40, 0, 0.8)
Mix_tau2 = rt.TH1F("Mix_tau2","Mix_tau2;#tau_{2}(j);Events",40, 0, 0.65)
Mix_tau3 = rt.TH1F("Mix_tau3","Mix_tau3;#tau_{3}(j);Events",40, 0, 0.35)
Mix_tau21 = rt.TH1F("Mix_tau21","Mix_tau21;#tau_{21}(j);Events",40, 0, 1.3)
Mix_tau32 = rt.TH1F("Mix_tau32","Mix_tau32;#tau_{32}(j);Events",40, 0, 1.3)
Mix_SDM = rt.TH1F("Mix_SDM","Mix_SDM; m_{SD}(j);Events",40,0,msd_max)
Mix_SDPt = rt.TH1F("Mix_SDPt","Mix_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
Mix_pmult = rt.TH1F("Mix_pmult","Mix_pmult;Number of Particles in a Jet;Events",40,0,400)

DarkerMix_Pt = rt.TH1F("DarkerMix_Pt","DarkerMix_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
DarkerMix_Eta = rt.TH1F("DarkerMix_Eta","DarkerMix_Eta;#eta(j);Events",15, -7, 7)
DarkerMix_Phi = rt.TH1F("DarkerMix_Phi","DarkerMix_Phi;#phi(j);Events",15, -3.2, 3.2)
DarkerMix_DPhi = rt.TH1F("DarkerMix_DPhi","DarkerMix_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
DarkerMix_M = rt.TH1F("DarkerMix_M","DarkerMix_M; M(j);Events",40,0,1200)
DarkerMix_mult = rt.TH1F("DarkerMix_mult","DarkerMix_mult;Number of Jets;Events",7,0,7)
DarkerMix_axisMajor = rt.TH1F("DarkerMix_axisMajor","DarkerMix_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
DarkerMix_axisMinor = rt.TH1F("DarkerMix_axisMinor","DarkerMix_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
DarkerMix_momentGirth = rt.TH1F("DarkerMix_momentGirth","DarkerMix_momentGirth;girth(j);Events",40, 0, 0.5)
DarkerMix_ptD = rt.TH1F("DarkerMix_ptD","DarkerMix_ptD;p_{T}D(j);Events",40, 0, 1.2)
DarkerMix_tau1 = rt.TH1F("DarkerMix_tau1","DarkerMix_tau1;#tau_{1}(j);Events",40, 0, 0.8)
DarkerMix_tau2 = rt.TH1F("DarkerMix_tau2","DarkerMix_tau2;#tau_{2}(j);Events",40, 0, 0.65)
DarkerMix_tau3 = rt.TH1F("DarkerMix_tau3","DarkerMix_tau3;#tau_{3}(j);Events",40, 0, 0.35)
DarkerMix_tau21 = rt.TH1F("DarkerMix_tau21","DarkerMix_tau21;#tau_{21}(j);Events",40, 0, 1.3)
DarkerMix_tau32 = rt.TH1F("DarkerMix_tau32","DarkerMix_tau32;#tau_{32}(j);Events",40, 0, 1.3)
DarkerMix_SDM = rt.TH1F("DarkerMix_SDM","DarkerMix_SDM; m_{SD}(j);Events",40,0,msd_max)
DarkerMix_SDPt = rt.TH1F("DarkerMix_SDPt","DarkerMix_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
DarkerMix_pmult = rt.TH1F("DarkerMix_pmult","DarkerMix_pmult;Number of Particles in a Jet;Events",40,0,400)

Dark_Pt = rt.TH1F("Dark_Pt","Dark_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
Dark_Eta = rt.TH1F("Dark_Eta","Dark_Eta;#eta(j);Events",15, -7, 7)
Dark_Phi = rt.TH1F("Dark_Phi","Dark_Phi;#phi(j);Events",15, -3.2, 3.2)
Dark_DPhi = rt.TH1F("Dark_DPhi","Dark_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Dark_DPhiMin = rt.TH1F("Dark_DPhiMin","Dark_DPhiMin; #Delta#phi_{min}(j,#slash{E}_{T});Events",30,0,3.5)
Dark_M = rt.TH1F("Dark_M","Dark_M; M(j);Events",40,0,1200)
Dark_mult = rt.TH1F("Dark_mult","Dark_mult;Number of Jets;Events",7,0,7)
Dark_axisMajor = rt.TH1F("Dark_axisMajor","Dark_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
Dark_axisMinor = rt.TH1F("Dark_axisMinor","Dark_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
Dark_momentGirth = rt.TH1F("Dark_momentGirth","Dark_momentGirth;girth(j);Events",40, 0, 0.5)
Dark_ptD = rt.TH1F("Dark_ptD","Dark_ptD;p_{T}D(j);Events",40, 0, 1.2)
Dark_tau1 = rt.TH1F("Dark_tau1","Dark_tau1;#tau_{1}(j);Events",40, 0, 0.8)
Dark_tau2 = rt.TH1F("Dark_tau2","Dark_tau2;#tau_{2}(j);Events",40, 0, 0.65)
Dark_tau3 = rt.TH1F("Dark_tau3","Dark_tau3;#tau_{3}(j);Events",40, 0, 0.35)
Dark_tau21 = rt.TH1F("Dark_tau21","Dark_tau21;#tau_{21}(j);Events",40, 0, 1.3)
Dark_tau32 = rt.TH1F("Dark_tau32","Dark_tau32;#tau_{32}(j);Events",40, 0, 1.3)
Dark_SDM = rt.TH1F("Dark_SDM","Dark_SDM; m_{SD}(j);Events",40,0,msd_max)
Dark_SDPt = rt.TH1F("Dark_SDPt","Dark_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
Dark_pmult = rt.TH1F("Dark_pmult","Dark_pmult;Number of Particles in a Jet;Events",40,0,400)

SoftDark_Pt = rt.TH1F("SoftDark_Pt","SoftDark_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
SoftDark_Eta = rt.TH1F("SoftDark_Eta","SoftDark_Eta;#eta(j);Events",15, -7, 7)
SoftDark_Phi = rt.TH1F("SoftDark_Phi","SoftDark_Phi;#phi(j);Events",15, -3.2, 3.2)
SoftDark_DPhi = rt.TH1F("SoftDark_DPhi","SoftDark_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
SoftDark_M = rt.TH1F("SoftDark_M","SoftDark_M; M(j);Events",40,0,1200)
SoftDark_mult = rt.TH1F("SoftDark_mult","SoftDark_mult;Number of Jets;Events",7,0,7)
SoftDark_axisMajor = rt.TH1F("SoftDark_axisMajor","SoftDark_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
SoftDark_axisMinor = rt.TH1F("SoftDark_axisMinor","SoftDark_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
SoftDark_momentGirth = rt.TH1F("SoftDark_momentGirth","SoftDark_momentGirth;girth(j);Events",40, 0, 0.5)
SoftDark_ptD = rt.TH1F("SoftDark_ptD","SoftDark_ptD;p_{T}D(j);Events",40, 0, 1.2)
SoftDark_tau1 = rt.TH1F("SoftDark_tau1","SoftDark_tau1;#tau_{1}(j);Events",40, 0, 0.8)
SoftDark_tau2 = rt.TH1F("SoftDark_tau2","SoftDark_tau2;#tau_{2}(j);Events",40, 0, 0.65)
SoftDark_tau3 = rt.TH1F("SoftDark_tau3","SoftDark_tau3;#tau_{3}(j);Events",40, 0, 0.35)
SoftDark_tau21 = rt.TH1F("SoftDark_tau21","SoftDark_tau21;#tau_{21}(j);Events",40, 0, 1.3)
SoftDark_tau32 = rt.TH1F("SoftDark_tau32","SoftDark_tau32;#tau_{32}(j);Events",40, 0, 1.3)
SoftDark_SDM = rt.TH1F("SoftDark_SDM","SoftDark_SDM; m_{SD}(j);Events",40,0,msd_max)
SoftDark_SDPt = rt.TH1F("SoftDark_SDPt","SoftDark_SDPt; pT_{SD}(j);Events",40,0,ptsd_max)
SoftDark_pmult = rt.TH1F("SoftDark_pmult","SoftDark_pmult;Number of Particles in a Jet;Events",40,0,400)

jet_SM_OnlydpTFrac = rt.TH1F("SM_OnlydpTFrac","SM_OnlydpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_Mix_dpTFrac = rt.TH1F("Mix_dpTFrac","Mix_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_DarkerMix_dpTFrac = rt.TH1F("DarkerMix_dpTFrac","DarkerMix_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_Dark_dpTFrac = rt.TH1F("Dark_dpTFrac","Dark_dpTFrac;Dark pT Fraction;Events",40, 0, 1)
jet_SoftDark_dpTFrac = rt.TH1F("SoftDark_dpTFrac","SoftDark_dpTFrac;Dark pT Fraction;Events",40, 0, 1)

# pT normalized version
SM_OnlyM_Norm = rt.TH1F("SM_OnlyM_Norm","SM_OnlyM_Norm; M(j);Events",40,0,1200)
SM_OnlyaxisMajor_Norm = rt.TH1F("SM_OnlyaxisMajor_Norm","SM_OnlyaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
SM_OnlyaxisMinor_Norm = rt.TH1F("SM_OnlyaxisMinor_Norm","SM_OnlyaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
SM_OnlymomentGirth_Norm = rt.TH1F("SM_OnlymomentGirth_Norm","SM_OnlymomentGirth_Norm;girth(j);Events",40, 0, 0.5)
SM_OnlyptD_Norm = rt.TH1F("SM_OnlyptD_Norm","SM_OnlyptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
SM_Onlytau1_Norm = rt.TH1F("SM_Onlytau1_Norm","SM_Onlytau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
SM_Onlytau2_Norm = rt.TH1F("SM_Onlytau2_Norm","SM_Onlytau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
SM_Onlytau3_Norm = rt.TH1F("SM_Onlytau3_Norm","SM_Onlytau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
SM_Onlytau21_Norm = rt.TH1F("SM_Onlytau21_Norm","SM_Onlytau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
SM_Onlytau32_Norm = rt.TH1F("SM_Onlytau32_Norm","SM_Onlytau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
SM_OnlySDM_Norm = rt.TH1F("SM_OnlySDM_Norm","SM_OnlySDM_Norm; m_{SD}(j);Events",40,0,msd_max)
SM_Onlypmult_Norm = rt.TH1F("SM_Onlypmult_Norm","SM_Onlypmult_Norm;Number of Particles in a Jet;Events",40,0,400)

Mix_M_Norm = rt.TH1F("Mix_M_Norm","Mix_M_Norm; M(j);Events",40,0,1200)
Mix_axisMajor_Norm = rt.TH1F("Mix_axisMajor_Norm","Mix_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
Mix_axisMinor_Norm = rt.TH1F("Mix_axisMinor_Norm","Mix_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
Mix_momentGirth_Norm = rt.TH1F("Mix_momentGirth_Norm","Mix_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
Mix_ptD_Norm = rt.TH1F("Mix_ptD_Norm","Mix_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
Mix_tau1_Norm = rt.TH1F("Mix_tau1_Norm","Mix_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
Mix_tau2_Norm = rt.TH1F("Mix_tau2_Norm","Mix_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
Mix_tau3_Norm = rt.TH1F("Mix_tau3_Norm","Mix_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
Mix_tau21_Norm = rt.TH1F("Mix_tau21_Norm","Mix_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
Mix_tau32_Norm = rt.TH1F("Mix_tau32_Norm","Mix_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
Mix_SDM_Norm = rt.TH1F("Mix_SDM_Norm","Mix_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
Mix_pmult_Norm = rt.TH1F("Mix_pmult_Norm","Mix_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

DarkerMix_M_Norm = rt.TH1F("DarkerMix_M_Norm","DarkerMix_M_Norm; M(j);Events",40,0,1200)
DarkerMix_axisMajor_Norm = rt.TH1F("DarkerMix_axisMajor_Norm","DarkerMix_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
DarkerMix_axisMinor_Norm = rt.TH1F("DarkerMix_axisMinor_Norm","DarkerMix_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
DarkerMix_momentGirth_Norm = rt.TH1F("DarkerMix_momentGirth_Norm","DarkerMix_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
DarkerMix_ptD_Norm = rt.TH1F("DarkerMix_ptD_Norm","DarkerMix_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
DarkerMix_tau1_Norm = rt.TH1F("DarkerMix_tau1_Norm","DarkerMix_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
DarkerMix_tau2_Norm = rt.TH1F("DarkerMix_tau2_Norm","DarkerMix_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
DarkerMix_tau3_Norm = rt.TH1F("DarkerMix_tau3_Norm","DarkerMix_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
DarkerMix_tau21_Norm = rt.TH1F("DarkerMix_tau21_Norm","DarkerMix_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
DarkerMix_tau32_Norm = rt.TH1F("DarkerMix_tau32_Norm","DarkerMix_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
DarkerMix_SDM_Norm = rt.TH1F("DarkerMix_SDM_Norm","DarkerMix_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
DarkerMix_pmult_Norm = rt.TH1F("DarkerMix_pmult_Norm","DarkerMix_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

Dark_M_Norm = rt.TH1F("Dark_M_Norm","Dark_M_Norm; M(j);Events",40,0,1200)
Dark_axisMajor_Norm = rt.TH1F("Dark_axisMajor_Norm","Dark_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
Dark_axisMinor_Norm = rt.TH1F("Dark_axisMinor_Norm","Dark_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
Dark_momentGirth_Norm = rt.TH1F("Dark_momentGirth_Norm","Dark_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
Dark_ptD_Norm = rt.TH1F("Dark_ptD_Norm","Dark_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
Dark_tau1_Norm = rt.TH1F("Dark_tau1_Norm","Dark_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
Dark_tau2_Norm = rt.TH1F("Dark_tau2_Norm","Dark_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
Dark_tau3_Norm = rt.TH1F("Dark_tau3_Norm","Dark_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
Dark_tau21_Norm = rt.TH1F("Dark_tau21_Norm","Dark_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
Dark_tau32_Norm = rt.TH1F("Dark_tau32_Norm","Dark_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
Dark_SDM_Norm = rt.TH1F("Dark_SDM_Norm","Dark_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
Dark_pmult_Norm = rt.TH1F("Dark_pmult_Norm","Dark_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

SoftDark_M_Norm = rt.TH1F("SoftDark_M_Norm","SoftDark_M_Norm; M(j);Events",40,0,1200)
SoftDark_axisMajor_Norm = rt.TH1F("SoftDark_axisMajor_Norm","SoftDark_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
SoftDark_axisMinor_Norm = rt.TH1F("SoftDark_axisMinor_Norm","SoftDark_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
SoftDark_momentGirth_Norm = rt.TH1F("SoftDark_momentGirth_Norm","SoftDark_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
SoftDark_ptD_Norm = rt.TH1F("SoftDark_ptD_Norm","SoftDark_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
SoftDark_tau1_Norm = rt.TH1F("SoftDark_tau1_Norm","SoftDark_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
SoftDark_tau2_Norm = rt.TH1F("SoftDark_tau2_Norm","SoftDark_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
SoftDark_tau3_Norm = rt.TH1F("SoftDark_tau3_Norm","SoftDark_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
SoftDark_tau21_Norm = rt.TH1F("SoftDark_tau21_Norm","SoftDark_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
SoftDark_tau32_Norm = rt.TH1F("SoftDark_tau32_Norm","SoftDark_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
SoftDark_SDM_Norm = rt.TH1F("SoftDark_SDM_Norm","SoftDark_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
SoftDark_pmult_Norm = rt.TH1F("SoftDark_pmult_Norm","SoftDark_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

jet_SM_OnlydpTFrac_Norm = rt.TH1F("SM_OnlydpTFrac_Norm","SM_OnlydpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_Mix_dpTFrac_Norm = rt.TH1F("Mix_dpTFrac_Norm","Mix_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_DarkerMix_dpTFrac_Norm = rt.TH1F("DarkerMix_dpTFrac_Norm","DarkerMix_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_Dark_dpTFrac_Norm = rt.TH1F("Dark_dpTFrac_Norm","Dark_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_SoftDark_dpTFrac_Norm = rt.TH1F("SoftDark_dpTFrac_Norm","SoftDark_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)

# pT normalized version for the stricter categories
SM_M_Norm = rt.TH1F("SM_M_Norm","SM_M_Norm; M(j);Events",40,0,1200)
SM_axisMajor_Norm = rt.TH1F("SM_axisMajor_Norm","SM_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
SM_axisMinor_Norm = rt.TH1F("SM_axisMinor_Norm","SM_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
SM_momentGirth_Norm = rt.TH1F("SM_momentGirth_Norm","SM_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
SM_ptD_Norm = rt.TH1F("SM_ptD_Norm","SM_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
SM_tau1_Norm = rt.TH1F("SM_tau1_Norm","SM_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
SM_tau2_Norm = rt.TH1F("SM_tau2_Norm","SM_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
SM_tau3_Norm = rt.TH1F("SM_tau3_Norm","SM_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
SM_tau21_Norm = rt.TH1F("SM_tau21_Norm","SM_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
SM_tau32_Norm = rt.TH1F("SM_tau32_Norm","SM_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
SM_SDM_Norm = rt.TH1F("SM_SDM_Norm","SM_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
SM_pmult_Norm = rt.TH1F("SM_pmult_Norm","SM_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

SMM_M_Norm = rt.TH1F("SMM_M_Norm","SMM_M_Norm; M(j);Events",40,0,1200)
SMM_axisMajor_Norm = rt.TH1F("SMM_axisMajor_Norm","SMM_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
SMM_axisMinor_Norm = rt.TH1F("SMM_axisMinor_Norm","SMM_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
SMM_momentGirth_Norm = rt.TH1F("SMM_momentGirth_Norm","SMM_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
SMM_ptD_Norm = rt.TH1F("SMM_ptD_Norm","SMM_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
SMM_tau1_Norm = rt.TH1F("SMM_tau1_Norm","SMM_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
SMM_tau2_Norm = rt.TH1F("SMM_tau2_Norm","SMM_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
SMM_tau3_Norm = rt.TH1F("SMM_tau3_Norm","SMM_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
SMM_tau21_Norm = rt.TH1F("SMM_tau21_Norm","SMM_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
SMM_tau32_Norm = rt.TH1F("SMM_tau32_Norm","SMM_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
SMM_SDM_Norm = rt.TH1F("SMM_SDM_Norm","SMM_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
SMM_pmult_Norm = rt.TH1F("SMM_pmult_Norm","SMM_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

G_M_Norm = rt.TH1F("G_M_Norm","G_M_Norm; M(j);Events",40,0,1200)
G_axisMajor_Norm = rt.TH1F("G_axisMajor_Norm","G_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
G_axisMinor_Norm = rt.TH1F("G_axisMinor_Norm","G_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
G_momentGirth_Norm = rt.TH1F("G_momentGirth_Norm","G_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
G_ptD_Norm = rt.TH1F("G_ptD_Norm","G_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
G_tau1_Norm = rt.TH1F("G_tau1_Norm","G_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
G_tau2_Norm = rt.TH1F("G_tau2_Norm","G_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
G_tau3_Norm = rt.TH1F("G_tau3_Norm","G_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
G_tau21_Norm = rt.TH1F("G_tau21_Norm","G_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
G_tau32_Norm = rt.TH1F("G_tau32_Norm","G_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
G_SDM_Norm = rt.TH1F("G_SDM_Norm","G_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
G_pmult_Norm = rt.TH1F("G_pmult_Norm","G_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

QM_M_Norm = rt.TH1F("QM_M_Norm","QM_M_Norm; M(j);Events",40,0,1200)
QM_axisMajor_Norm = rt.TH1F("QM_axisMajor_Norm","QM_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
QM_axisMinor_Norm = rt.TH1F("QM_axisMinor_Norm","QM_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_momentGirth_Norm = rt.TH1F("QM_momentGirth_Norm","QM_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
QM_ptD_Norm = rt.TH1F("QM_ptD_Norm","QM_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
QM_tau1_Norm = rt.TH1F("QM_tau1_Norm","QM_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
QM_tau2_Norm = rt.TH1F("QM_tau2_Norm","QM_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
QM_tau3_Norm = rt.TH1F("QM_tau3_Norm","QM_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
QM_tau21_Norm = rt.TH1F("QM_tau21_Norm","QM_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
QM_tau32_Norm = rt.TH1F("QM_tau32_Norm","QM_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
QM_SDM_Norm = rt.TH1F("QM_SDM_Norm","QM_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
QM_pmult_Norm = rt.TH1F("QM_pmult_Norm","QM_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

Q_M_Norm = rt.TH1F("Q_M_Norm","Q_M_Norm; M(j);Events",40,0,1200)
Q_axisMajor_Norm = rt.TH1F("Q_axisMajor_Norm","Q_axisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
Q_axisMinor_Norm = rt.TH1F("Q_axisMinor_Norm","Q_axisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
Q_momentGirth_Norm = rt.TH1F("Q_momentGirth_Norm","Q_momentGirth_Norm;girth(j);Events",40, 0, 0.5)
Q_ptD_Norm = rt.TH1F("Q_ptD_Norm","Q_ptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
Q_tau1_Norm = rt.TH1F("Q_tau1_Norm","Q_tau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
Q_tau2_Norm = rt.TH1F("Q_tau2_Norm","Q_tau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
Q_tau3_Norm = rt.TH1F("Q_tau3_Norm","Q_tau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
Q_tau21_Norm = rt.TH1F("Q_tau21_Norm","Q_tau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
Q_tau32_Norm = rt.TH1F("Q_tau32_Norm","Q_tau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
Q_SDM_Norm = rt.TH1F("Q_SDM_Norm","Q_SDM_Norm; m_{SD}(j);Events",40,0,msd_max)
Q_pmult_Norm = rt.TH1F("Q_pmult_Norm","Q_pmult_Norm;Number of Particles in a Jet;Events",40,0,400)

QM_QM_Norm = rt.TH1F("QM_QM_Norm","QM_QM_Norm; M(j);Events",40,0,1200)
QM_QaxisMajor_Norm = rt.TH1F("QM_QaxisMajor_Norm","QM_QaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
QM_QaxisMinor_Norm = rt.TH1F("QM_QaxisMinor_Norm","QM_QaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_QmomentGirth_Norm = rt.TH1F("QM_QmomentGirth_Norm","QM_QmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
QM_QptD_Norm = rt.TH1F("QM_QptD_Norm","QM_QptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
QM_Qtau1_Norm = rt.TH1F("QM_Qtau1_Norm","QM_Qtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
QM_Qtau2_Norm = rt.TH1F("QM_Qtau2_Norm","QM_Qtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
QM_Qtau3_Norm = rt.TH1F("QM_Qtau3_Norm","QM_Qtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
QM_Qtau21_Norm = rt.TH1F("QM_Qtau21_Norm","QM_Qtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
QM_Qtau32_Norm = rt.TH1F("QM_Qtau32_Norm","QM_Qtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
QM_QSDM_Norm = rt.TH1F("QM_QSDM_Norm","QM_QSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
QM_Qpmult_Norm = rt.TH1F("QM_Qpmult_Norm","QM_Qpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

QM_GM_Norm = rt.TH1F("QM_GM_Norm","QM_GM_Norm; M(j);Events",40,0,1200)
QM_GaxisMajor_Norm = rt.TH1F("QM_GaxisMajor_Norm","QM_GaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
QM_GaxisMinor_Norm = rt.TH1F("QM_GaxisMinor_Norm","QM_GaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_GmomentGirth_Norm = rt.TH1F("QM_GmomentGirth_Norm","QM_GmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
QM_GptD_Norm = rt.TH1F("QM_GptD_Norm","QM_GptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
QM_Gtau1_Norm = rt.TH1F("QM_Gtau1_Norm","QM_Gtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
QM_Gtau2_Norm = rt.TH1F("QM_Gtau2_Norm","QM_Gtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
QM_Gtau3_Norm = rt.TH1F("QM_Gtau3_Norm","QM_Gtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
QM_Gtau21_Norm = rt.TH1F("QM_Gtau21_Norm","QM_Gtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
QM_Gtau32_Norm = rt.TH1F("QM_Gtau32_Norm","QM_Gtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
QM_GSDM_Norm = rt.TH1F("QM_GSDM_Norm","QM_GSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
QM_Gpmult_Norm = rt.TH1F("QM_Gpmult_Norm","QM_Gpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

Q_GM_Norm = rt.TH1F("Q_GM_Norm","Q_GM_Norm; M(j);Events",40,0,1200)
Q_GaxisMajor_Norm = rt.TH1F("Q_GaxisMajor_Norm","Q_GaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
Q_GaxisMinor_Norm = rt.TH1F("Q_GaxisMinor_Norm","Q_GaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
Q_GmomentGirth_Norm = rt.TH1F("Q_GmomentGirth_Norm","Q_GmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
Q_GptD_Norm = rt.TH1F("Q_GptD_Norm","Q_GptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
Q_Gtau1_Norm = rt.TH1F("Q_Gtau1_Norm","Q_Gtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
Q_Gtau2_Norm = rt.TH1F("Q_Gtau2_Norm","Q_Gtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
Q_Gtau3_Norm = rt.TH1F("Q_Gtau3_Norm","Q_Gtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
Q_Gtau21_Norm = rt.TH1F("Q_Gtau21_Norm","Q_Gtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
Q_Gtau32_Norm = rt.TH1F("Q_Gtau32_Norm","Q_Gtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
Q_GSDM_Norm = rt.TH1F("Q_GSDM_Norm","Q_GSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
Q_Gpmult_Norm = rt.TH1F("Q_Gpmult_Norm","Q_Gpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

G_SMM_Norm = rt.TH1F("G_SMM_Norm","G_SMM_Norm; M(j);Events",40,0,1200)
G_SMaxisMajor_Norm = rt.TH1F("G_SMaxisMajor_Norm","G_SMaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
G_SMaxisMinor_Norm = rt.TH1F("G_SMaxisMinor_Norm","G_SMaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
G_SMmomentGirth_Norm = rt.TH1F("G_SMmomentGirth_Norm","G_SMmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
G_SMptD_Norm = rt.TH1F("G_SMptD_Norm","G_SMptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
G_SMtau1_Norm = rt.TH1F("G_SMtau1_Norm","G_SMtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
G_SMtau2_Norm = rt.TH1F("G_SMtau2_Norm","G_SMtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
G_SMtau3_Norm = rt.TH1F("G_SMtau3_Norm","G_SMtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
G_SMtau21_Norm = rt.TH1F("G_SMtau21_Norm","G_SMtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
G_SMtau32_Norm = rt.TH1F("G_SMtau32_Norm","G_SMtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
G_SMSDM_Norm = rt.TH1F("G_SMSDM_Norm","G_SMSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
G_SMpmult_Norm = rt.TH1F("G_SMpmult_Norm","G_SMpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

QM_SMM_Norm = rt.TH1F("QM_SMM_Norm","QM_SMM_Norm; M(j);Events",40,0,1200)
QM_SMaxisMajor_Norm = rt.TH1F("QM_SMaxisMajor_Norm","QM_SMaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
QM_SMaxisMinor_Norm = rt.TH1F("QM_SMaxisMinor_Norm","QM_SMaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
QM_SMmomentGirth_Norm = rt.TH1F("QM_SMmomentGirth_Norm","QM_SMmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
QM_SMptD_Norm = rt.TH1F("QM_SMptD_Norm","QM_SMptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
QM_SMtau1_Norm = rt.TH1F("QM_SMtau1_Norm","QM_SMtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
QM_SMtau2_Norm = rt.TH1F("QM_SMtau2_Norm","QM_SMtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
QM_SMtau3_Norm = rt.TH1F("QM_SMtau3_Norm","QM_SMtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
QM_SMtau21_Norm = rt.TH1F("QM_SMtau21_Norm","QM_SMtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
QM_SMtau32_Norm = rt.TH1F("QM_SMtau32_Norm","QM_SMtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
QM_SMSDM_Norm = rt.TH1F("QM_SMSDM_Norm","QM_SMSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
QM_SMpmult_Norm = rt.TH1F("QM_SMpmult_Norm","QM_SMpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

Q_SMM_Norm = rt.TH1F("Q_SMM_Norm","Q_SMM_Norm; M(j);Events",40,0,1200)
Q_SMaxisMajor_Norm = rt.TH1F("Q_SMaxisMajor_Norm","Q_SMaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
Q_SMaxisMinor_Norm = rt.TH1F("Q_SMaxisMinor_Norm","Q_SMaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
Q_SMmomentGirth_Norm = rt.TH1F("Q_SMmomentGirth_Norm","Q_SMmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
Q_SMptD_Norm = rt.TH1F("Q_SMptD_Norm","Q_SMptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
Q_SMtau1_Norm = rt.TH1F("Q_SMtau1_Norm","Q_SMtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
Q_SMtau2_Norm = rt.TH1F("Q_SMtau2_Norm","Q_SMtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
Q_SMtau3_Norm = rt.TH1F("Q_SMtau3_Norm","Q_SMtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
Q_SMtau21_Norm = rt.TH1F("Q_SMtau21_Norm","Q_SMtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
Q_SMtau32_Norm = rt.TH1F("Q_SMtau32_Norm","Q_SMtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
Q_SMSDM_Norm = rt.TH1F("Q_SMSDM_Norm","Q_SMSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
Q_SMpmult_Norm = rt.TH1F("Q_SMpmult_Norm","Q_SMpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

LD_highDFM_Norm = rt.TH1F("LD_highDFM_Norm","LD_highDFM_Norm; M(j);Events",40,0,1200)
LD_highDFaxisMajor_Norm = rt.TH1F("LD_highDFaxisMajor_Norm","LD_highDFaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
LD_highDFaxisMinor_Norm = rt.TH1F("LD_highDFaxisMinor_Norm","LD_highDFaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
LD_highDFmomentGirth_Norm = rt.TH1F("LD_highDFmomentGirth_Norm","LD_highDFmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
LD_highDFptD_Norm = rt.TH1F("LD_highDFptD_Norm","LD_highDFptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
LD_highDFtau1_Norm = rt.TH1F("LD_highDFtau1_Norm","LD_highDFtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
LD_highDFtau2_Norm = rt.TH1F("LD_highDFtau2_Norm","LD_highDFtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
LD_highDFtau3_Norm = rt.TH1F("LD_highDFtau3_Norm","LD_highDFtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
LD_highDFtau21_Norm = rt.TH1F("LD_highDFtau21_Norm","LD_highDFtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
LD_highDFtau32_Norm = rt.TH1F("LD_highDFtau32_Norm","LD_highDFtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
LD_highDFSDM_Norm = rt.TH1F("LD_highDFSDM_Norm","LD_highDFSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
LD_highDFpmult_Norm = rt.TH1F("LD_highDFpmult_Norm","LD_highDFpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

LD_lowDFM_Norm = rt.TH1F("LD_lowDFM_Norm","LD_lowDFM_Norm; M(j);Events",40,0,1200)
LD_lowDFaxisMajor_Norm = rt.TH1F("LD_lowDFaxisMajor_Norm","LD_lowDFaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
LD_lowDFaxisMinor_Norm = rt.TH1F("LD_lowDFaxisMinor_Norm","LD_lowDFaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
LD_lowDFmomentGirth_Norm = rt.TH1F("LD_lowDFmomentGirth_Norm","LD_lowDFmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
LD_lowDFptD_Norm = rt.TH1F("LD_lowDFptD_Norm","LD_lowDFptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
LD_lowDFtau1_Norm = rt.TH1F("LD_lowDFtau1_Norm","LD_lowDFtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
LD_lowDFtau2_Norm = rt.TH1F("LD_lowDFtau2_Norm","LD_lowDFtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
LD_lowDFtau3_Norm = rt.TH1F("LD_lowDFtau3_Norm","LD_lowDFtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
LD_lowDFtau21_Norm = rt.TH1F("LD_lowDFtau21_Norm","LD_lowDFtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
LD_lowDFtau32_Norm = rt.TH1F("LD_lowDFtau32_Norm","LD_lowDFtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
LD_lowDFSDM_Norm = rt.TH1F("LD_lowDFSDM_Norm","LD_lowDFSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
LD_lowDFpmult_Norm = rt.TH1F("LD_lowDFpmult_Norm","LD_lowDFpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

LD_SMM_Norm = rt.TH1F("LD_SMM_Norm","LD_SMM_Norm; M(j);Events",40,0,1200)
LD_SMaxisMajor_Norm = rt.TH1F("LD_SMaxisMajor_Norm","LD_SMaxisMajor_Norm;#sigma_{major}(j);Events",40, 0, 0.5)
LD_SMaxisMinor_Norm = rt.TH1F("LD_SMaxisMinor_Norm","LD_SMaxisMinor_Norm;#sigma_{minor}(j);Events",40, 0, 0.3)
LD_SMmomentGirth_Norm = rt.TH1F("LD_SMmomentGirth_Norm","LD_SMmomentGirth_Norm;girth(j);Events",40, 0, 0.5)
LD_SMptD_Norm = rt.TH1F("LD_SMptD_Norm","LD_SMptD_Norm;p_{T}D(j);Events",40, 0, 1.2)
LD_SMtau1_Norm = rt.TH1F("LD_SMtau1_Norm","LD_SMtau1_Norm;#tau_{1}(j);Events",40, 0, 0.8)
LD_SMtau2_Norm = rt.TH1F("LD_SMtau2_Norm","LD_SMtau2_Norm;#tau_{2}(j);Events",40, 0, 0.65)
LD_SMtau3_Norm = rt.TH1F("LD_SMtau3_Norm","LD_SMtau3_Norm;#tau_{3}(j);Events",40, 0, 0.35)
LD_SMtau21_Norm = rt.TH1F("LD_SMtau21_Norm","LD_SMtau21_Norm;#tau_{21}(j);Events",40, 0, 1.3)
LD_SMtau32_Norm = rt.TH1F("LD_SMtau32_Norm","LD_SMtau32_Norm;#tau_{32}(j);Events",40, 0, 1.3)
LD_SMSDM_Norm = rt.TH1F("LD_SMSDM_Norm","LD_SMSDM_Norm; m_{SD}(j);Events",40,0,msd_max)
LD_SMpmult_Norm = rt.TH1F("LD_SMpmult_Norm","LD_SMpmult_Norm;Number of Particles in a Jet;Events",40,0,400)

jet_SM_dpTFrac_Norm = rt.TH1F("SM_dpTFrac_Norm","SM_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_SMM_dpTFrac_Norm = rt.TH1F("SMM_dpTFrac_Norm","SMM_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_G_dpTFrac_Norm = rt.TH1F("G_dpTFrac_Norm","G_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_QM_dpTFrac_Norm = rt.TH1F("QM_dpTFrac_Norm","QM_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_Q_dpTFrac_Norm = rt.TH1F("Q_dpTFrac_Norm","Q_dpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_QM_QdpTFrac_Norm = rt.TH1F("QM_QdpTFrac_Norm","QM_QdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_QM_GdpTFrac_Norm = rt.TH1F("QM_GdpTFrac_Norm","QM_GdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_Q_GdpTFrac_Norm = rt.TH1F("Q_GdpTFrac_Norm","Q_GdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_G_SMdpTFrac_Norm = rt.TH1F("G_SMdpTFrac_Norm","G_SMdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_QM_SMdpTFrac_Norm = rt.TH1F("QM_SMdpTFrac_Norm","QM_SMdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_Q_SMdpTFrac_Norm = rt.TH1F("Q_SMdpTFrac_Norm","Q_SMdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_LD_highDFdpTFrac_Norm = rt.TH1F("LD_highDFdpTFrac_Norm","LD_highDFdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_LD_lowDFdpTFrac_Norm = rt.TH1F("LD_lowDFdpTFrac_Norm","LD_lowDFdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)
jet_LD_SMdpTFrac_Norm = rt.TH1F("LD_SMdpTFrac_Norm","LD_SMdpTFrac_Norm;Dark pT Fraction;Events",40, 0, 1)

# jet closest to MET (divided into categories)
SM_OnlyPt_Cl = rt.TH1F("SM_OnlyPt_Cl","SM_OnlyPt_Cl;p_{T}(j) [GeV];Events",40, 0, 3500)
SM_OnlyEta_Cl = rt.TH1F("SM_OnlyEta_Cl","SM_OnlyEta_Cl;#eta(j);Events",15, -7, 7)
SM_OnlyPhi_Cl = rt.TH1F("SM_OnlyPhi_Cl","SM_OnlyPhi_Cl;#phi(j);Events",15, -3.2, 3.2)
SM_OnlyDPhi_Cl = rt.TH1F("SM_OnlyDPhi_Cl","SM_OnlyDPhi_Cl; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
SM_OnlyM_Cl = rt.TH1F("SM_OnlyM_Cl","SM_OnlyM_Cl; M(j);Events",40,0,1200)
SM_Onlymult_Cl = rt.TH1F("SM_Onlymult_Cl","SM_Onlymult_Cl;Number of Jets;Events",7,0,7)
SM_OnlyaxisMajor_Cl = rt.TH1F("SM_OnlyaxisMajor_Cl","SM_OnlyaxisMajor_Cl;#sigma_{major}(j);Events",40, 0, 0.5)
SM_OnlyaxisMinor_Cl = rt.TH1F("SM_OnlyaxisMinor_Cl","SM_OnlyaxisMinor_Cl;#sigma_{minor}(j);Events",40, 0, 0.3)
SM_OnlymomentGirth_Cl = rt.TH1F("SM_OnlymomentGirth_Cl","SM_OnlymomentGirth_Cl;girth(j);Events",40, 0, 0.5)
SM_OnlyptD_Cl = rt.TH1F("SM_OnlyptD_Cl","SM_OnlyptD_Cl;p_{T}D(j);Events",40, 0, 1.2)
SM_Onlytau1_Cl = rt.TH1F("SM_Onlytau1_Cl","SM_Onlytau1_Cl;#tau_{1}(j);Events",40, 0, 0.8)
SM_Onlytau2_Cl = rt.TH1F("SM_Onlytau2_Cl","SM_Onlytau2_Cl;#tau_{2}(j);Events",40, 0, 0.65)
SM_Onlytau3_Cl = rt.TH1F("SM_Onlytau3_Cl","SM_Onlytau3_Cl;#tau_{3}(j);Events",40, 0, 0.35)
SM_Onlytau21_Cl = rt.TH1F("SM_Onlytau21_Cl","SM_Onlytau21_Cl;#tau_{21}(j);Events",40, 0, 1.3)
SM_Onlytau32_Cl = rt.TH1F("SM_Onlytau32_Cl","SM_Onlytau32_Cl;#tau_{32}(j);Events",40, 0, 1.3)
SM_OnlySDM_Cl = rt.TH1F("SM_OnlySDM_Cl","SM_OnlySDM_Cl; m_{SD}(j);Events",40,0,msd_max)
SM_Onlypmult_Cl = rt.TH1F("SM_Onlypmult_Cl","SM_Onlypmult_Cl;Number of Particles in a Jet;Events",40,0,400)

Mix_Pt_Cl = rt.TH1F("Mix_Pt_Cl","Mix_Pt_Cl;p_{T}(j) [GeV];Events",40, 0, 3500)
Mix_Eta_Cl = rt.TH1F("Mix_Eta_Cl","Mix_Eta_Cl;#eta(j);Events",15, -7, 7)
Mix_Phi_Cl = rt.TH1F("Mix_Phi_Cl","Mix_Phi_Cl;#phi(j);Events",15, -3.2, 3.2)
Mix_DPhi_Cl = rt.TH1F("Mix_DPhi_Cl","Mix_DPhi_Cl; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Mix_M_Cl = rt.TH1F("Mix_M_Cl","Mix_M_Cl; M(j);Events",40,0,1200)
Mix_mult_Cl = rt.TH1F("Mix_mult_Cl","Mix_mult_Cl;Number of Jets;Events",7,0,7)
Mix_axisMajor_Cl = rt.TH1F("Mix_axisMajor_Cl","Mix_axisMajor_Cl;#sigma_{major}(j);Events",40, 0, 0.5)
Mix_axisMinor_Cl = rt.TH1F("Mix_axisMinor_Cl","Mix_axisMinor_Cl;#sigma_{minor}(j);Events",40, 0, 0.3)
Mix_momentGirth_Cl = rt.TH1F("Mix_momentGirth_Cl","Mix_momentGirth_Cl;girth(j);Events",40, 0, 0.5)
Mix_ptD_Cl = rt.TH1F("Mix_ptD_Cl","Mix_ptD_Cl;p_{T}D(j);Events",40, 0, 1.2)
Mix_tau1_Cl = rt.TH1F("Mix_tau1_Cl","Mix_tau1_Cl;#tau_{1}(j);Events",40, 0, 0.8)
Mix_tau2_Cl = rt.TH1F("Mix_tau2_Cl","Mix_tau2_Cl;#tau_{2}(j);Events",40, 0, 0.65)
Mix_tau3_Cl = rt.TH1F("Mix_tau3_Cl","Mix_tau3_Cl;#tau_{3}(j);Events",40, 0, 0.35)
Mix_tau21_Cl = rt.TH1F("Mix_tau21_Cl","Mix_tau21_Cl;#tau_{21}(j);Events",40, 0, 1.3)
Mix_tau32_Cl = rt.TH1F("Mix_tau32_Cl","Mix_tau32_Cl;#tau_{32}(j);Events",40, 0, 1.3)
Mix_SDM_Cl = rt.TH1F("Mix_SDM_Cl","Mix_SDM_Cl; m_{SD}(j);Events",40,0,msd_max)
Mix_pmult_Cl = rt.TH1F("Mix_pmult_Cl","Mix_pmult_Cl;Number of Particles in a Jet;Events",40,0,400)

Dark_Pt_Cl = rt.TH1F("Dark_Pt_Cl","Dark_Pt_Cl;p_{T}(j) [GeV];Events",40, 0, 3500)
Dark_Eta_Cl = rt.TH1F("Dark_Eta_Cl","Dark_Eta_Cl;#eta(j);Events",15, -7, 7)
Dark_Phi_Cl = rt.TH1F("Dark_Phi_Cl","Dark_Phi_Cl;#phi(j);Events",15, -3.2, 3.2)
Dark_DPhi_Cl = rt.TH1F("Dark_DPhi_Cl","Dark_DPhi_Cl; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
Dark_M_Cl = rt.TH1F("Dark_M_Cl","Dark_M_Cl; M(j);Events",40,0,1200)
Dark_mult_Cl = rt.TH1F("Dark_mult_Cl","Dark_mult_Cl;Number of Jets;Events",7,0,7)
Dark_axisMajor_Cl = rt.TH1F("Dark_axisMajor_Cl","Dark_axisMajor_Cl;#sigma_{major}(j);Events",40, 0, 0.5)
Dark_axisMinor_Cl = rt.TH1F("Dark_axisMinor_Cl","Dark_axisMinor_Cl;#sigma_{minor}(j);Events",40, 0, 0.3)
Dark_momentGirth_Cl = rt.TH1F("Dark_momentGirth_Cl","Dark_momentGirth_Cl;girth(j);Events",40, 0, 0.5)
Dark_ptD_Cl = rt.TH1F("Dark_ptD_Cl","Dark_ptD_Cl;p_{T}D(j);Events",40, 0, 1.2)
Dark_tau1_Cl = rt.TH1F("Dark_tau1_Cl","Dark_tau1_Cl;#tau_{1}(j);Events",40, 0, 0.8)
Dark_tau2_Cl = rt.TH1F("Dark_tau2_Cl","Dark_tau2_Cl;#tau_{2}(j);Events",40, 0, 0.65)
Dark_tau3_Cl = rt.TH1F("Dark_tau3_Cl","Dark_tau3_Cl;#tau_{3}(j);Events",40, 0, 0.35)
Dark_tau21_Cl = rt.TH1F("Dark_tau21_Cl","Dark_tau21_Cl;#tau_{21}(j);Events",40, 0, 1.3)
Dark_tau32_Cl = rt.TH1F("Dark_tau32_Cl","Dark_tau32_Cl;#tau_{32}(j);Events",40, 0, 1.3)
Dark_SDM_Cl = rt.TH1F("Dark_SDM_Cl","Dark_SDM_Cl; m_{SD}(j);Events",40,0,msd_max)
Dark_pmult_Cl = rt.TH1F("Dark_pmult_Cl","Dark_pmult_Cl;Number of Particles in a Jet;Events",40,0,400)

# jet closest to MET (divided into jet number)
All_DPhiMin = rt.TH1F("All_DPhiMin","All_DPhiMin; #Delta#phi_{min}(j,#slash{E}_{T});Events",30,0,3.5)
All_Pt_Cl = rt.TH1F("All_Pt_Cl","All_Pt_Cl;p_{T}(j) [GeV];Events",40, 0, 3500)
All_Eta_Cl = rt.TH1F("All_Eta_Cl","All_Eta_Cl;#eta(j);Events",15, -7, 7)
All_Phi_Cl = rt.TH1F("All_Phi_Cl","All_Phi_Cl;#phi(j);Events",15, -3.2, 3.2)
All_DPhi_Cl = rt.TH1F("All_DPhi_Cl","All_DPhi_Cl; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
All_M_Cl = rt.TH1F("All_M_Cl","All_M_Cl; M(j);Events",40,0,1200)
All_mult_Cl = rt.TH1F("All_mult_Cl","All_mult_Cl;Number of Jets;Events",7,0,7)
All_axisMajor_Cl = rt.TH1F("All_axisMajor_Cl","All_axisMajor_Cl;#sigma_{major}(j);Events",40, 0, 0.5)
All_axisMinor_Cl = rt.TH1F("All_axisMinor_Cl","All_axisMinor_Cl;#sigma_{minor}(j);Events",40, 0, 0.3)
All_momentGirth_Cl = rt.TH1F("All_momentGirth_Cl","All_momentGirth_Cl;girth(j);Events",40, 0, 0.5)
All_ptD_Cl = rt.TH1F("All_ptD_Cl","All_ptD_Cl;p_{T}D(j);Events",40, 0, 1.2)
All_tau1_Cl = rt.TH1F("All_tau1_Cl","All_tau1_Cl;#tau_{1}(j);Events",40, 0, 0.8)
All_tau2_Cl = rt.TH1F("All_tau2_Cl","All_tau2_Cl;#tau_{2}(j);Events",40, 0, 0.65)
All_tau3_Cl = rt.TH1F("All_tau3_Cl","All_tau3_Cl;#tau_{3}(j);Events",40, 0, 0.35)
All_tau21_Cl = rt.TH1F("All_tau21_Cl","All_tau21_Cl;#tau_{21}(j);Events",40, 0, 1.3)
All_tau32_Cl = rt.TH1F("All_tau32_Cl","All_tau32_Cl;#tau_{32}(j);Events",40, 0, 1.3)
All_SDM_Cl = rt.TH1F("All_SDM_Cl","All_SDM_Cl; m_{SD}(j);Events",40,0,msd_max)
All_pmult_Cl = rt.TH1F("All_pmult_Cl","All_pmult_Cl;Number of Particles in a Jet;Events",40,0,400)

All_Pt_j1 = rt.TH1F("All_Pt_j1","All_Pt_j1;p_{T}(j) [GeV];Events",40, 0, 3500)
All_Eta_j1 = rt.TH1F("All_Eta_j1","All_Eta_j1;#eta(j);Events",15, -7, 7)
All_Phi_j1 = rt.TH1F("All_Phi_j1","All_Phi_j1;#phi(j);Events",15, -3.2, 3.2)
All_DPhi_j1 = rt.TH1F("All_DPhi_j1","All_DPhi_j1; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
All_M_j1 = rt.TH1F("All_M_j1","All_M_j1; M(j);Events",40,0,1200)
All_mult_j1 = rt.TH1F("All_mult_j1","All_mult_j1;Number of Jets;Events",7,0,7)
All_axisMajor_j1 = rt.TH1F("All_axisMajor_j1","All_axisMajor_j1;#sigma_{major}(j);Events",40, 0, 0.5)
All_axisMinor_j1 = rt.TH1F("All_axisMinor_j1","All_axisMinor_j1;#sigma_{minor}(j);Events",40, 0, 0.3)
All_momentGirth_j1 = rt.TH1F("All_momentGirth_j1","All_momentGirth_j1;girth(j);Events",40, 0, 0.5)
All_ptD_j1 = rt.TH1F("All_ptD_j1","All_ptD_j1;p_{T}D(j);Events",40, 0, 1.2)
All_tau1_j1 = rt.TH1F("All_tau1_j1","All_tau1_j1;#tau_{1}(j);Events",40, 0, 0.8)
All_tau2_j1 = rt.TH1F("All_tau2_j1","All_tau2_j1;#tau_{2}(j);Events",40, 0, 0.65)
All_tau3_j1 = rt.TH1F("All_tau3_j1","All_tau3_j1;#tau_{3}(j);Events",40, 0, 0.35)
All_tau21_j1 = rt.TH1F("All_tau21_j1","All_tau21_j1;#tau_{21}(j);Events",40, 0, 1.3)
All_tau32_j1 = rt.TH1F("All_tau32_j1","All_tau32_j1;#tau_{32}(j);Events",40, 0, 1.3)
All_SDM_j1 = rt.TH1F("All_SDM_j1","All_SDM_j1; m_{SD}(j);Events",40,0,msd_max)
All_pmult_j1 = rt.TH1F("All_pmult_j1","All_pmult_j1;Number of Particles in a Jet;Events",40,0,400)

All_Pt_j2 = rt.TH1F("All_Pt_j2","All_Pt_j2;p_{T}(j) [GeV];Events",40, 0, 3500)
All_Eta_j2 = rt.TH1F("All_Eta_j2","All_Eta_j2;#eta(j);Events",15, -7, 7)
All_Phi_j2 = rt.TH1F("All_Phi_j2","All_Phi_j2;#phi(j);Events",15, -3.2, 3.2)
All_DPhi_j2 = rt.TH1F("All_DPhi_j2","All_DPhi_j2; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
All_M_j2 = rt.TH1F("All_M_j2","All_M_j2; M(j);Events",40,0,1200)
All_mult_j2 = rt.TH1F("All_mult_j2","All_mult_j2;Number of Jets;Events",7,0,7)
All_axisMajor_j2 = rt.TH1F("All_axisMajor_j2","All_axisMajor_j2;#sigma_{major}(j);Events",40, 0, 0.5)
All_axisMinor_j2 = rt.TH1F("All_axisMinor_j2","All_axisMinor_j2;#sigma_{minor}(j);Events",40, 0, 0.3)
All_momentGirth_j2 = rt.TH1F("All_momentGirth_j2","All_momentGirth_j2;girth(j);Events",40, 0, 0.5)
All_ptD_j2 = rt.TH1F("All_ptD_j2","All_ptD_j2;p_{T}D(j);Events",40, 0, 1.2)
All_tau1_j2 = rt.TH1F("All_tau1_j2","All_tau1_j2;#tau_{1}(j);Events",40, 0, 0.8)
All_tau2_j2 = rt.TH1F("All_tau2_j2","All_tau2_j2;#tau_{2}(j);Events",40, 0, 0.65)
All_tau3_j2 = rt.TH1F("All_tau3_j2","All_tau3_j2;#tau_{3}(j);Events",40, 0, 0.35)
All_tau21_j2 = rt.TH1F("All_tau21_j2","All_tau21_j2;#tau_{21}(j);Events",40, 0, 1.3)
All_tau32_j2 = rt.TH1F("All_tau32_j2","All_tau32_j2;#tau_{32}(j);Events",40, 0, 1.3)
All_SDM_j2 = rt.TH1F("All_SDM_j2","All_SDM_j2; m_{SD}(j);Events",40,0,msd_max)
All_pmult_j2 = rt.TH1F("All_pmult_j2","All_pmult_j2;Number of Particles in a Jet;Events",40,0,400)

All_Pt_j3 = rt.TH1F("All_Pt_j3","All_Pt_j3;p_{T}(j) [GeV];Events",40, 0, 3500)
All_Eta_j3 = rt.TH1F("All_Eta_j3","All_Eta_j3;#eta(j);Events",15, -7, 7)
All_Phi_j3 = rt.TH1F("All_Phi_j3","All_Phi_j3;#phi(j);Events",15, -3.2, 3.2)
All_DPhi_j3 = rt.TH1F("All_DPhi_j3","All_DPhi_j3; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
All_M_j3 = rt.TH1F("All_M_j3","All_M_j3; M(j);Events",40,0,1200)
All_mult_j3 = rt.TH1F("All_mult_j3","All_mult_j3;Number of Jets;Events",7,0,7)
All_axisMajor_j3 = rt.TH1F("All_axisMajor_j3","All_axisMajor_j3;#sigma_{major}(j);Events",40, 0, 0.5)
All_axisMinor_j3 = rt.TH1F("All_axisMinor_j3","All_axisMinor_j3;#sigma_{minor}(j);Events",40, 0, 0.3)
All_momentGirth_j3 = rt.TH1F("All_momentGirth_j3","All_momentGirth_j3;girth(j);Events",40, 0, 0.5)
All_ptD_j3 = rt.TH1F("All_ptD_j3","All_ptD_j3;p_{T}D(j);Events",40, 0, 1.2)
All_tau1_j3 = rt.TH1F("All_tau1_j3","All_tau1_j3;#tau_{1}(j);Events",40, 0, 0.8)
All_tau2_j3 = rt.TH1F("All_tau2_j3","All_tau2_j3;#tau_{2}(j);Events",40, 0, 0.65)
All_tau3_j3 = rt.TH1F("All_tau3_j3","All_tau3_j3;#tau_{3}(j);Events",40, 0, 0.35)
All_tau21_j3 = rt.TH1F("All_tau21_j3","All_tau21_j3;#tau_{21}(j);Events",40, 0, 1.3)
All_tau32_j3 = rt.TH1F("All_tau32_j3","All_tau32_j3;#tau_{32}(j);Events",40, 0, 1.3)
All_SDM_j3 = rt.TH1F("All_SDM_j3","All_SDM_j3; m_{SD}(j);Events",40,0,msd_max)
All_pmult_j3 = rt.TH1F("All_pmult_j3","All_pmult_j3;Number of Particles in a Jet;Events",40,0,400)

All_Pt_j4 = rt.TH1F("All_Pt_j4","All_Pt_j4;p_{T}(j) [GeV];Events",40, 0, 3500)
All_Eta_j4 = rt.TH1F("All_Eta_j4","All_Eta_j4;#eta(j);Events",15, -7, 7)
All_Phi_j4 = rt.TH1F("All_Phi_j4","All_Phi_j4;#phi(j);Events",15, -3.2, 3.2)
All_DPhi_j4 = rt.TH1F("All_DPhi_j4","All_DPhi_j4; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
All_M_j4 = rt.TH1F("All_M_j4","All_M_j4; M(j);Events",40,0,1200)
All_mult_j4 = rt.TH1F("All_mult_j4","All_mult_j4;Number of Jets;Events",7,0,7)
All_axisMajor_j4 = rt.TH1F("All_axisMajor_j4","All_axisMajor_j4;#sigma_{major}(j);Events",40, 0, 0.5)
All_axisMinor_j4 = rt.TH1F("All_axisMinor_j4","All_axisMinor_j4;#sigma_{minor}(j);Events",40, 0, 0.3)
All_momentGirth_j4 = rt.TH1F("All_momentGirth_j4","All_momentGirth_j4;girth(j);Events",40, 0, 0.5)
All_ptD_j4 = rt.TH1F("All_ptD_j4","All_ptD_j4;p_{T}D(j);Events",40, 0, 1.2)
All_tau1_j4 = rt.TH1F("All_tau1_j4","All_tau1_j4;#tau_{1}(j);Events",40, 0, 0.8)
All_tau2_j4 = rt.TH1F("All_tau2_j4","All_tau2_j4;#tau_{2}(j);Events",40, 0, 0.65)
All_tau3_j4 = rt.TH1F("All_tau3_j4","All_tau3_j4;#tau_{3}(j);Events",40, 0, 0.35)
All_tau21_j4 = rt.TH1F("All_tau21_j4","All_tau21_j4;#tau_{21}(j);Events",40, 0, 1.3)
All_tau32_j4 = rt.TH1F("All_tau32_j4","All_tau32_j4;#tau_{32}(j);Events",40, 0, 1.3)
All_SDM_j4 = rt.TH1F("All_SDM_j4","All_SDM_j4; m_{SD}(j);Events",40,0,msd_max)
All_pmult_j4 = rt.TH1F("All_pmult_j4","All_pmult_j4;Number of Particles in a Jet;Events",40,0,400)

All_Pt_j4p = rt.TH1F("All_Pt_j4p","All_Pt_j4p;p_{T}(j) [GeV];Events",40, 0, 3500)
All_Eta_j4p = rt.TH1F("All_Eta_j4p","All_Eta_j4p;#eta(j);Events",15, -7, 7)
All_Phi_j4p = rt.TH1F("All_Phi_j4p","All_Phi_j4p;#phi(j);Events",15, -3.2, 3.2)
All_DPhi_j4p = rt.TH1F("All_DPhi_j4p","All_DPhi_j4p; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
All_M_j4p = rt.TH1F("All_M_j4p","All_M_j4p; M(j);Events",40,0,1200)
All_mult_j4p = rt.TH1F("All_mult_j4p","All_mult_j4p;Number of Jets;Events",7,0,7)
All_axisMajor_j4p = rt.TH1F("All_axisMajor_j4p","All_axisMajor_j4p;#sigma_{major}(j);Events",40, 0, 0.5)
All_axisMinor_j4p = rt.TH1F("All_axisMinor_j4p","All_axisMinor_j4p;#sigma_{minor}(j);Events",40, 0, 0.3)
All_momentGirth_j4p = rt.TH1F("All_momentGirth_j4p","All_momentGirth_j4p;girth(j);Events",40, 0, 0.5)
All_ptD_j4p = rt.TH1F("All_ptD_j4p","All_ptD_j4p;p_{T}D(j);Events",40, 0, 1.2)
All_tau1_j4p = rt.TH1F("All_tau1_j4p","All_tau1_j4p;#tau_{1}(j);Events",40, 0, 0.8)
All_tau2_j4p = rt.TH1F("All_tau2_j4p","All_tau2_j4p;#tau_{2}(j);Events",40, 0, 0.65)
All_tau3_j4p = rt.TH1F("All_tau3_j4p","All_tau3_j4p;#tau_{3}(j);Events",40, 0, 0.35)
All_tau21_j4p = rt.TH1F("All_tau21_j4p","All_tau21_j4p;#tau_{21}(j);Events",40, 0, 1.3)
All_tau32_j4p = rt.TH1F("All_tau32_j4p","All_tau32_j4p;#tau_{32}(j);Events",40, 0, 1.3)
All_SDM_j4p = rt.TH1F("All_SDM_j4p","All_SDM_j4p; m_{SD}(j);Events",40,0,msd_max)
All_pmult_j4p = rt.TH1F("All_pmult_j4p","All_pmult_j4p;Number of Particles in a Jet;Events",40,0,400)

# ST variables
ST = rt.TH1F("ST","ST;S_{T} [GeV];Events",50, 0, 10000)
MET = rt.TH1F("MET","MET;#slash{E}_{T} [GeV];Events",50, 0, 2500)
ST_2jet = rt.TH1F("ST_2jet","ST_2jet;S_{T} [GeV];Events",50, 0, 10000)
ST_3jet = rt.TH1F("ST_3jet","ST_3jet;S_{T} [GeV];Events",50, 0, 10000)
ST_4jet = rt.TH1F("ST_4jet","ST_4jet;S_{T} [GeV];Events",50, 0, 10000)
ST_4plusjet = rt.TH1F("ST_4plusjet","ST_4plusjet;S_{T} [GeV];Events",50, 0, 10000)

# RT variable
RT = rt.TH1F("RT","RT;R_{T};Events",50, 0,metRx*1.01)

# Mass variables
## MT
MT_2j = rt.TH1F("MT_2j","MT_2j;m_{T} (GeV);Events",40, 0, 9000)
MT_3j = rt.TH1F("MT_3j","MT_3j;m_{T} (GeV);Events",40, 0, 9000)
MT_4j = rt.TH1F("MT_4j","MT_4j;m_{T} (GeV);Events",40, 0, 9000)
MT_4pj = rt.TH1F("MT_4pj","MT_4pj;m_{T} (GeV);Events",40, 0, 9000)

MT_m2000_DR_f2 = rt.TH1F("MT_m2000_DR_f2","MT_m2000_DR_f2;#Delta R(j_{1},j_{2});Events",40, 0, 8) # delta R of first 2 jets
MT_l2000_DR_f2 = rt.TH1F("MT_l2000_DR_f2","MT_l2000_DR_f2;#Delta R(j_{1},j_{2});Events",40, 0, 8)
MT_m2000_DP_f2 = rt.TH1F("MT_m2000_DP_f2","MT_m2000_DP_f2;#Delta #phi(j_{1},j_{2});Events",40, 0, 3.5) # delta phi of first 2 jets
MT_l2000_DP_f2 = rt.TH1F("MT_l2000_DP_f2","MT_l2000_DP_f2;#Delta #phi(j_{1},j_{2});Events",40, 0, 3.5)
MT_m2000_DE_f2 = rt.TH1F("MT_m2000_DE_f2","MT_m2000_DE_f2;#Delta #eta(j_{1},j_{2});Events",40, 0, 5) # delta eta of first 2 jets
MT_l2000_DE_f2 = rt.TH1F("MT_l2000_DE_f2","MT_l2000_DE_f2;#Delta #eta(j_{1},j_{2});Events",40, 0, 5)

MT_R = rt.TH1F("MT_R","MT_R;m_{T} (GeV);Events",40, 0, 9000)
MT_W = rt.TH1F("MT_W","MT_W;m_{T} (GeV);Events",40, 0, 9000)
MT_R_DR_f2 = rt.TH1F("MT_R_DR_f2","MT_R_DR_f2;#Delta R(j_{1},j_{2});Events",40, 0, 8) # delta R of first 2 jets
MT_W_DR_f2 = rt.TH1F("MT_W_DR_f2","MT_W_DR_f2;#Delta R(j_{1},j_{2});Events",40, 0, 8)
MT_R_DP_f2 = rt.TH1F("MT_R_DP_f2","MT_R_DP_f2;#Delta #phi(j_{1},j_{2});Events",40, 0, 3.5) # delta phi of first 2 jets
MT_W_DP_f2 = rt.TH1F("MT_W_DP_f2","MT_W_DP_f2;#Delta #phi(j_{1},j_{2});Events",40, 0, 3.5)
MT_R_DE_f2 = rt.TH1F("MT_R_DE_f2","MT_R_DE_f2;#Delta #eta(j_{1},j_{2});Events",40, 0, 5) # delta eta of first 2 jets
MT_W_DE_f2 = rt.TH1F("MT_W_DE_f2","MT_W_DE_f2;#Delta #eta(j_{1},j_{2});Events",40, 0, 5)

## MT2
MT2_R = rt.TH1F("MT2_R","MT2_R; M_{T2} (GeV);Events",40, 0, 7000)
MT2_W = rt.TH1F("MT2_W","MT2_W; M_{T2} (GeV);Events",40, 0, 7000)
MT2_DeltaR_R = rt.TH1F("MT2_DeltaR_R","MT2_DeltaR_R; #Delta R(j_{1},j_{2});Events",40, 0, 8)
MT2_DeltaR_W = rt.TH1F("MT2_DeltaR_W","MT2_DeltaR_W; #Delta R(j_{1},j_{2});Events",40, 0, 8)
MT2_DeltaP_R = rt.TH1F("MT2_DeltaP_R","MT2_DeltaP_R; #Delta #phi(j_{1},j_{2});Events",40, 0, 3.5)
MT2_DeltaP_W = rt.TH1F("MT2_DeltaP_W","MT2_DeltaP_W; #Delta #phi(j_{1},j_{2});Events",40, 0, 3.5)
MT2_DeltaE_R = rt.TH1F("MT2_DeltaE_R","MT2_DeltaE_R; #Delta #eta(j_{1},j_{2});Events",40, 0, 5)
MT2_DeltaE_W = rt.TH1F("MT2_DeltaE_W","MT2_DeltaE_W; #Delta #eta(j_{1},j_{2});Events",40, 0, 5)
MT2_mp = rt.TH1F("MT2_mp","MT2_mp; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_DPhiPair = rt.TH1F("MT2_mp_DPhiPair","MT2_mp_DPhiPair; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_DPhiPair_Pt = rt.TH1F("MT2_mp_DPhiPair_Pt","MT2_mp_DPhiPair_Pt; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_PhiReq = rt.TH1F("MT2_mp_PhiReq","MT2_mp_PhiReq; M_{T2} (GeV);Events",40, 0, 7000)
MT2_Pt_mp = rt.TH1F("MT2_Pt_mp","MT2_Pt_mp; M_{T2} (GeV);Events",40, 0, 7000)
MT2_Pt_mp_DPhiPair = rt.TH1F("MT2_Pt_mp_DPhiPair","MT2_Pt_mp_DPhiPair; M_{T2} (GeV);Events",40, 0, 7000)


MT2_mp_f4 = rt.TH1F("MT2_mp_f4","MT2_mp_f4; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_DPhiPair_f4 = rt.TH1F("MT2_mp_DPhiPair_f4","MT2_mp_DPhiPair_f4; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_PhiReq_f4 = rt.TH1F("MT2_mp_PhiReq_f4","MT2_mp_PhiReq_f4; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_DPhiPair_Pt_f4 = rt.TH1F("MT2_mp_DPhiPair_Pt_f4","MT2_mp_DPhiPair_Pt_f4; M_{T2} (GeV);Events",40, 0, 7000)

MT2_mp_MDiff = rt.TH1F("MT2_mp_MDiff","MT2_mp_MDiff; #Delta M(j_{MT2_mp}) (GeV);Events",40, 0, 2000)
MT2_mp_0M = rt.TH1F("MT2_mp_0M","MT2_mp_0M; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_1M = rt.TH1F("MT2_mp_1M","MT2_mp_1M; M_{T2} (GeV);Events",40, 0, 7000)
MT2_mp_2M = rt.TH1F("MT2_mp_2M","MT2_mp_2M; M_{T2} (GeV);Events",40, 0, 7000)
MT2_R_MDiff = rt.TH1F("MT2_R_MDiff","MT2_R_MDiff; #Delta M(j_{MT2_R}) (GeV);Events",40, 0, 2000)
MT2_W_MDiff = rt.TH1F("MT2_W_MDiff","MT2_W_MDiff; #Delta M(j_{MT2_W}) (GeV);Events",40, 0, 2000)
MT2_R_M = rt.TH1F("MT2_R_M","MT2_R_M; M(j_{MT2_R}) (GeV);Events",40, 0, 5000)
MT2_W_M = rt.TH1F("MT2_W_M","MT2_W_M; M(j_{MT2_W}) (GeV);Events",40, 0, 5000)
MT2_R_Pt = rt.TH1F("MT2_R_Pt","MT2_R_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
MT2_mp_Pt = rt.TH1F("MT2_mp_Pt","MT2_mp_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
MT2_R_DPhiMET = rt.TH1F("MT2_R_DPhiMET","MT2_R_DPhiMET; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
MT2_mp_DPhiMET = rt.TH1F("MT2_mp_DPhiMET","MT2_mp_DPhiMET; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
MT2_R_DPhiMET_d0r = rt.TH1F("MT2_R_DPhiMET_d0r","MT2_R_DPhiMET_d0r; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
MT2_R_DPhiMET_s0r = rt.TH1F("MT2_R_DPhiMET_s0r","MT2_R_DPhiMET_s0r; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
MT2_R_DPhiMET_d1r = rt.TH1F("MT2_R_DPhiMET_d1r","MT2_R_DPhiMET_d1r; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
MT2_R_DPhiMET_s1r = rt.TH1F("MT2_R_DPhiMET_s1r","MT2_R_DPhiMET_s1r; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
MT2_R_DPhi_pair = rt.TH1F("MT2_R_DPhi_pair","MT2_R_DPhi_pair; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)

MT2_mp_MT2_R = rt.TH2F("MT2_mp_MT2_R","MT2_mp_MT2_R;M_{T2};Events",40,0, 7000, 40, 0, 7000)

# Just a filler function that doesn't do anything more than allowing
# implementation of functions that require 2D histograms as arguments
filler2D = rt.TH2F("LD_SMpmult_2D","LD_SMpmult_2D;Number of Particles in a Jet;Events",40, 0, 7000, 40, 0, 7000)

tf = rt.TFile.Open("Skims/s1_mMed-{}_mDark-{}_rinv-{}_alpha-{}.root".format(mMed,mDark,rInv,Alpha))
tr = tf.Get('GenMassAnalyzer/tree')
tf2 = rt.TFile.Open("Skims/skims_softDrop_tau.root")
tr2 = tf2.Get('GenMassAnalyzer/tree')

nEvents = tr.GetEntries()

# used for scatterplots
pT_2D = []
dphi_2D = []
eta_2D = []
phi_2D = []
M_2D = []
gir_2D = []
axisMj_2D = []
axisMn_2D = []
ptD_2D = []
np_2D = []
tau1_2D = []
tau2_2D = []
tau3_2D = []
tau21_2D = []
tau32_2D = []
SDM_2D = []

isMT2R = 0
is4lead = 0
MT2_mp_f4_Events = []
MT2_mp_DPhiPair_f4_Events = []
MT2_mp_DPhiPair_Pt_f4_Events = []
MT2_mp_PhiReq_f4_Events = []

MT2_mp_Events = []
MT2_mp_DPhiPair_Events = []
MT2_mp_DPhiPair_Pt_Events = []
MT2_mp_PhiReq_Events = []

MT2_Pt_mp_Events = []
MT2_Pt_mp_DPhiPair_Events = []

MT2_R_Events = []

for count in range(nEvents):
	tr.GetEntry(count)
	tr2.GetEntry(count)
	if count % 3000 == 0:
		print "Event number " + str(count)

	jets = VecToList(tr.AK8Jets)
	met4p = tr.MET[0]
	nfdMPart = tr.nfdMPart
	SM_Jets = VecToList(tr.AK8_SM_Jets)
	SMM_Jets = VecToList(tr.AK8_SMM_Jets)
	G_Jets = VecToList(tr.AK8_G_Jets)
	QM_Jets = VecToList(tr.AK8_QM_Jets)
	Q_Jets = VecToList(tr.AK8_Q_Jets)
	QM_QJets = VecToList(tr.AK8_QM_QJets)
	QM_GJets = VecToList(tr.AK8_QM_GJets)
	Q_GJets = VecToList(tr.AK8_Q_GJets)
	G_SMJets = VecToList(tr.AK8_G_SMJets)
	QM_SMJets = VecToList(tr.AK8_QM_SMJets)
	Q_SMJets = VecToList(tr.AK8_Q_SMJets)
	LD_lowDFJets = VecToList(tr.AK8_LD_lowDFJets)
	LD_highDFJets = VecToList(tr.AK8_LD_highDFJets)
	LD_SMJets = VecToList(tr.AK8_LD_SMJets)

	G_dpTFrac = VecToList(tr.G_dpTFrac)
	QM_dpTFrac = VecToList(tr.QM_dpTFrac)
	Q_dpTFrac = VecToList(tr.Q_dpTFrac)
	QM_QdpTFrac = VecToList(tr.QM_QdpTFrac)
	QM_GdpTFrac = VecToList(tr.QM_GdpTFrac)
	Q_GdpTFrac = VecToList(tr.Q_GdpTFrac)
	G_SMdpTFrac = VecToList(tr.G_SMdpTFrac)
	QM_SMdpTFrac = VecToList(tr.QM_SMdpTFrac)
	Q_SMdpTFrac = VecToList(tr.Q_SMdpTFrac)
	LD_lowDFdpTFrac = VecToList(tr.LD_lowDFdpTFrac)
	LD_highDFdpTFrac = VecToList(tr.LD_highDFdpTFrac)
	LD_SMdpTFrac = VecToList(tr.LD_SMdpTFrac)

	dgir = VecToList(tr.Dark_girth)
	daMj = VecToList(tr.Dark_axisMj)
	daMn = VecToList(tr.Dark_axisMn)
	dpD = VecToList(tr.Dark_ptD)

	All_gir = VecToList(tr.AK8_girth)
	All_aMj = VecToList(tr.AK8_axisMj)
	All_aMn = VecToList(tr.AK8_axisMn)
	All_pD = VecToList(tr.AK8_ptD)

	SM_gir = VecToList(tr.SM_girth)
	SM_aMj = VecToList(tr.SM_axisMj)
	SM_aMn = VecToList(tr.SM_axisMn)
	SM_pD = VecToList(tr.SM_ptD)

	SMM_gir = VecToList(tr.SMM_girth)
	SMM_aMj = VecToList(tr.SMM_axisMj)
	SMM_aMn = VecToList(tr.SMM_axisMn)
	SMM_pD = VecToList(tr.SMM_ptD)

	G_gir = VecToList(tr.G_girth)
	G_aMj = VecToList(tr.G_axisMj)
	G_aMn = VecToList(tr.G_axisMn)
	G_pD = VecToList(tr.G_ptD)

	QM_gir = VecToList(tr.QM_girth)
	QM_aMj = VecToList(tr.QM_axisMj)
	QM_aMn = VecToList(tr.QM_axisMn)
	QM_pD = VecToList(tr.QM_ptD)

	Q_gir = VecToList(tr.Q_girth)
	Q_aMj = VecToList(tr.Q_axisMj)
	Q_aMn = VecToList(tr.Q_axisMn)
	Q_pD = VecToList(tr.Q_ptD)

	QM_Qgir = VecToList(tr.QM_Qgirth)
	QM_QaMj = VecToList(tr.QM_QaxisMj)
	QM_QaMn = VecToList(tr.QM_QaxisMn)
	QM_QpD = VecToList(tr.QM_QptD)

	QM_Ggir = VecToList(tr.QM_Ggirth)
	QM_GaMj = VecToList(tr.QM_GaxisMj)
	QM_GaMn = VecToList(tr.QM_GaxisMn)
	QM_GpD = VecToList(tr.QM_GptD)

	Q_Ggir = VecToList(tr.Q_Ggirth)
	Q_GaMj = VecToList(tr.Q_GaxisMj)
	Q_GaMn = VecToList(tr.Q_GaxisMn)
	Q_GpD = VecToList(tr.Q_GptD)

	G_SMgir = VecToList(tr.G_SMgirth)
	G_SMaMj = VecToList(tr.G_SMaxisMj)
	G_SMaMn = VecToList(tr.G_SMaxisMn)
	G_SMpD = VecToList(tr.G_SMptD)

	QM_SMgir = VecToList(tr.QM_SMgirth)
	QM_SMaMj = VecToList(tr.QM_SMaxisMj)
	QM_SMaMn = VecToList(tr.QM_SMaxisMn)
	QM_SMpD = VecToList(tr.QM_SMptD)

	Q_SMgir = VecToList(tr.Q_SMgirth)
	Q_SMaMj = VecToList(tr.Q_SMaxisMj)
	Q_SMaMn = VecToList(tr.Q_SMaxisMn)
	Q_SMpD = VecToList(tr.Q_SMptD)

	LD_lowDFgir = VecToList(tr.LD_lowDFgirth)
	LD_lowDFaMj = VecToList(tr.LD_lowDFaxisMj)
	LD_lowDFaMn = VecToList(tr.LD_lowDFaxisMn)
	LD_lowDFpD = VecToList(tr.LD_lowDFptD)

	LD_highDFgir = VecToList(tr.LD_highDFgirth)
	LD_highDFaMj = VecToList(tr.LD_highDFaxisMj)
	LD_highDFaMn = VecToList(tr.LD_highDFaxisMn)
	LD_highDFpD = VecToList(tr.LD_highDFptD)

	LD_SMgir = VecToList(tr.LD_SMgirth)
	LD_SMaMj = VecToList(tr.LD_SMaxisMj)
	LD_SMaMn = VecToList(tr.LD_SMaxisMn)
	LD_SMpD = VecToList(tr.LD_SMptD)

	npAll_Jets = VecToList(tr.npAK8Jets)
	npDJets = VecToList(tr.npDJets)
	npSM_Jets = VecToList(tr.npSM_Jets)
	npSMM_Jets = VecToList(tr.npSMM_Jets)
	npG_Jets = VecToList(tr.npG_Jets)
	npQM_Jets = VecToList(tr.npQM_Jets)
	npQ_Jets = VecToList(tr.npQ_Jets)
	npQM_QJets = VecToList(tr.npQM_QJets)
	npQM_GJets = VecToList(tr.npQM_GJets)
	npQ_GJets = VecToList(tr.npQ_GJets)
	npG_SMJets = VecToList(tr.npG_SMJets)
	npQM_SMJets = VecToList(tr.npQM_SMJets)
	npQ_SMJets = VecToList(tr.npQ_SMJets)
	npLD_lowDFJets = VecToList(tr.npLD_lowDFJets)
	npLD_highDFJets = VecToList(tr.npLD_highDFJets)
	npLD_SMJets = VecToList(tr.npLD_SMJets)

	# taus
	All_Tau1 = VecToList(tr2.AK8_Tau1)
	All_Tau2 = VecToList(tr2.AK8_Tau2)
	All_Tau3 = VecToList(tr2.AK8_Tau3)
	SM_Tau1 = VecToList(tr2.SM_Tau1)
	SM_Tau2 = VecToList(tr2.SM_Tau2)
	SM_Tau3 = VecToList(tr2.SM_Tau3)
	SMM_Tau1 = VecToList(tr2.SMM_Tau1)
	SMM_Tau2 = VecToList(tr2.SMM_Tau2)
	SMM_Tau3 = VecToList(tr2.SMM_Tau3)
	G_Tau1 = VecToList(tr2.G_Tau1)
	G_Tau2 = VecToList(tr2.G_Tau2)
	G_Tau3 = VecToList(tr2.G_Tau3)
	QM_Tau1 = VecToList(tr2.QM_Tau1)
	QM_Tau2 = VecToList(tr2.QM_Tau2)
	QM_Tau3 = VecToList(tr2.QM_Tau3)
	Q_Tau1 = VecToList(tr2.Q_Tau1)
	Q_Tau2 = VecToList(tr2.Q_Tau2)
	Q_Tau3 = VecToList(tr2.Q_Tau3)
	QM_QTau1 = VecToList(tr2.QM_QTau1)
	QM_QTau2 = VecToList(tr2.QM_QTau2)
	QM_QTau3 = VecToList(tr2.QM_QTau3)
	QM_GTau1 = VecToList(tr2.QM_GTau1)
	QM_GTau2 = VecToList(tr2.QM_GTau2)
	QM_GTau3 = VecToList(tr2.QM_GTau3)
	Q_GTau1 = VecToList(tr2.Q_GTau1)
	Q_GTau2 = VecToList(tr2.Q_GTau2)
	Q_GTau3 = VecToList(tr2.Q_GTau3)
	G_SMTau1 = VecToList(tr2.G_SMTau1)
	G_SMTau2 = VecToList(tr2.G_SMTau2)
	G_SMTau3 = VecToList(tr2.G_SMTau3)
	QM_SMTau1 = VecToList(tr2.QM_SMTau1)
	QM_SMTau2 = VecToList(tr2.QM_SMTau2)
	QM_SMTau3 = VecToList(tr2.QM_SMTau3)
	Q_SMTau1 = VecToList(tr2.Q_SMTau1)
	Q_SMTau2 = VecToList(tr2.Q_SMTau2)
	Q_SMTau3 = VecToList(tr2.Q_SMTau3)
	LD_lowDFTau1 = VecToList(tr2.LD_lowDFTau1)
	LD_lowDFTau2 = VecToList(tr2.LD_lowDFTau2)
	LD_lowDFTau3 = VecToList(tr2.LD_lowDFTau3)
	LD_highDFTau1 = VecToList(tr2.LD_highDFTau1)
	LD_highDFTau2 = VecToList(tr2.LD_highDFTau2)
	LD_highDFTau3 = VecToList(tr2.LD_highDFTau3)
	LD_SMTau1 = VecToList(tr2.LD_SMTau1)
	LD_SMTau2 = VecToList(tr2.LD_SMTau2)
	LD_SMTau3 = VecToList(tr2.LD_SMTau3)

	# soft drop jets and corresponding dark pT fraction
	All_SDJets = VecToList(tr2.AK8_SDJets)
	SM_SDJets = VecToList(tr2.SM_SDJets)
	SMM_SDJets = VecToList(tr2.SMM_SDJets)
	G_SDJets = VecToList(tr2.G_SDJets)
	QM_SDJets = VecToList(tr2.QM_SDJets)
	Q_SDJets = VecToList(tr2.Q_SDJets)
	QM_QSDJets = VecToList(tr2.QM_QSDJets)
	QM_GSDJets = VecToList(tr2.QM_GSDJets)
	Q_GSDJets = VecToList(tr2.Q_GSDJets)
	G_SMSDJets = VecToList(tr2.G_SMSDJets)
	QM_SMSDJets = VecToList(tr2.QM_SMSDJets)
	Q_SMSDJets = VecToList(tr2.Q_SMSDJets)
	LD_lowDFSDJets = VecToList(tr2.LD_lowDFSDJets)
	LD_highDFSDJets = VecToList(tr2.LD_highDFSDJets)
	LD_SMSDJets = VecToList(tr2.LD_SMSDJets)

	G_SDdptFrac = VecToList(tr2.G_SDdptFrac)
	QM_SDdptFrac = VecToList(tr2.QM_SDdptFrac)
	Q_SDdptFrac = VecToList(tr2.Q_SDdptFrac)
	QM_QSDdptFrac = VecToList(tr2.QM_QSDdptFrac)
	QM_GSDdptFrac = VecToList(tr2.QM_GSDdptFrac)
	Q_GSDdptFrac = VecToList(tr2.Q_GSDdptFrac)
	G_SMSDdptFrac = VecToList(tr2.G_SMSDdptFrac)
	QM_SMSDdptFrac = VecToList(tr2.QM_SMSDdptFrac)
	Q_SMSDdptFrac = VecToList(tr2.Q_SMSDdptFrac)
	LD_lowDFSDdptFrac = VecToList(tr2.LD_lowDFSDdptFrac)
	LD_highDFSDdptFrac = VecToList(tr2.LD_highDFSDdptFrac)
	LD_SMSDdptFrac = VecToList(tr2.LD_SMSDdptFrac)
	SM_SDdptFrac = VecToList(tr2.SM_SDdptFrac)
	SMM_SDdptFrac = VecToList(tr2.SMM_SDdptFrac)

	# calculating mT
	MedPair1 = tr.AK8_MedPair1
	MedPair2 = tr.AK8_MedPair2
	MTPair = tr.AK8_MTPair

	METv = met4p.Pt()
	METPhiv = met4p.Phi()

	SM_OnlyJets = SM_Jets + SMM_Jets
	SM_OnlydpTFrac = [0]*len(SM_OnlyJets)
	SM_Onlygir = SM_gir + SMM_gir
	SM_OnlyaMj = SM_aMj + SMM_aMj
	SM_OnlyaMn = SM_aMn + SMM_aMn
	SM_OnlypD = SM_pD + SMM_pD
	npSM_OnlyJets = npSM_Jets + npSMM_Jets
	SM_OnlyTau1 = SM_Tau1 + SMM_Tau1
	SM_OnlyTau2 = SM_Tau2 + SMM_Tau2
	SM_OnlyTau3 = SM_Tau3 + SMM_Tau3
	SM_OnlySDJets = SM_SDJets + SMM_SDJets
	SM_OnlySDdpTFrac = [0]*len(SM_OnlySDJets)

	Mix_Jets = G_SMJets + LD_SMJets + LD_lowDFJets
	Mix_dpTFrac = G_SMdpTFrac + LD_SMdpTFrac + LD_lowDFdpTFrac
	Mix_gir = G_SMgir + LD_SMgir + LD_lowDFgir
	Mix_aMj = G_SMaMj + LD_SMaMj + LD_lowDFaMj
	Mix_aMn = G_SMaMn + LD_SMaMn + LD_lowDFaMn
	Mix_pD = G_SMpD + LD_SMpD + LD_lowDFpD
	npMix_Jets = npG_SMJets + npLD_SMJets + npLD_lowDFJets
	Mix_Tau1 = G_SMTau1 + LD_SMTau1 + LD_lowDFTau1
	Mix_Tau2 = G_SMTau2 + LD_SMTau2 + LD_lowDFTau2
	Mix_Tau3 = G_SMTau3 + LD_SMTau3 + LD_lowDFTau3
	Mix_SDJets = G_SMSDJets + LD_SMSDJets + LD_lowDFSDJets
	Mix_SDdpTFrac = G_SMSDdptFrac + LD_SMSDdptFrac + LD_lowDFSDdptFrac

	DarkerMix_Jets = QM_SMJets + Q_SMJets + QM_QJets
	DarkerMix_dpTFrac = QM_SMdpTFrac + Q_SMdpTFrac + QM_QdpTFrac
	DarkerMix_gir = QM_SMgir + Q_SMgir + QM_Qgir
	DarkerMix_aMj = QM_SMaMj + Q_SMaMj + QM_QaMj
	DarkerMix_aMn = QM_SMaMn + Q_SMaMn + QM_QaMn
	DarkerMix_pD = QM_SMpD + Q_SMpD + QM_QpD
	npDarkerMix_Jets = npQM_SMJets + npQ_SMJets + npQM_QJets
	DarkerMix_Tau1 = QM_SMTau1 + Q_SMTau1 + QM_QTau1
	DarkerMix_Tau2 = QM_SMTau2 + Q_SMTau2 + QM_QTau2
	DarkerMix_Tau3 = QM_SMTau3 + Q_SMTau3 + QM_QTau3
	DarkerMix_SDJets = QM_SMSDJets + Q_SMSDJets + QM_QSDJets
	DarkerMix_SDdpTFrac = QM_SMSDdptFrac + Q_SMSDdptFrac + QM_QSDdptFrac

	Dark_Jets = QM_Jets + QM_GJets + Q_GJets + Q_Jets
	Dark_dpTFrac = QM_dpTFrac + QM_GdpTFrac + Q_GdpTFrac + Q_dpTFrac
	Dark_gir = QM_gir + QM_Ggir + Q_Ggir + Q_gir
	Dark_aMj = QM_aMj + QM_GaMj + Q_GaMj + Q_aMj
	Dark_aMn = QM_aMn + QM_GaMn + Q_GaMn + Q_aMn
	Dark_pD = QM_pD + QM_GpD + Q_GpD + Q_pD
	npDark_Jets = npQM_Jets + npQM_GJets + npQ_GJets + npQ_Jets
	Dark_Tau1 = QM_Tau1 + QM_GTau1 + Q_GTau1 + Q_Tau1
	Dark_Tau2 = QM_Tau2 + QM_GTau2 + Q_GTau2 + Q_Tau2
	Dark_Tau3 = QM_Tau3 + QM_GTau3 + Q_GTau3 + Q_Tau3
	Dark_SDJets = QM_SDJets + QM_GSDJets + Q_GSDJets + Q_SDJets
	Dark_SDdpTFrac = QM_SDdptFrac + QM_GSDdptFrac + Q_GSDdptFrac + Q_SDdptFrac

	SoftDark_Jets = G_Jets + LD_highDFJets
	SoftDark_dpTFrac = G_dpTFrac + LD_highDFdpTFrac
	SoftDark_gir = G_gir + LD_highDFgir
	SoftDark_aMj = G_aMj + LD_highDFaMj
	SoftDark_aMn = G_aMn + LD_highDFaMn
	SoftDark_pD = G_pD + LD_highDFpD
	npSoftDark_Jets = npG_Jets + npLD_highDFJets
	SoftDark_Tau1 = G_Tau1 + LD_highDFTau1
	SoftDark_Tau2 = G_Tau2 + LD_highDFTau2
	SoftDark_Tau3 = G_Tau3 + LD_highDFTau3
	SoftDark_SDJets = G_SDJets + LD_highDFSDJets
	SoftDark_SDdpTFrac = G_SDdptFrac + LD_highDFSDdptFrac

	# Delta Phi Min for all jets
	# and kinematic variables for the jet closest to MET
	closestToMET(SM_OnlyJets,SM_OnlySDJets,SM_OnlyDPhiMin,SM_OnlyPt_Cl,SM_OnlyEta_Cl,
	SM_OnlyPhi_Cl,SM_OnlyDPhi_Cl,SM_OnlyM_Cl,
	SM_OnlyaxisMajor_Cl,SM_OnlyaxisMinor_Cl,SM_OnlymomentGirth_Cl,SM_OnlyptD_Cl,
	SM_Onlytau1_Cl,SM_Onlytau2_Cl,SM_Onlytau3_Cl,
	SM_Onlytau21_Cl,SM_Onlytau32_Cl,SM_Onlymult_Cl,SM_Onlypmult_Cl,SM_OnlySDM_Cl,
	SM_OnlyaMj,SM_OnlyaMn,SM_Onlygir,SM_OnlypD,
	SM_OnlyTau1,SM_OnlyTau2,SM_OnlyTau3,npSM_OnlyJets)
	closestToMET(Mix_Jets,Mix_SDJets,Mix_DPhiMin,Mix_Pt_Cl,Mix_Eta_Cl,
	Mix_Phi_Cl,Mix_DPhi_Cl,Mix_M_Cl,
	Mix_axisMajor_Cl,Mix_axisMinor_Cl,Mix_momentGirth_Cl,Mix_ptD_Cl,
	Mix_tau1_Cl,Mix_tau2_Cl,Mix_tau3_Cl,
	Mix_tau21_Cl,Mix_tau32_Cl,Mix_mult_Cl,Mix_pmult_Cl,Mix_SDM_Cl,
	Mix_aMj,Mix_aMn,Mix_gir,Mix_pD,
	Mix_Tau1,Mix_Tau2,Mix_Tau3,npMix_Jets)
	closestToMET(Dark_Jets,Dark_SDJets,Dark_DPhiMin,Dark_Pt_Cl,Dark_Eta_Cl,
	Dark_Phi_Cl,Dark_DPhi_Cl,Dark_M_Cl,
	Dark_axisMajor_Cl,Dark_axisMinor_Cl,Dark_momentGirth_Cl,Dark_ptD_Cl,
	Dark_tau1_Cl,Dark_tau2_Cl,Dark_tau3_Cl,
	Dark_tau21_Cl,Dark_tau32_Cl,Dark_mult_Cl,Dark_pmult_Cl,Dark_SDM_Cl,
	Dark_aMj,Dark_aMn,Dark_gir,Dark_pD,
	Dark_Tau1,Dark_Tau2,Dark_Tau3,npDark_Jets)
	closestToMET(jets,All_SDJets,All_DPhiMin,All_Pt_Cl,All_Eta_Cl,
	All_Phi_Cl,All_DPhi_Cl,All_M_Cl,
	All_axisMajor_Cl,All_axisMinor_Cl,All_momentGirth_Cl,All_ptD_Cl,
	All_tau1_Cl,All_tau2_Cl,All_tau3_Cl,
	All_tau21_Cl,All_tau32_Cl,All_mult_Cl,All_pmult_Cl,All_SDM_Cl,
	All_aMj,All_aMn,All_gir,All_pD,
	All_Tau1,All_Tau2,All_Tau3,npAll_Jets)

	if len(jets) > 0:
		jNFill(jets,All_Pt_j1,All_Eta_j1,All_Phi_j1,All_DPhi_j1,
		All_M_j1,All_axisMajor_j1,All_axisMinor_j1,All_momentGirth_j1,All_ptD_j1,
		All_tau1_j1,All_tau2_j1,All_tau3_j1,All_tau21_j1,All_tau32_j1,All_mult_j1,
		All_pmult_j1,All_aMj,All_aMn,All_gir,All_pD,
		All_Tau1,All_Tau2,All_Tau3,npAll_Jets,[0])
		if len(jets) > 1:
			jNFill(jets,All_Pt_j2,All_Eta_j2,All_Phi_j2,All_DPhi_j2,
			All_M_j2,All_axisMajor_j2,All_axisMinor_j2,All_momentGirth_j2,All_ptD_j2,
			All_tau1_j2,All_tau2_j2,All_tau3_j2,All_tau21_j2,All_tau32_j2,All_mult_j2,
			All_pmult_j2,All_aMj,All_aMn,All_gir,All_pD,
			All_Tau1,All_Tau2,All_Tau3,npAll_Jets,[1])
			if len(jets) > 2:
				jNFill(jets,All_Pt_j3,All_Eta_j3,All_Phi_j3,All_DPhi_j3,
				All_M_j3,All_axisMajor_j3,All_axisMinor_j3,All_momentGirth_j3,All_ptD_j3,
				All_tau1_j3,All_tau2_j3,All_tau3_j3,All_tau21_j3,All_tau32_j3,All_mult_j3,
				All_pmult_j3,All_aMj,All_aMn,All_gir,All_pD,
				All_Tau1,All_Tau2,All_Tau3,npAll_Jets,[2])
				if len(jets) > 3:
					jNFill(jets,All_Pt_j4,All_Eta_j4,All_Phi_j4,All_DPhi_j4,
					All_M_j4,All_axisMajor_j4,All_axisMinor_j4,All_momentGirth_j4,All_ptD_j4,
					All_tau1_j4,All_tau2_j4,All_tau3_j4,All_tau21_j4,All_tau32_j4,All_mult_j4,
					All_pmult_j4,All_aMj,All_aMn,All_gir,All_pD,
					All_Tau1,All_Tau2,All_Tau3,npAll_Jets,[3])
					if len(jets) > 4:
						jNFill(jets,All_Pt_j4p,All_Eta_j4p,All_Phi_j4p,All_DPhi_j4p,
						All_M_j4p,All_axisMajor_j4p,All_axisMinor_j4p,All_momentGirth_j4p,All_ptD_j4p,
						All_tau1_j4p,All_tau2_j4p,All_tau3_j4p,All_tau21_j4p,All_tau32_j4p,All_mult_j4p,
						All_pmult_j4p,All_aMj,All_aMn,All_gir,All_pD,
						All_Tau1,All_Tau2,All_Tau3,npAll_Jets,range(4,len(jets)))

	if len(All_SDJets) > 0:
		jNSDFill(All_SDJets,All_SDM_j1,[0])
		if len(All_SDJets) > 1:
			jNSDFill(All_SDJets,All_SDM_j2,[1])
			if len(All_SDJets) > 2:
				jNSDFill(All_SDJets,All_SDM_j3,[2])
				if len(All_SDJets) > 3:
					jNSDFill(All_SDJets,All_SDM_j4,[3])
					if len(All_SDJets) > 4:
						jNSDFill(All_SDJets,All_SDM_j4p,range(4,len(All_SDJets)))

	fillDPT([0]*len(SM_Jets),jet_SM_dpTFrac,SM_Jets)
	fillDPT([0]*len(SMM_Jets),jet_SMM_dpTFrac,SMM_Jets)
	fillDPT(G_dpTFrac,jet_G_dpTFrac,G_Jets)
	fillDPT(QM_dpTFrac,jet_QM_dpTFrac,QM_Jets)
	fillDPT(Q_dpTFrac,jet_Q_dpTFrac,Q_Jets)
	fillDPT(QM_QdpTFrac,jet_QM_QdpTFrac,QM_QJets)
	fillDPT(QM_GdpTFrac,jet_QM_GdpTFrac,QM_GJets)
	fillDPT(Q_GdpTFrac,jet_Q_GdpTFrac,Q_GJets)
	fillDPT(G_SMdpTFrac,jet_G_SMdpTFrac,G_SMJets)
	fillDPT(QM_SMdpTFrac,jet_QM_SMdpTFrac,QM_SMJets)
	fillDPT(Q_SMdpTFrac,jet_Q_SMdpTFrac,Q_SMJets)
	fillDPT(LD_lowDFdpTFrac,jet_LD_lowDFdpTFrac,LD_lowDFJets)
	fillDPT(LD_highDFdpTFrac,jet_LD_highDFdpTFrac,LD_highDFJets)
	fillDPT(LD_SMdpTFrac,jet_LD_SMdpTFrac,LD_SMJets)

	fillDPT(SM_OnlydpTFrac,jet_SM_OnlydpTFrac,SM_OnlyJets)
	fillDPT(Mix_dpTFrac,jet_Mix_dpTFrac,Mix_Jets)
	fillDPT(DarkerMix_dpTFrac,jet_DarkerMix_dpTFrac,DarkerMix_Jets)
	fillDPT(Dark_dpTFrac,jet_Dark_dpTFrac,Dark_Jets)
	fillDPT(SoftDark_dpTFrac,jet_SoftDark_dpTFrac,SoftDark_Jets)

	fillM(SM_Jets,SM_M,SM_M_2D,[0]*len(SM_Jets))
	fillM(SMM_Jets,SMM_M,SMM_M_2D,[0]*len(SMM_Jets))
	fillM(G_Jets,G_M,G_M_2D,G_dpTFrac)
	fillM(QM_Jets,QM_M,QM_M_2D,QM_dpTFrac)
	fillM(Q_Jets,Q_M,Q_M_2D,Q_dpTFrac)
	fillM(QM_QJets,QM_QM,QM_QM_2D,QM_QdpTFrac)
	fillM(QM_GJets,QM_GM,QM_GM_2D,QM_GdpTFrac)
	fillM(Q_GJets,Q_GM,Q_GM_2D,Q_GdpTFrac)
	fillM(G_SMJets,G_SMM,G_SMM_2D,G_SMdpTFrac)
	fillM(QM_SMJets,QM_SMM,QM_SMM_2D,QM_SMdpTFrac)
	fillM(Q_SMJets,Q_SMM,Q_SMM_2D,Q_SMdpTFrac)
	fillM(LD_lowDFJets,LD_lowDFM,LD_lowDFM_2D,LD_lowDFdpTFrac)
	fillM(LD_highDFJets,LD_highDFM,LD_highDFM_2D,LD_highDFdpTFrac)
	fillM(LD_SMJets,LD_SMM,LD_SMM_2D,LD_SMdpTFrac)

	fillM(SM_OnlyJets,SM_OnlyM,filler2D,[0]*len(SM_OnlyJets))
	fillM(Mix_Jets,Mix_M,filler2D,Mix_dpTFrac)
	fillM(DarkerMix_Jets,DarkerMix_M,filler2D,DarkerMix_dpTFrac)
	fillM(Dark_Jets,Dark_M,filler2D,Dark_dpTFrac)
	fillM(SoftDark_Jets,SoftDark_M,filler2D,SoftDark_dpTFrac)

	makejsubplot(SM_Jets,SM_aMj,SM_aMn,SM_gir,SM_pD,npSM_Jets,
	SM_axisMajor,SM_axisMinor,SM_momentGirth,SM_ptD,SM_pmult,
	SM_axisMajor_2D,SM_axisMinor_2D,SM_momentGirth_2D,SM_ptD_2D,SM_pmult_2D,[0]*len(SM_Jets))
	makejsubplot(SMM_Jets,SMM_aMj,SMM_aMn,SMM_gir,SMM_pD,npSMM_Jets,
	SMM_axisMajor,SMM_axisMinor,SMM_momentGirth,SMM_ptD,SMM_pmult,
	SMM_axisMajor_2D,SMM_axisMinor_2D,SMM_momentGirth_2D,SMM_ptD_2D,SMM_pmult_2D,[0]*len(SMM_Jets))
	makejsubplot(G_Jets,G_aMj,G_aMn,G_gir,G_pD,npG_Jets,
	G_axisMajor,G_axisMinor,G_momentGirth,G_ptD,G_pmult,
	G_axisMajor_2D,G_axisMinor_2D,G_momentGirth_2D,G_ptD_2D,G_pmult_2D,G_dpTFrac)
	makejsubplot(QM_Jets,QM_aMj,QM_aMn,QM_gir,QM_pD,npQM_Jets,
	QM_axisMajor,QM_axisMinor,QM_momentGirth,QM_ptD,QM_pmult,
	QM_axisMajor_2D,QM_axisMinor_2D,QM_momentGirth_2D,QM_ptD_2D,QM_pmult_2D,QM_dpTFrac)
	makejsubplot(Q_Jets,Q_aMj,Q_aMn,Q_gir,Q_pD,npQ_Jets,
	Q_axisMajor,Q_axisMinor,Q_momentGirth,Q_ptD,Q_pmult,
	Q_axisMajor_2D,Q_axisMinor_2D,Q_momentGirth_2D,Q_ptD_2D,Q_pmult_2D,Q_dpTFrac)
	makejsubplot(QM_QJets,QM_QaMj,QM_QaMn,QM_Qgir,QM_QpD,npQM_QJets,
	QM_QaxisMajor,QM_QaxisMinor,QM_QmomentGirth,QM_QptD,QM_Qpmult,
	QM_QaxisMajor_2D,QM_QaxisMinor_2D,QM_QmomentGirth_2D,QM_QptD_2D,QM_Qpmult_2D,QM_QdpTFrac)
	makejsubplot(QM_GJets,QM_GaMj,QM_GaMn,QM_Ggir,QM_GpD,npQM_GJets,
	QM_GaxisMajor,QM_GaxisMinor,QM_GmomentGirth,QM_GptD,QM_Gpmult,
	QM_GaxisMajor_2D,QM_GaxisMinor_2D,QM_GmomentGirth_2D,QM_GptD_2D,QM_Gpmult_2D,QM_GdpTFrac)
	makejsubplot(Q_GJets,Q_GaMj,Q_GaMn,Q_Ggir,Q_GpD,npQ_GJets,
	Q_GaxisMajor,Q_GaxisMinor,Q_GmomentGirth,Q_GptD,Q_Gpmult,
	Q_GaxisMajor_2D,Q_GaxisMinor_2D,Q_GmomentGirth_2D,Q_GptD_2D,Q_Gpmult_2D,Q_GdpTFrac)
	makejsubplot(G_SMJets,G_SMaMj,G_SMaMn,G_SMgir,G_SMpD,npG_SMJets,
	G_SMaxisMajor,G_SMaxisMinor,G_SMmomentGirth,G_SMptD,G_SMpmult,
	G_SMaxisMajor_2D,G_SMaxisMinor_2D,G_SMmomentGirth_2D,G_SMptD_2D,G_SMpmult_2D,G_SMdpTFrac)
	makejsubplot(QM_SMJets,QM_SMaMj,QM_SMaMn,QM_SMgir,QM_SMpD,npQM_SMJets,
	QM_SMaxisMajor,QM_SMaxisMinor,QM_SMmomentGirth,QM_SMptD,QM_SMpmult,
	QM_SMaxisMajor_2D,QM_SMaxisMinor_2D,QM_SMmomentGirth_2D,QM_SMptD_2D,QM_SMpmult_2D,QM_SMdpTFrac)
	makejsubplot(Q_SMJets,Q_SMaMj,Q_SMaMn,Q_SMgir,Q_SMpD,npQ_SMJets,
	Q_SMaxisMajor,Q_SMaxisMinor,Q_SMmomentGirth,Q_SMptD,Q_SMpmult,
	Q_SMaxisMajor_2D,Q_SMaxisMinor_2D,Q_SMmomentGirth_2D,Q_SMptD_2D,Q_SMpmult_2D,Q_SMdpTFrac)
	makejsubplot(LD_lowDFJets,LD_lowDFaMj,LD_lowDFaMn,LD_lowDFgir,LD_lowDFpD,npLD_lowDFJets,
	LD_lowDFaxisMajor,LD_lowDFaxisMinor,LD_lowDFmomentGirth,LD_lowDFptD,LD_lowDFpmult,
	LD_lowDFaxisMajor_2D,LD_lowDFaxisMinor_2D,LD_lowDFmomentGirth_2D,LD_lowDFptD_2D,LD_lowDFpmult_2D,LD_lowDFdpTFrac)
	makejsubplot(LD_highDFJets,LD_highDFaMj,LD_highDFaMn,LD_highDFgir,LD_highDFpD,npLD_highDFJets,
	LD_highDFaxisMajor,LD_highDFaxisMinor,LD_highDFmomentGirth,LD_highDFptD,LD_highDFpmult,
	LD_highDFaxisMajor_2D,LD_highDFaxisMinor_2D,LD_highDFmomentGirth_2D,LD_highDFptD_2D,LD_highDFpmult_2D,LD_highDFdpTFrac)
	makejsubplot(LD_SMJets,LD_SMaMj,LD_SMaMn,LD_SMgir,LD_SMpD,npLD_SMJets,
	LD_SMaxisMajor,LD_SMaxisMinor,LD_SMmomentGirth,LD_SMptD,LD_SMpmult,
	LD_SMaxisMajor_2D,LD_SMaxisMinor_2D,LD_SMmomentGirth_2D,LD_SMptD_2D,LD_SMpmult_2D,LD_SMdpTFrac)

	makejsubplot(SM_OnlyJets,SM_OnlyaMj,SM_OnlyaMn,SM_Onlygir,SM_OnlypD,npSM_OnlyJets,
	SM_OnlyaxisMajor,SM_OnlyaxisMinor,SM_OnlymomentGirth,SM_OnlyptD,SM_Onlypmult,
	filler2D,filler2D,filler2D,filler2D,filler2D,[0]*len(SM_OnlyJets))
	makejsubplot(Mix_Jets,Mix_aMj,Mix_aMn,Mix_gir,Mix_pD,npMix_Jets,
	Mix_axisMajor,Mix_axisMinor,Mix_momentGirth,Mix_ptD,Mix_pmult,
	filler2D,filler2D,filler2D,filler2D,filler2D,Mix_dpTFrac)
	makejsubplot(DarkerMix_Jets,DarkerMix_aMj,DarkerMix_aMn,DarkerMix_gir,DarkerMix_pD,npDarkerMix_Jets,
	DarkerMix_axisMajor,DarkerMix_axisMinor,DarkerMix_momentGirth,DarkerMix_ptD,DarkerMix_pmult,
	filler2D,filler2D,filler2D,filler2D,filler2D,DarkerMix_dpTFrac)
	makejsubplot(Dark_Jets,Dark_aMj,Dark_aMn,Dark_gir,Dark_pD,npDark_Jets,
	Dark_axisMajor,Dark_axisMinor,Dark_momentGirth,Dark_ptD,Dark_pmult,
	filler2D,filler2D,filler2D,filler2D,filler2D,Dark_dpTFrac)
	makejsubplot(SoftDark_Jets,SoftDark_aMj,SoftDark_aMn,SoftDark_gir,SoftDark_pD,npSoftDark_Jets,
	SoftDark_axisMajor,SoftDark_axisMinor,SoftDark_momentGirth,SoftDark_ptD,SoftDark_pmult,
	filler2D,filler2D,filler2D,filler2D,filler2D,SoftDark_dpTFrac)

	## tau plots
	tauplot(SM_Jets,SM_Tau1,SM_Tau2,SM_Tau3,SM_tau1,
	SM_tau2,SM_tau3,SM_tau21,SM_tau32,SM_tau1_2D,
	SM_tau2_2D,SM_tau3_2D,SM_tau21_2D,SM_tau32_2D,[0]*len(SM_Jets))
	tauplot(SMM_Jets,SMM_Tau1,SMM_Tau2,SMM_Tau3,SMM_tau1,
	SMM_tau2,SMM_tau3,SMM_tau21,SMM_tau32,SMM_tau1_2D,
	SMM_tau2_2D,SMM_tau3_2D,SMM_tau21_2D,SMM_tau32_2D,[0]*len(SMM_Jets))
	tauplot(G_Jets,G_Tau1,G_Tau2,G_Tau3,G_tau1,
	G_tau2,G_tau3,G_tau21,G_tau32,G_tau1_2D,
	G_tau2_2D,G_tau3_2D,G_tau21_2D,G_tau32_2D,G_dpTFrac)
	tauplot(QM_Jets,QM_Tau1,QM_Tau2,QM_Tau3,QM_tau1,
	QM_tau2,QM_tau3,QM_tau21,QM_tau32,QM_tau1_2D,
	QM_tau2_2D,QM_tau3_2D,QM_tau21_2D,QM_tau32_2D,QM_dpTFrac)
	tauplot(Q_Jets,Q_Tau1,Q_Tau2,Q_Tau3,Q_tau1,
	Q_tau2,Q_tau3,Q_tau21,Q_tau32,Q_tau1_2D,
	Q_tau2_2D,Q_tau3_2D,Q_tau21_2D,Q_tau32_2D,Q_dpTFrac)
	tauplot(QM_QJets,QM_QTau1,QM_QTau2,QM_QTau3,QM_Qtau1,
	QM_Qtau2,QM_Qtau3,QM_Qtau21,QM_Qtau32,QM_Qtau1_2D,
	QM_Qtau2_2D,QM_Qtau3_2D,QM_Qtau21_2D,QM_Qtau32_2D,QM_QdpTFrac)
	tauplot(QM_GJets,QM_GTau1,QM_GTau2,QM_GTau3,QM_Gtau1,
	QM_Gtau2,QM_Gtau3,QM_Gtau21,QM_Gtau32,QM_Gtau1_2D,
	QM_Gtau2_2D,QM_Gtau3_2D,QM_Gtau21_2D,QM_Gtau32_2D,QM_GdpTFrac)
	tauplot(Q_GJets,Q_GTau1,Q_GTau2,Q_GTau3,Q_Gtau1,
	Q_Gtau2,Q_Gtau3,Q_Gtau21,Q_Gtau32,Q_Gtau1_2D,
	Q_Gtau2_2D,Q_Gtau3_2D,Q_Gtau21_2D,Q_Gtau32_2D,Q_GdpTFrac)
	tauplot(G_SMJets,G_SMTau1,G_SMTau2,G_SMTau3,G_SMtau1,
	G_SMtau2,G_SMtau3,G_SMtau21,G_SMtau32,G_SMtau1_2D,
	G_SMtau2_2D,G_SMtau3_2D,G_SMtau21_2D,G_SMtau32_2D,G_SMdpTFrac)
	tauplot(QM_SMJets,QM_SMTau1,QM_SMTau2,QM_SMTau3,QM_SMtau1,
	QM_SMtau2,QM_SMtau3,QM_SMtau21,QM_SMtau32,QM_SMtau1_2D,
	QM_SMtau2_2D,QM_SMtau3_2D,QM_SMtau21_2D,QM_SMtau32_2D,QM_SMdpTFrac)
	tauplot(Q_SMJets,Q_SMTau1,Q_SMTau2,Q_SMTau3,Q_SMtau1,
	Q_SMtau2,Q_SMtau3,Q_SMtau21,Q_SMtau32,Q_SMtau1_2D,
	Q_SMtau2_2D,Q_SMtau3_2D,Q_SMtau21_2D,Q_SMtau32_2D,Q_SMdpTFrac)
	tauplot(LD_lowDFJets,LD_lowDFTau1,LD_lowDFTau2,LD_lowDFTau3,LD_lowDFtau1,
	LD_lowDFtau2,LD_lowDFtau3,LD_lowDFtau21,LD_lowDFtau32,LD_lowDFtau1_2D,
	LD_lowDFtau2_2D,LD_lowDFtau3_2D,LD_lowDFtau21_2D,LD_lowDFtau32_2D,LD_lowDFdpTFrac)
	tauplot(LD_highDFJets,LD_highDFTau1,LD_highDFTau2,LD_highDFTau3,LD_highDFtau1,
	LD_highDFtau2,LD_highDFtau3,LD_highDFtau21,LD_highDFtau32,LD_highDFtau1_2D,
	LD_highDFtau2_2D,LD_highDFtau3_2D,LD_highDFtau21_2D,LD_highDFtau32_2D,LD_highDFdpTFrac)
	tauplot(LD_SMJets,LD_SMTau1,LD_SMTau2,LD_SMTau3,LD_SMtau1,
	LD_SMtau2,LD_SMtau3,LD_SMtau21,LD_SMtau32,LD_SMtau1_2D,
	LD_SMtau2_2D,LD_SMtau3_2D,LD_SMtau21_2D,LD_SMtau32_2D,LD_SMdpTFrac)

	tauplot(SM_OnlyJets,SM_OnlyTau1,SM_OnlyTau2,SM_OnlyTau3,SM_Onlytau1,
	SM_Onlytau2,SM_Onlytau3,SM_Onlytau21,SM_Onlytau32,filler2D,
	filler2D,filler2D,filler2D,filler2D,[0]*len(SM_OnlyJets))
	tauplot(Mix_Jets,Mix_Tau1,Mix_Tau2,Mix_Tau3,Mix_tau1,
	Mix_tau2,Mix_tau3,Mix_tau21,Mix_tau32,filler2D,
	filler2D,filler2D,filler2D,filler2D,Mix_dpTFrac)
	tauplot(DarkerMix_Jets,DarkerMix_Tau1,DarkerMix_Tau2,DarkerMix_Tau3,DarkerMix_tau1,
	DarkerMix_tau2,DarkerMix_tau3,DarkerMix_tau21,DarkerMix_tau32,filler2D,
	filler2D,filler2D,filler2D,filler2D,DarkerMix_dpTFrac)
	tauplot(Dark_Jets,Dark_Tau1,Dark_Tau2,Dark_Tau3,Dark_tau1,
	Dark_tau2,Dark_tau3,Dark_tau21,Dark_tau32,filler2D,
	filler2D,filler2D,filler2D,filler2D,Dark_dpTFrac)
	tauplot(SoftDark_Jets,SoftDark_Tau1,SoftDark_Tau2,SoftDark_Tau3,SoftDark_tau1,
	SoftDark_tau2,SoftDark_tau3,SoftDark_tau21,SoftDark_tau32,filler2D,
	filler2D,filler2D,filler2D,filler2D,SoftDark_dpTFrac)

	## SD Mass plots
	SDFill(SM_SDJets,SM_SDM,SM_SDM_2D,[0]*len(SM_SDJets))
	SDFill(SMM_SDJets,SMM_SDM,SMM_SDM_2D,[0]*len(SMM_SDJets))
	SDFill(G_SDJets,G_SDM,G_SDM_2D,G_dpTFrac)
	SDFill(QM_SDJets,QM_SDM,QM_SDM_2D,QM_dpTFrac)
	SDFill(Q_SDJets,Q_SDM,Q_SDM_2D,Q_dpTFrac)
	SDFill(QM_QSDJets,QM_QSDM,QM_QSDM_2D,QM_QdpTFrac)
	SDFill(QM_GSDJets,QM_GSDM,QM_GSDM_2D,QM_GdpTFrac)
	SDFill(Q_GSDJets,Q_GSDM,Q_GSDM_2D,Q_GdpTFrac)
	SDFill(G_SMSDJets,G_SMSDM,G_SMSDM_2D,G_SMdpTFrac)
	SDFill(QM_SMSDJets,QM_SMSDM,QM_SMSDM_2D,QM_SMdpTFrac)
	SDFill(Q_SMSDJets,Q_SMSDM,Q_SMSDM_2D,Q_SMdpTFrac)
	SDFill(LD_lowDFSDJets,LD_lowDFSDM,LD_lowDFSDM_2D,LD_lowDFdpTFrac)
	SDFill(LD_highDFSDJets,LD_highDFSDM,LD_highDFSDM_2D,LD_highDFdpTFrac)
	SDFill(LD_SMSDJets,LD_SMSDM,LD_SMSDM_2D,LD_SMdpTFrac)

	SDFill(SM_OnlySDJets,SM_OnlySDM,filler2D,[0]*len(SM_OnlySDJets))
	SDFill(Mix_SDJets,Mix_SDM,filler2D,Mix_dpTFrac)
	SDFill(DarkerMix_SDJets,DarkerMix_SDM,filler2D,DarkerMix_dpTFrac)
	SDFill(Dark_SDJets,Dark_SDM,filler2D,Dark_dpTFrac)
	SDFill(SoftDark_SDJets,SoftDark_SDM,filler2D,SoftDark_dpTFrac)

	dSfillJetVar(SM_Jets,SM_mult,SM_Pt,SM_Eta,SM_Phi,SM_DPhi,METPhiv,
	SM_Pt_2D,SM_Eta_2D,SM_Phi_2D,SM_DPhi_2D,[0]*len(SM_Jets))
	dSfillJetVar(SMM_Jets,SMM_mult,SMM_Pt,SMM_Eta,SMM_Phi,SMM_DPhi,METPhiv,
	SMM_Pt_2D,SMM_Eta_2D,SMM_Phi_2D,SMM_DPhi_2D,[0]*len(SMM_Jets))
	dSfillJetVar(G_Jets,G_mult,G_Pt,G_Eta,G_Phi,G_DPhi,METPhiv,
	G_Pt_2D,G_Eta_2D,G_Phi_2D,G_DPhi_2D,G_dpTFrac)
	dSfillJetVar(QM_Jets,QM_mult,QM_Pt,QM_Eta,QM_Phi,QM_DPhi,METPhiv,
	QM_Pt_2D,QM_Eta_2D,QM_Phi_2D,QM_DPhi_2D,QM_dpTFrac)
	dSfillJetVar(Q_Jets,Q_mult,Q_Pt,Q_Eta,Q_Phi,Q_DPhi,METPhiv,
	Q_Pt_2D,Q_Eta_2D,Q_Phi_2D,Q_DPhi_2D,Q_dpTFrac)
	dSfillJetVar(QM_QJets,QM_Qmult,QM_QPt,QM_QEta,QM_QPhi,QM_QDPhi,METPhiv,
	QM_QPt_2D,QM_QEta_2D,QM_QPhi_2D,QM_QDPhi_2D,QM_QdpTFrac)
	dSfillJetVar(QM_GJets,QM_Gmult,QM_GPt,QM_GEta,QM_GPhi,QM_GDPhi,METPhiv,
	QM_GPt_2D,QM_GEta_2D,QM_GPhi_2D,QM_GDPhi_2D,QM_GdpTFrac)
	dSfillJetVar(Q_GJets,Q_Gmult,Q_GPt,Q_GEta,Q_GPhi,Q_GDPhi,METPhiv,
	Q_GPt_2D,Q_GEta_2D,Q_GPhi_2D,Q_GDPhi_2D,Q_GdpTFrac)
	dSfillJetVar(G_SMJets,G_SMmult,G_SMPt,G_SMEta,G_SMPhi,G_SMDPhi,METPhiv,
	G_SMPt_2D,G_SMEta_2D,G_SMPhi_2D,G_SMDPhi_2D,G_SMdpTFrac)
	dSfillJetVar(QM_SMJets,QM_SMmult,QM_SMPt,QM_SMEta,QM_SMPhi,QM_SMDPhi,METPhiv,
	QM_SMPt_2D,QM_SMEta_2D,QM_SMPhi_2D,QM_SMDPhi_2D,QM_SMdpTFrac)
	dSfillJetVar(Q_SMJets,Q_SMmult,Q_SMPt,Q_SMEta,Q_SMPhi,Q_SMDPhi,METPhiv,
	Q_SMPt_2D,Q_SMEta_2D,Q_SMPhi_2D,Q_SMDPhi_2D,Q_SMdpTFrac)
	dSfillJetVar(LD_lowDFJets,LD_lowDFmult,LD_lowDFPt,LD_lowDFEta,LD_lowDFPhi,LD_lowDFDPhi,METPhiv,
	LD_lowDFPt_2D,LD_lowDFEta_2D,LD_lowDFPhi_2D,LD_lowDFDPhi_2D,LD_lowDFdpTFrac)
	dSfillJetVar(LD_highDFJets,LD_highDFmult,LD_highDFPt,LD_highDFEta,LD_highDFPhi,LD_highDFDPhi,METPhiv,
	LD_highDFPt_2D,LD_highDFEta_2D,LD_highDFPhi_2D,LD_highDFDPhi_2D,LD_highDFdpTFrac)
	dSfillJetVar(LD_SMJets,LD_SMmult,LD_SMPt,LD_SMEta,LD_SMPhi,LD_SMDPhi,METPhiv,
	LD_SMPt_2D,LD_SMEta_2D,LD_SMPhi_2D,LD_SMDPhi_2D,LD_SMdpTFrac)

	dSfillJetVar_no2D(SM_OnlyJets,SM_Onlymult,SM_OnlyPt,
	SM_OnlyEta,SM_OnlyPhi,SM_OnlyDPhi,METPhiv,SM_OnlydpTFrac)
	dSfillJetVar_no2D(Mix_Jets,Mix_mult,Mix_Pt,
	Mix_Eta,Mix_Phi,Mix_DPhi,METPhiv,Mix_dpTFrac)
	dSfillJetVar_no2D(DarkerMix_Jets,DarkerMix_mult,DarkerMix_Pt,
	DarkerMix_Eta,DarkerMix_Phi,DarkerMix_DPhi,METPhiv,DarkerMix_dpTFrac)
	dSfillJetVar_no2D(Dark_Jets,Dark_mult,Dark_Pt,
	Dark_Eta,Dark_Phi,Dark_DPhi,METPhiv,Dark_dpTFrac)
	dSfillJetVar_no2D(SoftDark_Jets,SoftDark_mult,SoftDark_Pt,
	SoftDark_Eta,SoftDark_Phi,SoftDark_DPhi,METPhiv,SoftDark_dpTFrac)

	fillSDPt(SM_SDJets,SM_SDPt)
	fillSDPt(SMM_SDJets,SMM_SDPt)
	fillSDPt(G_SDJets,G_SDPt)
	fillSDPt(QM_SDJets,QM_SDPt)
	fillSDPt(Q_SDJets,Q_SDPt)
	fillSDPt(QM_QSDJets,QM_QSDPt)
	fillSDPt(QM_GSDJets,QM_GSDPt)
	fillSDPt(Q_GSDJets,Q_GSDPt)
	fillSDPt(G_SMSDJets,G_SMSDPt)
	fillSDPt(QM_SMSDJets,QM_SMSDPt)
	fillSDPt(Q_SMSDJets,Q_SMSDPt)
	fillSDPt(LD_lowDFSDJets,LD_lowDFSDPt)
	fillSDPt(LD_highDFSDJets,LD_highDFSDPt)
	fillSDPt(LD_SMSDJets,LD_SMSDPt)

	fillSDPt(SM_OnlySDJets,SM_OnlySDPt)
	fillSDPt(Mix_SDJets,Mix_SDPt)
	fillSDPt(DarkerMix_SDJets,DarkerMix_SDPt)
	fillSDPt(Dark_SDJets,Dark_SDPt)
	fillSDPt(SoftDark_SDJets,SoftDark_SDPt)

	# Calculating the ST variable
	if len(jets) >= 1:
		ST.Fill(ST_Val(METv,jets))
		MET.Fill(METv)
	if len(jets) == 2:
		ST_2jet.Fill(ST_Val(METv,jets))
	if len(jets) == 3:
		ST_3jet.Fill(ST_Val(METv,jets))
	if len(jets) == 4:
		ST_4jet.Fill(ST_Val(METv,jets))
	if len(jets) > 4:
		ST_4plusjet.Fill(ST_Val(METv,jets))

	if len(jets) >= 2:
		MTv = trans_mass_Njet(jets[0],jets[1],METv,METPhiv)
		RT.Fill(METv/MTv[0])
		if len(jets) == 2 and jets[0].Pt() > 200 and jets[1].Pt() > 200:
			MT_2j.Fill(MTv[0])
		if len(jets) == 3 and jets[0].Pt() > 200 and jets[1].Pt() > 200:
			MT_3j.Fill(MTv[0])
		if len(jets) == 4 and jets[0].Pt() > 200 and jets[1].Pt() > 200:
			MT_4j.Fill(MTv[0])
		if len(jets) > 4 and jets[0].Pt() > 200 and jets[1].Pt() > 200:
			MT_4pj.Fill(MTv[0])

		if MTv[0] >= 2000:
			MT_m2000_DR_f2.Fill(MTv[1])
			MT_m2000_DP_f2.Fill(MTv[2])
			MT_m2000_DE_f2.Fill(MTv[3])
		else:
			MT_l2000_DR_f2.Fill(MTv[1])
			MT_l2000_DP_f2.Fill(MTv[2])
			MT_l2000_DE_f2.Fill(MTv[3])

	if len(MTPair) > 0:
		MTvR = trans_mass_Njet(MTPair[0],MTPair[1],METv,METPhiv)
		MT_R.Fill(MTvR[0])
		MT_R_DR_f2.Fill(MTvR[1])
		MT_R_DP_f2.Fill(MTvR[2])
		MT_R_DE_f2.Fill(MTvR[3])

		jlist = []
		for jla in jets:
			jlist.append(jla)

		rind0 = jlist.index(MTPair[0])
		rind1 = jlist.index(MTPair[1])
		rComb = tuple(sorted([rind0,rind1])) # tuple of indices of the right combination

		com_obj = itertools.combinations(range(len(jlist)),2) # getting all combinations of 2 indices for calculating mT
		allComb = list(com_obj)

		wComb = allComb[:]
		wComb.remove(rComb)

		for wind in wComb:
			MTvW = trans_mass_Njet(jlist[wind[0]],jlist[wind[1]],METv,METPhiv)
			MT_W.Fill(MTvW[0])
			MT_W_DR_f2.Fill(MTvW[1])
			MT_W_DP_f2.Fill(MTvW[2])
			MT_W_DE_f2.Fill(MTvW[3])

	# calculating mT2 most probable: for different number of mediators in an event
	if len(jets) >= 4:
		## get all the combinations of 2 pairs of jets from all the jets in an event
		List4jets_3Com = allJetCom()
		mpCJ = MT_mpCalc(List4jets_3Com)[1]
		MT2_mp.Fill(MT_mpCalc(List4jets_3Com)[0][0])
		MT2_mp_Events.append(count)
		if deltaPhi(mpCJ[0].Phi(),mpCJ[1].Phi()) > 1.5 and deltaPhi(mpCJ[2].Phi(),mpCJ[3].Phi()) > 1.5:
			MT2_mp_DPhiPair.Fill(MT_mpCalc(List4jets_3Com)[0][0])
			MT2_mp_DPhiPair_Events.append(count)
			if (mpCJ[0].Pt() > 500 and mpCJ[1].Pt() > 500 and mpCJ[2].Pt() > 500 and mpCJ[3].Pt() > 500):
				MT2_mp_DPhiPair_Pt.Fill(MT_mpCalc(List4jets_3Com)[0][0])
				MT2_mp_DPhiPair_Pt_Events.append(count)
		if ( (mpCJ[0].Phi() > 2.0 and mpCJ[1].Phi() < 0.5) or (mpCJ[0].Phi() < 0.5 and mpCJ[1].Phi() > 2.0) and
		(mpCJ[2].Phi() > 2.0 and mpCJ[3].Phi() < 0.5) or (mpCJ[2].Phi() < 0.5 and mpCJ[3].Phi() > 2.0) ):
			MT2_mp_PhiReq.Fill(MT_mpCalc(List4jets_3Com)[0][0])
			MT2_mp_PhiReq_Events.append(count)

		## get all the combinations with pT > 500
		List4jets_3Com = allJetCom(1)
		if len(List4jets_3Com) > 0:
			mpCJ = MT_mpCalc(List4jets_3Com)[1]
			MT2_Pt_mp.Fill(MT_mpCalc(List4jets_3Com)[0][0])
			MT2_Pt_mp_Events.append(count)
			if deltaPhi(mpCJ[0].Phi(),mpCJ[1].Phi()) > 1.5 and deltaPhi(mpCJ[2].Phi(),mpCJ[3].Phi()) > 1.5:
				MT2_Pt_mp_DPhiPair.Fill(MT_mpCalc(List4jets_3Com)[0][0])
				MT2_Pt_mp_DPhiPair_Events.append(count)

		## only use first 4 leading jets for MT2_mp calculation
		List4jets_3Com = [[0,1,2,3],[0,2,1,3],[0,3,1,2]]
		mpCJ = MT_mpCalc(List4jets_3Com)[1]
		MT2_mp_f4.Fill(MT_mpCalc(List4jets_3Com)[0][0])
		MT2_mp_f4_Events.append(count)
		if deltaPhi(mpCJ[0].Phi(),mpCJ[1].Phi()) > 1.5 and deltaPhi(mpCJ[2].Phi(),mpCJ[3].Phi()) > 1.5:
			MT2_mp_DPhiPair_f4.Fill(MT_mpCalc(List4jets_3Com)[0][0])
			MT2_mp_DPhiPair_f4_Events.append(count)
			if (mpCJ[0].Pt() > 500 and mpCJ[1].Pt() > 500 and mpCJ[2].Pt() > 500 and mpCJ[3].Pt() > 500):
				MT2_mp_DPhiPair_Pt_f4.Fill(MT_mpCalc(List4jets_3Com)[0][0])
				MT2_mp_DPhiPair_Pt_f4_Events.append(count)
		if ( (mpCJ[0].Phi() > 2.0 and mpCJ[1].Phi() < 0.5) or (mpCJ[0].Phi() < 0.5 and mpCJ[1].Phi() > 2.0) and
		(mpCJ[2].Phi() > 2.0 and mpCJ[3].Phi() < 0.5) or (mpCJ[2].Phi() < 0.5 and mpCJ[3].Phi() > 2.0) ):
			MT2_mp_PhiReq_f4.Fill(MT_mpCalc(List4jets_3Com)[0][0])
			MT2_mp_PhiReq_f4_Events.append(count)

		MT2_mp_MDiff.Fill(abs(M_2J(mpCJ[0],mpCJ[1])-M_2J(mpCJ[2],mpCJ[3])))
		MT2_mp_DPhiMET.Fill(deltaPhi(SumJet(mpCJ[0],mpCJ[1]).Phi(),METPhiv))
		MT2_mp_DPhiMET.Fill(deltaPhi(SumJet(mpCJ[2],mpCJ[3]).Phi(),METPhiv))
		MT2_mp_Pt.Fill(mpCJ[0].Pt())
		MT2_mp_Pt.Fill(mpCJ[1].Pt())
		MT2_mp_Pt.Fill(mpCJ[2].Pt())
		MT2_mp_Pt.Fill(mpCJ[3].Pt())
		if nfdMPart == 0:
			MT2_mp_0M.Fill(MT_mpCalc(List4jets_3Com)[0][0])
		if nfdMPart == 1:
			MT2_mp_1M.Fill(MT_mpCalc(List4jets_3Com)[0][0])
		if nfdMPart == 2:
			MT2_mp_2M.Fill(MT_mpCalc(List4jets_3Com)[0][0])

		if len(MedPair1) > 0 and len(MedPair2) > 0:
			MT2_R_Events.append(count)
		# right permutations
			d0r = MedPair1[0]
			s0r = MedPair1[1]
			d1r = MedPair2[0]
			s1r = MedPair2[1]

			# test how often the right combination comes from the first 4 leading jets
			isMT2R += 1
			MT2_RJets = [9,9,9,9]
			for i in range(len(jets)):
				if d0r == jets[i]:
					MT2_RJets[0] = i
				elif s0r == jets[i]:
					MT2_RJets[1] = i
				elif d1r == jets[i]:
					MT2_RJets[2] = i
				elif s1r == jets[i]:
					MT2_RJets[3] = i

			comval = np.intersect1d([0,1,2,3], MT2_RJets)
			if len(comval) == 4:
				is4lead += 1

			MT2_R_MDiff.Fill(abs(M_2J(d0r,s0r)-M_2J(d1r,s1r)))
			MT2_R_M.Fill(M_2J(d0r,s0r))
			MT2_R_M.Fill(M_2J(d1r,s1r))
			MT2_R_Pt.Fill(d0r.Pt())
			MT2_R_Pt.Fill(s0r.Pt())
			MT2_R_Pt.Fill(d1r.Pt())
			MT2_R_Pt.Fill(s1r.Pt())
			MT2_R_DPhiMET.Fill(deltaPhi(SumJet(d0r,s0r).Phi(),METPhiv))
			MT2_R_DPhiMET.Fill(deltaPhi(SumJet(d1r,s1r).Phi(),METPhiv))
			MT2_R_DPhiMET_d0r.Fill(deltaPhi(d0r.Phi(),METPhiv))
			MT2_R_DPhiMET_s0r.Fill(deltaPhi(s0r.Phi(),METPhiv))
			MT2_R_DPhiMET_d1r.Fill(deltaPhi(d1r.Phi(),METPhiv))
			MT2_R_DPhiMET_s1r.Fill(deltaPhi(s1r.Phi(),METPhiv))
			MT2_R_DPhi_pair.Fill(deltaPhi(d0r.Phi(),s0r.Phi()))
			MT2_R_DPhi_pair.Fill(deltaPhi(d1r.Phi(),s1r.Phi()))

			MT2DRRight = MT2DRCal(d0r,s0r,d1r,s1r,met4p)
			MT2_R.Fill(MT2DRRight[0])
			MT2_DeltaR_R.Fill(MT2DRRight[1])
			MT2_DeltaP_R.Fill(MT2DRRight[2])
			MT2_DeltaE_R.Fill(MT2DRRight[3])
			# wrong permutations
			allMT2Wrong = allWrongCom(List4jets_3Com,MT2_RJets)
			MT2_W_MDiff.Fill(abs(M_2J(d0r,s1r)-M_2J(d1r,s0r)))
			MT2_W_MDiff.Fill(abs(M_2J(d0r,d1r)-M_2J(s0r,s1r)))
			MT2_W_M.Fill(M_2J(d0r,s1r))
			MT2_W_M.Fill(M_2J(d1r,s0r))
			MT2_W_M.Fill(M_2J(d0r,d1r))
			MT2_W_M.Fill(M_2J(s0r,s1r))
			MT2DRWrong1 = MT2DRCal(d0r,s1r,d1r,s0r,met4p)
			MT2_W.Fill(MT2DRWrong1[0])
			MT2_DeltaR_W.Fill(MT2DRWrong1[1])
			MT2_DeltaP_W.Fill(MT2DRWrong1[2])
			MT2_DeltaE_W.Fill(MT2DRWrong1[3])
			MT2DRWrong2 = MT2DRCal(d0r,d1r,s0r,s1r,met4p)
			MT2_W.Fill(MT2DRWrong2[0])
			MT2_DeltaR_W.Fill(MT2DRWrong2[1])
			MT2_DeltaP_W.Fill(MT2DRWrong2[2])
			MT2_DeltaE_W.Fill(MT2DRWrong2[3])

			MT2_mp_MT2_R.Fill(MT_mpCalc(List4jets_3Com)[0][0],MT2DRRight[0])

	# delta phi for different jet categories
		fill_DPhi_PP(SM_Jets,SM_DPhi_PP,METPhiv)
		fill_DPhi_PP(SMM_Jets,SMM_DPhi_PP,METPhiv)
		fill_DPhi_PP(G_Jets,G_DPhi_PP,METPhiv)
		fill_DPhi_PP(QM_Jets,QM_DPhi_PP,METPhiv,QM_Dhi_JJ_PP,QM_Dhi_JMET_Near,QM_Dhi_JMET_Far)
		fill_DPhi_PP(Q_Jets,Q_DPhi_PP,METPhiv)
		fill_DPhi_PP(QM_QJets,QM_QDPhi_PP,METPhiv)
		fill_DPhi_PP(QM_GJets,QM_GDPhi_PP,METPhiv)
		fill_DPhi_PP(Q_GJets,Q_GDPhi_PP,METPhiv)
		fill_DPhi_PP(G_SMJets,G_SMDPhi_PP,METPhiv)
		fill_DPhi_PP(QM_SMJets,QM_SMDPhi_PP,METPhiv)
		fill_DPhi_PP(Q_SMJets,Q_SMDPhi_PP,METPhiv)
		fill_DPhi_PP(LD_lowDFJets,LD_lowDFDPhi_PP,METPhiv)
		fill_DPhi_PP(LD_highDFJets,LD_highDFDPhi_PP,METPhiv)
		fill_DPhi_PP(LD_SMJets,LD_SMDPhi_PP,METPhiv)
	# if count == 1000: # which event number to stop the code
	# 	break

print "Is from 4 Lead/Number of MT2_R Event: " + str(is4lead/float(isMT2R))

remainMT2R_mp_f4 = np.intersect1d(MT2_mp_f4_Events, MT2_R_Events)
remainMT2R_DPhiPair_f4 = np.intersect1d(MT2_mp_DPhiPair_f4_Events, MT2_R_Events)
remainMT2R_DPhiPair_Pt_f4 = np.intersect1d(MT2_mp_DPhiPair_Pt_f4_Events, MT2_R_Events)
remainMT2R_PhiReq_f4 = np.intersect1d(MT2_mp_PhiReq_f4_Events, MT2_R_Events)

remainMT2R_mp = np.intersect1d(MT2_mp_Events, MT2_R_Events)
remainMT2R_DPhiPair = np.intersect1d(MT2_mp_DPhiPair_Events, MT2_R_Events)
remainMT2R_DPhiPair_Pt = np.intersect1d(MT2_mp_DPhiPair_Pt_Events, MT2_R_Events)
remainMT2R_PhiReq = np.intersect1d(MT2_mp_PhiReq_Events, MT2_R_Events)

remainMT2_Pt_mp = np.intersect1d(MT2_Pt_mp_Events, MT2_R_Events)
remainMT2_Pt_mp_DPhiPair = np.intersect1d(MT2_Pt_mp_DPhiPair_Events, MT2_R_Events)

MT2_R_mp_f4_sigEff = len(remainMT2R_mp_f4)/float(len(MT2_R_Events))
MT2_R_DPhiPair_f4_sigEff = len(remainMT2R_DPhiPair_f4)/float(len(MT2_R_Events))
MT2_R_DPhiPair_Pt_f4_sigEff = len(remainMT2R_DPhiPair_Pt_f4)/float(len(MT2_R_Events))
MT2_R_PhiReq_f4_sigEff = len(remainMT2R_PhiReq_f4)/float(len(MT2_R_Events))

MT2_R_mp_sigEff = len(remainMT2R_mp)/float(len(MT2_R_Events))
MT2_R_DPhiPair_sigEff = len(remainMT2R_DPhiPair)/float(len(MT2_R_Events))
MT2_R_DPhiPair_Pt_sigEff = len(remainMT2R_DPhiPair_Pt)/float(len(MT2_R_Events))
MT2_R_PhiReq_sigEff = len(remainMT2R_PhiReq)/float(len(MT2_R_Events))

MT2_R_Pt_mp_sigEff = len(remainMT2_Pt_mp)/float(len(MT2_R_Events))
MT2_R_Pt_mp_DPhiPair_sigEff = len(remainMT2_Pt_mp_DPhiPair)/float(len(MT2_R_Events))

MT2_R_mp_f4_SB = len(remainMT2R_mp_f4)/float(len(MT2_mp_f4_Events))
MT2_R_DPhiPair_f4_SB = len(remainMT2R_DPhiPair_f4)/float(len(MT2_mp_DPhiPair_f4_Events))
MT2_R_DPhiPair_Pt_f4_SB = len(remainMT2R_DPhiPair_Pt_f4)/float(len(MT2_mp_DPhiPair_Pt_f4_Events))
MT2_R_PhiReq_f4_SB = len(remainMT2R_PhiReq_f4)/float(len(MT2_mp_PhiReq_f4_Events))

MT2_R_mp_SB = len(remainMT2R_mp)/float(len(MT2_mp_Events))
MT2_R_DPhiPair_SB = len(remainMT2R_DPhiPair)/float(len(MT2_mp_DPhiPair_Events))
MT2_R_DPhiPair_Pt_SB = len(remainMT2R_DPhiPair_Pt)/float(len(MT2_mp_DPhiPair_Pt_Events))
MT2_R_PhiReq_SB = len(remainMT2R_PhiReq)/float(len(MT2_mp_PhiReq_Events))

MT2_R_Pt_mp_SB = len(remainMT2_Pt_mp)/float(len(MT2_Pt_mp_Events))
MT2_R_Pt_mp_DPhiPair_SB = len(remainMT2_Pt_mp_DPhiPair)/float(len(MT2_Pt_mp_DPhiPair_Events))

print "MT2_R_mp_f4_sigEff: " + str(MT2_R_mp_f4_sigEff)
print "MT2_R_DPhiPair_f4_sigEff: " + str(MT2_R_DPhiPair_f4_sigEff)
print "MT2_R_DPhiPair_Pt_f4_sigEff: " + str(MT2_R_DPhiPair_Pt_f4_sigEff)
print "MT2_R_PhiReq_f4_sigEff: " + str(MT2_R_PhiReq_f4_sigEff)

print "MT2_R_mp_sigEff: " + str(MT2_R_mp_sigEff)
print "MT2_R_DPhiPair_sigEff: " + str(MT2_R_DPhiPair_sigEff)
print "MT2_R_DPhiPair_Pt_sigEff: " + str(MT2_R_DPhiPair_Pt_sigEff)
print "MT2_R_PhiReq_sigEff: " + str(MT2_R_PhiReq_sigEff)

print "MT2_R_Pt_mp_sigEff: " + str(MT2_R_Pt_mp_sigEff)
print "MT2_R_Pt_mp_DPhiPair_sigEff: " + str(MT2_R_Pt_mp_DPhiPair_sigEff)

print "MT2_R_mp_f4_SB: " + str(MT2_R_mp_f4_SB)
print "MT2_R_DPhiPair_f4_SB: " + str(MT2_R_DPhiPair_f4_SB)
print "MT2_R_DPhiPair_Pt_f4_SB: " + str(MT2_R_DPhiPair_Pt_f4_SB)
print "MT2_R_PhiReq_f4_SB: " + str(MT2_R_PhiReq_f4_SB)

print "MT2_R_mp_SB: " + str(MT2_R_mp_SB)
print "MT2_R_DPhiPair_SB: " + str(MT2_R_DPhiPair_SB)
print "MT2_R_DPhiPair_Pt_SB: " + str(MT2_R_DPhiPair_Pt_SB)
print "MT2_R_PhiReq_SB: " + str(MT2_R_PhiReq_SB)

print "MT2_R_Pt_mp_SB: " + str(MT2_R_Pt_mp_SB)
print "MT2_R_Pt_mp_DPhiPair_SB: " + str(MT2_R_Pt_mp_DPhiPair_SB)

# Setting up a way to normalize the pT distribution to QM pT distribution.
# The idea is to weight an event based on the jet pT.
SM_Pt_Norm = Dark_Pt.Clone("SM_Pt_Norm")
SMM_Pt_Norm = Dark_Pt.Clone("SMM_Pt_Norm")
G_Pt_Norm = Dark_Pt.Clone("G_Pt_Norm")
QM_Pt_Norm = Dark_Pt.Clone("QM_Pt_Norm")
Q_Pt_Norm = Dark_Pt.Clone("Q_Pt_Norm")
QM_QPt_Norm = Dark_Pt.Clone("QM_QPt_Norm")
QM_GPt_Norm = Dark_Pt.Clone("QM_GPt_Norm")
Q_GPt_Norm = Dark_Pt.Clone("Q_GPt_Norm")
G_SMPt_Norm = Dark_Pt.Clone("G_SMPt_Norm")
QM_SMPt_Norm = Dark_Pt.Clone("QM_SMPt_Norm")
Q_SMPt_Norm = Dark_Pt.Clone("Q_SMPt_Norm")
LD_lowDFPt_Norm = Dark_Pt.Clone("LD_lowDFPt_Norm")
LD_highDFPt_Norm = Dark_Pt.Clone("LD_highDFPt_Norm")
LD_SMPt_Norm = Dark_Pt.Clone("LD_SMPt_Norm")

SM_Pt_Norm.Divide(SM_Pt)
SMM_Pt_Norm.Divide(SMM_Pt)
G_Pt_Norm.Divide(G_Pt)
QM_Pt_Norm.Divide(QM_Pt)
Q_Pt_Norm.Divide(Q_Pt)
QM_QPt_Norm.Divide(QM_QPt)
QM_GPt_Norm.Divide(QM_GPt)
Q_GPt_Norm.Divide(Q_GPt)
G_SMPt_Norm.Divide(G_SMPt)
QM_SMPt_Norm.Divide(QM_SMPt)
Q_SMPt_Norm.Divide(Q_SMPt)
LD_lowDFPt_Norm.Divide(LD_lowDFPt)
LD_highDFPt_Norm.Divide(LD_highDFPt)
LD_SMPt_Norm.Divide(LD_SMPt)

SM_SDPt_Norm = Dark_SDPt.Clone("SM_SDPt_Norm")
SMM_SDPt_Norm = Dark_SDPt.Clone("SMM_SDPt_Norm")
G_SDPt_Norm = Dark_SDPt.Clone("G_SDPt_Norm")
QM_SDPt_Norm = Dark_SDPt.Clone("QM_SDPt_Norm")
Q_SDPt_Norm = Dark_SDPt.Clone("Q_SDPt_Norm")
QM_QSDPt_Norm = Dark_SDPt.Clone("QM_QSDPt_Norm")
QM_GSDPt_Norm = Dark_SDPt.Clone("QM_GSDPt_Norm")
Q_GSDPt_Norm = Dark_SDPt.Clone("Q_GSDPt_Norm")
G_SMSDPt_Norm = Dark_SDPt.Clone("G_SMSDPt_Norm")
QM_SMSDPt_Norm = Dark_SDPt.Clone("QM_SMSDPt_Norm")
Q_SMSDPt_Norm = Dark_SDPt.Clone("Q_SMSDPt_Norm")
LD_lowDFSDPt_Norm = Dark_SDPt.Clone("LD_lowDFSDPt_Norm")
LD_highDFSDPt_Norm = Dark_SDPt.Clone("LD_highDFSDPt_Norm")
LD_SMSDPt_Norm = Dark_SDPt.Clone("LD_SMSDPt_Norm")
#
# print ("Normal pT weight")
# getPtWeight(Dark_Pt)
# print ("SD pT weight")
# getPtWeight(Dark_SDPt)

SM_SDPt_Norm.Divide(SM_SDPt)
SMM_SDPt_Norm.Divide(SMM_SDPt)
G_SDPt_Norm.Divide(G_SDPt)
QM_SDPt_Norm.Divide(QM_SDPt)
Q_SDPt_Norm.Divide(Q_SDPt)
QM_QSDPt_Norm.Divide(QM_QSDPt)
QM_GSDPt_Norm.Divide(QM_GSDPt)
Q_GSDPt_Norm.Divide(Q_GSDPt)
G_SMSDPt_Norm.Divide(G_SMSDPt)
QM_SMSDPt_Norm.Divide(QM_SMSDPt)
Q_SMSDPt_Norm.Divide(Q_SMSDPt)
LD_lowDFSDPt_Norm.Divide(LD_lowDFSDPt)
LD_highDFSDPt_Norm.Divide(LD_highDFSDPt)
LD_SMSDPt_Norm.Divide(LD_SMSDPt)

# Normalizing the looser categories to the pT of Dark_
SM_OnlyPt_Norm = Dark_Pt.Clone("SM_OnlyPt_Norm")
Mix_Pt_Norm = Dark_Pt.Clone("Mix_Pt_Norm")
DarkerMix_Pt_Norm = Dark_Pt.Clone("DarkerMix_Pt_Norm")
Dark_Pt_Norm = Dark_Pt.Clone("Dark_Pt_Norm")
SoftDark_Pt_Norm = Dark_Pt.Clone("SoftDark_Pt_Norm")

SM_OnlyPt_Norm.Divide(SM_OnlyPt)
Mix_Pt_Norm.Divide(Mix_Pt)
DarkerMix_Pt_Norm.Divide(DarkerMix_Pt)
Dark_Pt_Norm.Divide(Dark_Pt)
SoftDark_Pt_Norm.Divide(SoftDark_Pt)

SM_OnlySDPt_Norm = Dark_SDPt.Clone("SM_OnlySDPt_Norm")
Mix_SDPt_Norm = Dark_SDPt.Clone("Mix_SDPt_Norm")
DarkerMix_SDPt_Norm = Dark_SDPt.Clone("DarkerMix_SDPt_Norm")
Dark_SDPt_Norm = Dark_SDPt.Clone("Dark_SDPt_Norm")
SoftDark_SDPt_Norm = Dark_SDPt.Clone("SoftDark_SDPt_Norm")

SM_OnlySDPt_Norm.Divide(SM_OnlySDPt)
Mix_SDPt_Norm.Divide(Mix_SDPt)
DarkerMix_SDPt_Norm.Divide(DarkerMix_SDPt)
Dark_SDPt_Norm.Divide(Dark_SDPt)
SoftDark_SDPt_Norm.Divide(SoftDark_SDPt)

for count in range(nEvents):
	tr.GetEntry(count)
	tr2.GetEntry(count)
	# print "Event number " + str(count)

	jets = VecToList(tr.AK8Jets)
	met4p = tr.MET[0]
	SM_Jets = VecToList(tr.AK8_SM_Jets)
	SMM_Jets = VecToList(tr.AK8_SMM_Jets)
	G_Jets = VecToList(tr.AK8_G_Jets)
	QM_Jets = VecToList(tr.AK8_QM_Jets)
	Q_Jets = VecToList(tr.AK8_Q_Jets)
	QM_QJets = VecToList(tr.AK8_QM_QJets)
	QM_GJets = VecToList(tr.AK8_QM_GJets)
	Q_GJets = VecToList(tr.AK8_Q_GJets)
	G_SMJets = VecToList(tr.AK8_G_SMJets)
	QM_SMJets = VecToList(tr.AK8_QM_SMJets)
	Q_SMJets = VecToList(tr.AK8_Q_SMJets)
	LD_lowDFJets = VecToList(tr.AK8_LD_lowDFJets)
	LD_highDFJets = VecToList(tr.AK8_LD_highDFJets)
	LD_SMJets = VecToList(tr.AK8_LD_SMJets)

	G_dpTFrac = VecToList(tr.G_dpTFrac)
	QM_dpTFrac = VecToList(tr.QM_dpTFrac)
	Q_dpTFrac = VecToList(tr.Q_dpTFrac)
	QM_QdpTFrac = VecToList(tr.QM_QdpTFrac)
	QM_GdpTFrac = VecToList(tr.QM_GdpTFrac)
	Q_GdpTFrac = VecToList(tr.Q_GdpTFrac)
	G_SMdpTFrac = VecToList(tr.G_SMdpTFrac)
	QM_SMdpTFrac = VecToList(tr.QM_SMdpTFrac)
	Q_SMdpTFrac = VecToList(tr.Q_SMdpTFrac)
	LD_lowDFdpTFrac = VecToList(tr.LD_lowDFdpTFrac)
	LD_highDFdpTFrac = VecToList(tr.LD_highDFdpTFrac)
	LD_SMdpTFrac = VecToList(tr.LD_SMdpTFrac)

	dgir = VecToList(tr.Dark_girth)
	daMj = VecToList(tr.Dark_axisMj)
	daMn = VecToList(tr.Dark_axisMn)
	dpD = VecToList(tr.Dark_ptD)

	SM_gir = VecToList(tr.SM_girth)
	SM_aMj = VecToList(tr.SM_axisMj)
	SM_aMn = VecToList(tr.SM_axisMn)
	SM_pD = VecToList(tr.SM_ptD)

	SMM_gir = VecToList(tr.SMM_girth)
	SMM_aMj = VecToList(tr.SMM_axisMj)
	SMM_aMn = VecToList(tr.SMM_axisMn)
	SMM_pD = VecToList(tr.SMM_ptD)

	G_gir = VecToList(tr.G_girth)
	G_aMj = VecToList(tr.G_axisMj)
	G_aMn = VecToList(tr.G_axisMn)
	G_pD = VecToList(tr.G_ptD)

	QM_gir = VecToList(tr.QM_girth)
	QM_aMj = VecToList(tr.QM_axisMj)
	QM_aMn = VecToList(tr.QM_axisMn)
	QM_pD = VecToList(tr.QM_ptD)

	Q_gir = VecToList(tr.Q_girth)
	Q_aMj = VecToList(tr.Q_axisMj)
	Q_aMn = VecToList(tr.Q_axisMn)
	Q_pD = VecToList(tr.Q_ptD)

	QM_Qgir = VecToList(tr.QM_Qgirth)
	QM_QaMj = VecToList(tr.QM_QaxisMj)
	QM_QaMn = VecToList(tr.QM_QaxisMn)
	QM_QpD = VecToList(tr.QM_QptD)

	QM_Ggir = VecToList(tr.QM_Ggirth)
	QM_GaMj = VecToList(tr.QM_GaxisMj)
	QM_GaMn = VecToList(tr.QM_GaxisMn)
	QM_GpD = VecToList(tr.QM_GptD)

	Q_Ggir = VecToList(tr.Q_Ggirth)
	Q_GaMj = VecToList(tr.Q_GaxisMj)
	Q_GaMn = VecToList(tr.Q_GaxisMn)
	Q_GpD = VecToList(tr.Q_GptD)

	G_SMgir = VecToList(tr.G_SMgirth)
	G_SMaMj = VecToList(tr.G_SMaxisMj)
	G_SMaMn = VecToList(tr.G_SMaxisMn)
	G_SMpD = VecToList(tr.G_SMptD)

	QM_SMgir = VecToList(tr.QM_SMgirth)
	QM_SMaMj = VecToList(tr.QM_SMaxisMj)
	QM_SMaMn = VecToList(tr.QM_SMaxisMn)
	QM_SMpD = VecToList(tr.QM_SMptD)

	Q_SMgir = VecToList(tr.Q_SMgirth)
	Q_SMaMj = VecToList(tr.Q_SMaxisMj)
	Q_SMaMn = VecToList(tr.Q_SMaxisMn)
	Q_SMpD = VecToList(tr.Q_SMptD)

	LD_lowDFgir = VecToList(tr.LD_lowDFgirth)
	LD_lowDFaMj = VecToList(tr.LD_lowDFaxisMj)
	LD_lowDFaMn = VecToList(tr.LD_lowDFaxisMn)
	LD_lowDFpD = VecToList(tr.LD_lowDFptD)

	LD_highDFgir = VecToList(tr.LD_highDFgirth)
	LD_highDFaMj = VecToList(tr.LD_highDFaxisMj)
	LD_highDFaMn = VecToList(tr.LD_highDFaxisMn)
	LD_highDFpD = VecToList(tr.LD_highDFptD)

	LD_SMgir = VecToList(tr.LD_SMgirth)
	LD_SMaMj = VecToList(tr.LD_SMaxisMj)
	LD_SMaMn = VecToList(tr.LD_SMaxisMn)
	LD_SMpD = VecToList(tr.LD_SMptD)

	npDJets = VecToList(tr.npDJets)
	npSM_Jets = VecToList(tr.npSM_Jets)
	npSMM_Jets = VecToList(tr.npSMM_Jets)
	npG_Jets = VecToList(tr.npG_Jets)
	npQM_Jets = VecToList(tr.npQM_Jets)
	npQ_Jets = VecToList(tr.npQ_Jets)
	npQM_QJets = VecToList(tr.npQM_QJets)
	npQM_GJets = VecToList(tr.npQM_GJets)
	npQ_GJets = VecToList(tr.npQ_GJets)
	npG_SMJets = VecToList(tr.npG_SMJets)
	npQM_SMJets = VecToList(tr.npQM_SMJets)
	npQ_SMJets = VecToList(tr.npQ_SMJets)
	npLD_lowDFJets = VecToList(tr.npLD_lowDFJets)
	npLD_highDFJets = VecToList(tr.npLD_highDFJets)
	npLD_SMJets = VecToList(tr.npLD_SMJets)

	# taus
	SM_Tau1 = VecToList(tr2.SM_Tau1)
	SM_Tau2 = VecToList(tr2.SM_Tau2)
	SM_Tau3 = VecToList(tr2.SM_Tau3)
	SMM_Tau1 = VecToList(tr2.SMM_Tau1)
	SMM_Tau2 = VecToList(tr2.SMM_Tau2)
	SMM_Tau3 = VecToList(tr2.SMM_Tau3)
	G_Tau1 = VecToList(tr2.G_Tau1)
	G_Tau2 = VecToList(tr2.G_Tau2)
	G_Tau3 = VecToList(tr2.G_Tau3)
	QM_Tau1 = VecToList(tr2.QM_Tau1)
	QM_Tau2 = VecToList(tr2.QM_Tau2)
	QM_Tau3 = VecToList(tr2.QM_Tau3)
	Q_Tau1 = VecToList(tr2.Q_Tau1)
	Q_Tau2 = VecToList(tr2.Q_Tau2)
	Q_Tau3 = VecToList(tr2.Q_Tau3)
	QM_QTau1 = VecToList(tr2.QM_QTau1)
	QM_QTau2 = VecToList(tr2.QM_QTau2)
	QM_QTau3 = VecToList(tr2.QM_QTau3)
	QM_GTau1 = VecToList(tr2.QM_GTau1)
	QM_GTau2 = VecToList(tr2.QM_GTau2)
	QM_GTau3 = VecToList(tr2.QM_GTau3)
	Q_GTau1 = VecToList(tr2.Q_GTau1)
	Q_GTau2 = VecToList(tr2.Q_GTau2)
	Q_GTau3 = VecToList(tr2.Q_GTau3)
	G_SMTau1 = VecToList(tr2.G_SMTau1)
	G_SMTau2 = VecToList(tr2.G_SMTau2)
	G_SMTau3 = VecToList(tr2.G_SMTau3)
	QM_SMTau1 = VecToList(tr2.QM_SMTau1)
	QM_SMTau2 = VecToList(tr2.QM_SMTau2)
	QM_SMTau3 = VecToList(tr2.QM_SMTau3)
	Q_SMTau1 = VecToList(tr2.Q_SMTau1)
	Q_SMTau2 = VecToList(tr2.Q_SMTau2)
	Q_SMTau3 = VecToList(tr2.Q_SMTau3)
	LD_lowDFTau1 = VecToList(tr2.LD_lowDFTau1)
	LD_lowDFTau2 = VecToList(tr2.LD_lowDFTau2)
	LD_lowDFTau3 = VecToList(tr2.LD_lowDFTau3)
	LD_highDFTau1 = VecToList(tr2.LD_highDFTau1)
	LD_highDFTau2 = VecToList(tr2.LD_highDFTau2)
	LD_highDFTau3 = VecToList(tr2.LD_highDFTau3)
	LD_SMTau1 = VecToList(tr2.LD_SMTau1)
	LD_SMTau2 = VecToList(tr2.LD_SMTau2)
	LD_SMTau3 = VecToList(tr2.LD_SMTau3)

	# soft drop jets and corresponding dark pT fraction
	SM_SDJets = VecToList(tr2.SM_SDJets)
	SMM_SDJets = VecToList(tr2.SMM_SDJets)
	G_SDJets = VecToList(tr2.G_SDJets)
	QM_SDJets = VecToList(tr2.QM_SDJets)
	Q_SDJets = VecToList(tr2.Q_SDJets)
	QM_QSDJets = VecToList(tr2.QM_QSDJets)
	QM_GSDJets = VecToList(tr2.QM_GSDJets)
	Q_GSDJets = VecToList(tr2.Q_GSDJets)
	G_SMSDJets = VecToList(tr2.G_SMSDJets)
	QM_SMSDJets = VecToList(tr2.QM_SMSDJets)
	Q_SMSDJets = VecToList(tr2.Q_SMSDJets)
	LD_lowDFSDJets = VecToList(tr2.LD_lowDFSDJets)
	LD_highDFSDJets = VecToList(tr2.LD_highDFSDJets)
	LD_SMSDJets = VecToList(tr2.LD_SMSDJets)

	G_SDdptFrac = VecToList(tr2.G_SDdptFrac)
	QM_SDdptFrac = VecToList(tr2.QM_SDdptFrac)
	Q_SDdptFrac = VecToList(tr2.Q_SDdptFrac)
	QM_QSDdptFrac = VecToList(tr2.QM_QSDdptFrac)
	QM_GSDdptFrac = VecToList(tr2.QM_GSDdptFrac)
	Q_GSDdptFrac = VecToList(tr2.Q_GSDdptFrac)
	G_SMSDdptFrac = VecToList(tr2.G_SMSDdptFrac)
	QM_SMSDdptFrac = VecToList(tr2.QM_SMSDdptFrac)
	Q_SMSDdptFrac = VecToList(tr2.Q_SMSDdptFrac)
	LD_lowDFSDdptFrac = VecToList(tr2.LD_lowDFSDdptFrac)
	LD_highDFSDdptFrac = VecToList(tr2.LD_highDFSDdptFrac)
	LD_SMSDdptFrac = VecToList(tr2.LD_SMSDdptFrac)
	SM_SDdptFrac = VecToList(tr2.SM_SDdptFrac)
	SMM_SDdptFrac = VecToList(tr2.SMM_SDdptFrac)

	METv = met4p.Pt()
	METPhiv = met4p.Phi()

	SM_OnlyJets = SM_Jets + SMM_Jets
	SM_OnlydpTFrac = [0]*len(SM_OnlyJets)
	SM_Onlygir = SM_gir + SMM_gir
	SM_OnlyaMj = SM_aMj + SMM_aMj
	SM_OnlyaMn = SM_aMn + SMM_aMn
	SM_OnlypD = SM_pD + SMM_pD
	npSM_OnlyJets = npSM_Jets + npSMM_Jets
	SM_OnlyTau1 = SM_Tau1 + SMM_Tau1
	SM_OnlyTau2 = SM_Tau2 + SMM_Tau2
	SM_OnlyTau3 = SM_Tau3 + SMM_Tau3
	SM_OnlySDJets = SM_SDJets + SMM_SDJets
	SM_OnlySDdpTFrac = [0]*len(SM_OnlySDJets)

	Mix_Jets = G_SMJets + LD_SMJets + LD_lowDFJets
	Mix_dpTFrac = G_SMdpTFrac + LD_SMdpTFrac + LD_lowDFdpTFrac
	Mix_gir = G_SMgir + LD_SMgir + LD_lowDFgir
	Mix_aMj = G_SMaMj + LD_SMaMj + LD_lowDFaMj
	Mix_aMn = G_SMaMn + LD_SMaMn + LD_lowDFaMn
	Mix_pD = G_SMpD + LD_SMpD + LD_lowDFpD
	npMix_Jets = npG_SMJets + npLD_SMJets + npLD_lowDFJets
	Mix_Tau1 = G_SMTau1 + LD_SMTau1 + LD_lowDFTau1
	Mix_Tau2 = G_SMTau2 + LD_SMTau2 + LD_lowDFTau2
	Mix_Tau3 = G_SMTau3 + LD_SMTau3 + LD_lowDFTau3
	Mix_SDJets = G_SMSDJets + LD_SMSDJets + LD_lowDFSDJets
	Mix_SDdpTFrac = G_SMSDdptFrac + LD_SMSDdptFrac + LD_lowDFSDdptFrac

	DarkerMix_Jets = QM_SMJets + Q_SMJets + QM_QJets
	DarkerMix_dpTFrac = QM_SMdpTFrac + Q_SMdpTFrac + QM_QdpTFrac
	DarkerMix_gir = QM_SMgir + Q_SMgir + QM_Qgir
	DarkerMix_aMj = QM_SMaMj + Q_SMaMj + QM_QaMj
	DarkerMix_aMn = QM_SMaMn + Q_SMaMn + QM_QaMn
	DarkerMix_pD = QM_SMpD + Q_SMpD + QM_QpD
	npDarkerMix_Jets = npQM_SMJets + npQ_SMJets + npQM_QJets
	DarkerMix_Tau1 = QM_SMTau1 + Q_SMTau1 + QM_QTau1
	DarkerMix_Tau2 = QM_SMTau2 + Q_SMTau2 + QM_QTau2
	DarkerMix_Tau3 = QM_SMTau3 + Q_SMTau3 + QM_QTau3
	DarkerMix_SDJets = QM_SMSDJets + Q_SMSDJets + QM_QSDJets
	DarkerMix_SDdpTFrac = QM_SMSDdptFrac + Q_SMSDdptFrac + QM_QSDdptFrac

	Dark_Jets = QM_Jets + QM_GJets + Q_GJets + Q_Jets
	Dark_dpTFrac = QM_dpTFrac + QM_GdpTFrac + Q_GdpTFrac + Q_dpTFrac
	Dark_gir = QM_gir + QM_Ggir + Q_Ggir + Q_gir
	Dark_aMj = QM_aMj + QM_GaMj + Q_GaMj + Q_aMj
	Dark_aMn = QM_aMn + QM_GaMn + Q_GaMn + Q_aMn
	Dark_pD = QM_pD + QM_GpD + Q_GpD + Q_pD
	npDark_Jets = npQM_Jets + npQM_GJets + npQ_GJets + npQ_Jets
	Dark_Tau1 = QM_Tau1 + QM_GTau1 + Q_GTau1 + Q_Tau1
	Dark_Tau2 = QM_Tau2 + QM_GTau2 + Q_GTau2 + Q_Tau2
	Dark_Tau3 = QM_Tau3 + QM_GTau3 + Q_GTau3 + Q_Tau3
	Dark_SDJets = QM_SDJets + QM_GSDJets + Q_GSDJets + Q_SDJets
	Dark_SDdpTFrac = QM_SDdptFrac + QM_GSDdptFrac + Q_GSDdptFrac + Q_SDdptFrac

	SoftDark_Jets = G_Jets + LD_highDFJets
	SoftDark_dpTFrac = G_dpTFrac + LD_highDFdpTFrac
	SoftDark_gir = G_gir + LD_highDFgir
	SoftDark_aMj = G_aMj + LD_highDFaMj
	SoftDark_aMn = G_aMn + LD_highDFaMn
	SoftDark_pD = G_pD + LD_highDFpD
	npSoftDark_Jets = npG_Jets + npLD_highDFJets
	SoftDark_Tau1 = G_Tau1 + LD_highDFTau1
	SoftDark_Tau2 = G_Tau2 + LD_highDFTau2
	SoftDark_Tau3 = G_Tau3 + LD_highDFTau3
	SoftDark_SDJets = G_SDJets + LD_highDFSDJets
	SoftDark_SDdpTFrac = G_SDdptFrac + LD_highDFSDdptFrac

	fillDPT([0]*len(SM_Jets),jet_SM_dpTFrac_Norm,SM_Jets,SM_Pt_Norm)
	fillDPT([0]*len(SMM_Jets),jet_SMM_dpTFrac_Norm,SMM_Jets,SMM_Pt_Norm)
	fillDPT(G_dpTFrac,jet_G_dpTFrac_Norm,G_Jets,G_Pt_Norm)
	fillDPT(QM_dpTFrac,jet_QM_dpTFrac_Norm,QM_Jets,QM_Pt_Norm)
	fillDPT(Q_dpTFrac,jet_Q_dpTFrac_Norm,Q_Jets,Q_Pt_Norm)
	fillDPT(QM_QdpTFrac,jet_QM_QdpTFrac_Norm,QM_QJets,QM_QPt_Norm)
	fillDPT(QM_GdpTFrac,jet_QM_GdpTFrac_Norm,QM_GJets,QM_GPt_Norm)
	fillDPT(Q_GdpTFrac,jet_Q_GdpTFrac_Norm,Q_GJets,Q_GPt_Norm)
	fillDPT(G_SMdpTFrac,jet_G_SMdpTFrac_Norm,G_SMJets,G_SMPt_Norm)
	fillDPT(QM_SMdpTFrac,jet_QM_SMdpTFrac_Norm,QM_SMJets,QM_SMPt_Norm)
	fillDPT(Q_SMdpTFrac,jet_Q_SMdpTFrac_Norm,Q_SMJets,Q_SMPt_Norm)
	fillDPT(LD_highDFdpTFrac,jet_LD_highDFdpTFrac_Norm,LD_highDFJets,LD_highDFPt_Norm)
	fillDPT(LD_lowDFdpTFrac,jet_LD_lowDFdpTFrac_Norm,LD_lowDFJets,LD_lowDFPt_Norm)
	fillDPT(LD_SMdpTFrac,jet_LD_SMdpTFrac_Norm,LD_SMJets,LD_SMPt_Norm)

	fillDPT(SM_OnlydpTFrac,jet_SM_OnlydpTFrac_Norm,SM_OnlyJets,SM_OnlyPt_Norm)
	fillDPT(Mix_dpTFrac,jet_Mix_dpTFrac_Norm,Mix_Jets,Mix_Pt_Norm)
	fillDPT(DarkerMix_dpTFrac,jet_DarkerMix_dpTFrac_Norm,DarkerMix_Jets,DarkerMix_Pt_Norm)
	fillDPT(Dark_dpTFrac,jet_Dark_dpTFrac_Norm,Dark_Jets,Dark_Pt_Norm)
	fillDPT(SoftDark_dpTFrac,jet_SoftDark_dpTFrac_Norm,SoftDark_Jets,SoftDark_Pt_Norm)

	# jet mass plots
	fillM(SM_Jets,SM_M_Norm,SM_M_2D,[0]*len(SM_Jets),SM_Pt_Norm)
	fillM(SMM_Jets,SMM_M_Norm,SMM_M_2D,[0]*len(SMM_Jets),SMM_Pt_Norm)
	fillM(G_Jets,G_M_Norm,G_M_2D,G_dpTFrac,G_Pt_Norm)
	fillM(QM_Jets,QM_M_Norm,QM_M_2D,QM_dpTFrac,QM_Pt_Norm)
	fillM(Q_Jets,Q_M_Norm,Q_M_2D,Q_dpTFrac,Q_Pt_Norm)
	fillM(QM_QJets,QM_QM_Norm,QM_QM_2D,QM_QdpTFrac,QM_QPt_Norm)
	fillM(QM_GJets,QM_GM_Norm,QM_GM_2D,QM_GdpTFrac,QM_GPt_Norm)
	fillM(Q_GJets,Q_GM_Norm,Q_GM_2D,Q_GdpTFrac,Q_GPt_Norm)
	fillM(G_SMJets,G_SMM_Norm,G_SMM_2D,G_SMdpTFrac,G_SMPt_Norm)
	fillM(QM_SMJets,QM_SMM_Norm,QM_SMM_2D,QM_SMdpTFrac,QM_SMPt_Norm)
	fillM(Q_SMJets,Q_SMM_Norm,Q_SMM_2D,Q_SMdpTFrac,Q_SMPt_Norm)
	fillM(LD_highDFJets,LD_highDFM_Norm,LD_highDFM_2D,LD_highDFdpTFrac,LD_highDFPt_Norm)
	fillM(LD_lowDFJets,LD_lowDFM_Norm,LD_lowDFM_2D,LD_lowDFdpTFrac,LD_lowDFPt_Norm)
	fillM(LD_SMJets,LD_SMM_Norm,LD_SMM_2D,LD_SMdpTFrac,LD_SMPt_Norm)

	fillM(SM_OnlyJets,SM_OnlyM_Norm,filler2D,[0]*len(SM_OnlyJets),SM_OnlyPt_Norm)
	fillM(Mix_Jets,Mix_M_Norm,filler2D,Mix_dpTFrac,Mix_Pt_Norm)
	fillM(DarkerMix_Jets,DarkerMix_M_Norm,filler2D,DarkerMix_dpTFrac,DarkerMix_Pt_Norm)
	fillM(Dark_Jets,Dark_M_Norm,filler2D,Dark_dpTFrac,Dark_Pt_Norm)
	fillM(SoftDark_Jets,SoftDark_M_Norm,filler2D,SoftDark_dpTFrac,SoftDark_Pt_Norm)

	# jet substructure plots
	makejsubplot(SM_Jets,SM_aMj,SM_aMn,SM_gir,SM_pD,npSM_Jets,
	SM_axisMajor_Norm,SM_axisMinor_Norm,SM_momentGirth_Norm,SM_ptD_Norm,SM_pmult_Norm,SM_axisMajor_2D,
	SM_axisMinor_2D,SM_momentGirth_2D,SM_ptD_2D,SM_pmult_2D,[0]*len(SM_Jets),SM_Pt_Norm)
	makejsubplot(SMM_Jets,SMM_aMj,SMM_aMn,SMM_gir,SMM_pD,npSMM_Jets,
	SMM_axisMajor_Norm,SMM_axisMinor_Norm,SMM_momentGirth_Norm,SMM_ptD_Norm,SMM_pmult_Norm,SMM_axisMajor_2D,
	SMM_axisMinor_2D,SMM_momentGirth_2D,SMM_ptD_2D,SMM_pmult_2D,[0]*len(SMM_Jets),SMM_Pt_Norm)
	makejsubplot(G_Jets,G_aMj,G_aMn,G_gir,G_pD,npG_Jets,
	G_axisMajor_Norm,G_axisMinor_Norm,G_momentGirth_Norm,G_ptD_Norm,G_pmult_Norm,G_axisMajor_2D,
	G_axisMinor_2D,G_momentGirth_2D,G_ptD_2D,G_pmult_2D,G_dpTFrac,G_Pt_Norm)
	makejsubplot(QM_Jets,QM_aMj,QM_aMn,QM_gir,QM_pD,npQM_Jets,
	QM_axisMajor_Norm,QM_axisMinor_Norm,QM_momentGirth_Norm,QM_ptD_Norm,QM_pmult_Norm,QM_axisMajor_2D,
	QM_axisMinor_2D,QM_momentGirth_2D,QM_ptD_2D,QM_pmult_2D,QM_dpTFrac,QM_Pt_Norm)
	makejsubplot(Q_Jets,Q_aMj,Q_aMn,Q_gir,Q_pD,npQ_Jets,
	Q_axisMajor_Norm,Q_axisMinor_Norm,Q_momentGirth_Norm,Q_ptD_Norm,Q_pmult_Norm,Q_axisMajor_2D,
	Q_axisMinor_2D,Q_momentGirth_2D,Q_ptD_2D,Q_pmult_2D,Q_dpTFrac,Q_Pt_Norm)
	makejsubplot(QM_QJets,QM_QaMj,QM_QaMn,QM_Qgir,QM_QpD,npQM_QJets,
	QM_QaxisMajor_Norm,QM_QaxisMinor_Norm,QM_QmomentGirth_Norm,QM_QptD_Norm,QM_Qpmult_Norm,QM_QaxisMajor_2D,
	QM_QaxisMinor_2D,QM_QmomentGirth_2D,QM_QptD_2D,QM_Qpmult_2D,QM_QdpTFrac,QM_QPt_Norm)
	makejsubplot(QM_GJets,QM_GaMj,QM_GaMn,QM_Ggir,QM_GpD,npQM_GJets,
	QM_GaxisMajor_Norm,QM_GaxisMinor_Norm,QM_GmomentGirth_Norm,QM_GptD_Norm,QM_Gpmult_Norm,QM_GaxisMajor_2D,
	QM_GaxisMinor_2D,QM_GmomentGirth_2D,QM_GptD_2D,QM_Gpmult_2D,QM_GdpTFrac,QM_GPt_Norm)
	makejsubplot(Q_GJets,Q_GaMj,Q_GaMn,Q_Ggir,Q_GpD,npQ_GJets,
	Q_GaxisMajor_Norm,Q_GaxisMinor_Norm,Q_GmomentGirth_Norm,Q_GptD_Norm,Q_Gpmult_Norm,Q_GaxisMajor_2D,
	Q_GaxisMinor_2D,Q_GmomentGirth_2D,Q_GptD_2D,Q_Gpmult_2D,Q_GdpTFrac,Q_GPt_Norm)
	makejsubplot(G_SMJets,G_SMaMj,G_SMaMn,G_SMgir,G_SMpD,npG_SMJets,
	G_SMaxisMajor_Norm,G_SMaxisMinor_Norm,G_SMmomentGirth_Norm,G_SMptD_Norm,G_SMpmult_Norm,G_SMaxisMajor_2D,
	G_SMaxisMinor_2D,G_SMmomentGirth_2D,G_SMptD_2D,G_SMpmult_2D,G_SMdpTFrac,G_SMPt_Norm)
	makejsubplot(QM_SMJets,QM_SMaMj,QM_SMaMn,QM_SMgir,QM_SMpD,npQM_SMJets,
	QM_SMaxisMajor_Norm,QM_SMaxisMinor_Norm,QM_SMmomentGirth_Norm,QM_SMptD_Norm,QM_SMpmult_Norm,QM_SMaxisMajor_2D,
	QM_SMaxisMinor_2D,QM_SMmomentGirth_2D,QM_SMptD_2D,QM_SMpmult_2D,QM_SMdpTFrac,QM_SMPt_Norm)
	makejsubplot(Q_SMJets,Q_SMaMj,Q_SMaMn,Q_SMgir,Q_SMpD,npQ_SMJets,
	Q_SMaxisMajor_Norm,Q_SMaxisMinor_Norm,Q_SMmomentGirth_Norm,Q_SMptD_Norm,Q_SMpmult_Norm,Q_SMaxisMajor_2D,
	Q_SMaxisMinor_2D,Q_SMmomentGirth_2D,Q_SMptD_2D,Q_SMpmult_2D,Q_SMdpTFrac,Q_SMPt_Norm)
	makejsubplot(LD_highDFJets,LD_highDFaMj,LD_highDFaMn,LD_highDFgir,LD_highDFpD,npLD_highDFJets,
	LD_highDFaxisMajor_Norm,LD_highDFaxisMinor_Norm,LD_highDFmomentGirth_Norm,LD_highDFptD_Norm,LD_highDFpmult_Norm,LD_highDFaxisMajor_2D,
	LD_highDFaxisMinor_2D,LD_highDFmomentGirth_2D,LD_highDFptD_2D,LD_highDFpmult_2D,LD_highDFdpTFrac,LD_highDFPt_Norm)
	makejsubplot(LD_lowDFJets,LD_lowDFaMj,LD_lowDFaMn,LD_lowDFgir,LD_lowDFpD,npLD_lowDFJets,
	LD_lowDFaxisMajor_Norm,LD_lowDFaxisMinor_Norm,LD_lowDFmomentGirth_Norm,LD_lowDFptD_Norm,LD_lowDFpmult_Norm,LD_lowDFaxisMajor_2D,
	LD_lowDFaxisMinor_2D,LD_lowDFmomentGirth_2D,LD_lowDFptD_2D,LD_lowDFpmult_2D,LD_lowDFdpTFrac,LD_lowDFPt_Norm)
	makejsubplot(LD_SMJets,LD_SMaMj,LD_SMaMn,LD_SMgir,LD_SMpD,npLD_SMJets,
	LD_SMaxisMajor_Norm,LD_SMaxisMinor_Norm,LD_SMmomentGirth_Norm,LD_SMptD_Norm,LD_SMpmult_Norm,LD_SMaxisMajor_2D,
	LD_SMaxisMinor_2D,LD_SMmomentGirth_2D,LD_SMptD_2D,LD_SMpmult_2D,LD_SMdpTFrac,LD_SMPt_Norm)

	makejsubplot(SM_OnlyJets,SM_OnlyaMj,SM_OnlyaMn,SM_Onlygir,SM_OnlypD,npSM_OnlyJets,
	SM_OnlyaxisMajor_Norm,SM_OnlyaxisMinor_Norm,SM_OnlymomentGirth_Norm,SM_OnlyptD_Norm,SM_Onlypmult_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,[0]*len(SM_OnlyJets),SM_OnlyPt_Norm)
	makejsubplot(Mix_Jets,Mix_aMj,Mix_aMn,Mix_gir,Mix_pD,npMix_Jets,
	Mix_axisMajor_Norm,Mix_axisMinor_Norm,Mix_momentGirth_Norm,Mix_ptD_Norm,Mix_pmult_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,Mix_dpTFrac,Mix_Pt_Norm)
	makejsubplot(DarkerMix_Jets,DarkerMix_aMj,DarkerMix_aMn,DarkerMix_gir,DarkerMix_pD,npDarkerMix_Jets,
	DarkerMix_axisMajor_Norm,DarkerMix_axisMinor_Norm,DarkerMix_momentGirth_Norm,DarkerMix_ptD_Norm,DarkerMix_pmult_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,DarkerMix_dpTFrac,DarkerMix_Pt_Norm)
	makejsubplot(Dark_Jets,Dark_aMj,Dark_aMn,Dark_gir,Dark_pD,npDark_Jets,
	Dark_axisMajor_Norm,Dark_axisMinor_Norm,Dark_momentGirth_Norm,Dark_ptD_Norm,Dark_pmult_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,Dark_dpTFrac,Dark_Pt_Norm)
	makejsubplot(SoftDark_Jets,SoftDark_aMj,SoftDark_aMn,SoftDark_gir,SoftDark_pD,npSoftDark_Jets,
	SoftDark_axisMajor_Norm,SoftDark_axisMinor_Norm,SoftDark_momentGirth_Norm,SoftDark_ptD_Norm,SoftDark_pmult_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,SoftDark_dpTFrac,SoftDark_Pt_Norm)

	## tau plots

	tauplot(SM_Jets,SM_Tau1,SM_Tau2,SM_Tau3,SM_tau1_Norm,
	SM_tau2_Norm,SM_tau3_Norm,SM_tau21_Norm,SM_tau32_Norm,SM_tau1_2D,
	SM_tau2_2D,SM_tau3_2D,SM_tau21_2D,SM_tau32_2D,[0]*len(SM_Jets),SM_Pt_Norm)
	tauplot(SMM_Jets,SMM_Tau1,SMM_Tau2,SMM_Tau3,SMM_tau1_Norm,
	SMM_tau2_Norm,SMM_tau3_Norm,SMM_tau21_Norm,SMM_tau32_Norm,SMM_tau1_2D,
	SMM_tau2_2D,SMM_tau3_2D,SMM_tau21_2D,SMM_tau32_2D,[0]*len(SMM_Jets),SMM_Pt_Norm)
	tauplot(G_Jets,G_Tau1,G_Tau2,G_Tau3,G_tau1_Norm,
	G_tau2_Norm,G_tau3_Norm,G_tau21_Norm,G_tau32_Norm,G_tau1_2D,
	G_tau2_2D,G_tau3_2D,G_tau21_2D,G_tau32_2D,G_dpTFrac,G_Pt_Norm)
	tauplot(QM_Jets,QM_Tau1,QM_Tau2,QM_Tau3,QM_tau1_Norm,
	QM_tau2_Norm,QM_tau3_Norm,QM_tau21_Norm,QM_tau32_Norm,QM_tau1_2D,
	QM_tau2_2D,QM_tau3_2D,QM_tau21_2D,QM_tau32_2D,QM_dpTFrac,QM_Pt_Norm)
	tauplot(Q_Jets,Q_Tau1,Q_Tau2,Q_Tau3,Q_tau1_Norm,
	Q_tau2_Norm,Q_tau3_Norm,Q_tau21_Norm,Q_tau32_Norm,Q_tau1_2D,
	Q_tau2_2D,Q_tau3_2D,Q_tau21_2D,Q_tau32_2D,Q_dpTFrac,Q_Pt_Norm)
	tauplot(QM_QJets,QM_QTau1,QM_QTau2,QM_QTau3,QM_Qtau1_Norm,
	QM_Qtau2_Norm,QM_Qtau3_Norm,QM_Qtau21_Norm,QM_Qtau32_Norm,QM_Qtau1_2D,
	QM_Qtau2_2D,QM_Qtau3_2D,QM_Qtau21_2D,QM_Qtau32_2D,QM_QdpTFrac,QM_QPt_Norm)
	tauplot(QM_GJets,QM_GTau1,QM_GTau2,QM_GTau3,QM_Gtau1_Norm,
	QM_Gtau2_Norm,QM_Gtau3_Norm,QM_Gtau21_Norm,QM_Gtau32_Norm,QM_Gtau1_2D,
	QM_Gtau2_2D,QM_Gtau3_2D,QM_Gtau21_2D,QM_Gtau32_2D,QM_GdpTFrac,QM_GPt_Norm)
	tauplot(Q_GJets,Q_GTau1,Q_GTau2,Q_GTau3,Q_Gtau1_Norm,
	Q_Gtau2_Norm,Q_Gtau3_Norm,Q_Gtau21_Norm,Q_Gtau32_Norm,Q_Gtau1_2D,
	Q_Gtau2_2D,Q_Gtau3_2D,Q_Gtau21_2D,Q_Gtau32_2D,Q_GdpTFrac,Q_GPt_Norm)
	tauplot(G_SMJets,G_SMTau1,G_SMTau2,G_SMTau3,G_SMtau1_Norm,
	G_SMtau2_Norm,G_SMtau3_Norm,G_SMtau21_Norm,G_SMtau32_Norm,G_SMtau1_2D,
	G_SMtau2_2D,G_SMtau3_2D,G_SMtau21_2D,G_SMtau32_2D,G_SMdpTFrac,G_SMPt_Norm)
	tauplot(QM_SMJets,QM_SMTau1,QM_SMTau2,QM_SMTau3,QM_SMtau1_Norm,
	QM_SMtau2_Norm,QM_SMtau3_Norm,QM_SMtau21_Norm,QM_SMtau32_Norm,QM_SMtau1_2D,
	QM_SMtau2_2D,QM_SMtau3_2D,QM_SMtau21_2D,QM_SMtau32_2D,QM_SMdpTFrac,QM_SMPt_Norm)
	tauplot(Q_SMJets,Q_SMTau1,Q_SMTau2,Q_SMTau3,Q_SMtau1_Norm,
	Q_SMtau2_Norm,Q_SMtau3_Norm,Q_SMtau21_Norm,Q_SMtau32_Norm,Q_SMtau1_2D,
	Q_SMtau2_2D,Q_SMtau3_2D,Q_SMtau21_2D,Q_SMtau32_2D,Q_SMdpTFrac,Q_SMPt_Norm)
	tauplot(LD_highDFJets,LD_highDFTau1,LD_highDFTau2,LD_highDFTau3,LD_highDFtau1_Norm,
	LD_highDFtau2_Norm,LD_highDFtau3_Norm,LD_highDFtau21_Norm,LD_highDFtau32_Norm,LD_highDFtau1_2D,
	LD_highDFtau2_2D,LD_highDFtau3_2D,LD_highDFtau21_2D,LD_highDFtau32_2D,LD_highDFdpTFrac,LD_highDFPt_Norm)
	tauplot(LD_lowDFJets,LD_lowDFTau1,LD_lowDFTau2,LD_lowDFTau3,LD_lowDFtau1_Norm,
	LD_lowDFtau2_Norm,LD_lowDFtau3_Norm,LD_lowDFtau21_Norm,LD_lowDFtau32_Norm,LD_lowDFtau1_2D,
	LD_lowDFtau2_2D,LD_lowDFtau3_2D,LD_lowDFtau21_2D,LD_lowDFtau32_2D,LD_lowDFdpTFrac,LD_lowDFPt_Norm)
	tauplot(LD_SMJets,LD_SMTau1,LD_SMTau2,LD_SMTau3,LD_SMtau1_Norm,
	LD_SMtau2_Norm,LD_SMtau3_Norm,LD_SMtau21_Norm,LD_SMtau32_Norm,LD_SMtau1_2D,
	LD_SMtau2_2D,LD_SMtau3_2D,LD_SMtau21_2D,LD_SMtau32_2D,LD_SMdpTFrac,LD_SMPt_Norm)

	tauplot(SM_OnlyJets,SM_OnlyTau1,SM_OnlyTau2,SM_OnlyTau3,SM_Onlytau1_Norm,
	SM_Onlytau2_Norm,SM_Onlytau3_Norm,SM_Onlytau21_Norm,SM_Onlytau32_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,[0]*len(SM_OnlyJets),SM_OnlyPt_Norm)
	tauplot(Mix_Jets,Mix_Tau1,Mix_Tau2,Mix_Tau3,Mix_tau1_Norm,
	Mix_tau2_Norm,Mix_tau3_Norm,Mix_tau21_Norm,Mix_tau32_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,Mix_dpTFrac,Mix_Pt_Norm)
	tauplot(DarkerMix_Jets,DarkerMix_Tau1,DarkerMix_Tau2,DarkerMix_Tau3,DarkerMix_tau1_Norm,
	DarkerMix_tau2_Norm,DarkerMix_tau3_Norm,DarkerMix_tau21_Norm,DarkerMix_tau32_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,DarkerMix_dpTFrac,DarkerMix_Pt_Norm)
	tauplot(Dark_Jets,Dark_Tau1,Dark_Tau2,Dark_Tau3,Dark_tau1_Norm,
	Dark_tau2_Norm,Dark_tau3_Norm,Dark_tau21_Norm,Dark_tau32_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,Dark_dpTFrac,Dark_Pt_Norm)
	tauplot(SoftDark_Jets,SoftDark_Tau1,SoftDark_Tau2,SoftDark_Tau3,SoftDark_tau1_Norm,
	SoftDark_tau2_Norm,SoftDark_tau3_Norm,SoftDark_tau21_Norm,SoftDark_tau32_Norm,filler2D,
	filler2D,filler2D,filler2D,filler2D,SoftDark_dpTFrac,SoftDark_Pt_Norm)

	## SD Mass plots
	SDFill(SM_SDJets,SM_SDM_Norm,SM_SDM_2D,[0]*len(SM_SDJets),SM_SDPt_Norm)
	SDFill(SMM_SDJets,SMM_SDM_Norm,SMM_SDM_2D,[0]*len(SMM_SDJets),SMM_SDPt_Norm)
	SDFill(G_SDJets,G_SDM_Norm,G_SDM_2D,G_dpTFrac,G_SDPt_Norm)
	SDFill(QM_SDJets,QM_SDM_Norm,QM_SDM_2D,QM_dpTFrac,QM_SDPt_Norm)
	SDFill(Q_SDJets,Q_SDM_Norm,Q_SDM_2D,Q_dpTFrac,Q_SDPt_Norm)
	SDFill(QM_QSDJets,QM_QSDM_Norm,QM_QSDM_2D,QM_QdpTFrac,QM_QSDPt_Norm)
	SDFill(QM_GSDJets,QM_GSDM_Norm,QM_GSDM_2D,QM_GdpTFrac,QM_GSDPt_Norm)
	SDFill(Q_GSDJets,Q_GSDM_Norm,Q_GSDM_2D,Q_GdpTFrac,Q_GSDPt_Norm)
	SDFill(G_SMSDJets,G_SMSDM_Norm,G_SMSDM_2D,G_SMdpTFrac,G_SMSDPt_Norm)
	SDFill(QM_SMSDJets,QM_SMSDM_Norm,QM_SMSDM_2D,QM_SMdpTFrac,QM_SMSDPt_Norm)
	SDFill(Q_SMSDJets,Q_SMSDM_Norm,Q_SMSDM_2D,Q_SMdpTFrac,Q_SMSDPt_Norm)
	SDFill(LD_highDFSDJets,LD_highDFSDM_Norm,LD_highDFSDM_2D,LD_highDFdpTFrac,LD_highDFSDPt_Norm)
	SDFill(LD_lowDFSDJets,LD_lowDFSDM_Norm,LD_lowDFSDM_2D,LD_lowDFdpTFrac,LD_lowDFSDPt_Norm)
	SDFill(LD_SMSDJets,LD_SMSDM_Norm,LD_SMSDM_2D,LD_SMdpTFrac,LD_SMSDPt_Norm)

	SDFill(SM_OnlySDJets,SM_OnlySDM_Norm,filler2D,[0]*len(SM_OnlySDJets),SM_OnlySDPt_Norm)
	SDFill(Mix_SDJets,Mix_SDM_Norm,filler2D,Mix_dpTFrac,Mix_SDPt_Norm)
	SDFill(DarkerMix_SDJets,DarkerMix_SDM_Norm,filler2D,DarkerMix_dpTFrac,DarkerMix_SDPt_Norm)
	SDFill(Dark_SDJets,Dark_SDM_Norm,filler2D,Dark_dpTFrac,Dark_SDPt_Norm)
	SDFill(SoftDark_SDJets,SoftDark_SDM_Norm,filler2D,SoftDark_dpTFrac,SoftDark_SDPt_Norm)

	if count == 10: # which event number to stop the code
		break

olist = [
SM_Pt,SM_Eta,SM_Phi,SM_DPhi,SM_M,
SMM_Pt,SMM_Eta,SMM_Phi,SMM_DPhi,SMM_M,
G_Pt,G_Eta,G_Phi,G_DPhi,G_M,
QM_Pt,QM_Eta,QM_Phi,QM_DPhi,QM_M,
Q_Pt,Q_Eta,Q_Phi,Q_DPhi,Q_M,
QM_QPt,QM_QEta,QM_QPhi,QM_QDPhi,QM_QM,
QM_GPt,QM_GEta,QM_GPhi,QM_GDPhi,QM_GM,
Q_GPt,Q_GEta,Q_GPhi,Q_GDPhi,Q_GM,
G_SMPt,G_SMEta,G_SMPhi,G_SMDPhi,G_SMM,
QM_SMPt,QM_SMEta,QM_SMPhi,QM_SMDPhi,QM_SMM,
Q_SMPt,Q_SMEta,Q_SMPhi,Q_SMDPhi,Q_SMM,
LD_lowDFPt,LD_lowDFEta,LD_lowDFPhi,LD_lowDFDPhi,LD_lowDFM,
LD_highDFPt,LD_highDFEta,LD_highDFPhi,LD_highDFDPhi,LD_highDFM,
LD_SMPt,LD_SMEta,LD_SMPhi,LD_SMDPhi,LD_SMM,
jet_SM_dpTFrac,jet_SMM_dpTFrac,jet_G_dpTFrac,jet_QM_dpTFrac,jet_Q_dpTFrac,
jet_QM_QdpTFrac,jet_QM_GdpTFrac,jet_Q_GdpTFrac,jet_G_SMdpTFrac,jet_QM_SMdpTFrac,
jet_Q_SMdpTFrac,jet_LD_lowDFdpTFrac,jet_LD_highDFdpTFrac,jet_LD_SMdpTFrac,SM_mult,
SMM_mult,G_mult,QM_mult,Q_mult,QM_Qmult,
QM_Gmult,Q_Gmult,G_SMmult,QM_SMmult,Q_SMmult,
LD_lowDFmult,LD_highDFmult,LD_SMmult,SM_axisMajor,SM_axisMinor,
SM_momentGirth,SM_ptD,SMM_axisMajor,SMM_axisMinor,SMM_momentGirth,
SMM_ptD,G_axisMajor,G_axisMinor,G_momentGirth,G_ptD,
QM_axisMajor,QM_axisMinor,QM_momentGirth,QM_ptD,Q_axisMajor,
Q_axisMinor,Q_momentGirth,Q_ptD,QM_QaxisMajor,QM_QaxisMinor,
QM_QmomentGirth,QM_QptD,QM_GaxisMajor,QM_GaxisMinor,QM_GmomentGirth,
QM_GptD,Q_GaxisMajor,Q_GaxisMinor,Q_GmomentGirth,Q_GptD,
G_SMaxisMajor,G_SMaxisMinor,G_SMmomentGirth,G_SMptD,QM_SMaxisMajor,
QM_SMaxisMinor,QM_SMmomentGirth,QM_SMptD,Q_SMaxisMajor,Q_SMaxisMinor,
Q_SMmomentGirth,Q_SMptD,LD_lowDFaxisMajor,LD_lowDFaxisMinor,LD_lowDFmomentGirth,
LD_lowDFptD,LD_highDFaxisMajor,LD_highDFaxisMinor,LD_highDFmomentGirth,LD_highDFptD,
LD_SMaxisMajor,LD_SMaxisMinor,LD_SMmomentGirth,LD_SMptD,SM_tau1,
SM_tau2,SM_tau3,SM_tau21,SM_tau32,SMM_tau1,
SMM_tau2,SMM_tau3,SMM_tau21,SMM_tau32,G_tau1,
G_tau2,G_tau3,G_tau21,G_tau32,QM_tau1,
QM_tau2,QM_tau3,QM_tau21,QM_tau32,Q_tau1,
Q_tau2,Q_tau3,Q_tau21,Q_tau32,QM_Qtau1,
QM_Qtau2,QM_Qtau3,QM_Qtau21,QM_Qtau32,QM_Gtau1,
QM_Gtau2,QM_Gtau3,QM_Gtau21,QM_Gtau32,Q_Gtau1,
Q_Gtau2,Q_Gtau3,Q_Gtau21,Q_Gtau32,G_SMtau1,
G_SMtau2,G_SMtau3,G_SMtau21,G_SMtau32,QM_SMtau1,
QM_SMtau2,QM_SMtau3,QM_SMtau21,QM_SMtau32,Q_SMtau1,
Q_SMtau2,Q_SMtau3,Q_SMtau21,Q_SMtau32,LD_lowDFtau1,
LD_lowDFtau2,LD_lowDFtau3,LD_lowDFtau21,LD_lowDFtau32,LD_highDFtau1,
LD_highDFtau2,LD_highDFtau3,LD_highDFtau21,LD_highDFtau32,LD_SMtau1,
LD_SMtau2,LD_SMtau3,LD_SMtau21,LD_SMtau32,SM_SDM,
SMM_SDM,G_SDM,QM_SDM,Q_SDM,QM_QSDM,
QM_GSDM,Q_GSDM,G_SMSDM,QM_SMSDM,Q_SMSDM,
LD_lowDFSDM,LD_highDFSDM,LD_SMSDM,SM_SDPt,SMM_SDPt,
G_SDPt,QM_SDPt,Q_SDPt,QM_QSDPt,QM_GSDPt,
Q_GSDPt,G_SMSDPt,QM_SMSDPt,Q_SMSDPt,LD_lowDFSDPt,
LD_highDFSDPt,LD_SMSDPt,SM_pmult,SMM_pmult,G_pmult,
QM_pmult,Q_pmult,QM_Qpmult,QM_Gpmult,Q_Gpmult,
G_SMpmult,QM_SMpmult,Q_SMpmult,LD_lowDFpmult,LD_highDFpmult,
LD_SMpmult,SM_Pt_2D,SM_Eta_2D,SM_Phi_2D,SM_DPhi_2D,
SMM_Pt_2D,SMM_Eta_2D,SMM_Phi_2D,SMM_DPhi_2D,G_Pt_2D,
G_Eta_2D,G_Phi_2D,G_DPhi_2D,QM_Pt_2D,QM_Eta_2D,
QM_Phi_2D,QM_DPhi_2D,Q_Pt_2D,Q_Eta_2D,Q_Phi_2D,
Q_DPhi_2D,QM_QPt_2D,QM_QEta_2D,QM_QPhi_2D,QM_QDPhi_2D,
QM_GPt_2D,QM_GEta_2D,QM_GPhi_2D,QM_GDPhi_2D,Q_GPt_2D,
Q_GEta_2D,Q_GPhi_2D,Q_GDPhi_2D,G_SMPt_2D,G_SMEta_2D,
G_SMPhi_2D,G_SMDPhi_2D,QM_SMPt_2D,QM_SMEta_2D,QM_SMPhi_2D,
QM_SMDPhi_2D,Q_SMPt_2D,Q_SMEta_2D,Q_SMPhi_2D,Q_SMDPhi_2D,
LD_lowDFPt_2D,LD_lowDFEta_2D,LD_lowDFPhi_2D,LD_lowDFDPhi_2D,LD_highDFPt_2D,
LD_highDFEta_2D,LD_highDFPhi_2D,LD_highDFDPhi_2D,LD_SMPt_2D,LD_SMEta_2D,
LD_SMPhi_2D,LD_SMDPhi_2D,SM_axisMajor_2D,SM_axisMinor_2D,SM_momentGirth_2D,
SM_ptD_2D,SMM_axisMajor_2D,SMM_axisMinor_2D,SMM_momentGirth_2D,SMM_ptD_2D,
G_axisMajor_2D,G_axisMinor_2D,G_momentGirth_2D,G_ptD_2D,QM_axisMajor_2D,
QM_axisMinor_2D,QM_momentGirth_2D,QM_ptD_2D,Q_axisMajor_2D,Q_axisMinor_2D,
Q_momentGirth_2D,Q_ptD_2D,QM_QaxisMajor_2D,QM_QaxisMinor_2D,QM_QmomentGirth_2D,
QM_QptD_2D,QM_GaxisMajor_2D,QM_GaxisMinor_2D,QM_GmomentGirth_2D,QM_GptD_2D,
Q_GaxisMajor_2D,Q_GaxisMinor_2D,Q_GmomentGirth_2D,Q_GptD_2D,G_SMaxisMajor_2D,
G_SMaxisMinor_2D,G_SMmomentGirth_2D,G_SMptD_2D,QM_SMaxisMajor_2D,QM_SMaxisMinor_2D,
QM_SMmomentGirth_2D,QM_SMptD_2D,Q_SMaxisMajor_2D,Q_SMaxisMinor_2D,Q_SMmomentGirth_2D,
Q_SMptD_2D,LD_lowDFaxisMajor_2D,LD_lowDFaxisMinor_2D,LD_lowDFmomentGirth_2D,LD_lowDFptD_2D,
LD_highDFaxisMajor_2D,LD_highDFaxisMinor_2D,LD_highDFmomentGirth_2D,LD_highDFptD_2D,LD_SMaxisMajor_2D,
LD_SMaxisMinor_2D,LD_SMmomentGirth_2D,LD_SMptD_2D,SM_M_2D,SMM_M_2D,
G_M_2D,QM_M_2D,Q_M_2D,QM_QM_2D,QM_GM_2D,
Q_GM_2D,G_SMM_2D,QM_SMM_2D,Q_SMM_2D,LD_lowDFM_2D,
LD_highDFM_2D,LD_SMM_2D,SM_tau1_2D,SM_tau2_2D,SM_tau3_2D,
SM_tau21_2D,SM_tau32_2D,SM_SDM_2D,SMM_tau1_2D,SMM_tau2_2D,
SMM_tau3_2D,SMM_tau21_2D,SMM_tau32_2D,SMM_SDM_2D,G_tau1_2D,
G_tau2_2D,G_tau3_2D,G_tau21_2D,G_tau32_2D,G_SDM_2D,
QM_tau1_2D,QM_tau2_2D,QM_tau3_2D,QM_tau21_2D,QM_tau32_2D,
QM_SDM_2D,Q_tau1_2D,Q_tau2_2D,Q_tau3_2D,Q_tau21_2D,
Q_tau32_2D,Q_SDM_2D,QM_Qtau1_2D,QM_Qtau2_2D,QM_Qtau3_2D,
QM_Qtau21_2D,QM_Qtau32_2D,QM_QSDM_2D,QM_Gtau1_2D,QM_Gtau2_2D,
QM_Gtau3_2D,QM_Gtau21_2D,QM_Gtau32_2D,QM_GSDM_2D,Q_Gtau1_2D,
Q_Gtau2_2D,Q_Gtau3_2D,Q_Gtau21_2D,Q_Gtau32_2D,Q_GSDM_2D,
G_SMtau1_2D,G_SMtau2_2D,G_SMtau3_2D,G_SMtau21_2D,G_SMtau32_2D,
G_SMSDM_2D,QM_SMtau1_2D,QM_SMtau2_2D,QM_SMtau3_2D,QM_SMtau21_2D,
QM_SMtau32_2D,QM_SMSDM_2D,Q_SMtau1_2D,Q_SMtau2_2D,Q_SMtau3_2D,
Q_SMtau21_2D,Q_SMtau32_2D,Q_SMSDM_2D,LD_lowDFtau1_2D,LD_lowDFtau2_2D,
LD_lowDFtau3_2D,LD_lowDFtau21_2D,LD_lowDFtau32_2D,LD_lowDFSDM_2D,LD_highDFtau1_2D,
LD_highDFtau2_2D,LD_highDFtau3_2D,LD_highDFtau21_2D,LD_highDFtau32_2D,LD_highDFSDM_2D,
LD_SMtau1_2D,LD_SMtau2_2D,LD_SMtau3_2D,LD_SMtau21_2D,LD_SMtau32_2D,
LD_SMSDM_2D,SM_pmult_2D,SMM_pmult_2D,G_pmult_2D,QM_pmult_2D,
Q_pmult_2D,QM_Qpmult_2D,QM_Gpmult_2D,Q_Gpmult_2D,G_SMpmult_2D,
QM_SMpmult_2D,Q_SMpmult_2D,LD_lowDFpmult_2D,LD_highDFpmult_2D,LD_SMpmult_2D,
SM_DPhi_PP,SMM_DPhi_PP,G_DPhi_PP,QM_DPhi_PP,
Q_DPhi_PP,QM_QDPhi_PP,QM_GDPhi_PP,Q_GDPhi_PP,G_SMDPhi_PP,
QM_SMDPhi_PP,Q_SMDPhi_PP,LD_lowDFDPhi_PP,LD_highDFDPhi_PP,LD_SMDPhi_PP,
QM_Dhi_JJ_PP,QM_Dhi_JMET_Near,QM_Dhi_JMET_Far,SM_OnlyPt,SM_OnlyEta,
SM_OnlyPhi,SM_OnlyDPhi,SM_OnlyDPhiMin,SM_OnlyM,SM_Onlymult,
SM_OnlyaxisMajor,SM_OnlyaxisMinor,SM_OnlymomentGirth,SM_OnlyptD,SM_Onlytau1,
SM_Onlytau2,SM_Onlytau3,SM_Onlytau21,SM_Onlytau32,SM_OnlySDM,
SM_OnlySDPt,SM_Onlypmult,Mix_Pt,Mix_Eta,Mix_Phi,
Mix_DPhi,Mix_DPhiMin,Mix_M,Mix_mult,Mix_axisMajor,
Mix_axisMinor,Mix_momentGirth,Mix_ptD,Mix_tau1,Mix_tau2,
Mix_tau3,Mix_tau21,Mix_tau32,Mix_SDM,Mix_SDPt,
Mix_pmult,DarkerMix_Pt,DarkerMix_Eta,DarkerMix_Phi,DarkerMix_DPhi,
DarkerMix_M,DarkerMix_mult,DarkerMix_axisMajor,DarkerMix_axisMinor,DarkerMix_momentGirth,
DarkerMix_ptD,DarkerMix_tau1,DarkerMix_tau2,DarkerMix_tau3,DarkerMix_tau21,
DarkerMix_tau32,DarkerMix_SDM,DarkerMix_SDPt,DarkerMix_pmult,Dark_Pt,
Dark_Eta,Dark_Phi,Dark_DPhi,Dark_DPhiMin,Dark_M,
Dark_mult,Dark_axisMajor,Dark_axisMinor,Dark_momentGirth,Dark_ptD,
Dark_tau1,Dark_tau2,Dark_tau3,Dark_tau21,Dark_tau32,
Dark_SDM,Dark_SDPt,Dark_pmult,SoftDark_Pt,SoftDark_Eta,
SoftDark_Phi,SoftDark_DPhi,SoftDark_M,SoftDark_mult,SoftDark_axisMajor,
SoftDark_axisMinor,SoftDark_momentGirth,SoftDark_ptD,SoftDark_tau1,SoftDark_tau2,
SoftDark_tau3,SoftDark_tau21,SoftDark_tau32,SoftDark_SDM,SoftDark_SDPt,
SoftDark_pmult,jet_SM_OnlydpTFrac,jet_Mix_dpTFrac,jet_DarkerMix_dpTFrac,jet_Dark_dpTFrac,
jet_SoftDark_dpTFrac,SM_OnlyM_Norm,SM_OnlyaxisMajor_Norm,SM_OnlyaxisMinor_Norm,SM_OnlymomentGirth_Norm,
SM_OnlyptD_Norm,SM_Onlytau1_Norm,SM_Onlytau2_Norm,SM_Onlytau3_Norm,SM_Onlytau21_Norm,
SM_Onlytau32_Norm,SM_OnlySDM_Norm,SM_Onlypmult_Norm,Mix_M_Norm,Mix_axisMajor_Norm,
Mix_axisMinor_Norm,Mix_momentGirth_Norm,Mix_ptD_Norm,Mix_tau1_Norm,Mix_tau2_Norm,
Mix_tau3_Norm,Mix_tau21_Norm,Mix_tau32_Norm,Mix_SDM_Norm,Mix_pmult_Norm,
DarkerMix_M_Norm,DarkerMix_axisMajor_Norm,DarkerMix_axisMinor_Norm,DarkerMix_momentGirth_Norm,DarkerMix_ptD_Norm,
DarkerMix_tau1_Norm,DarkerMix_tau2_Norm,DarkerMix_tau3_Norm,DarkerMix_tau21_Norm,DarkerMix_tau32_Norm,
DarkerMix_SDM_Norm,DarkerMix_pmult_Norm,Dark_M_Norm,Dark_axisMajor_Norm,Dark_axisMinor_Norm,
Dark_momentGirth_Norm,Dark_ptD_Norm,Dark_tau1_Norm,Dark_tau2_Norm,Dark_tau3_Norm,
Dark_tau21_Norm,Dark_tau32_Norm,Dark_SDM_Norm,Dark_pmult_Norm,SoftDark_M_Norm,
SoftDark_axisMajor_Norm,SoftDark_axisMinor_Norm,SoftDark_momentGirth_Norm,SoftDark_ptD_Norm,SoftDark_tau1_Norm,
SoftDark_tau2_Norm,SoftDark_tau3_Norm,SoftDark_tau21_Norm,SoftDark_tau32_Norm,SoftDark_SDM_Norm,
SoftDark_pmult_Norm,jet_SM_OnlydpTFrac_Norm,jet_Mix_dpTFrac_Norm,jet_DarkerMix_dpTFrac_Norm,jet_Dark_dpTFrac_Norm,
jet_SoftDark_dpTFrac_Norm,SM_M_Norm,SM_axisMajor_Norm,SM_axisMinor_Norm,SM_momentGirth_Norm,
SM_ptD_Norm,SM_tau1_Norm,SM_tau2_Norm,SM_tau3_Norm,SM_tau21_Norm,
SM_tau32_Norm,SM_SDM_Norm,SM_pmult_Norm,SMM_M_Norm,SMM_axisMajor_Norm,
SMM_axisMinor_Norm,SMM_momentGirth_Norm,SMM_ptD_Norm,SMM_tau1_Norm,SMM_tau2_Norm,
SMM_tau3_Norm,SMM_tau21_Norm,SMM_tau32_Norm,SMM_SDM_Norm,SMM_pmult_Norm,
G_M_Norm,G_axisMajor_Norm,G_axisMinor_Norm,G_momentGirth_Norm,G_ptD_Norm,
G_tau1_Norm,G_tau2_Norm,G_tau3_Norm,G_tau21_Norm,G_tau32_Norm,
G_SDM_Norm,G_pmult_Norm,QM_M_Norm,QM_axisMajor_Norm,QM_axisMinor_Norm,
QM_momentGirth_Norm,QM_ptD_Norm,QM_tau1_Norm,QM_tau2_Norm,QM_tau3_Norm,
QM_tau21_Norm,QM_tau32_Norm,QM_SDM_Norm,QM_pmult_Norm,Q_M_Norm,
Q_axisMajor_Norm,Q_axisMinor_Norm,Q_momentGirth_Norm,Q_ptD_Norm,Q_tau1_Norm,
Q_tau2_Norm,Q_tau3_Norm,Q_tau21_Norm,Q_tau32_Norm,Q_SDM_Norm,
Q_pmult_Norm,QM_QM_Norm,QM_QaxisMajor_Norm,QM_QaxisMinor_Norm,QM_QmomentGirth_Norm,
QM_QptD_Norm,QM_Qtau1_Norm,QM_Qtau2_Norm,QM_Qtau3_Norm,QM_Qtau21_Norm,
QM_Qtau32_Norm,QM_QSDM_Norm,QM_Qpmult_Norm,QM_GM_Norm,QM_GaxisMajor_Norm,
QM_GaxisMinor_Norm,QM_GmomentGirth_Norm,QM_GptD_Norm,QM_Gtau1_Norm,QM_Gtau2_Norm,
QM_Gtau3_Norm,QM_Gtau21_Norm,QM_Gtau32_Norm,QM_GSDM_Norm,QM_Gpmult_Norm,
Q_GM_Norm,Q_GaxisMajor_Norm,Q_GaxisMinor_Norm,Q_GmomentGirth_Norm,Q_GptD_Norm,
Q_Gtau1_Norm,Q_Gtau2_Norm,Q_Gtau3_Norm,Q_Gtau21_Norm,Q_Gtau32_Norm,
Q_GSDM_Norm,Q_Gpmult_Norm,G_SMM_Norm,G_SMaxisMajor_Norm,G_SMaxisMinor_Norm,
G_SMmomentGirth_Norm,G_SMptD_Norm,G_SMtau1_Norm,G_SMtau2_Norm,G_SMtau3_Norm,
G_SMtau21_Norm,G_SMtau32_Norm,G_SMSDM_Norm,G_SMpmult_Norm,QM_SMM_Norm,
QM_SMaxisMajor_Norm,QM_SMaxisMinor_Norm,QM_SMmomentGirth_Norm,QM_SMptD_Norm,QM_SMtau1_Norm,
QM_SMtau2_Norm,QM_SMtau3_Norm,QM_SMtau21_Norm,QM_SMtau32_Norm,QM_SMSDM_Norm,
QM_SMpmult_Norm,Q_SMM_Norm,Q_SMaxisMajor_Norm,Q_SMaxisMinor_Norm,Q_SMmomentGirth_Norm,
Q_SMptD_Norm,Q_SMtau1_Norm,Q_SMtau2_Norm,Q_SMtau3_Norm,Q_SMtau21_Norm,
Q_SMtau32_Norm,Q_SMSDM_Norm,Q_SMpmult_Norm,LD_highDFM_Norm,LD_highDFaxisMajor_Norm,
LD_highDFaxisMinor_Norm,LD_highDFmomentGirth_Norm,LD_highDFptD_Norm,LD_highDFtau1_Norm,LD_highDFtau2_Norm,
LD_highDFtau3_Norm,LD_highDFtau21_Norm,LD_highDFtau32_Norm,LD_highDFSDM_Norm,LD_highDFpmult_Norm,
LD_lowDFM_Norm,LD_lowDFaxisMajor_Norm,LD_lowDFaxisMinor_Norm,LD_lowDFmomentGirth_Norm,LD_lowDFptD_Norm,
LD_lowDFtau1_Norm,LD_lowDFtau2_Norm,LD_lowDFtau3_Norm,LD_lowDFtau21_Norm,LD_lowDFtau32_Norm,
LD_lowDFSDM_Norm,LD_lowDFpmult_Norm,LD_SMM_Norm,LD_SMaxisMajor_Norm,LD_SMaxisMinor_Norm,
LD_SMmomentGirth_Norm,LD_SMptD_Norm,LD_SMtau1_Norm,LD_SMtau2_Norm,LD_SMtau3_Norm,
LD_SMtau21_Norm,LD_SMtau32_Norm,LD_SMSDM_Norm,LD_SMpmult_Norm,jet_SM_dpTFrac_Norm,
jet_SMM_dpTFrac_Norm,jet_G_dpTFrac_Norm,jet_QM_dpTFrac_Norm,jet_Q_dpTFrac_Norm,jet_QM_QdpTFrac_Norm,
jet_QM_GdpTFrac_Norm,jet_Q_GdpTFrac_Norm,jet_G_SMdpTFrac_Norm,jet_QM_SMdpTFrac_Norm,jet_Q_SMdpTFrac_Norm,
jet_LD_highDFdpTFrac_Norm,jet_LD_lowDFdpTFrac_Norm,jet_LD_SMdpTFrac_Norm,SM_OnlyPt_Cl,SM_OnlyEta_Cl,
SM_OnlyPhi_Cl,SM_OnlyDPhi_Cl,SM_OnlyM_Cl,SM_Onlymult_Cl,SM_OnlyaxisMajor_Cl,
SM_OnlyaxisMinor_Cl,SM_OnlymomentGirth_Cl,SM_OnlyptD_Cl,SM_Onlytau1_Cl,SM_Onlytau2_Cl,
SM_Onlytau3_Cl,SM_Onlytau21_Cl,SM_Onlytau32_Cl,SM_OnlySDM_Cl,SM_Onlypmult_Cl,
Mix_Pt_Cl,Mix_Eta_Cl,Mix_Phi_Cl,Mix_DPhi_Cl,Mix_M_Cl,
Mix_mult_Cl,Mix_axisMajor_Cl,Mix_axisMinor_Cl,Mix_momentGirth_Cl,Mix_ptD_Cl,
Mix_tau1_Cl,Mix_tau2_Cl,Mix_tau3_Cl,Mix_tau21_Cl,Mix_tau32_Cl,
Mix_SDM_Cl,Mix_pmult_Cl,Dark_Pt_Cl,Dark_Eta_Cl,Dark_Phi_Cl,
Dark_DPhi_Cl,Dark_M_Cl,Dark_mult_Cl,Dark_axisMajor_Cl,Dark_axisMinor_Cl,
Dark_momentGirth_Cl,Dark_ptD_Cl,Dark_tau1_Cl,Dark_tau2_Cl,Dark_tau3_Cl,
Dark_tau21_Cl,Dark_tau32_Cl,Dark_SDM_Cl,Dark_pmult_Cl,All_DPhiMin,
All_Pt_Cl,All_Eta_Cl,All_Phi_Cl,All_DPhi_Cl,All_M_Cl,
All_mult_Cl,All_axisMajor_Cl,All_axisMinor_Cl,All_momentGirth_Cl,All_ptD_Cl,
All_tau1_Cl,All_tau2_Cl,All_tau3_Cl,All_tau21_Cl,All_tau32_Cl,
All_SDM_Cl,All_pmult_Cl,All_Pt_j1,All_Eta_j1,All_Phi_j1,
All_DPhi_j1,All_M_j1,All_mult_j1,All_axisMajor_j1,All_axisMinor_j1,
All_momentGirth_j1,All_ptD_j1,All_tau1_j1,All_tau2_j1,All_tau3_j1,
All_tau21_j1,All_tau32_j1,All_SDM_j1,All_pmult_j1,All_Pt_j2,
All_Eta_j2,All_Phi_j2,All_DPhi_j2,All_M_j2,All_mult_j2,
All_axisMajor_j2,All_axisMinor_j2,All_momentGirth_j2,All_ptD_j2,All_tau1_j2,
All_tau2_j2,All_tau3_j2,All_tau21_j2,All_tau32_j2,All_SDM_j2,
All_pmult_j2,All_Pt_j3,All_Eta_j3,All_Phi_j3,All_DPhi_j3,
All_M_j3,All_mult_j3,All_axisMajor_j3,All_axisMinor_j3,All_momentGirth_j3,
All_ptD_j3,All_tau1_j3,All_tau2_j3,All_tau3_j3,All_tau21_j3,
All_tau32_j3,All_SDM_j3,All_pmult_j3,All_Pt_j4,All_Eta_j4,
All_Phi_j4,All_DPhi_j4,All_M_j4,All_mult_j4,All_axisMajor_j4,
All_axisMinor_j4,All_momentGirth_j4,All_ptD_j4,All_tau1_j4,All_tau2_j4,
All_tau3_j4,All_tau21_j4,All_tau32_j4,All_SDM_j4,All_pmult_j4,
All_Pt_j4p,All_Eta_j4p,All_Phi_j4p,All_DPhi_j4p,All_M_j4p,
All_mult_j4p,All_axisMajor_j4p,All_axisMinor_j4p,All_momentGirth_j4p,All_ptD_j4p,
All_tau1_j4p,All_tau2_j4p,All_tau3_j4p,All_tau21_j4p,All_tau32_j4p,
All_SDM_j4p,All_pmult_j4p,ST,MET,ST_2jet,ST_3jet,ST_4jet,
ST_4plusjet,RT,MT_2j,MT_3j,MT_4j,
MT_4pj,MT_m2000_DR_f2,MT_l2000_DR_f2,MT_m2000_DP_f2,MT_l2000_DP_f2,
MT_m2000_DE_f2,MT_l2000_DE_f2,MT_R,MT_W,MT_R_DR_f2,
MT_W_DR_f2,MT_R_DP_f2,MT_W_DP_f2,MT_R_DE_f2,MT_W_DE_f2,
MT2_R,MT2_W,MT2_DeltaR_R,MT2_DeltaR_W,MT2_DeltaP_R,
MT2_DeltaP_W,MT2_DeltaE_R,MT2_DeltaE_W,MT2_mp,MT2_mp_0M,
MT2_mp_1M,MT2_mp_2M,MT2_R_MDiff,MT2_W_MDiff,MT2_R_M,
MT2_W_M,MT2_mp_MT2_R,MT2_mp_MDiff,MT2_R_Pt,MT2_mp_Pt,MT2_R_DPhiMET,MT2_mp_DPhiMET,
MT2_R_DPhiMET_d0r,MT2_R_DPhiMET_s0r,MT2_R_DPhiMET_d1r,MT2_R_DPhiMET_s1r,MT2_R_DPhi_pair,
MT2_mp_DPhiPair,MT2_mp_PhiReq,MT2_mp_f4,MT2_mp_DPhiPair_f4,MT2_mp_PhiReq_f4,MT2_mp_DPhiPair_Pt_f4,
MT2_Pt_mp,MT2_Pt_mp_DPhiPair,MT2_mp_DPhiPair_Pt
]

rootOutFile = rt.TFile.Open("Skims/s2_mMed-{}_mDark-{}_rinv-{}_alpha-{}.root".format(mMed,mDark,rInv,Alpha), "recreate")
rootOutFile.cd()
for histo in olist:
	histo.Write()
rootOutFile.Close()
