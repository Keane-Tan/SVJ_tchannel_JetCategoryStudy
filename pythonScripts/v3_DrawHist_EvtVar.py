import ROOT as rt
import sys
import numpy as np
import CMS_lumi, tdrstyle

tdrstyle.setTDRStyle()

#rt.gStyle.SetPalette(54,0)
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

def norm(hist):
	if hist.Integral(0,hist.GetNbinsX()) > 0:
		hist.Scale(1.0/hist.Integral(0,hist.GetNbinsX())) # normalizing the histograms

def cleanStr(ast):
	dump = ["s2","mMed-","mDark-","rinv-","alpha-"]
	s = ast
	for ch in dump:
	    s = s.replace(ch,"")
	print s
	return s

def grabHist(allCat):
	for categ in allCat:
		_filet.GetObject(categ[1],categ[0])
		categ[0] = categ[0].Clone(categ[0].GetName()+"_")
		categ[0].SetDirectory(0)
		norm(categ[0])

def grabHist_File(categ,_filet):
	_filet.GetObject("ST",categ[0])
	categ[0] = categ[0].Clone(categ[0].GetName()+"_")
	categ[0].SetDirectory(0)
	norm(categ[0])

def drawHist(allC,alabel):

	H=700
	W=700

	H_ref = 700
	W_ref = 700

	T = 0.08*H_ref
	B = 0.12*H_ref
	L = 0.12*W_ref
	R = 0.08*W_ref

	c = rt.TCanvas("c","canvas",0,0,W,H)
	c.SetFillColor(0);
	c.SetBorderMode(0);
	c.SetFrameFillStyle(0);
	c.SetFrameBorderMode(0);
	#c.SetLeftMargin( 0.15 );#L/W
	c.SetRightMargin(0.04);
	#c.SetRightMargin( R/W );
	c.SetTopMargin( T/H );
	#c.SetBottomMargin(0.2);
	c.SetTickx(0);
	c.SetTicky(0);
	c.Draw()
	c.SetLogy()

	pad1 = rt.TPad("pad1", "pad1", 0, 0, 1, 0.95)

	#pad1.SetGridx()         # Vertical grid
	#pad1.SetGridy()         # Horizontal grid
	pad1.SetLogy()
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad

	firstHist = allC[0][0]
	firstHist.Draw("hist")
	for ca in allC:
		ca[0].Draw("hist SAME")

	leg = rt.TLegend(0.73,0.65,0.93,0.9)
	if "DP" in alabel or "DeltaP" in alabel:
		leg = rt.TLegend(0.2,0.65,0.35,0.9)

	leg.SetTextFont(42)

	ymaxList = []

	for ci in range(len(allC)):
		chist = allC[ci]
		chist[0].SetLineColor(chist[3])
		chist[0].SetLineStyle(chist[2])
		chist[0].SetLineWidth(2)
		ymaxList.append(chist[0].GetMaximum())
		leg.AddEntry(chist[0],chist[1],"l")

	ymax = max(ymaxList)
	leg.SetTextSize(.03)
	leg.Draw("SAME")

	firstHist.GetXaxis().SetTitleSize(30)
	firstHist.GetXaxis().SetTitleFont(43)
	firstHist.GetXaxis().SetLabelFont(43)
	firstHist.GetXaxis().SetLabelSize(25)
	firstHist.GetYaxis().SetTitleSize(30)
	firstHist.GetYaxis().SetTitleFont(43)
	firstHist.GetYaxis().SetTitleOffset(1.35)
	firstHist.GetYaxis().SetLabelFont(43)
	firstHist.GetYaxis().SetLabelSize(25)
	# firstHist.SetMinimum(10)
	# firstHist.SetMaximum(ymax*1.1)
	firstHist.SetMinimum(0.001) # for normalized plots
	firstHist.SetMaximum(ymax*1.1)

	c.Update()

	# CMS style
	CMS_lumi.lumi_sqrtS = "(13 TeV)"
	CMS_lumi.extraText   = "  Simulation"

	iPeriod = 0
	iPos = 0

	CMS_lumi.CMS_lumi(c, iPeriod, iPos)
	c.cd()
	c.Update();
	c.RedrawAxis()

	c.SaveAs(saveDir+"EvtVar_"+cleanStr(bid)+"_"+alabel+"_DS.png")

# # benchmark
bid = "s2_mMed-3000_mDark-20_rinv-0.3_alpha-peak"
idM4000 = "s2_mMed-4000_mDark-20_rinv-0.3_alpha-peak"
idM2000 = "s2_mMed-2000_mDark-20_rinv-0.3_alpha-peak"
saveDir = "../plots/Draw/Test/"

# 1D plots

# ST variables
ST = [rt.TH1F(),"ST",1,1]
MET = [rt.TH1F(),"MET",1,rt.kCyan+2]
ST_2jet = [rt.TH1F(),"ST_2jet",1,2]
ST_3jet = [rt.TH1F(),"ST_3jet",1,3]
ST_4jet = [rt.TH1F(),"ST_4jet",1,rt.kOrange+7]
ST_4plusjet = [rt.TH1F(),"ST_4plusjet",1,6]

ST_2000_20_0p3_peak = [rt.TH1F(),"ST_2000_20_0p3_peak",1,2]
ST_3000_20_0p3_peak = [rt.TH1F(),"ST_3000_20_0p3_peak",1,1]
ST_4000_20_0p3_peak = [rt.TH1F(),"ST_4000_20_0p3_peak",1,3]

ST_MET = [ST,MET]
ST_nJet = [ST_2jet,ST_3jet,ST_4jet,ST_4plusjet,MET]
ST_Mass = [ST_2000_20_0p3_peak,ST_3000_20_0p3_peak,ST_4000_20_0p3_peak,MET]

# MT variables
MT_2j = [rt.TH1F(),"MT_2j",1,1]
MT_3j = [rt.TH1F(),"MT_3j",1,rt.kCyan+2]
MT_4j = [rt.TH1F(),"MT_4j",1,2]
MT_4pj = [rt.TH1F(),"MT_4pj",1,3]

MT_R = [rt.TH1F(),"MT_R",1,1]
MT_W = [rt.TH1F(),"MT_W",1,rt.kCyan+2]

MT_m2000_DR_f2 = [rt.TH1F(),"MT_m2000_DR_f2",1,1]
MT_l2000_DR_f2 = [rt.TH1F(),"MT_l2000_DR_f2",1,rt.kCyan+2]

MT_m2000_DP_f2 = [rt.TH1F(),"MT_m2000_DP_f2",1,1]
MT_l2000_DP_f2 = [rt.TH1F(),"MT_l2000_DP_f2",1,rt.kCyan+2]

MT_m2000_DE_f2 = [rt.TH1F(),"MT_m2000_DE_f2",1,1]
MT_l2000_DE_f2 = [rt.TH1F(),"MT_l2000_DE_f2",1,rt.kCyan+2]

MT_R_DR_f2 = [rt.TH1F(),"MT_R_DR_f2",1,1]
MT_W_DR_f2 = [rt.TH1F(),"MT_W_DR_f2",1,rt.kCyan+2]

MT_R_DP_f2 = [rt.TH1F(),"MT_R_DP_f2",1,1]
MT_W_DP_f2 = [rt.TH1F(),"MT_W_DP_f2",1,rt.kCyan+2]

MT_R_DE_f2 = [rt.TH1F(),"MT_R_DE_f2",1,1]
MT_W_DE_f2 = [rt.TH1F(),"MT_W_DE_f2",1,rt.kCyan+2]

MT_nj = [MT_2j,MT_3j,MT_4j,MT_4pj]
MT_RW = [MT_R,MT_W]
MT_2000_DR = [MT_m2000_DR_f2,MT_l2000_DR_f2]
MT_2000_DP = [MT_m2000_DP_f2,MT_l2000_DP_f2]
MT_2000_DE = [MT_m2000_DE_f2,MT_l2000_DE_f2]
MT_RW_DR = [MT_R_DR_f2,MT_W_DR_f2]
MT_RW_DP = [MT_R_DP_f2,MT_W_DP_f2]
MT_RW_DE = [MT_R_DE_f2,MT_W_DE_f2]

# MT2 variables
MT2_R = [rt.TH1F(),"MT2_R",1,1]
MT2_W = [rt.TH1F(),"MT2_W",1,rt.kCyan+2]

MT2_DeltaR_R = [rt.TH1F(),"MT2_DeltaR_R",1,1]
MT2_DeltaR_W = [rt.TH1F(),"MT2_DeltaR_W",1,rt.kCyan+2]

MT2_DeltaP_R = [rt.TH1F(),"MT2_DeltaP_R",1,1]
MT2_DeltaP_W = [rt.TH1F(),"MT2_DeltaP_W",1,rt.kCyan+2]

MT2_DeltaE_R = [rt.TH1F(),"MT2_DeltaE_R",1,1]
MT2_DeltaE_W = [rt.TH1F(),"MT2_DeltaE_W",1,rt.kCyan+2]

MT2_RW = [MT2_R,MT2_W]
MT2_DeltaR_RW = [MT2_DeltaR_R,MT2_DeltaR_W]
MT2_DeltaP_RW = [MT2_DeltaP_R,MT2_DeltaP_W]
MT2_DeltaE_RW = [MT2_DeltaE_R,MT2_DeltaE_W]

_filet = rt.TFile.Open("Skims/"+bid+".root","read")
grabHist(ST_MET)
grabHist(ST_nJet)
grabHist(MT_nj)
grabHist(MT_RW)
grabHist(MT_2000_DR)
grabHist(MT_2000_DP)
grabHist(MT_2000_DE)
grabHist(MT_RW_DR)
grabHist(MT_RW_DP)
grabHist(MT_RW_DE)
grabHist(MT2_RW)
grabHist(MT2_DeltaR_RW)
grabHist(MT2_DeltaP_RW)
grabHist(MT2_DeltaE_RW)
grabHist_File(ST_3000_20_0p3_peak,_filet)
_filet.Close()

_filet_m2000 = rt.TFile.Open("Skims/"+idM2000+".root","read")
grabHist_File(ST_2000_20_0p3_peak,_filet_m2000)
_filet_m2000.Close()

_filet_m4000 = rt.TFile.Open("Skims/"+idM4000+".root","read")
grabHist_File(ST_4000_20_0p3_peak,_filet_m4000)
_filet_m4000.Close()

drawHist(ST_MET,"ST_MET")
drawHist(ST_nJet,"ST_nJet")
drawHist(MT_nj,"MT_nj")
drawHist(MT_RW,"MT_RW")
drawHist(MT_2000_DR,"MT_2000_DR")
drawHist(MT_2000_DP,"MT_2000_DP")
drawHist(MT_2000_DE,"MT_2000_DE")
drawHist(MT_RW_DR,"MT_RW_DR")
drawHist(MT_RW_DP,"MT_RW_DP")
drawHist(MT_RW_DE,"MT_RW_DE")
drawHist(MT2_RW,"MT2_RW")
drawHist(MT2_DeltaR_RW,"MT2_DeltaR_RW")
drawHist(MT2_DeltaP_RW,"MT2_DeltaP_RW")
drawHist(MT2_DeltaE_RW,"MT2_DeltaE_RW")
drawHist(ST_Mass,"ST_Mass")
