# include plotting for the 5 general categories
# plotting stuff made by v5_02_MakeHist.py
# calculate root mean squared distance between general categories and specific categories

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
		if var == "dpTFrac" and (categ[1] == "SM_" or categ[1] == "SMM_" or categ[1] == "SM_Only"):
			continue
		_filet.GetObject(categ[1]+var,categ[0])
		categ[0] = categ[0].Clone(categ[0].GetName()+"_")
		categ[0].SetDirectory(0)
		# norm(categ[0])

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
		# if var == "dpTFrac" and (ca[1] == "SM_" or ca[1] == "SMM_" or ca[1] == "SM_Only"):
		# 	continue
		# norm(ca[0])
		ca[0].Draw("hist SAME")

	leg = rt.TLegend(0.73,0.45,0.93,0.9)
	leg.SetTextFont(42)

	# Mean Squared Distance
	SpecHist = allC[1][0]
	GenHist = allC[0][0]
	DiffHist = GenHist.Clone(GenHist.GetName()+"_Diff")
	DiffHist.Add(SpecHist,-1)
	numBins = DiffHist.GetNbinsX()

	sumOfDiff = 0
	sumOfBins = 0
	for i in range(numBins):
		# if CombinedHist.GetBinContent(i) > 10:
		sumOfDiff += DiffHist.GetBinContent(i)**2
		sumOfBins += 1
	rmsd = np.sqrt(sumOfDiff)/sumOfBins
	rmsVals.append(rmsd)

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
	# leg.Draw("SAME")

	if var == "mult":
		firstHist.GetXaxis().SetTitle("Number of Jets")
	firstHist.GetXaxis().SetTitleSize(30)
	firstHist.GetXaxis().SetTitleFont(43)
	firstHist.GetXaxis().SetLabelFont(43)
	firstHist.GetXaxis().SetLabelSize(25)
	firstHist.GetYaxis().SetTitleSize(30)
	firstHist.GetYaxis().SetTitleFont(43)
	firstHist.GetYaxis().SetTitleOffset(1.35)
	firstHist.GetYaxis().SetLabelFont(43)
	firstHist.GetYaxis().SetLabelSize(25)
	firstHist.SetMinimum(10)
	firstHist.SetMaximum(ymax*1.1)
	# firstHist.SetMinimum(0.001) # for normalized plots
	# firstHist.SetMaximum(ymax*1.1) # for normalized plots

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

	c.SaveAs(saveDir+alabel+"_"+cleanStr(bid)+"_"+var+"_DS.png")

# # benchmark
bid = "s2_mMed-3000_mDark-20_rinv-0.3_alpha-peak"

# dShistLimsDict = [
# "dpTFrac",
# "Pt",
# "Eta",
# "Phi",
# "DPhi",
# "mult"
# ]

dShistLimsDict = [
# "dpTFrac",
# "Pt",
# "Eta",
# "Phi",
# "DPhi",
# "mult",
# "M",
# 'axisMajor',
# 'axisMinor',
# 'momentGirth',
# 'ptD',
# 'pmult',
# 'tau1',
# 'tau2',
# 'tau3',
# 'tau21',
# 'tau32',
# 'SDM',
'dpTFrac_Norm',
"M_Norm",
'pmult_Norm',
'momentGirth_Norm',
'axisMajor_Norm',
'axisMinor_Norm',
'ptD_Norm',
'tau1_Norm',
'tau2_Norm',
'tau3_Norm',
'tau21_Norm',
'tau32_Norm',
'SDM_Norm'
]

saveDir = "../plots/Draw/Test/"

rmsVals = []

######## No LD (make it easier to compare the graphs) ##########
# 1D plots
for var in dShistLimsDict:
	# Dark = [rt.TH1F(),"Dark",1,7]
	SM_ = [rt.TH1F(),"SM_",1,1]
	SMM_ = [rt.TH1F(),"SMM_",1,rt.kCyan+2]
	G_ = [rt.TH1F(),"G_",2,1]
	QM_ = [rt.TH1F(),"QM_",1,2]
	Q_ = [rt.TH1F(),"Q_",1,3] #
	QM_Q = [rt.TH1F(),"QM_Q",1,rt.kCyan+2] #
	QM_G = [rt.TH1F(),"QM_G",2,rt.kOrange+7]
	Q_G = [rt.TH1F(),"Q_G",2,6]
	G_SM = [rt.TH1F(),"G_SM",2,2]
	QM_SM = [rt.TH1F(),"QM_SM",1,3] #
	Q_SM = [rt.TH1F(),"Q_SM",2,4]
	LD_lowDF = [rt.TH1F(),"LD_lowDF",1,8] #
	LD_highDF = [rt.TH1F(),"LD_highDF",1,4]
	LD_SM = [rt.TH1F(),"LD_SM",1,6]

	SM_Only = [rt.TH1F(),"SM_Only",1,1]
	Mix_ = [rt.TH1F(),"Mix_",1,2]
	DarkerMix_ = [rt.TH1F(),"DarkerMix_",1,2]
	Dark_ = [rt.TH1F(),"Dark_",1,6]
	SoftDark_ = [rt.TH1F(),"SoftDark_",1,4]

	allCat = [SM_, SMM_, G_, QM_, Q_, QM_Q,
	QM_G, Q_G, G_SM, QM_SM, Q_SM, LD_lowDF, LD_highDF, LD_SM,
	SM_Only,Mix_,DarkerMix_,Dark_,SoftDark_]

	Dark_SM_ = [Dark_, SM_]
	Mix_SM_ = [Mix_, SM_]
	DarkerMix_SM_ = [DarkerMix_, SM_]
	SM_OnlySM_ = [SM_Only, SM_]
	SoftDark_SM_ = [SoftDark_, SM_]

	Dark_SMM_ = [Dark_, SMM_]
	Mix_SMM_ = [Mix_, SMM_]
	DarkerMix_SMM_ = [DarkerMix_, SMM_]
	SM_OnlySMM_ = [SM_Only, SMM_]
	SoftDark_SMM_ = [SoftDark_, SMM_]

	Dark_G_ = [Dark_, G_]
	Mix_G_ = [Mix_, G_]
	DarkerMix_G_ = [DarkerMix_, G_]
	SM_OnlyG_ = [SM_Only, G_]
	SoftDark_G_ = [SoftDark_, G_]

	Dark_QM_ = [Dark_, QM_]
	Mix_QM_ = [Mix_, QM_]
	DarkerMix_QM_ = [DarkerMix_, QM_]
	SM_OnlyQM_ = [SM_Only, QM_]
	SoftDark_QM_ = [SoftDark_, QM_]

	Dark_Q_ = [Dark_, Q_]
	Mix_Q_ = [Mix_, Q_]
	DarkerMix_Q_ = [DarkerMix_, Q_]
	SM_OnlyQ_ = [SM_Only, Q_]
	SoftDark_Q_ = [SoftDark_, Q_]

	Dark_QM_Q = [Dark_, QM_Q]
	Mix_QM_Q = [Mix_, QM_Q]
	DarkerMix_QM_Q = [DarkerMix_, QM_Q]
	SM_OnlyQM_Q = [SM_Only, QM_Q]
	SoftDark_QM_Q = [SoftDark_, QM_Q]

	Dark_QM_G = [Dark_, QM_G]
	Mix_QM_G = [Mix_, QM_G]
	DarkerMix_QM_G = [DarkerMix_, QM_G]
	SM_OnlyQM_G = [SM_Only, QM_G]
	SoftDark_QM_G = [SoftDark_, QM_G]

	Dark_Q_G = [Dark_, Q_G]
	Mix_Q_G = [Mix_, Q_G]
	DarkerMix_Q_G = [DarkerMix_, Q_G]
	SM_OnlyQ_G = [SM_Only, Q_G]
	SoftDark_Q_G = [SoftDark_, Q_G]

	Dark_G_SM = [Dark_, G_SM]
	Mix_G_SM = [Mix_, G_SM]
	DarkerMix_G_SM = [DarkerMix_, G_SM]
	SM_OnlyG_SM = [SM_Only, G_SM]
	SoftDark_G_SM = [SoftDark_, G_SM]

	Dark_QM_SM = [Dark_, QM_SM]
	Mix_QM_SM = [Mix_, QM_SM]
	DarkerMix_QM_SM = [DarkerMix_, QM_SM]
	SM_OnlyQM_SM = [SM_Only, QM_SM]
	SoftDark_QM_SM = [SoftDark_, QM_SM]

	Dark_Q_SM = [Dark_, Q_SM]
	Mix_Q_SM = [Mix_, Q_SM]
	DarkerMix_Q_SM = [DarkerMix_, Q_SM]
	SM_OnlyQ_SM = [SM_Only, Q_SM]
	SoftDark_Q_SM = [SoftDark_, Q_SM]

	Dark_LD_lowDF = [Dark_, LD_lowDF]
	Mix_LD_lowDF = [Mix_, LD_lowDF]
	DarkerMix_LD_lowDF = [DarkerMix_, LD_lowDF]
	SM_OnlyLD_lowDF = [SM_Only, LD_lowDF]
	SoftDark_LD_lowDF = [SoftDark_, LD_lowDF]

	Dark_LD_highDF = [Dark_, LD_highDF]
	Mix_LD_highDF = [Mix_, LD_highDF]
	DarkerMix_LD_highDF = [DarkerMix_, LD_highDF]
	SM_OnlyLD_highDF = [SM_Only, LD_highDF]
	SoftDark_LD_highDF = [SoftDark_, LD_highDF]

	Dark_LD_SM = [Dark_, LD_SM]
	Mix_LD_SM = [Mix_, LD_SM]
	DarkerMix_LD_SM = [DarkerMix_, LD_SM]
	SM_OnlyLD_SM = [SM_Only, LD_SM]
	SoftDark_LD_SM = [SoftDark_, LD_SM]

	_filet = rt.TFile.Open("Skims/"+bid+".root","read")
	grabHist(allCat)
	_filet.Close()

	Mix_Q_hist = Mix_[0].Clone(Mix_[0].GetName()+"_Mix")
	Mix_Q_hist.Add(Q_[0])
	Mix_Q = [Mix_Q_hist,"Mix_Q",1,1]

	Mix_Q_Com = [Mix_Q,Mix_,Q_]

	### Can only draw so many plots at a time
	## Only uncomment one line at a time
	## Otherwise, will get break segmentation error if run too many at a time.

	# drawHist(Dark_SM_,"Dark_SM_")
	# drawHist(Mix_SM_,"Mix_SM_")
	# drawHist(DarkerMix_SM_,"DarkerMix_SM_")
	# drawHist(SM_OnlySM_,"SM_OnlySM_")
	# drawHist(SoftDark_SM_,"SoftDark_SM_")

	# drawHist(Dark_SMM_,"Dark_SMM_")
	# drawHist(Mix_SMM_,"Mix_SMM_")
	# drawHist(DarkerMix_SMM_,"DarkerMix_SMM_")
	# drawHist(SM_OnlySMM_,"SM_OnlySMM_")
	# drawHist(SoftDark_SMM_,"SoftDark_SMM_")

	# drawHist(Dark_G_,"Dark_G_")
	# drawHist(Mix_G_,"Mix_G_")
	# drawHist(DarkerMix_G_,"DarkerMix_G_")
	# drawHist(SM_OnlyG_,"SM_OnlyG_")
	# drawHist(SoftDark_G_,"SoftDark_G_")

	# drawHist(Dark_QM_,"Dark_QM_")
	# drawHist(Mix_QM_,"Mix_QM_")
	# drawHist(DarkerMix_QM_,"DarkerMix_QM_")
	# drawHist(SM_OnlyQM_,"SM_OnlyQM_")
	# drawHist(SoftDark_QM_,"SoftDark_QM_")

	# drawHist(Dark_Q_,"Dark_Q_")
	# drawHist(Mix_Q_,"Mix_Q_")
	# drawHist(DarkerMix_Q_,"DarkerMix_Q_")
	# drawHist(SM_OnlyQ_,"SM_OnlyQ_")
	# drawHist(SoftDark_Q_,"SoftDark_Q_")

	# drawHist(Dark_QM_Q,"Dark_QM_Q")
	# drawHist(Mix_QM_Q,"Mix_QM_Q")
	# drawHist(DarkerMix_QM_Q,"DarkerMix_QM_Q")
	# drawHist(SM_OnlyQM_Q,"SM_OnlyQM_Q")
	# drawHist(SoftDark_QM_Q,"SoftDark_QM_Q")

	# drawHist(Dark_QM_G,"Dark_QM_G")
	# drawHist(Mix_QM_G,"Mix_QM_G")
	# drawHist(DarkerMix_QM_G,"DarkerMix_QM_G")
	# drawHist(SM_OnlyQM_G,"SM_OnlyQM_G")
	# drawHist(SoftDark_QM_G,"SoftDark_QM_G")

	# drawHist(Dark_Q_G,"Dark_Q_G")
	# drawHist(Mix_Q_G,"Mix_Q_G")
	# drawHist(DarkerMix_Q_G,"DarkerMix_Q_G")
	# drawHist(SM_OnlyQ_G,"SM_OnlyQ_G")
	# drawHist(SoftDark_Q_G,"SoftDark_Q_G")

	# drawHist(Dark_G_SM,"Dark_G_SM")
	# drawHist(Mix_G_SM,"Mix_G_SM")
	# drawHist(DarkerMix_G_SM,"DarkerMix_G_SM")
	# drawHist(SM_OnlyG_SM,"SM_OnlyG_SM")
	# drawHist(SoftDark_G_SM,"SoftDark_G_SM")

	# drawHist(Dark_QM_SM,"Dark_QM_SM")
	# drawHist(Mix_QM_SM,"Mix_QM_SM")
	# drawHist(DarkerMix_QM_SM,"DarkerMix_QM_SM")
	# drawHist(SM_OnlyQM_SM,"SM_OnlyQM_SM")
	# drawHist(SoftDark_QM_SM,"SoftDark_QM_SM")

	# drawHist(Dark_Q_SM,"Dark_Q_SM")
	# drawHist(Mix_Q_SM,"Mix_Q_SM")
	# drawHist(DarkerMix_Q_SM,"DarkerMix_Q_SM")
	# drawHist(SM_OnlyQ_SM,"SM_OnlyQ_SM")
	# drawHist(SoftDark_Q_SM,"SoftDark_Q_SM")

	# drawHist(Dark_LD_lowDF,"Dark_LD_lowDF")
	# drawHist(Mix_LD_lowDF,"Mix_LD_lowDF")
	# drawHist(DarkerMix_LD_lowDF,"DarkerMix_LD_lowDF")
	# drawHist(SM_OnlyLD_lowDF,"SM_OnlyLD_lowDF")
	# drawHist(SoftDark_LD_lowDF,"SoftDark_LD_lowDF")

	# drawHist(Dark_LD_highDF,"Dark_LD_highDF")
	# drawHist(Mix_LD_highDF,"Mix_LD_highDF")
	# drawHist(DarkerMix_LD_highDF,"DarkerMix_LD_highDF")
	# drawHist(SM_OnlyLD_highDF,"SM_OnlyLD_highDF")
	# drawHist(SoftDark_LD_highDF,"SoftDark_LD_highDF")

	drawHist(Dark_LD_SM,"Dark_LD_SM")
	drawHist(Mix_LD_SM,"Mix_LD_SM")
	drawHist(DarkerMix_LD_SM,"DarkerMix_LD_SM")
	drawHist(SM_OnlyLD_SM,"SM_OnlyLD_SM")
	drawHist(SoftDark_LD_SM,"SoftDark_LD_SM")

print "LD_SM"
cc = 1
sen = ""
for val in rmsVals:
	if cc % 5 == 0:
		sen += str(val) + "\n"
	else:
		sen += str(val) + "\t"
	cc += 1

print sen
