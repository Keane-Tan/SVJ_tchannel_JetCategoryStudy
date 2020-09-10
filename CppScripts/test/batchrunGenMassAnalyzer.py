import os

for ith in range(1,101):
    print "Running File " + str(ith)
    cmsRunCom = 'cmsRun runSVJ.py config=genmassanalyzer_cfg output=TFileService outpre=genmassanalysis inpre=root://cmseos.fnal.gov//store/user/keanet/tchannel/GEN/step1_LHE-GEN mMediator=3000.0 mDark=20.0 rinv=0.3 alpha=peak part='+ str(ith) +' channel="t" maxEvents=1000 madgraph=True'
    os.system(cmsRunCom)

outfile = 'genmassanalysis_t-channel_mMed-3000_mDark-20_rinv-0.3_alpha-peak_13TeV-madgraphMLM-pythia8_n-1000_part-'
resultfile = "skims_softDrop_tau"

startlist = []

def hadder(start,end):
    startlist.append(start)
    haddCom = "hadd -f " + resultfile + str(start) + ".root "
    for ith in range(start,end):
        haddCom += outfile + str(ith) + ".root "
    print haddCom
    print (" ")
    os.system(haddCom)

hadder(1,11)
hadder(11,21)
hadder(21,31)
hadder(31,41)
hadder(41,51)
hadder(51,61)
hadder(61,71)
hadder(71,81)
hadder(81,91)
hadder(91,101)

haddCom = "hadd -f " + resultfile + ".root "
for start in startlist:
    haddCom += resultfile + str(start) + ".root "

os.system(haddCom)
