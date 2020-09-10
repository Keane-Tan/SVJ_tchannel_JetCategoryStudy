cd ../plugins
scram b
cd ../test
cmsRun runSVJ.py config=genmassanalyzer_cfg output=TFileService outpre=genmassanalysis inpre=root://cmseos.fnal.gov//store/user/keanet/tchannel/GEN/step1_LHE-GEN mMediator=3000.0 mDark=20.0 rinv=0.3 alpha=peak part=1 channel="t" maxEvents=1000 madgraph=True > output.txt 
