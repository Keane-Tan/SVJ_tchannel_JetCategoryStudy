# SVJ_tchannel_JetCategoryStudy

## v2_darkIdentifier_newJetCat.py
This python script uses the jet categorization algorithm implemented in https://github.com/cms-svj/SVJProduction/pull/9.
This script also stores jets for calculating MT (L395-418) and MT2 (L423-462).
In each case, the jet chosen is only matched to one and only one particle (a DQM or SMM) from the mediator.
DQM = dark quark from mediator; SMM = SM quark from mediator.
If a jet contains more than one particle from the mediator, or if a particle from the mediator is found in more than one jet, then the event is ignored when calculating MT or MT2.
`storTLV` is a function that stores a list of jets as a vector of TLorentzVector in the root file.
