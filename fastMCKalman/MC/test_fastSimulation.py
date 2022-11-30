# include fastMCKalman to the path
import sys,os;
sys.path.insert(1, os.environ[f"fastMCKalman"]+'/fastMCKalman/MC/');
import ROOT
from ROOT import TFile, gSystem
from fastSimulation import *
from RootInteractive.Tools.aliTreePlayer import *
from RootInteractive.MLpipeline.MIForestErrPDF import *


# load tree
inputData = "fastParticle.list"
ROOT.initTreeFast(inputData)
tree=ROOT.treeFast

# define aliases
def makeAliases(tree):
    for pType in ["In","Out","Refit","MC"]:
        tree.SetAlias(f"statusMaskFull{pType}",f"partFull.fStatusMask{pType}")
        for iPar in [0,1,2,3,4]:
            tree.SetAlias(f"paramFull{pType}{iPar}",f"partFull.fParam{pType}[].fP[{iPar}]")
            iCov=ROOT.AliExternalTrackParam.GetIndex(iPar, iPar)
            tree.SetAlias(f"pullFull{pType}{iPar}",
                          f"(partFull.fParam{pType}[].fP[{iPar}]-partFull.fParamMC[].fP[{iPar}])/sqrt(partFull.fParam{pType}[].fC[{iCov}])")


def loadPanda(tree):
    variables=[".*pullFull.*",".*statusMaskFull.*",".*paramFull.*",
               ".*gx.*",".*gy.*",".*gz.*",
               "ptMC","tglMC","fPdgCodeMC",
               # to add  lever arm at given layer
               ]
    df=tree2Panda(tree,variables, "partFull.fLengthIn>5",columnMask=[["_fElements",""]],nEntries=10000)
    df["statusMaskFullRefit"]=df["statusMaskFullRefit"].astype("int")
    return df

def makeRegression(df):
    varIn=["fPdgCodeMC","paramFullMC2","paramFullMC3","paramFullMC4"]
    target="pullFullRefit4"
    n_estimators=200
    n_jobs=100
    nPoints=500000
    max_depthBase=14
    max_samples=0.05

    regressor = RandomForestRegressor(n_estimators =n_estimators,n_jobs=n_jobs,max_depth=max_depthBase,max_samples=max_samples)
    #dfFit=df.query("(statusMaskFullRefit&0x2000)>0")
    dfFit0=df[((df["statusMaskFullRefit"]&0x1000)>0)]
    dfFit1=dfFit0.query(f"abs({target})<5")
    #
    nAll=dfFit1.shape[0]
    regressor.fit(dfFit1[varIn][:nAll//2], dfFit1[target][:nAll//2])
    dfFit1[f"{target}Pred0"]=regressor.predict(dfFit1[varIn])
    #
    ((dfFit1[f"{target}"]-dfFit1[f"{target}Pred0"])[nAll//2:]).std()


