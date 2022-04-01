#include <iostream>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "utility.h"
#include "matrix.h"
#include "fieldMap.h"
#include "planeHit.h"
#include "plane.h"
#include "target.h"
#include "sft.h"
#include "tof1.h"
#include "ac.h"
#include "mwpc.h"
#include "tof2.h"
#include "trackState.h"
#include "trackSite.h"
#include "trackSystem.h"
#include "track.h"
#include "trackFinder.h"


using namespace std;


int main(){
  TFile *file=new TFile("mwpcCookedDataV21Run3994.root");
  TTree *tree=(TTree*)file->Get("MWPCCookedEventV21");
  UInt_t run,event;
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);
  UInt_t nHitsC2X[12],nHitsC2Y[12],nHitsC2Z[12],nHitsC3X[12],nHitsC3Y[12],nHitsC3Z[12],nHitsC4X[12],nHitsC4Y[12],nHitsC4Z[12];
  Double_t posC2X[12][28], posC2Y[12][8], posC2Z[12][28], posC3X[12][32], posC3Y[12][8], posC3Z[12][1], posC4X[12][36], posC4Y[12][8], posC4Z[12][1];
  tree->SetBranchAddress("nHitsC2X",nHitsC2X);
  tree->SetBranchAddress("nHitsC2Y",nHitsC2Y);
  tree->SetBranchAddress("nHitsC2Z",nHitsC2Z);
  tree->SetBranchAddress("nHitsC3X",nHitsC3X);
  tree->SetBranchAddress("nHitsC3Y",nHitsC3Y);
  tree->SetBranchAddress("nHitsC3Z",nHitsC3Z);
  tree->SetBranchAddress("nHitsC4X",nHitsC4X);
  tree->SetBranchAddress("nHitsC4Y",nHitsC4Y);
  tree->SetBranchAddress("nHitsC4Z",nHitsC4Z);


  tree->SetBranchAddress("posC2X",posC2X);
  tree->SetBranchAddress("posC2Y",posC2Y);
  tree->SetBranchAddress("posC2Z",posC2Z);
  tree->SetBranchAddress("posC3X",posC3X);
  tree->SetBranchAddress("posC3Y",posC3Y);
  tree->SetBranchAddress("posC3Z",posC3Z);
  tree->SetBranchAddress("posC4X",posC4X);
  tree->SetBranchAddress("posC4Y",posC4Y);
  tree->SetBranchAddress("posC4Z",posC4Z);

  TFile *fileToF=new TFile("tofRun3994.root");
  TTree *treeToF=(TTree*)fileToF->Get("cookedTof");
  Int_t gapNumToF1;
  treeToF->SetBranchAddress("gapNumToF1",&gapNumToF1);
  Int_t gapNumToF2;
  treeToF->SetBranchAddress("gapNumToF2",&gapNumToF2);
  Double_t timeUToF1, timeDToF1;
  Double_t timeAIToF2, timeAOToF2, timeBIToF2, timeBOToF2;
  treeToF->SetBranchAddress("timeUToF1",&timeUToF1);
  treeToF->SetBranchAddress("timeDToF1",&timeDToF1);
  treeToF->SetBranchAddress("timeAIToF2",&timeAIToF2);
  treeToF->SetBranchAddress("timeAOToF2",&timeAOToF2);
  treeToF->SetBranchAddress("timeBIToF2",&timeBIToF2);
  treeToF->SetBranchAddress("timeBOToF2",&timeBOToF2);

  TFile *fileTar=new TFile("targetRun3994.root");
  TTree *treeTar=(TTree*)fileTar->Get("targetRun3994");
  Int_t tag;
  Double_t vertX,vertY,vertPhi;
  Double_t vertXErr,vertYErr,vertPhiErr;
  treeTar->SetBranchAddress("tag",&tag);
  treeTar->SetBranchAddress("vertX",&vertX);
  treeTar->SetBranchAddress("vertY",&vertY);
  treeTar->SetBranchAddress("vertPhi",&vertPhi);
  treeTar->SetBranchAddress("vertXErr",&vertXErr);
  treeTar->SetBranchAddress("vertYErr",&vertYErr);
  treeTar->SetBranchAddress("vertPhiErr",&vertPhiErr);

  TFile* rootFile=new TFile("trackInformation.root","recreate");
  TTree* treeTrack=new TTree("tracks","");
  treeTrack->Branch("run",&run,"run/I");
  treeTrack->Branch("event",&event,"event/I");

  //----------------------------------- Muon ---------------------------------------//
  Int_t nTracks=MAXNTRACKS;
  Int_t nHits[nTracks];
  Double_t chi2PerNDF[nTracks];
  Int_t flagTargetToAC[nTracks];
  Bool_t flagTofMatch[nTracks];
  treeTrack->Branch("nTracks",&nTracks,"nTracks/I");
  treeTrack->Branch("nHits",nHits,"nHits[nTracks]/I");
  treeTrack->Branch("chi2PerNDF",chi2PerNDF,"chi2PerNDF[nTracks]/D");
  treeTrack->Branch("flagTargetToAC",flagTargetToAC,"flagTargetToAC[nTracks]/I");
  treeTrack->Branch("flagTofMatch",flagTofMatch,"flagTofMatch[nTracks]/O");

  Double_t sMwpc4[nTracks];
  Double_t tMwpc4[nTracks];
  Double_t eLossMwpc4[nTracks];
  Double_t mXMwpc4[nTracks];
  Double_t mYMwpc4[nTracks];
  Double_t mZMwpc4[nTracks];
  Double_t mXErrMwpc4[nTracks];
  Double_t mYErrMwpc4[nTracks];
  Double_t mZErrMwpc4[nTracks];
  Double_t xMwpc4[nTracks];
  Double_t yMwpc4[nTracks];
  Double_t zMwpc4[nTracks];
  Double_t nXMwpc4[nTracks];
  Double_t nYMwpc4[nTracks];
  Double_t nZMwpc4[nTracks];
  Double_t pMwpc4[nTracks];
  treeTrack->Branch("sMwpc4",sMwpc4,"sMwpc4[nTracks]/D");
  treeTrack->Branch("tMwpc4",tMwpc4,"tMwpc4[nTracks]/D");
  treeTrack->Branch("eLossMwpc4",eLossMwpc4,"eLossMwpc4[nTracks]/D");
  treeTrack->Branch("mXMwpc4",mXMwpc4,"mXMwpc4[nTracks]/D");
  treeTrack->Branch("mYMwpc4",mYMwpc4,"mYMwpc4[nTracks]/D");
  treeTrack->Branch("mZMwpc4",mZMwpc4,"mZMwpc4[nTracks]/D");
  treeTrack->Branch("mXErrMwpc4",mXErrMwpc4,"mXErrMwpc4[nTracks]/D");
  treeTrack->Branch("mYErrMwpc4",mYErrMwpc4,"mYErrMwpc4[nTracks]/D");
  treeTrack->Branch("mZErrMwpc4",mZErrMwpc4,"mZErrMwpc4[nTracks]/D");
  treeTrack->Branch("xMwpc4",xMwpc4,"xMwpc4[nTracks]/D");
  treeTrack->Branch("yMwpc4",yMwpc4,"yMwpc4[nTracks]/D");
  treeTrack->Branch("zMwpc4",zMwpc4,"zMwpc4[nTracks]/D");
  treeTrack->Branch("nXMwpc4",nXMwpc4,"nXMwpc4[nTracks]/D");
  treeTrack->Branch("nYMwpc4",nYMwpc4,"nYMwpc4[nTracks]/D");
  treeTrack->Branch("nZMwpc4",nZMwpc4,"nZMwpc4[nTracks]/D");
  treeTrack->Branch("pMwpc4",pMwpc4,"pMwpc4[nTracks]/D");

  Double_t sMwpc3[nTracks];
  Double_t tMwpc3[nTracks];
  Double_t eLossMwpc3[nTracks];
  Double_t mXMwpc3[nTracks];
  Double_t mYMwpc3[nTracks];
  Double_t mZMwpc3[nTracks];
  Double_t mXErrMwpc3[nTracks];
  Double_t mYErrMwpc3[nTracks];
  Double_t mZErrMwpc3[nTracks];
  Double_t xMwpc3[nTracks];
  Double_t yMwpc3[nTracks];
  Double_t zMwpc3[nTracks];
  Double_t nXMwpc3[nTracks];
  Double_t nYMwpc3[nTracks];
  Double_t nZMwpc3[nTracks];
  Double_t pMwpc3[nTracks];
  treeTrack->Branch("sMwpc3",sMwpc3,"sMwpc3[nTracks]/D");
  treeTrack->Branch("tMwpc3",tMwpc3,"tMwpc3[nTracks]/D");
  treeTrack->Branch("eLossMwpc3",eLossMwpc3,"eLossMwpc3[nTracks]/D");
  treeTrack->Branch("mXMwpc3",mXMwpc3,"mXMwpc3[nTracks]/D");
  treeTrack->Branch("mYMwpc3",mYMwpc3,"mYMwpc3[nTracks]/D");
  treeTrack->Branch("mZMwpc3",mZMwpc3,"mZMwpc3[nTracks]/D");
  treeTrack->Branch("mXErrMwpc3",mXErrMwpc3,"mXErrMwpc3[nTracks]/D");
  treeTrack->Branch("mYErrMwpc3",mYErrMwpc3,"mYErrMwpc3[nTracks]/D");
  treeTrack->Branch("mZErrMwpc3",mZErrMwpc3,"mZErrMwpc3[nTracks]/D");
  treeTrack->Branch("xMwpc3",xMwpc3,"xMwpc3[nTracks]/D");
  treeTrack->Branch("yMwpc3",yMwpc3,"yMwpc3[nTracks]/D");
  treeTrack->Branch("zMwpc3",zMwpc3,"zMwpc3[nTracks]/D");
  treeTrack->Branch("nXMwpc3",nXMwpc3,"nXMwpc3[nTracks]/D");
  treeTrack->Branch("nYMwpc3",nYMwpc3,"nYMwpc3[nTracks]/D");
  treeTrack->Branch("nZMwpc3",nZMwpc3,"nZMwpc3[nTracks]/D");
  treeTrack->Branch("pMwpc3",pMwpc3,"pMwpc3[nTracks]/D");

  Double_t sMwpc2[nTracks];
  Double_t tMwpc2[nTracks];
  Double_t eLossMwpc2[nTracks];
  Double_t mXMwpc2[nTracks];
  Double_t mYMwpc2[nTracks];
  Double_t mZMwpc2[nTracks];
  Double_t mXErrMwpc2[nTracks];
  Double_t mYErrMwpc2[nTracks];
  Double_t mZErrMwpc2[nTracks];
  Double_t xMwpc2[nTracks];
  Double_t yMwpc2[nTracks];
  Double_t zMwpc2[nTracks];
  Double_t nXMwpc2[nTracks];
  Double_t nYMwpc2[nTracks];
  Double_t nZMwpc2[nTracks];
  Double_t pMwpc2[nTracks];
  treeTrack->Branch("sMwpc2",sMwpc2,"sMwpc2[nTracks]/D");
  treeTrack->Branch("tMwpc2",tMwpc2,"tMwpc2[nTracks]/D");
  treeTrack->Branch("eLossMwpc2",eLossMwpc2,"eLossMwpc2[nTracks]/D");
  treeTrack->Branch("mXMwpc2",mXMwpc2,"mXMwpc2[nTracks]/D");
  treeTrack->Branch("mYMwpc2",mYMwpc2,"mYMwpc2[nTracks]/D");
  treeTrack->Branch("mZMwpc2",mZMwpc2,"mZMwpc2[nTracks]/D");
  treeTrack->Branch("mXErrMwpc2",mXErrMwpc2,"mXErrMwpc2[nTracks]/D");
  treeTrack->Branch("mYErrMwpc2",mYErrMwpc2,"mYErrMwpc2[nTracks]/D");
  treeTrack->Branch("mZErrMwpc2",mZErrMwpc2,"mZErrMwpc2[nTracks]/D");
  treeTrack->Branch("xMwpc2",xMwpc2,"xMwpc2[nTracks]/D");
  treeTrack->Branch("yMwpc2",yMwpc2,"yMwpc2[nTracks]/D");
  treeTrack->Branch("zMwpc2",zMwpc2,"zMwpc2[nTracks]/D");
  treeTrack->Branch("nXMwpc2",nXMwpc2,"nXMwpc2[nTracks]/D");
  treeTrack->Branch("nYMwpc2",nYMwpc2,"nYMwpc2[nTracks]/D");
  treeTrack->Branch("nZMwpc2",nZMwpc2,"nZMwpc2[nTracks]/D");
  treeTrack->Branch("pMwpc2",pMwpc2,"pMwpc2[nTracks]/D");

  Int_t gapNumTof2[nTracks];
  Char_t partTof2[nTracks];
  Double_t sTof2[nTracks];
  Double_t tTof2[nTracks];
  Double_t timeITof2[nTracks];
  Double_t timeOTof2[nTracks];
  Double_t timeMeanTof2[nTracks];
  Double_t xTof2[nTracks];
  Double_t yTof2[nTracks];
  Double_t zTof2[nTracks];
  Double_t nXTof2[nTracks];
  Double_t nYTof2[nTracks];
  Double_t nZTof2[nTracks];
  Double_t pTof2[nTracks];
  treeTrack->Branch("gapNumTof2",gapNumTof2,"gapNumTof2[nTracks]/I");
  treeTrack->Branch("partTof2",partTof2,"partTof2[nTracks]/C");
  treeTrack->Branch("sTof2",sTof2,"sTof2[nTracks]/D");
  treeTrack->Branch("tTof2",tTof2,"tTof2[nTracks]/D");
  treeTrack->Branch("timeITof2",timeITof2,"timeITof2[nTracks]/D");
  treeTrack->Branch("timeOTof2",timeOTof2,"timeOTof2[nTracks]/D");
  treeTrack->Branch("timeMeanTof2",timeMeanTof2,"timeMeanTof2[nTracks]/D");
  treeTrack->Branch("xTof2",xTof2,"xTof2[nTracks]/D");
  treeTrack->Branch("yTof2",yTof2,"yTof2[nTracks]/D");
  treeTrack->Branch("zTof2",zTof2,"zTof2[nTracks]/D");
  treeTrack->Branch("nXTof2",nXTof2,"nXTof2[nTracks]/D");
  treeTrack->Branch("nYTof2",nYTof2,"nYTof2[nTracks]/D");
  treeTrack->Branch("nZTof2",nZTof2,"nZTof2[nTracks]/D");
  treeTrack->Branch("pTof2",pTof2,"pTof2[nTracks]/D");

  Int_t gapNumTof1[nTracks];
  Double_t sTof1[nTracks];
  Double_t tTof1[nTracks];
  Double_t timeUTof1[nTracks];
  Double_t timeDTof1[nTracks];
  Double_t timeMeanTof1[nTracks];
  Double_t xTof1[nTracks];
  Double_t yTof1[nTracks];
  Double_t zTof1[nTracks];
  Double_t nXTof1[nTracks];
  Double_t nYTof1[nTracks];
  Double_t nZTof1[nTracks];
  Double_t pTof1[nTracks];
  treeTrack->Branch("gapNumTof1",gapNumTof1,"gapNumTof1[nTracks]/I");
  treeTrack->Branch("sTof1",sTof1,"sTof1[nTracks]/D");
  treeTrack->Branch("tTof1",tTof1,"tTof1[nTracks]/D");
  treeTrack->Branch("timeUTof1",timeUTof1,"timeUTof1[nTracks]/D");
  treeTrack->Branch("timeDTof1",timeDTof1,"timeDTof1[nTracks]/D");
  treeTrack->Branch("timeMeanTof1",timeMeanTof1,"timeMeanTof1[nTracks]/D");
  treeTrack->Branch("xTof1",xTof1,"xTof1[nTracks]/D");
  treeTrack->Branch("yTof1",yTof1,"yTof1[nTracks]/D");
  treeTrack->Branch("zTof1",zTof1,"zTof1[nTracks]/D");
  treeTrack->Branch("nXTof1",nXTof1,"nXTof1[nTracks]/D");
  treeTrack->Branch("nYTof1",nYTof1,"nYTof1[nTracks]/D");
  treeTrack->Branch("nZTof1",nZTof1,"nZTof1[nTracks]/D");
  treeTrack->Branch("pTof1",pTof1,"pTof1[nTracks]/D");

  Double_t xSft[nTracks];
  Double_t ySft[nTracks];
  Double_t zSft[nTracks];
  Double_t nXSft[nTracks];
  Double_t nYSft[nTracks];
  Double_t nZSft[nTracks];
  Double_t pSft[nTracks];
  treeTrack->Branch("xSft",xSft,"xSft[nTracks]/D");
  treeTrack->Branch("ySft",ySft,"ySft[nTracks]/D");
  treeTrack->Branch("zSft",zSft,"zSft[nTracks]/D");
  treeTrack->Branch("nXSft",nXSft,"nXSft[nTracks]/D");
  treeTrack->Branch("nYSft",nYSft,"nYSft[nTracks]/D");
  treeTrack->Branch("nZSft",nZSft,"nZSft[nTracks]/D");
  treeTrack->Branch("pSft",pSft,"pSft[nTracks]/D");

  Double_t sVert[nTracks];
  Double_t tVert[nTracks];
  Double_t mXVert[nTracks];
  Double_t mYVert[nTracks];
  Double_t mPhiVert[nTracks];
  Double_t mXErrVert[nTracks];
  Double_t mYErrVert[nTracks];
  Double_t mPhiErrVert[nTracks];
  Double_t xVert[nTracks];
  Double_t yVert[nTracks];
  Double_t zVert[nTracks];
  Double_t nXVert[nTracks];
  Double_t nYVert[nTracks];
  Double_t nZVert[nTracks];
  Double_t pVert[nTracks];
  Double_t phiVert[nTracks];
  treeTrack->Branch("sVert",sVert,"sVert[nTracks]/D");
  treeTrack->Branch("tVert",tVert,"tVert[nTracks]/D");
  treeTrack->Branch("mXVert",mXVert,"mXVert[nTracks]/D");
  treeTrack->Branch("mYVert",mYVert,"mYVert[nTracks]/D");
  treeTrack->Branch("mPhiVert",mPhiVert,"mPhiVert[nTracks]/D");
  treeTrack->Branch("mXErrVert",mXErrVert,"mXErrVert[nTracks]/D");
  treeTrack->Branch("mYErrVert",mYErrVert,"mYErrVert[nTracks]/D");
  treeTrack->Branch("mPhiErrVert",mPhiErrVert,"mPhiErrVert[nTracks]/D");
  treeTrack->Branch("xVert",xVert,"xVert[nTracks]/D");
  treeTrack->Branch("yVert",yVert,"yVert[nTracks]/D");
  treeTrack->Branch("zVert",zVert,"zVert[nTracks]/D");
  treeTrack->Branch("nXVert",nXVert,"nXVert[nTracks]/D");
  treeTrack->Branch("nYVert",nYVert,"nYVert[nTracks]/D");
  treeTrack->Branch("nZVert",nZVert,"nZVert[nTracks]/D");
  treeTrack->Branch("pVert",pVert,"pVert[nTracks]/D");
  treeTrack->Branch("phiVert",phiVert,"phiVert[nTracks]/D");

  Double_t xBLowEdge[nTracks];
  Double_t yBLowEdge[nTracks];
  Double_t zBLowEdge[nTracks];
  Double_t nXBLowEdge[nTracks];
  Double_t nYBLowEdge[nTracks];
  Double_t nZBLowEdge[nTracks];
  Double_t pBLowEdge[nTracks];
  treeTrack->Branch("xBLowEdge",xBLowEdge,"xBLowEdge[nTracks]/D");
  treeTrack->Branch("yBLowEdge",yBLowEdge,"yBLowEdge[nTracks]/D");
  treeTrack->Branch("zBLowEdge",zBLowEdge,"zBLowEdge[nTracks]/D");
  treeTrack->Branch("nXBLowEdge",nXBLowEdge,"nXBLowEdge[nTracks]/D");
  treeTrack->Branch("nYBLowEdge",nYBLowEdge,"nYBLowEdge[nTracks]/D");
  treeTrack->Branch("nZBLowEdge",nZBLowEdge,"nZBLowEdge[nTracks]/D");
  treeTrack->Branch("pBLowEdge",pBLowEdge,"pBLowEdge[nTracks]/D");

  //-----------------------------------PionPlus---------------------------------------//
  Double_t tTof1PionPlus[nTracks];
  Double_t pTof1PionPlus[nTracks];
  treeTrack->Branch("tTof1PionPlus",tTof1PionPlus,"tTof1PionPlus[nTracks]/D");
  treeTrack->Branch("pTof1PionPlus",pTof1PionPlus,"pTof1PionPlus[nTracks]/D");
  Double_t pSftPionPlus[nTracks];
  treeTrack->Branch("pSftPionPlus",pSftPionPlus,"pSftPionPlus[nTracks]/D");
  Double_t tVertPionPlus[nTracks];
  Double_t pVertPionPlus[nTracks];
  treeTrack->Branch("tVertPionPlus",tVertPionPlus,"tVertPionPlus[nTracks]/D");
  treeTrack->Branch("pVertPionPlus",pVertPionPlus,"pVertPionPlus[nTracks]/D");

  //-----------------------------------Positron---------------------------------------//
  Double_t tTof1Positron[nTracks];
  Double_t pTof1Positron[nTracks];
  treeTrack->Branch("tTof1Positron",tTof1Positron,"tTof1Positron[nTracks]/D");
  treeTrack->Branch("pTof1Positron",pTof1Positron,"pTof1Positron[nTracks]/D");
  Double_t pSftPositron[nTracks];
  treeTrack->Branch("pSftPositron",pSftPositron,"pSftPositron[nTracks]/D");
  Double_t tVertPositron[nTracks];
  Double_t pVertPositron[nTracks];
  treeTrack->Branch("tVertPositron",tVertPositron,"tVertPositron[nTracks]/D");
  treeTrack->Branch("pVertPositron",pVertPositron,"pVertPositron[nTracks]/D");

  sft *uSft=new sft(kSft, kPlaneDirX);
  vector<double> sftX, sftY, sftZ;
  vector<double> sftXErr, sftYErr, sftZErr;
  uSft->init(sftX, sftY, sftZ, sftXErr, sftYErr, sftZErr);
  ac *uAc=new ac(kAc, kPlaneDirX);

  ofstream output("output");

  int nentries=treeTar->GetEntries();
  cout<<"Total number of events: "<<nentries<<endl;

  particleId partId;
  for(int i=0;i<10000;i++){

    if(i%100==0) cout<<i<<" events are processed"<<endl;
    //cout<<i<<" events are processed"<<endl;
    //cout<<"----------------------------------------------------"<<endl;

    tree->GetEntry(i);
    treeToF->GetEntry(i);
    treeTar->GetEntry(i);

    mwpc *mwpc2=new mwpc(kMwpc2,kPlaneDirX);
    mwpc *mwpc3=new mwpc(kMwpc3,kPlaneDirZ);
    mwpc *mwpc4=new mwpc(kMwpc4,kPlaneDirZ);
    mwpc2->init(nHitsC2X, nHitsC2Y, nHitsC2Z, posC2X, posC2Y, posC2Z);
    mwpc3->init(nHitsC3X, nHitsC3Y, nHitsC3Z, posC3X, posC3Y, posC3Z);
    mwpc4->init(nHitsC4X, nHitsC4Y, nHitsC4Z, posC4X, posC4Y, posC4Z);

    tof2 *uTof2=new tof2(kTof2,kPlaneDirZ, gapNumToF2, timeAIToF2, timeAOToF2, timeBIToF2, timeBOToF2);

    tof1 *uTof1=new tof1(kTof1,kPlaneDirX, gapNumToF1, timeUToF1, timeDToF1);

    target *uTarget=new target(kTarget, kPlaneDirX, tag, vertX, vertY, vertPhi, vertXErr, vertYErr, vertPhiErr);
    uTarget->setNHits(0);

    output<<"Event: "<<i<<endl;

    partId=kMuonId;
    trackFinder *trackFind=new trackFinder(partId, kParticleMass[partId], kParticleCharge[partId], uTarget, uSft, uTof1, uAc, mwpc2, mwpc3, mwpc4, uTof2);
    TClonesArray* theTracks=new TClonesArray("track", MAXNTRACKS, kTRUE);
    theTracks->Clear();
    if(!trackFind->processHits(theTracks)){
      nTracks=0;
    }
    else   nTracks=theTracks->GetEntries();
  
    output<<"  Number of candidate tracks: "<<nTracks<<endl<<endl;

    for(int j=0; j<nTracks; j++){
      //-----------------------------------Muon---------------------------------------//
      nHits[j]=((track*)theTracks->At(j))->getNHits();
      chi2PerNDF[j]=((track*)theTracks->At(j))->getChi2PerNDF();
      flagTargetToAC[j]=((track*)theTracks->At(j))->getFlagTargetToAC();
      flagTofMatch[j]=((track*)theTracks->At(j))->getFlagToFMatch();

      sMwpc4[j]=((track*)theTracks->At(j))->getMwpc4PathLocation();
      tMwpc4[j]=((track*)theTracks->At(j))->getMwpc4PathTime();
      eLossMwpc4[j]=((track*)theTracks->At(j))->getMwpc4EnergyLoss();
      mXMwpc4[j]=((track*)theTracks->At(j))->getMwpc4MX();
      mYMwpc4[j]=((track*)theTracks->At(j))->getMwpc4MY();
      mZMwpc4[j]=((track*)theTracks->At(j))->getMwpc4MZ();
      mXErrMwpc4[j]=((track*)theTracks->At(j))->getMwpc4MXErr();
      mYErrMwpc4[j]=((track*)theTracks->At(j))->getMwpc4MYErr();
      mZErrMwpc4[j]=((track*)theTracks->At(j))->getMwpc4MZErr();
      xMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SX();
      yMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SY();
      zMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SZ();
      nXMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SNx();
      nYMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SNy();
      nZMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SNz();
      pMwpc4[j]=((track*)theTracks->At(j))->getMwpc4SP();

      sMwpc3[j]=((track*)theTracks->At(j))->getMwpc3PathLocation();
      tMwpc3[j]=((track*)theTracks->At(j))->getMwpc3PathTime();
      eLossMwpc3[j]=((track*)theTracks->At(j))->getMwpc3EnergyLoss();
      mXMwpc3[j]=((track*)theTracks->At(j))->getMwpc3MX();
      mYMwpc3[j]=((track*)theTracks->At(j))->getMwpc3MY();
      mZMwpc3[j]=((track*)theTracks->At(j))->getMwpc3MZ();
      mXErrMwpc3[j]=((track*)theTracks->At(j))->getMwpc3MXErr();
      mYErrMwpc3[j]=((track*)theTracks->At(j))->getMwpc3MYErr();
      mZErrMwpc3[j]=((track*)theTracks->At(j))->getMwpc3MZErr();
      xMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SX();
      yMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SY();
      zMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SZ();
      nXMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SNx();
      nYMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SNy();
      nZMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SNz();
      pMwpc3[j]=((track*)theTracks->At(j))->getMwpc3SP();

      sMwpc2[j]=((track*)theTracks->At(j))->getMwpc2PathLocation();
      tMwpc2[j]=((track*)theTracks->At(j))->getMwpc2PathTime();
      eLossMwpc2[j]=((track*)theTracks->At(j))->getMwpc2EnergyLoss();
      mXMwpc2[j]=((track*)theTracks->At(j))->getMwpc2MX();
      mYMwpc2[j]=((track*)theTracks->At(j))->getMwpc2MY();
      mZMwpc2[j]=((track*)theTracks->At(j))->getMwpc2MZ();
      mXErrMwpc2[j]=((track*)theTracks->At(j))->getMwpc2MXErr();
      mYErrMwpc2[j]=((track*)theTracks->At(j))->getMwpc2MYErr();
      mZErrMwpc2[j]=((track*)theTracks->At(j))->getMwpc2MZErr();
      xMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SX();
      yMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SY();
      zMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SZ();
      nXMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SNx();
      nYMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SNy();
      nZMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SNz();
      pMwpc2[j]=((track*)theTracks->At(j))->getMwpc2SP();

      gapNumTof2[j]=((track*)theTracks->At(j))->getTof2GapNum();
      partTof2[j]=((track*)theTracks->At(j))->getTof2Part();
      sTof2[j]=((track*)theTracks->At(j))->getTof2PathLocation();
      tTof2[j]=((track*)theTracks->At(j))->getTof2PathTime();
      timeITof2[j]=((track*)theTracks->At(j))->getTof2TimeI();
      timeOTof2[j]=((track*)theTracks->At(j))->getTof2TimeO();
      timeMeanTof2[j]=((track*)theTracks->At(j))->getTof2TimeMean();
      xTof2[j]=((track*)theTracks->At(j))->getTof2SX();
      yTof2[j]=((track*)theTracks->At(j))->getTof2SY();
      zTof2[j]=((track*)theTracks->At(j))->getTof2SZ();
      nXTof2[j]=((track*)theTracks->At(j))->getTof2SNx();
      nYTof2[j]=((track*)theTracks->At(j))->getTof2SNy();
      nZTof2[j]=((track*)theTracks->At(j))->getTof2SNz();
      pTof2[j]=((track*)theTracks->At(j))->getTof2SP();

      gapNumTof1[j]=((track*)theTracks->At(j))->getTof1GapNum();
      sTof1[j]=((track*)theTracks->At(j))->getTof1PathLocation();
      tTof1[j]=((track*)theTracks->At(j))->getTof1PathTime();
      timeUTof1[j]=((track*)theTracks->At(j))->getTof1TimeU();
      timeDTof1[j]=((track*)theTracks->At(j))->getTof1TimeD();
      timeMeanTof1[j]=((track*)theTracks->At(j))->getTof1TimeMean();
      xTof1[j]=((track*)theTracks->At(j))->getTof1SX();
      yTof1[j]=((track*)theTracks->At(j))->getTof1SY();
      zTof1[j]=((track*)theTracks->At(j))->getTof1SZ();
      nXTof1[j]=((track*)theTracks->At(j))->getTof1SNx();
      nYTof1[j]=((track*)theTracks->At(j))->getTof1SNy();
      nZTof1[j]=((track*)theTracks->At(j))->getTof1SNz();
      pTof1[j]=((track*)theTracks->At(j))->getTof1SP();

      xSft[j]=((track*)theTracks->At(j))->getSftSX();
      ySft[j]=((track*)theTracks->At(j))->getSftSY();
      zSft[j]=((track*)theTracks->At(j))->getSftSZ();
      nXSft[j]=((track*)theTracks->At(j))->getSftSNx();
      nYSft[j]=((track*)theTracks->At(j))->getSftSNy();
      nZSft[j]=((track*)theTracks->At(j))->getSftSNz();
      pSft[j]=((track*)theTracks->At(j))->getSftSP();

      sVert[j]=((track*)theTracks->At(j))->getVertPathLocation();
      tVert[j]=((track*)theTracks->At(j))->getVertPathTime();
      mXVert[j]=((track*)theTracks->At(j))->getVertMX();
      mYVert[j]=((track*)theTracks->At(j))->getVertMY();
      mPhiVert[j]=((track*)theTracks->At(j))->getVertMPhi();
      mXErrVert[j]=((track*)theTracks->At(j))->getVertMXErr();
      mYErrVert[j]=((track*)theTracks->At(j))->getVertMYErr();
      mPhiErrVert[j]=((track*)theTracks->At(j))->getVertMPhiErr();
      xVert[j]=((track*)theTracks->At(j))->getVertSX();
      yVert[j]=((track*)theTracks->At(j))->getVertSY();
      zVert[j]=((track*)theTracks->At(j))->getVertSZ();
      nXVert[j]=((track*)theTracks->At(j))->getVertSNx();
      nYVert[j]=((track*)theTracks->At(j))->getVertSNy();
      nZVert[j]=((track*)theTracks->At(j))->getVertSNz();
      pVert[j]=((track*)theTracks->At(j))->getVertSP();
      phiVert[j]=((track*)theTracks->At(j))->getVertSPhi();

      xBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSX();
      yBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSY();
      zBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSZ();
      nXBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSNx();
      nYBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSNy();
      nZBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSNz();
      pBLowEdge[j]=((track*)theTracks->At(j))->getBLowEdgeSP();

      //-----------------------------------PionPlus---------------------------------------//
      tTof1PionPlus[j]=((track*)theTracks->At(j))->getTof1PathTimePionPlus();
      pTof1PionPlus[j]=((track*)theTracks->At(j))->getTof1SPPionPlus();
      pSftPionPlus[j]=((track*)theTracks->At(j))->getSftSPPionPlus();
      tVertPionPlus[j]=((track*)theTracks->At(j))->getVertPathTimePionPlus();
      pVertPionPlus[j]=((track*)theTracks->At(j))->getVertSPPionPlus();

      //-----------------------------------Positron---------------------------------------//
      tTof1Positron[j]=((track*)theTracks->At(j))->getTof1PathTimePositron();
      pTof1Positron[j]=((track*)theTracks->At(j))->getTof1SPPositron();
      pSftPositron[j]=((track*)theTracks->At(j))->getSftSPPositron();
      tVertPositron[j]=((track*)theTracks->At(j))->getVertPathTimePositron();
      pVertPositron[j]=((track*)theTracks->At(j))->getVertSPPositron();

      //-----------------------------------Muon---------------------------------------//
      output<<"    track: "<<j<<endl;
      output<<"    Number of hits: "<<nHits[j]<<endl;
      output<<"    chi2 over NDF: "<<chi2PerNDF[j]<<endl;
      output<<"    flag from AC to Target: "<<flagTargetToAC[j]<<endl;
      output<<"    flag to denote if there is one or more gap number in the ToF1 list which matchup with gap number of ToF2: "<<flagTofMatch[j]<<endl;
      output<<"      MWPC4:"<<endl;
      output<<"        Measurement (mm): ("<<mXMwpc4[j]<<", "<<mYMwpc4[j]<<", "<<mZMwpc4[j]<<')'<<endl;
      output<<"        Error of Measurement (mm): ("<<mXErrMwpc4[j]<<", "<<mYErrMwpc4[j]<<", "<<mZErrMwpc4[j]<<')'<<endl;
      output<<"        Path location (mm): "<<sMwpc4[j]<<endl;
      output<<"        Path time (ns): "<<tMwpc4[j]<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xMwpc4[j]<<", "<<yMwpc4[j]<<", "<<zMwpc4[j]<<", "<<nXMwpc4[j]<<", "<<nYMwpc4[j]<<", "<<nZMwpc4[j]<<", "<<pMwpc4[j]<<')'<<endl;
      output<<"      MWPC3:"<<endl;
      output<<"        Measurement (mm): ("<<mXMwpc3[j]<<", "<<mYMwpc3[j]<<", "<<mZMwpc3[j]<<')'<<endl;
      output<<"        Error of Measurement (mm): ("<<mXErrMwpc3[j]<<", "<<mYErrMwpc3[j]<<", "<<mZErrMwpc3[j]<<')'<<endl;
      output<<"        Path location (mm): "<<sMwpc3[j]<<endl;
      output<<"        Path time (ns): "<<tMwpc3[j]<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xMwpc3[j]<<", "<<yMwpc3[j]<<", "<<zMwpc3[j]<<", "<<nXMwpc3[j]<<", "<<nYMwpc3[j]<<", "<<nZMwpc3[j]<<", "<<pMwpc3[j]<<')'<<endl;
      output<<"      MWPC2:"<<endl;
      output<<"        Measurement (mm): ("<<mXMwpc2[j]<<", "<<mYMwpc2[j]<<", "<<mZMwpc2[j]<<')'<<endl;
      output<<"        Error of Measurement (mm): ("<<mXErrMwpc2[j]<<", "<<mYErrMwpc2[j]<<", "<<mZErrMwpc2[j]<<')'<<endl;
      output<<"        Path location (mm): "<<sMwpc2[j]<<endl;
      output<<"        Path time (ns): "<<tMwpc2[j]<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xMwpc2[j]<<", "<<yMwpc2[j]<<", "<<zMwpc2[j]<<", "<<nXMwpc2[j]<<", "<<nYMwpc2[j]<<", "<<nZMwpc2[j]<<", "<<pMwpc2[j]<<')'<<endl;
      output<<"      ToF2:"<<endl;
      output<<"        GapNum staring from 0: "<<gapNumTof2[j]<<endl; 
      output<<"        Part: "<<partTof2[j]<<endl; 
      output<<"        Measurement (timeI, timeO, timeMean) (ns): ("<<timeITof2[j]<<", "<<timeOTof2[j]<<", "<<timeMeanTof2[j]<<')'<<endl;
      output<<"        Path location (mm): "<<sTof2[j]<<endl;
      output<<"        Path time (ns): "<<tTof2[j]<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xTof2[j]<<", "<<yTof2[j]<<", "<<zTof2[j]<<", "<<nXTof2[j]<<", "<<nYTof2[j]<<", "<<nZTof2[j]<<", "<<pTof2[j]<<')'<<endl;
      output<<"      ToF1:"<<endl;
      output<<"        GapNum staring from 0: "<<gapNumTof1[j]<<endl; 
      output<<"        Measurement (timeU, timeD, timeMean) (ns): ("<<timeUTof1[j]<<", "<<timeDTof1[j]<<", "<<timeMeanTof1[j]<<')'<<endl;
      output<<"        Path location (mm): "<<sTof1[j]<<endl;
      output<<"        Path time (ns): "<<tTof1[j]<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xTof1[j]<<", "<<yTof1[j]<<", "<<zTof1[j]<<", "<<nXTof1[j]<<", "<<nYTof1[j]<<", "<<nZTof1[j]<<", "<<pTof1[j]<<')'<<endl;

      output<<"      SFT:"<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xSft[j]<<", "<<ySft[j]<<", "<<zSft[j]<<", "<<nXSft[j]<<", "<<nYSft[j]<<", "<<nZSft[j]<<", "<<pSft[j]<<')'<<endl;

      output<<"      Vertex:"<<endl;
      output<<"        Measurement (mm): ("<<mXVert[j]<<", "<<mYVert[j]<<", "<<mPhiVert[j]<<')'<<endl;
      output<<"        Error of Measurement (mm): ("<<mXErrVert[j]<<", "<<mYErrVert[j]<<", "<<mPhiErrVert[j]<<')'<<endl;
      output<<"        Path location (mm): "<<sVert[j]<<endl;
      output<<"        Path time (ns): "<<tVert[j]<<endl;
      output<<"        State vector (mm, MeV/c): ("<<xVert[j]<<", "<<yVert[j]<<", "<<zVert[j]<<", "<<nXVert[j]<<", "<<nYVert[j]<<", "<<nZVert[j]<<", "<<pVert[j]<<", "<<phiVert[j]<<')'<<endl;
    }
    output<<endl<<endl;

    theTracks->SetOwner(kTRUE);
    delete theTracks;
    theTracks=NULL;
    delete trackFind;
    trackFind=NULL;

    treeTrack->Fill();

    delete mwpc2;
    mwpc2=NULL;
    delete mwpc3;
    mwpc3=NULL;
    delete mwpc4;
    mwpc4=NULL;
    delete uTof2;
    uTof2=NULL;
    delete uTof1;
    uTof1=NULL;

    delete uTarget;
    uTarget=NULL;
  }

  delete uSft;
  delete uAc;

  delete tree;
  delete file;

  output.close();
  treeTrack->Write();
  rootFile->Write();

  delete  treeTrack;
  delete  rootFile;

  return 0;
}

