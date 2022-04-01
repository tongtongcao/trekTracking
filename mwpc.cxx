#include "mwpc.h"
#include <vector>
using namespace std;

mwpc::mwpc(planeId id, planeDir dir):plane(id,dir){
  thickness=33.303; // mm

  for(Int_t i=0;i<12;i++){
    fNHits[i]=0;
    fHits[i] = new TClonesArray("planeHit");

    if(id == kMwpc2){
      frontPlaneHit[i]=new planeHit(id, dir, c2XCenter[i]-thickness/2);
      backPlaneHit[i]=new planeHit(id, dir, c2XCenter[i]+thickness/2);
    }
    if(id == kMwpc3){
      frontPlaneHit[i]=new planeHit(id, dir, c3ZCenter[i]-thickness/2);
      backPlaneHit[i]=new planeHit(id, dir, c3ZCenter[i]+thickness/2);
    }
    if(id == kMwpc4){
      frontPlaneHit[i]=new planeHit(id, dir, c4ZCenter[i]-thickness/2);
      backPlaneHit[i]=new planeHit(id, dir, c4ZCenter[i]+thickness/2);
    }
  }
}

mwpc::~mwpc(){
  for(Int_t i=0;i<12;i++){
    fHits[i]->SetOwner(kTRUE);
    fHits[i]->Clear("C");
    delete fHits[i];

    delete frontPlaneHit[i];
    delete backPlaneHit[i];
  }
}

void mwpc::init(UInt_t nHitsC2X[12], UInt_t nHitsC2Y[12], UInt_t nHitsC2Z[12], Double_t posC2X[12][28], Double_t posC2Y[12][8],Double_t posC2Z[12][28]){
  for(int i=0;i<12;i++){
    fNHits[i]=0;
    if(fHits[i]->GetEntries()!=0) fHits[i]->Clear();
    for(UInt_t j=0;j<nHitsC2X[i];j++){
      for(UInt_t k=0;k<nHitsC2Y[i];k++){
	new ( (*fHits[i])[fNHits[i]]) planeHit(fPlaneId, fPlaneDir, j, k, posC2X[i][j], posC2Y[i][k], posC2Z[i][j],0.01,1,0.2);  
	fNHits[i]++;
      }
    }
  }
}

void mwpc::init(UInt_t nHitsC3X[12], UInt_t nHitsC3Y[12], UInt_t nHitsC3Z[12], Double_t posC3X[12][32], Double_t posC3Y[12][8],Double_t posC3Z[12][1]){
  for(int i=0;i<12;i++){
    fNHits[i]=0;
    if(fHits[i]->GetEntries()!=0) fHits[i]->Clear();
    for(UInt_t j=0;j<nHitsC3X[i];j++){
      for(UInt_t k=0;k<nHitsC3Y[i];k++){
	new ( (*fHits[i])[fNHits[i]]) planeHit(fPlaneId, fPlaneDir, j, k, posC3X[i][j], posC3Y[i][k], posC3Z[i][0],0.2,1,0.01);  
	fNHits[i]++;
      }
    }
  }
}

void mwpc::init(UInt_t nHitsC4X[12], UInt_t nHitsC4Y[12], UInt_t nHitsC4Z[12], Double_t posC4X[12][36], Double_t posC4Y[12][8],Double_t posC4Z[12][1]){
  for(int i=0;i<12;i++){
    fNHits[i]=0;
    if(fHits[i]->GetEntries()!=0) fHits[i]->Clear();
    for(UInt_t j=0;j<nHitsC4X[i];j++){
      for(UInt_t k=0;k<nHitsC4Y[i];k++){
	new ( (*fHits[i])[fNHits[i]]) planeHit(fPlaneId, fPlaneDir, j, k, posC4X[i][j], posC4Y[i][k], posC4Z[i][0],0.2,1,0.01);  
	fNHits[i]++;
      }
    }
  }
}

ClassImp(mwpc);


