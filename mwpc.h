#ifndef ROOT_MWPC
#define ROOT_MWPC

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"

const Double_t c2XCenter[12]={629.6, 629.6, 629.2, 629.85, 629.45, 630.1, 629.95, 629.6, 629.2, 629.6, 629.6, 629.6};
const Double_t c3ZCenter[12]={480.5, 480.5, 480.5, 480.5, 480.5, 480.5, 480.5, 480.5, 480.5, 480.5, 480.5, 480.5};
const Double_t c4ZCenter[12]={699.5, 699.5, 699.5, 699.5, 699.5, 699.5, 699.5, 699.5, 699.5, 699.5, 699.5, 699.5};

class mwpc: public plane{
 public:
  mwpc(planeId id, planeDir dir);
  ~mwpc();
  void init(UInt_t nHitsC2X[12], UInt_t nHitsC2Y[12], UInt_t nHitsC2Z[12], Double_t posC2X[12][28], Double_t posC2Y[12][8],Double_t posC2Z[12][28]);
  void init(UInt_t nHitsC3X[12], UInt_t nHitsC3Y[12], UInt_t nHitsC3Z[12], Double_t posC3X[12][32], Double_t posC3Y[12][8],Double_t posC3Z[12][1]);
  void init(UInt_t nHitsC4X[12], UInt_t nHitsC4Y[12], UInt_t nHitsC4Z[12], Double_t posC4X[12][36], Double_t posC4Y[12][8],Double_t posC4Z[12][1]);

  Int_t getNHits(Int_t gapNumToF2) const{return fNHits[gapNumToF2];}
  TClonesArray* getHits(Int_t gapNumToF2) const{return fHits[gapNumToF2];}

  planeHit* getFrontPlaneHit(Int_t gapNumToF2) const{return frontPlaneHit[gapNumToF2];}
  planeHit* getBackPlaneHit(Int_t gapNumToF2) const{return backPlaneHit[gapNumToF2];}

 private:
  TClonesArray* fHits[12];
  Int_t fNHits[12];
  Double_t thickness;
  planeHit* frontPlaneHit[12]; // not real hit; only a coordinate (x/z) is set and valid before propagation; the other two coordinates can be determined after propagation; gap-dependent
  planeHit* backPlaneHit[12]; // not real hit; only a coordinate (x/z) is set and valid before propagation; the other two coordinates can be determined after propagation; gap-dependent
  ClassDef(mwpc,1)      

};

#endif
