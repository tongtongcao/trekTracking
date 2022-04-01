#ifndef ROOT_TOF2
#define ROOT_TOF2

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"

const Double_t tof2ZCenter[12]={1600, 1613, 1621, 1620, 1536, 1539, 1535, 1611, 1608, 1600, 1598, 1601};
const Double_t tof2XCenter[12]={1390, 1390, 1390, 1390, 1316, 1390, 1316, 1390, 1390, 1390, 1390, 1390};
const Double_t tof2YCenter[12]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

class tof2: public plane{
 public:
  tof2(planeId id, planeDir dir, Int_t gapNum, Double_t timeAI, Double_t timeAO, Double_t timeBI, Double_t timeBO);
  ~tof2();

  Int_t getGapNum() const{return fGapNum;}
  Char_t getPart() const{return fPart;}

  Double_t getTimeI() const{return fTimeI;}
  Double_t getTimeO() const{return fTimeO;}
  Double_t getTimeMean() const{return fTimeMean;}

  planeHit* getFrontPlaneHit(Int_t gapNumToF2) const{return frontPlaneHit[gapNumToF2];}
  planeHit* getMiddlePlaneHit(Int_t gapNumToF2) const{return middlePlaneHit[gapNumToF2];}
  planeHit* getBackPlaneHit(Int_t gapNumToF2) const{return backPlaneHit[gapNumToF2];}

 private:
  Int_t fGapNum;
  Char_t fPart;
  Double_t fTimeI;
  Double_t fTimeO;
  Double_t fTimeMean;

  Double_t thickness;
  Double_t xSize[12];
  Double_t ySize;

  planeHit* frontPlaneHit[12]; // not real hit; only a coordinate (x/z) is set and valid before propagation; the other two coordinates can be determined after propagation; gap-dependent
  planeHit* middlePlaneHit[12]; // not real hit; only a coordinate (x/z) is set and valid before propagation; the other two coordinates can be determined after propagation; gap-dependent
  planeHit* backPlaneHit[12]; // not real hit; only a coordinate (x/z) is set and valid before propagation; the other two coordinates can be determined after propagation; gap-dependent


  ClassDef(tof2,1)      
};

#endif
