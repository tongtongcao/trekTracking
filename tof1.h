#ifndef ROOT_TOF1
#define ROOT_TOF1

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"

class tof1: public plane{
 public:
  tof1(planeId id, planeDir dir, Int_t gapNum, Double_t timeU, Double_t timeD);
  ~tof1();

  Int_t getGapNum() const{return fGapNum;}
  Double_t getTimeU() const{return fTimeU;}
  Double_t getTimeD() const{return fTimeD;}
  Double_t getTimeMean() const{return fTimeMean;}

  Double_t getRadiusOuter() const{return radiusOuter;}
  Double_t getRadiusInner() const{return radiusInner;}


 private:
  Int_t fGapNum;
  Double_t fTimeU;
  Double_t fTimeD;
  Double_t fTimeMean;

  Double_t tof1XCenter;
  Double_t tof1YCenter;
  Double_t tof1ZCenter;

  Double_t thickness;
  Double_t ySizeTop;
  Double_t ySizeBottom;
  Double_t zSize;


  Double_t radiusOuter;
  Double_t radiusInner;

  ClassDef(tof1,1)      
};

#endif
