#ifndef ROOT_AC
#define ROOT_AC

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"

class ac: public plane{
 public:
  ac(planeId id, planeDir dir);
  ~ac();

  Double_t getRadiusTopHousingOuter() const{return radiusTopHousingOuter;}
  Double_t getRadiusTopHousingInner() const{return radiusTopHousingInner;}

  Double_t getRadiusOuter() const{return radiusOuter;}
  Double_t getRadiusInner() const{return radiusInner;}

  Double_t getRadiusBottomHousingOuter() const{return radiusBottomHousingOuter;}
  Double_t getRadiusBottomHousingInner() const{return radiusBottomHousingInner;}

 private:
  Double_t radiusTopHousingOuter;
  Double_t radiusTopHousingInner;
  Double_t radiusOuter;
  Double_t radiusInner;
  Double_t radiusBottomHousingOuter;
  Double_t radiusBottomHousingInner;

  ClassDef(ac,1)      
};

#endif
