#ifndef ROOT_TARGET
#define ROOT_TARGET

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"

class target: public plane{
 public:
  target(planeId id, planeDir dir, Int_t t, Double_t vertx, Double_t verty, Double_t vertphi, Double_t vertx_err, Double_t verty_err, Double_t vertphi_err);
  ~target();

  Double_t getTag() const{return tag;}

  Double_t getVertX() const{return vertX;}
  Double_t getVertY() const{return vertY;}
  Double_t getVertPhi() const{return vertPhi;}

  Double_t getVertXErr() const{return vertXErr;}
  Double_t getVertYErr() const{return vertYErr;}
  Double_t getVertPhiErr() const{return vertPhiErr;}

  void init(vector<double> x, vector<double> y, vector<double> z, vector<double> xErr, vector<double> yErr, vector<double> zErr);
  Int_t getNHits() const{return fNHits;}
  void setNHits(Int_t nhits) {fNHits=nhits;}
  TClonesArray* getHits() const{return fHits;}

  Double_t getRadiusHolderOuter() const{return radiusHolderOuter;}
  Double_t getRadiusHolderInner() const{return radiusHolderInner;}
  Double_t getRadius() const{return radius;}

 private:
  Int_t tag;
  Double_t vertX;
  Double_t vertY;
  Double_t vertPhi;

  Double_t vertXErr;
  Double_t vertYErr;
  Double_t vertPhiErr;

  TClonesArray* fHits;
  int fNHits;

  Double_t radius;
  Double_t radiusHolderInner;
  Double_t radiusHolderOuter;

  ClassDef(target,1)      
};

#endif
