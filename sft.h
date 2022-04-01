#ifndef ROOT_SFT
#define ROOT_SFT

#include "TClonesArray.h"
#include "TVector3.h"

#include "utility.h"
#include "planeHit.h"
#include "plane.h"

class sft: public plane{
 public:
  sft(planeId id, planeDir dir);
  ~sft();
  void init(vector<double> x, vector<double> y, vector<double> z, vector<double> xErr, vector<double> yErr, vector<double> zErr);
  int getNHits() const{return fNHits;}
  TClonesArray* getHits() const{return fHits;}

  Double_t getRadiusOuter() const{return radiusOuter;}
  Double_t getRadiusMiddle() const{return radiusMiddle;}
  Double_t getRadiusInner() const{return radiusInner;}

 private:
  TClonesArray* fHits;
  int fNHits;

  Double_t radiusOuter;
  Double_t radiusMiddle;
  Double_t radiusInner;

  ClassDef(sft,1)      
};

#endif
