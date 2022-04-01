#include "target.h"
#include <vector>
using namespace std;

target::target(planeId id, planeDir dir, Int_t t, Double_t vertx, Double_t verty, Double_t vertphi, Double_t vertx_err, Double_t verty_err, Double_t vertphi_err):plane(id,dir){
  tag=t;
  vertX=vertx;
  vertY=verty;
  vertPhi=vertphi;
  vertXErr=vertx_err;
  vertYErr=verty_err;
  vertPhiErr=vertphi_err;

  fNHits=0;
  fHits = new TClonesArray("planeHit");

  radiusHolderOuter=38.5;
  radiusHolderInner=37;
  radius=27.8;
}

target::~target(){
  fHits->SetOwner(kTRUE);
  fHits->Clear("C");
  delete fHits;
}

void target::init(vector<double> x, vector<double> y, vector<double> z, vector<double> xErr, vector<double> yErr, vector<double> zErr){
  fNHits=0;
  if(fHits->GetEntries()!=0) fHits->Clear();
}

ClassImp(target);


