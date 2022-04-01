#include "sft.h"
#include <vector>
using namespace std;

sft::sft(planeId id, planeDir dir):plane(id,dir){
  fNHits=0;
  fHits = new TClonesArray("planeHit");

  radiusOuter=43.95;
  radiusMiddle=41.6;
  radiusInner=40.25;
}

sft::~sft(){
  fHits->SetOwner(kTRUE);
  fHits->Clear("C");
  delete fHits;
}

void sft::init(vector<double> x, vector<double> y, vector<double> z, vector<double> xErr, vector<double> yErr, vector<double> zErr){
  fNHits=0;
  if(fHits->GetEntries()!=0) fHits->Clear();
}

ClassImp(sft);


