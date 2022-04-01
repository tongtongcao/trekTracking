#include "tof1.h"
#include <vector>
using namespace std;

tof1::tof1(planeId id, planeDir dir, Int_t gapNum, Double_t timeU, Double_t timeD):plane(id,dir){

  fGapNum=gapNum-1;
  fTimeU=timeU;
  fTimeD=timeD;
  fTimeMean=(timeU+timeD)/2;

  tof1XCenter=47.; //mm
  tof1YCenter=0.; //mm
  tof1ZCenter=0.; //mm

  thickness=5.; //mm
  ySizeTop=26.5; //mm
  ySizeBottom=23.8; //mm
  zSize=188; //mm

  //radiusOuter=49.5;
  //radiusInner=44.5;
  radiusOuter=50.078;
  radiusInner=45.0196;

}

tof1::~tof1(){}

ClassImp(tof1);


