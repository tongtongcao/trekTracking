#include "tof2.h"
#include <vector>
using namespace std;

tof2::tof2(planeId id, planeDir dir, Int_t gapNum, Double_t timeAI, Double_t timeAO, Double_t timeBI, Double_t timeBO):plane(id,dir){

  fGapNum=gapNum-1;

  if(timeAI>0){
    fPart='A';
    fTimeI=timeAI;
    fTimeO=timeAO;
    fTimeMean=(timeAI+timeAO)/2;
  }
  else{
    fPart='B';
    fTimeI=timeBI;
    fTimeO=timeBO;
    fTimeMean=(timeBI+timeBO)/2;
  }

  thickness=20; // mm
  ySize=2*150;
  for(Int_t i=0;i<12;i++){
    if(i==4 || i==6) xSize[i]=680;
    else xSize[i]=800;
    frontPlaneHit[i]=new planeHit(id, dir, tof2ZCenter[i]-thickness/2);
    middlePlaneHit[i]=new planeHit(id, dir, tof2ZCenter[i]);
    backPlaneHit[i]=new planeHit(id, dir, tof2ZCenter[i]+thickness/2);
  }
}

tof2::~tof2(){
  for(Int_t i=0;i<12;i++){
    delete frontPlaneHit[i];
    delete middlePlaneHit[i];
    delete backPlaneHit[i];
  }
}

ClassImp(tof2);


