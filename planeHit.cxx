#include "TVector3.h"
#include "planeHit.h"

planeHit::planeHit(planeId id, planeDir dir, Double_t xOrZ):fPlaneId(id), fPlaneDir(dir){
  if(dir==kPlaneDirX) fX=xOrZ;
  else if(dir==kPlaneDirZ) fZ=xOrZ;
}


TVector3 planeHit::getPosition () const{
  TVector3 pos(fX,fY,fZ);
  return pos;
}

TVector3 planeHit::getPositionError () const{
  TVector3 posErr(fXErr,fYErr,fZErr);
  return posErr;
}

ClassImp(planeHit);
