#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
using namespace std;

#include "trackSite.h"



trackSite::trackSite(Int_t m, Int_t p, Double_t chi2) : TObjArray(2), fCurStatePtr(0), fM(m,1), fV(m,m), fH(m,p), fHT(p,m), fDeltaChi2(0.), fMaxDeltaChi2(chi2){}
//_________________________________________________________________________
trackSite::trackSite(planeHit *hit, Int_t m, Int_t  p, Double_t chi2) : TObjArray(2), fCurStatePtr(0), fM(m,1), fV(m,m),fH(m,p), fHT(p,m), fDeltaChi2(0.), fMaxDeltaChi2(chi2), fPlaneHitPtr(hit) {

  fPlaneHit=*fPlaneHitPtr;

  fM(kIdxX, 0) = hit->getPosition().X(); 
  fM(kIdxY, 0) = hit->getPosition().Y();
  fM(kIdxZ, 0) = hit->getPosition().Z();

  fV(kIdxX, kIdxX) = pow(hit->getPositionError().X(), 2);
  fV(kIdxY, kIdxY) = pow(hit->getPositionError().Y(), 2);
  fV(kIdxZ, kIdxZ) = pow(hit->getPositionError().Z(), 2);

  fH.Zero();
  fH(0, 0) = 1.;
  fH(1, 1) = 1.;
  fH(2, 2) = 1.;

  fHT = matrix(matrix::kTransposed, fH);

  fPlaneDir=hit->getPlaneDir();
  if(fPlaneDir==kPlaneDirX) fPlaneLocation=hit->getPosition().X();
  else if(fPlaneDir==kPlaneDirZ) fPlaneLocation=hit->getPosition().Z();
  fPathLocation=0;
  fPathTime=0;
}
//_________________________________________________________________________
trackSite::~trackSite(){
  SetOwner(kTRUE);
  Delete();
}
//_________________________________________________________________________

void trackSite::Add(TObject *obj){
   TObjArray::Add(obj);
   fCurStatePtr = static_cast<trackState*>(obj);
   fCurStatePtr->setSitePtr(this);
   setPathLocation(fCurStatePtr->getPathLocation());
   setPathTime(fCurStatePtr->getPathTime());
}

//___________________________________________________________________________
Bool_t trackSite::filter(Double_t particleCharge){
   // prea and preC should be preset
   trackState *preSV = getStatePtr(trackSite::kPredicted);
   matrix preC  = getStatePtr(trackSite::kPredicted)->getCovMat();
   if(particleCharge/(*preSV)(kIdxQP,0)<kMomLowerLimit || particleCharge/(*preSV)(kIdxQP,0)>kMomUpperLimit){
     fDeltaChi2=1.e12;
     return kFALSE;
   }
   matrix preh(3, 1);
   preh(kIdxX, 0)=(*preSV)(kIdxX, 0);
   preh(kIdxY, 0)=(*preSV)(kIdxY, 0);
   preh(kIdxZ, 0)=(*preSV)(kIdxZ, 0);
   matrix pull  = fM - preh;
   matrix preR    = fV + fH * preC * fHT;
   try{
     if(preR.Determinant() == 0) throw 0.;
   }
   catch(Double_t & cond){
     cout<<"Don't worry. Here's error (matrix is singular) has been fixed. Go ahead."<<endl;
     fDeltaChi2=1.e12;
     return kFALSE;
   }
   matrix preRInv = matrix(matrix::kInverted, preR);
   matrix I(kSdim, kSdim);
   I.UnitMatrix();
   matrix K = preC * fHT *  preRInv;
   matrix curC = (I - K * fH) * preC;
   matrix curSV = *preSV + K * pull;
   if(particleCharge/curSV(kIdxQP,0)<kMomLowerLimit || particleCharge/curSV(kIdxQP,0)>kMomUpperLimit){
     fDeltaChi2=1.e12;
     return kFALSE;
   }
   trackState *newSV = new trackState(curSV, curC, kFiltered);
   newSV->setPathLocation(preSV->getPathLocation());
   newSV->setPathTime(preSV->getPathTime());
   newSV->setNoiseMat(preSV->getNoiseMat());
   newSV->setPropMat(preSV->getPropMat());
   newSV->setPropTranMat(preSV->getPropTranMat());
   Add(newSV);

   // Calculate chi2 increment
   matrix curh(3, 1);
   curh(kIdxX, 0)=curSV(kIdxX, 0);
   curh(kIdxY, 0)=curSV(kIdxY, 0);
   curh(kIdxZ, 0)=curSV(kIdxZ, 0);
   matrix resVec = fM - curh;
   matrix resVecT = matrix(matrix::kTransposed, resVec);
  
   try{
     if(fV.Determinant() == 0) throw 0.;
   }
   catch(Double_t & cond){
     cout<<"Don't worry. Here's error (matrix is singular) has been fixed. Go ahead."<<endl;
     fDeltaChi2=1.e12;
     return kFALSE;
   }
   matrix VInv = matrix(matrix::kInverted, fV);
   matrix resSV = curSV - *preSV;
   matrix resSVT = matrix(matrix::kTransposed, resSV);

   try{
     if(preC.Determinant() == 0) throw 0.;
   }
   catch(Double_t & cond){
     cout<<"Don't worry. Here's error (matrix is singular) has been fixed. Go ahead."<<endl;
     fDeltaChi2=1.e12;
     return kFALSE;
   }
   matrix preCInv = matrix(matrix::kInverted, preC);
   fDeltaChi2=(resSVT * preCInv * resSV + resVecT * VInv * resVec)(0,0);
 

   if(fDeltaChi2 < fMaxDeltaChi2) return kTRUE;
   else return kFALSE;
}
//______________________________________________________________________________
void trackSite::smooth(trackSite* pre){
   if (getStatePtr(trackSite::kSmoothed)) return;

   trackState* cura  = getStatePtr(trackSite::kFiltered);
   trackState* prea  = pre->getStatePtr(trackSite::kPredicted);
   trackState* sprea = pre->getStatePtr(trackSite::kSmoothed);

   matrix curC    = cura->getCovMat();
   matrix curFt   = cura->getPropTranMat();
   matrix preC    = prea->getCovMat();
   matrix spreC   = sprea->getCovMat();
   matrix preCinv = matrix(matrix::kInverted, preC);
   matrix curA    = curC * curFt * preCinv;
   matrix curAt   = matrix(matrix::kTransposed, curA);
   matrix scurC   = curC + curA * (spreC - preC) * curAt;

   matrix sv = *cura + curA * (*sprea - *prea);
   trackState *newSV = new trackState(sv, scurC, kSmoothed);
   newSV->setPathLocation(cura->getPathLocation());
   newSV->setPathTime(cura->getPathTime());
   Add(newSV);
   SetOwner();
}

trackState* trackSite::getStatePtr(trackSite::siteType t){
  trackState *ap = 0;
  if (t >= 0 && t < GetEntries()) {
    ap = static_cast<trackState*>(UncheckedAt(t));
  }
  return ap;
}



ClassImp(trackSite)

























