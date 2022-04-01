#include "math.h"
#include <vector>
#include <cassert>
#include "TF1.h"
using namespace std;

#include "TMath.h"
#include "TRandom3.h"

#include "trackFinder.h"

trackFinder::trackFinder(particleId id, double mass, double charge, target* uTarget, sft* uSft, tof1* uTof1, ac* uAc, mwpc* mwpc2, mwpc* mwpc3, mwpc* mwpc4, tof2* uTof2) : particle(id, mass, charge), fTarget(uTarget), fSft(uSft), fTof1(uTof1), fAc(uAc), fMwpc2(mwpc2), fMwpc3(mwpc3), fMwpc4(mwpc4), fTof2(uTof2), fNSeeds(0),fChi2PerMdimCut(kChi2PerMdimCut), fNGoodTracks(0){
  fCoarseTracks = new TClonesArray("trackSystem", MAXNTRACKS, kTRUE);
  vector<seed> targetMwpc34Seed;
  fSeedPool[kTargetMwpc34] = targetMwpc34Seed;
  vector<seed> sftMwpc34Seed;
  fSeedPool[kSftMwpc34] = sftMwpc34Seed;
  vector<seed> mwpc2Mwpc34Seed;
  fSeedPool[kMwpc2Mwpc34] = mwpc2Mwpc34Seed;
}


trackFinder::~trackFinder(){
  Clear();

  fCoarseTracks->SetOwner(kTRUE);
  delete fCoarseTracks;

}

void trackFinder::Clear( Option_t* opt ){
  map<seedType, vector<seed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
  fSeedPool.clear();

  map<Int_t, vector<planeHit> >::iterator it;
  for (it = fGoodHits.begin(); it != fGoodHits.end(); it++) { (it->second).clear(); }
  fGoodHits.clear();
}

Bool_t trackFinder::processHits(TClonesArray* theTracks){  
  if(!findSeed()) return kFALSE;
  if(!trackFilter()) return kFALSE;
  trackSmoothBack();
  if(!finalSelection(theTracks)) return kFALSE;

  return kTRUE;
}

//______________________________________________________________________________
Bool_t trackFinder::findSeed(){ // Three types of seeds: kTargetMwpc34, kSftMwpc34 and kTargetMwpc34; All seeds are stored in the pool
  if(fTarget->getTag()==0) return kFALSE;

  if(fTof2->getGapNum()<0 || fTof2->getGapNum()>11) return kFALSE;
  else{
    if(fMwpc4->getNHits(fTof2->getGapNum()) == 0 || fMwpc3->getNHits(fTof2->getGapNum()) == 0)  return kFALSE;
  }

  map<seedType, vector<seed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
  fSeedPool.clear();

  seedType type;
  vector<planeHit> vectHit;
  Double_t initCov[5];

  for(int i=0; i<fMwpc4->getNHits(fTof2->getGapNum()); i++){
    for(int j=0; j<fMwpc3->getNHits(fTof2->getGapNum()); j++){
      planeHit *hitC=(planeHit*)fMwpc4->getHits(fTof2->getGapNum())->At(i);
      planeHit *hitB=(planeHit*)fMwpc3->getHits(fTof2->getGapNum())->At(j); 
      calcInitCoV(hitB, hitC,initCov);

      //--------------------- Begin: seeds with type of kMwpc2Mwpc34 --------------------------
      type=kMwpc2Mwpc34;
      for(int k=0; k<fMwpc2->getNHits(fTof2->getGapNum());k++){
	planeHit* hitA=(planeHit*)fMwpc2->getHits(fTof2->getGapNum())->At(k);
	vectHit.clear();
	vectHit.push_back(*hitA);
	vectHit.push_back(*hitB);
	vectHit.push_back(*hitC);


	TVector3 initPos(hitC->getPosition().X(),hitC->getPosition().Y(),hitC->getPosition().Z());
	TVector3 initDir(hitC->getPosition().X() - hitB->getPosition().X(), hitC->getPosition().Y() - hitB->getPosition().Y(), hitC->getPosition().Z() - hitB->getPosition().Z());
	initDir=initDir.Unit();
	Double_t initTheta=initDir.Theta() -  (0.0002519+0.01303*atan(0.534*(hitC->getPosition().X() - hitB->getPosition().X())-1.353)) - (0.0008226 -2.224e-6*(hitC->getPosition().Y() - hitB->getPosition().Y()));
	Double_t initPhi=initDir.Phi() - (0.0001189 - 0.001142*(hitC->getPosition().Y() - hitB->getPosition().Y()));
	Double_t initMom;
	if(!determineInitMom(fTof2->getGapNum(), type, vectHit, initPos, initTheta, initPhi, initCov, initMom)) continue;	  
	map<seedType, vector<seed> >::iterator it=fSeedPool.find(type);
	if(it!=fSeedPool.end()){
	  (it->second).push_back(seed(fTof2->getGapNum(), type, vectHit, initPos, initTheta, initPhi, initMom, fParticleCharge, initCov));
	}
	else{
	  vector<seed> thisVector;
	  thisVector.push_back(seed(fTof2->getGapNum(), type, vectHit, initPos, initTheta, initPhi, initMom, fParticleCharge, initCov));
	  fSeedPool.insert(pair<seedType,vector<seed> >(type,thisVector));
	}	
      } 
      //--------------------- End: seeds with type of kMwpc2Mwpc34 --------------------------

    }
  }

  map<seedType, vector<seed> >::iterator it;
  Bool_t isEmpty=kTRUE;
  for(it=fSeedPool.begin();it!=fSeedPool.end();it++){
    if(it->second.size()!=0){
      isEmpty=kFALSE;
      break;
    }
  }
  if(isEmpty == kTRUE) return kFALSE;
  else return kTRUE;
}

//______________________________________________________________________________
Double_t trackFinder::fastTracker(vector<planeHit> hits){
  const Double_t kappa = 0.299792458;
  Double_t a1, a2, b2;
  Double_t z0, x0, r;
  Double_t c2x,c2z,c3x,c3z,c4x,c4z;
  Double_t mom;

  c2x=hits[0].getPosition().X();
  c2z=hits[0].getPosition().Z();
  c3x=hits[1].getPosition().X();
  c3z=hits[1].getPosition().Z();
  c4x=hits[2].getPosition().X();
  c4z=hits[2].getPosition().Z();

  Double_t c3y=hits[1].getPosition().Y();
  Double_t c4y=hits[2].getPosition().Y();

  a1 = (c4x-c3x)/(c4z-c3z);
  a2 = -1/a1;
  b2 = c3x - a2*c3z;
  z0 = ((c3z*c3z + c3x*c3x) - (c2z*c2z + c2x*c2x) + 2*b2*(c2x-c3x)) / 2 / ((c3z-c2z) + a2*(c3x-c2x));
  x0 = a2*z0 + b2;
  r = sqrt((z0-c3z)*(z0-c3z) + (x0-c3x)*(x0-c3x));

  TVector3 dir(c4x-c3x, c4y-c3y, c4z-c3z);
  dir=dir.Unit();

  TF1 *func=new TF1("func","1.313-0.4891*exp(-158.4/(x-518.4))",525,10000);
  if(r>525) mom=kappa*r*func->Eval(r)/sqrt(dir.X()*dir.X()+dir.Z()*dir.Z());
  else mom=kappa*r*func->Eval(525)/sqrt(dir.X()*dir.X()+dir.Z()*dir.Z());
  delete func;

  return mom;
}


//______________________________________________________________________________
Bool_t trackFinder::determineInitMom(Int_t gapNumToF2, seedType type, vector<planeHit> hits, TVector3 initPos, Double_t initTheta, Double_t initPhi, Double_t* initCov, Double_t &mom){

  Double_t initMom[5], chi2PerNDF[5];
  seed* gSeed[5];
  Double_t fInitMom[5],fChi2PerNDF[5];

  Double_t temp=0;
  Int_t index=0;

  Double_t momWithMinimChi2;
  Bool_t flagFastTracker; // fast tracker just can be applied for seedType::kMwpc2Mwpc34
  if(type == kMwpc2Mwpc34){
    flagFastTracker = kTRUE;
    momWithMinimChi2 = fastTracker(hits);
  }
  else{
    momWithMinimChi2=220;
    flagFastTracker = kFALSE;
  }

  if(flagFastTracker == kFALSE){
    for(Int_t j=0;j<5;j++){
      initMom[j]=momWithMinimChi2+(j-2)*pow(2, 5);
      gSeed[j]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[j], fParticleCharge, initCov);
      if(!seedTrackTry(gSeed[j], chi2PerNDF[j])) chi2PerNDF[j]=1e12;
      delete gSeed[j];
      gSeed[j]=NULL;
    }
    temp=chi2PerNDF[0];
    index=0;
    for(int j=1;j<5;j++){
      if(temp>chi2PerNDF[j]){
	temp=chi2PerNDF[j];
	index=j;
      }
    }
    for(int j=0;j<5;j++){
      //cout<<initMom[j]<<"  "<<chi2PerNDF[j]<<endl;
      fInitMom[j]=initMom[j];
      fChi2PerNDF[j]=chi2PerNDF[j];
    }

    for(Int_t i=6;i>=0;i--){
      momWithMinimChi2=fInitMom[index];
      if(index==0){
	initMom[0]=momWithMinimChi2-2*pow(2, i)*0.25;
	gSeed[0]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[0], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[0], chi2PerNDF[0])) chi2PerNDF[0]=1e12;
	delete gSeed[0];
	gSeed[0]=NULL;

	initMom[2]=momWithMinimChi2;
	chi2PerNDF[2]=fChi2PerNDF[index];

	initMom[4]=fInitMom[index+1];
	chi2PerNDF[4]=fChi2PerNDF[index+1];

	initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	delete gSeed[1];
	gSeed[1]=NULL;

	initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	delete gSeed[3];
	gSeed[3]=NULL;
      }

      else if(index==4){
	initMom[4]=momWithMinimChi2+2*pow(2, i)*0.25;
	gSeed[4]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[4], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[4], chi2PerNDF[4])) chi2PerNDF[4]=1e12;
	delete gSeed[4];
	gSeed[4]=NULL;

	initMom[2]=momWithMinimChi2;
	chi2PerNDF[2]=fChi2PerNDF[index];

	initMom[0]=fInitMom[index-1];
	chi2PerNDF[0]=fChi2PerNDF[index-1];

	initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	delete gSeed[1];
	gSeed[1]=NULL;

	initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	delete gSeed[3];
	gSeed[3]=NULL;
      }

      else{
	initMom[2]=momWithMinimChi2;
	chi2PerNDF[2]=fChi2PerNDF[index];

	initMom[0]=fInitMom[index-1];
	chi2PerNDF[0]=fChi2PerNDF[index-1];

	initMom[4]=fInitMom[index+1];
	chi2PerNDF[4]=fChi2PerNDF[index+1];

	initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	delete gSeed[1];
	gSeed[1]=NULL;

	initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	delete gSeed[3];
	gSeed[3]=NULL;
      }
      temp=chi2PerNDF[0];
      index=0;
      for(int j=1;j<5;j++){
	if(temp>chi2PerNDF[j]){
	  temp=chi2PerNDF[j];
	  index=j;
	}
      }
      for(int j=0;j<5;j++){
	//cout<<initMom[j]<<"  "<<chi2PerNDF[j]<<endl;
	fInitMom[j]=initMom[j];
	fChi2PerNDF[j]=chi2PerNDF[j];
      }
    }

    if(temp>1e8 || temp < 0) return kFALSE;
    else{
      mom=fInitMom[index];   
      return kTRUE;
    }
  }

  else{
    // Test if the best momentum candidate is between result from fastTracker +/- 4
    for(Int_t j=0;j<5;j++){
      initMom[j]=momWithMinimChi2+(j-2)*4;
      gSeed[j]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[j], fParticleCharge, initCov);
      if(!seedTrackTry(gSeed[j], chi2PerNDF[j])) chi2PerNDF[j]=1e12;
      delete gSeed[j];
      gSeed[j]=NULL;
    }
    temp=chi2PerNDF[0];
    index=0;
    for(int j=1;j<5;j++){
      if(temp>chi2PerNDF[j]){
	temp=chi2PerNDF[j];
	index=j;
      }
    }
    for(int j=0;j<5;j++){
      fInitMom[j]=initMom[j];
      fChi2PerNDF[j]=chi2PerNDF[j];
    }

    // If yes, apply the new momentum determination way
    if((index != 0) && (index != 4)){
      for(Int_t i=3;i>=0;i--){
	momWithMinimChi2=fInitMom[index];
	if(index==0){
	  initMom[0]=momWithMinimChi2-2*pow(2, i)*0.25;
	  gSeed[0]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[0], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[0], chi2PerNDF[0])) chi2PerNDF[0]=1e12;
	  delete gSeed[0];
	  gSeed[0]=NULL;

	  initMom[2]=momWithMinimChi2;
	  chi2PerNDF[2]=fChi2PerNDF[index];

	  initMom[4]=fInitMom[index+1];
	  chi2PerNDF[4]=fChi2PerNDF[index+1];

	  initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	  gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	  delete gSeed[1];
	  gSeed[1]=NULL;

	  initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	  gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	  delete gSeed[3];
	  gSeed[3]=NULL;
	}

	else if(index==4){
	  initMom[4]=momWithMinimChi2+2*pow(2, i)*0.25;
	  gSeed[4]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[4], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[4], chi2PerNDF[4])) chi2PerNDF[4]=1e12;
	  delete gSeed[4];
	  gSeed[4]=NULL;

	  initMom[2]=momWithMinimChi2;
	  chi2PerNDF[2]=fChi2PerNDF[index];

	  initMom[0]=fInitMom[index-1];
	  chi2PerNDF[0]=fChi2PerNDF[index-1];

	  initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	  gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	  delete gSeed[1];
	  gSeed[1]=NULL;

	  initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	  gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	  delete gSeed[3];
	  gSeed[3]=NULL;
	}

	else{
	  initMom[2]=momWithMinimChi2;
	  chi2PerNDF[2]=fChi2PerNDF[index];

	  initMom[0]=fInitMom[index-1];
	  chi2PerNDF[0]=fChi2PerNDF[index-1];

	  initMom[4]=fInitMom[index+1];
	  chi2PerNDF[4]=fChi2PerNDF[index+1];

	  initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	  gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	  delete gSeed[1];
	  gSeed[1]=NULL;

	  initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	  gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	  delete gSeed[3];
	  gSeed[3]=NULL;
	}
	temp=chi2PerNDF[0];
	index=0;
	for(int j=1;j<5;j++){
	  if(temp>chi2PerNDF[j]){
	    temp=chi2PerNDF[j];
	    index=j;
	  }
	}
	for(int j=0;j<5;j++){
	  //cout<<initMom[j]<<"  "<<chi2PerNDF[j]<<endl;
	  fInitMom[j]=initMom[j];
	  fChi2PerNDF[j]=chi2PerNDF[j];
	}
      }

      if(temp>1e8 || temp < 0) return kFALSE;
      else{
	mom=fInitMom[index];   
	return kTRUE;
      }
    }

    // If no, apply the old momentum determination way
    else{
      momWithMinimChi2=220;
      for(Int_t j=0;j<5;j++){
	initMom[j]=momWithMinimChi2+(j-2)*pow(2, 5);
	gSeed[j]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[j], fParticleCharge, initCov);
	if(!seedTrackTry(gSeed[j], chi2PerNDF[j])) chi2PerNDF[j]=1e12;
	delete gSeed[j];
	gSeed[j]=NULL;
      }
      temp=chi2PerNDF[0];
      index=0;
      for(int j=1;j<5;j++){
	if(temp>chi2PerNDF[j]){
	  temp=chi2PerNDF[j];
	  index=j;
	}
      }
      for(int j=0;j<5;j++){
	//cout<<initMom[j]<<"  "<<chi2PerNDF[j]<<endl;
	fInitMom[j]=initMom[j];
	fChi2PerNDF[j]=chi2PerNDF[j];
      }

      for(Int_t i=6;i>=0;i--){
	momWithMinimChi2=fInitMom[index];
	if(index==0){
	  initMom[0]=momWithMinimChi2-2*pow(2, i)*0.25;
	  gSeed[0]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[0], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[0], chi2PerNDF[0])) chi2PerNDF[0]=1e12;
	  delete gSeed[0];
	  gSeed[0]=NULL;

	  initMom[2]=momWithMinimChi2;
	  chi2PerNDF[2]=fChi2PerNDF[index];

	  initMom[4]=fInitMom[index+1];
	  chi2PerNDF[4]=fChi2PerNDF[index+1];

	  initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	  gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	  delete gSeed[1];
	  gSeed[1]=NULL;

	  initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	  gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	  delete gSeed[3];
	  gSeed[3]=NULL;
	}

	else if(index==4){
	  initMom[4]=momWithMinimChi2+2*pow(2, i)*0.25;
	  gSeed[4]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[4], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[4], chi2PerNDF[4])) chi2PerNDF[4]=1e12;
	  delete gSeed[4];
	  gSeed[4]=NULL;

	  initMom[2]=momWithMinimChi2;
	  chi2PerNDF[2]=fChi2PerNDF[index];

	  initMom[0]=fInitMom[index-1];
	  chi2PerNDF[0]=fChi2PerNDF[index-1];

	  initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	  gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	  delete gSeed[1];
	  gSeed[1]=NULL;

	  initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	  gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	  delete gSeed[3];
	  gSeed[3]=NULL;
	}

	else{
	  initMom[2]=momWithMinimChi2;
	  chi2PerNDF[2]=fChi2PerNDF[index];

	  initMom[0]=fInitMom[index-1];
	  chi2PerNDF[0]=fChi2PerNDF[index-1];

	  initMom[4]=fInitMom[index+1];
	  chi2PerNDF[4]=fChi2PerNDF[index+1];

	  initMom[1]=momWithMinimChi2-pow(2, i)*0.25;
	  gSeed[1]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[1], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[1], chi2PerNDF[1])) chi2PerNDF[1]=1e12;
	  delete gSeed[1];
	  gSeed[1]=NULL;

	  initMom[3]=momWithMinimChi2+pow(2, i)*0.25;
	  gSeed[3]=new seed(gapNumToF2, type, hits, initPos, initTheta, initPhi, initMom[3], fParticleCharge, initCov);
	  if(!seedTrackTry(gSeed[3], chi2PerNDF[3])) chi2PerNDF[3]=1e12;
	  delete gSeed[3];
	  gSeed[3]=NULL;
	}
	temp=chi2PerNDF[0];
	index=0;
	for(int j=1;j<5;j++){
	  if(temp>chi2PerNDF[j]){
	    temp=chi2PerNDF[j];
	    index=j;
	  }
	}
	for(int j=0;j<5;j++){
	  //cout<<initMom[j]<<"  "<<chi2PerNDF[j]<<endl;
	  fInitMom[j]=initMom[j];
	  fChi2PerNDF[j]=chi2PerNDF[j];
	}
      }

      if(temp>1e8 || temp < 0) return kFALSE;
      else{
	mom=fInitMom[index];   
	return kTRUE;
      }
    }
  }

}
//______________________________________________________________________________
Bool_t trackFinder::seedTrackTry(seed* gSeed, Double_t &chi2PerNDF){

  Bool_t magFlag;

  trackSystem *thisSystem=new trackSystem();
  trackSite* initSite=siteInit(gSeed);
  thisSystem->SetOwner();
  thisSystem->Add(initSite);
	
  trackSite* currentSite;
  trackState* currentState;
  trackState* predictState;
  trackSite* newSite[gSeed->hits.size()];

  currentSite=thisSystem->getCurSite();
  currentState=currentSite->getCurState();
  predictState=currentState->predictSVatPlane(&gSeed->hits[2], fParticleId, fParticleMass, fParticleCharge, magFlag);

  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  newSite[2] =new trackSite(&gSeed->hits[2], kMdim, kSdim, kChi2PerMdimCut);
  newSite[2]->Add(predictState);
  if(!newSite[2]->filter(fParticleCharge)){
    delete newSite[2];
    delete thisSystem;
    return kFALSE;
  }
  thisSystem->Add(newSite[2]);
  thisSystem->increaseChi2(newSite[2]->getDeltaChi2());

  currentSite=thisSystem->getCurSite();
  currentState=currentSite->getCurState();
  predictState=currentState->predictSVatFirstPlane(fMwpc4->getFrontPlaneHit(gSeed->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc4Mat);
  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  predictState->predictSVatNextPlane(fMwpc3->getBackPlaneHit(gSeed->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kAir);
  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  predictState->predictSVatLastPlane(&gSeed->hits[1], fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc3Mat);
  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  newSite[1] =new trackSite(&gSeed->hits[1], kMdim, kSdim, kChi2PerMdimCut);
  newSite[1]->Add(predictState);
  if(!newSite[1]->filter(fParticleCharge)){
    delete newSite[1];
    delete thisSystem;
    return kFALSE;
  }
  thisSystem->Add(newSite[1]);
  thisSystem->increaseChi2(newSite[1]->getDeltaChi2());

  currentSite=thisSystem->getCurSite();
  currentState=currentSite->getCurState();
  predictState=currentState->predictSVatFirstPlane(fMwpc3->getFrontPlaneHit(gSeed->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc3Mat);
  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  predictState->predictSVatNextPlane(fMwpc2->getBackPlaneHit(gSeed->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kAir);
  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  predictState->predictSVatLastPlane(&gSeed->hits[0], fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc2Mat);
  if(magFlag == kFALSE){
    delete predictState;
    delete thisSystem;
    return kFALSE;
  }
  newSite[0] =new trackSite( &gSeed->hits[0], kMdim, kSdim, kChi2PerMdimCut);
  newSite[0]->Add(predictState);
  if(!newSite[0]->filter(fParticleCharge)){
    delete newSite[0];
    delete thisSystem;
    return kFALSE;
  }
  thisSystem->Add(newSite[0]);
  thisSystem->increaseChi2(newSite[0]->getDeltaChi2());

  planeHit *planeBLowEdge=new planeHit(kBLowEdge,kPlaneDirX,300);
  currentSite=thisSystem->getCurSite();
  currentState=currentSite->getCurState();
  trackState* predictStateBLowEdge=currentState->predictSVatPlane(planeBLowEdge,fParticleId, fParticleMass, fParticleCharge, magFlag, kAir);
  Double_t chi2Vert;

  if(!calcVertChi2(gSeed->fGapNumTof2, predictStateBLowEdge, chi2Vert)){
    delete planeBLowEdge;
    delete predictStateBLowEdge;
    delete thisSystem;
    return kFALSE;
  }

  chi2PerNDF=(thisSystem->getChi2()+chi2Vert)/(thisSystem->getNDF()+3);
  delete planeBLowEdge;
  delete predictStateBLowEdge;
  delete thisSystem;
 
  if(chi2PerNDF<0) return kFALSE;
  else return kTRUE;
}
//______________________________________________________________________________
trackSite* trackFinder::siteInit(seed* thisSeed){
  assert(thisSeed->hits.size()>2 && thisSeed->hits[thisSeed->hits.size()-1].getPlaneId() == kMwpc4 && thisSeed->hits[thisSeed->hits.size()-2].getPlaneId() == kMwpc3);
  //-----------prepare seeds for Kalman Filter track finding------------//
  matrix svd(kSdim,1);
  svd.Zero();
  svd(kIdxX,0) = (thisSeed->initPos).X();
  svd(kIdxY,0) = (thisSeed->initPos).Y();
  svd(kIdxZ,0) = (thisSeed->initPos).Z();
  svd(kIdxTheta,0) = thisSeed->initTheta;
  svd(kIdxPhi,0) = thisSeed->initPhi;
  svd(kIdxQP,0) = thisSeed->charge/thisSeed->initMom;
  
  matrix C(kSdim,kSdim);
  C.Zero();
  C(kIdxX, kIdxX) = thisSeed->initCov[0]; 
  C(kIdxY, kIdxY) = thisSeed->initCov[1];
  C(kIdxZ, kIdxZ) = thisSeed->initCov[2];
  C(kIdxTheta, kIdxTheta) = thisSeed->initCov[3];
  C(kIdxPhi, kIdxPhi) = thisSeed->initCov[4];
  C(kIdxQP, kIdxQP) = pow(1/thisSeed->initMom/thisSeed->initMom*0.25, 2);

  trackSite* initSite = new trackSite(&thisSeed->hits[thisSeed->hits.size()-1], kMdim, kSdim, kMdim*fChi2PerMdimCut);
  initSite->SetOwner();
  initSite->setPathLocation(0);
  initSite->setPathTime(0);
  trackState* trackStatePredicted=new trackState(svd, C, *initSite, trackSite::kPredicted);
  trackStatePredicted->setPathLocation(0);
  trackStatePredicted->setPathTime(0);
  initSite->Add(trackStatePredicted);
  trackState* trackStateFiltered=new trackState(svd, C, *initSite, trackSite::kFiltered);
  trackStateFiltered->setPathLocation(0);
  trackStateFiltered->setPathTime(0);
  initSite->Add(trackStateFiltered);
  
  return initSite; //Initial site is a virtual site, and its location is the same as MWPC4
}
//______________________________________________________________________________
void trackFinder::calcInitCoV(planeHit* hitMwpc3, planeHit* hitMwpc4, Double_t initCov[]){
  assert(hitMwpc3->getPlaneId() == kMwpc3 && hitMwpc4->getPlaneId() == kMwpc4);
  TVector3 posMwpc3=hitMwpc3->getPosition();
  TVector3 posErrMwpc3=hitMwpc3->getPositionError();
  Double_t xMwpc3=posMwpc3.X();
  Double_t xErrMwpc3=posErrMwpc3.X();
  Double_t yMwpc3=posMwpc3.Y();
  Double_t yErrMwpc3=posErrMwpc3.Y();
  Double_t zMwpc3=posMwpc3.Z();
  Double_t zErrMwpc3=posErrMwpc3.Z();

  TVector3 posMwpc4=hitMwpc4->getPosition();
  TVector3 posErrMwpc4=hitMwpc4->getPositionError();
  Double_t xMwpc4=posMwpc4.X();
  Double_t xErrMwpc4=posErrMwpc4.X();
  Double_t yMwpc4=posMwpc4.Y();
  Double_t yErrMwpc4=posErrMwpc4.Y();
  Double_t zMwpc4=posMwpc4.Z();
  Double_t zErrMwpc4=posErrMwpc4.Z();

  Double_t deltaX=xMwpc4-xMwpc3, deltaY=yMwpc4-yMwpc3, deltaZ=zMwpc4-zMwpc3;
  Double_t deltaXErr2=pow(xErrMwpc4,2)+pow(xErrMwpc3,2), deltaYErr2=pow(yErrMwpc4,2)+pow(yErrMwpc3,2), deltaZErr2=pow(zErrMwpc4,2)+pow(zErrMwpc3,2);

  Double_t r2=pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ,2);
  Double_t r=sqrt(r2);
  Double_t cosTheta2=pow(deltaZ,2)/r2;
  Double_t partThetaDeltaX=1/sqrt(1-cosTheta2) * deltaX * deltaZ * pow(r,-3);
  Double_t partThetaDeltaY=1/sqrt(1-cosTheta2) * deltaY * deltaZ * pow(r,-3);
  Double_t partThetaDeltaZ=-1/sqrt(1-cosTheta2) * (1/r - pow(deltaZ,2) * pow(r,-3));
  Double_t varTheta=pow(partThetaDeltaX,2)*deltaXErr2 + pow(partThetaDeltaY,2)*deltaYErr2 + pow(partThetaDeltaZ,2)*deltaZErr2;

  Double_t tanPhi2=pow(deltaY,2)/pow(deltaX,2);
  Double_t partPhiDeltaX=1/(1+tanPhi2)*(-deltaY/pow(deltaX,2));
  Double_t partPhiDeltaY=1/(1+tanPhi2)/deltaX;
  Double_t varPhi=pow(partPhiDeltaX,2)*deltaXErr2 + pow(partPhiDeltaY,2)*deltaYErr2;

  initCov[0]=pow(xErrMwpc4,2);
  initCov[1]=pow(yErrMwpc4,2);
  initCov[2]=pow(zErrMwpc4,2);
  initCov[3]=varTheta+1e-4;
  initCov[4]=varPhi+1e-2;
}
//______________________________________________________________________________
Bool_t trackFinder::trackFilter(){   
  //Construct tracks from MWPC4 to Target
  //For a seed,  a track is attempted, and the best hits for detectors excluding the seed is choosen. So basically, one seed correponds one track.
  //All candidate tracks are stored in fCoarseTracks
  Bool_t magFlag;
  if(fCoarseTracks->GetEntries()!=0) fCoarseTracks->Clear();
  fNSeeds=0;
  map<seedType, vector<seed> >::iterator it;
  for(it=fSeedPool.begin();it!=fSeedPool.end();it++){

    //--------------------- Begin: seeds with type of kTargetMwpc34 --------------------------
    if(it->first == kTargetMwpc34){
    }
    //--------------------- End: seeds with type of kTargetMwpc34 --------------------------

    //--------------------- Begin: seeds with type of kSftMwpc34 --------------------------
    if(it->first == kTargetMwpc34){
    }
    //--------------------- End: seeds with type of kSftMwpc34 --------------------------

    //--------------------- Begin: seeds with type of kMwpc2Mwpc34 --------------------------
    if(it->first == kMwpc2Mwpc34){
      vector<seed> thisVector=(it->second);
      for(unsigned i=0;i<thisVector.size();i++){
	if(!thisVector[i].isActive) continue;
	Bool_t flagM=0, flagS=0, flagT=0;
	trackSystem *thisSystem=new ((*fCoarseTracks)[fNSeeds++]) trackSystem();
	thisSystem->setGapNumTof2(thisVector[i].getGapNumTof2());
	fCoarseTracks->SetOwner(kTRUE);
	thisSystem->setSeedType(thisVector[i].type);
	thisSystem->SetOwner();
	trackSite* initSite=siteInit(&(thisVector[i]));
	thisSystem->Add(initSite);

	trackSite* currentSite;
	trackState* currentState;
	trackState* predictState;
	trackSite* newSite;

	thisSystem->checkTrackStatus();
	if(!thisSystem->getTrackStatus()) continue;
	currentSite=thisSystem->getCurSite();
	currentState=currentSite->getCurState();
	predictState=currentState->predictSVatPlane(&thisVector[i].hits[2], fParticleId, fParticleMass, fParticleCharge, magFlag);
	newSite =new trackSite( &thisVector[i].hits[2], kMdim, kSdim, kMdim*fChi2PerMdimCut);
	newSite->Add(predictState);
	if(newSite->filter(fParticleCharge)){
	  thisSystem->Add(newSite);
	  thisSystem->increaseChi2(newSite->getDeltaChi2());
	}
	else{      
	  thisSystem->addMissingHits();
	  delete newSite;  
	}

	thisSystem->checkTrackStatus();
	if(!thisSystem->getTrackStatus()) continue;
	currentSite=thisSystem->getCurSite();
	currentState=currentSite->getCurState();
	predictState=currentState->predictSVatFirstPlane(fMwpc4->getFrontPlaneHit(thisVector[i].getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc4Mat);
	predictState->predictSVatNextPlane(fMwpc3->getBackPlaneHit(thisVector[i].getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kAir);
	predictState->predictSVatLastPlane(&thisVector[i].hits[1], fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc3Mat);
	newSite =new trackSite( &thisVector[i].hits[1], kMdim, kSdim, kMdim*fChi2PerMdimCut);
	newSite->Add(predictState);
	if(newSite->filter(fParticleCharge)){
	  thisSystem->Add(newSite);
	  thisSystem->increaseChi2(newSite->getDeltaChi2());
	}
	else{    
	  thisSystem->addMissingHits();
	  delete newSite;  
	}

	thisSystem->checkTrackStatus();
	if(!thisSystem->getTrackStatus()) continue;
	currentSite=thisSystem->getCurSite();
	currentState=currentSite->getCurState();
	predictState=currentState->predictSVatFirstPlane(fMwpc3->getFrontPlaneHit(thisVector[i].getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc3Mat);
	predictState->predictSVatNextPlane(fMwpc2->getBackPlaneHit(thisVector[i].getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kAir);
	predictState->predictSVatLastPlane(&thisVector[i].hits[0], fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc2Mat);
	newSite =new trackSite(&thisVector[i].hits[0], kMdim, kSdim, kMdim*fChi2PerMdimCut);
	newSite->Add(predictState);
	if(newSite->filter(fParticleCharge)){
	  thisSystem->Add(newSite);
	  thisSystem->increaseChi2(newSite->getDeltaChi2());
	}
	else{    
	  flagM=1;
	  thisSystem->addMissingHits();
	  delete newSite;  
	}

	thisSystem->checkTrackStatus();
	if(!thisSystem->getTrackStatus()) continue;
	if(fSft->getNHits() == 0) {
	  flagS=1;
	  thisSystem->addMissingHits();
	}
	else{
	  currentSite=thisSystem->getCurSite();
	  currentState=currentSite->getCurState();
	  Int_t bestHitIndex;
	  if(!chooseBestHit(currentState, fSft->getHits(), bestHitIndex)){//have problem!!!!
	    flagS=1;
	    thisSystem->addMissingHits();
	  }
	  else{
	    predictState=currentState->predictSVatPlane((planeHit*)fSft->getHits()->At(bestHitIndex),fParticleId, fParticleMass, fParticleCharge, magFlag);
	    newSite =new trackSite((planeHit*)fSft->getHits()->At(bestHitIndex), kMdim, kSdim, kMdim*fChi2PerMdimCut);
	    newSite->Add(predictState);
	    if(newSite->filter(fParticleCharge)){
	      thisSystem->Add(newSite);
	      thisSystem->increaseChi2(newSite->getDeltaChi2());
	    }
	    else{     
	      flagS=1;
	      thisSystem->addMissingHits();
	      delete newSite;  
	    }
	  }
	}
	
	thisSystem->checkTrackStatus();
	if(!thisSystem->getTrackStatus()) continue;
	if(fTarget->getNHits() == 0){
	  flagT=1;
	  thisSystem->addMissingHits();
	}
	else{
	  currentSite=thisSystem->getCurSite();
	  currentState=currentSite->getCurState();
	  Int_t bestHitIndex;
	  if(!chooseBestHit(currentState, fTarget->getHits(), bestHitIndex)){ //have problem!!!!
	    flagT=1;
	    thisSystem->addMissingHits();
	  }
	  else{
	    predictState=currentState->predictSVatPlane((planeHit*)fTarget->getHits()->At(bestHitIndex), fParticleId, fParticleMass, fParticleCharge, magFlag);
	    newSite =new trackSite((planeHit*)fTarget->getHits()->At(bestHitIndex), kMdim, kSdim, kMdim*fChi2PerMdimCut);
	    newSite->Add(predictState);
	    if(newSite->filter(fParticleCharge)){
	      thisSystem->Add(newSite);
	      thisSystem->increaseChi2(newSite->getDeltaChi2());
	    }
	    else{  
	      flagT=1;
	      thisSystem->addMissingHits();
	      delete newSite;  
	    }
	  }
	}
	thisSystem->checkTrackStatus();
	if(flagM ==0 && flagS==0 && flagT == 0) thisSystem->setTrackId(kMST);
	else if(flagM ==0 && flagS==0 && flagT == 1) thisSystem->setTrackId(kMS);
	else if(flagM ==0 && flagS==1 && flagT == 0) thisSystem->setTrackId(kMT);
	else if(flagM ==1 && flagS==0 && flagT == 0) thisSystem->setTrackId(kST);
	else if(flagM ==0 && flagS==1 && flagT == 1) thisSystem->setTrackId(kM);
	else if(flagM ==1 && flagS==0 && flagT == 1) thisSystem->setTrackId(kS);
	else if(flagM ==1 && flagS==1 && flagT == 0) thisSystem->setTrackId(kT);
      }
    }

    //--------------------- End: seeds with type of kMwpc2Mwpc34 --------------------------
  }

  if(fCoarseTracks->GetEntries()==0) return kFALSE;
  else return kTRUE;
}    
//______________________________________________________________________________
Bool_t trackFinder::chooseBestHit(trackState* currentState, TClonesArray* hits, Int_t &bestHitIndex){
  Bool_t magFlag;
  assert(currentState->getSiteType()!=trackSite::kPredicted);

  Double_t deltaChi2=kMdim*fChi2PerMdimCut;
  for(int i=0;i<hits->GetEntries();i++){
    trackState *predictState;
    predictState=currentState->predictSVatPlane((planeHit*)hits->At(i), fParticleId, fParticleMass, fParticleCharge, magFlag);
    trackSite *newSite =new trackSite((planeHit*)hits->At(i), kMdim, kSdim, kMdim*fChi2PerMdimCut);
    newSite->Add(predictState);
    if(newSite->filter(fParticleCharge)){
      if(deltaChi2 > newSite->getDeltaChi2()){
	deltaChi2=newSite->getDeltaChi2();
	bestHitIndex = i;
      }
    }
  }

  if(deltaChi2 == kMdim*fChi2PerMdimCut) return kFALSE;
  else return kTRUE;
}

void trackFinder::trackSmoothBack(){
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    trackSystem *thisSystem = (trackSystem*)(fCoarseTracks->At(i));
    if ( !(thisSystem->getTrackStatus()) ) continue; 
    thisSystem->smoothBackTo(1);
  }
}

Bool_t trackFinder::finalSelection(TClonesArray* theTracks){
  Bool_t magFlag;
  fNGoodTracks=0;
  fCoarseTracks->Sort();
  for(int i=0;i<fCoarseTracks->GetEntries();i++){
    trackSystem *thisSystem = (trackSystem*)(fCoarseTracks->At(i));
    if ( !(thisSystem->getTrackStatus()) ) continue; 
    trackId uTrackId=thisSystem->getTrackId();
    if(uTrackId == kM){
      trackSite* siteMwpc2=(trackSite*)thisSystem->At(3);
      trackState* stateMwpc2=siteMwpc2->getStatePtr(trackSite::kSmoothed);

      planeHit *planeHitBLowEdge=new planeHit(kMwpc2, kPlaneDirX, 300);
      trackState *predictStateBLowEdge=new trackState(trackSite::kFiltered);
      stateMwpc2->predictSVatPlane(predictStateBLowEdge,planeHitBLowEdge,fParticleId, fParticleMass, fParticleCharge, magFlag);

      trackState *predictStateVert=new trackState(trackSite::kFiltered);

      TVector3 measVert, measErrVert;
      Double_t phiTrack;
      Double_t chi2Vert=calcVertChi2(thisSystem->getGapNumTof2(), predictStateBLowEdge, predictStateVert, measVert, measErrVert, phiTrack);
      //cout<<chi2Vert<<"!!!!"<<endl;
      Int_t flag = 0;
      //start from 1 because the 0th is the dummy site that we used to initialize Kalman Filter
      //One hit cannot be shared by two or more tracks
      for (Int_t j=1;j!=thisSystem->GetLast()+1;j++){
	planeHit thisHit = ((trackSite*)thisSystem->At(j))->getPlaneHit();
	Int_t planeId = thisHit.getPlaneId();
	map< Int_t, vector<planeHit> >::iterator it = fGoodHits.find(planeId);
      
	if (it != fGoodHits.end()){
	  for (UInt_t n = 0; n<(it->second).size(); n++){
	    if (thisHit.getPlaneDir()==kPlaneDirX){
	      if ((thisHit.getY() == ((it->second).at(n)).getY()) || 
		  (thisHit.getZ() == ((it->second).at(n)).getZ())) { flag = 1; }
	    }
	    else{
	      if ((thisHit.getX() == ((it->second).at(n)).getX()) || 
		  (thisHit.getY() == ((it->second).at(n)).getY())) { flag = 1; }
	    }
	  }
	}
      }

      if (flag == 0){
	track* newTrack = new ((*theTracks)[fNGoodTracks++]) track();
	copyTrack(newTrack, thisSystem, predictStateBLowEdge, chi2Vert, predictStateVert, measVert, measErrVert, phiTrack);

	for (Int_t j=1; j!=thisSystem->GetLast()+1;j++){
	  planeHit thisHit = ((trackSite*)thisSystem->At(j))->getPlaneHit();
	  Int_t planeId = thisHit.getPlaneId();
	  map< Int_t, vector<planeHit> >::iterator it = fGoodHits.find(planeId);
	  if (it != fGoodHits.end()){
	    (it->second).push_back(thisHit);
	  }
	  else{
	    vector<planeHit> thisVector;
	    thisVector.push_back(thisHit);
	    fGoodHits.insert(std::pair<Int_t, vector<planeHit> >(planeId, thisVector));
	  }
	}

      }

      delete planeHitBLowEdge;
      delete predictStateBLowEdge;
      delete predictStateVert;
    }
  }

  if(theTracks->GetEntries()==0) return kFALSE;
  else return kTRUE;
}
//______________________________________________________________________________
void trackFinder::copyTrack(track* newTrack, trackSystem* thisSystem, trackState *predictStateBLowEdge, Double_t chi2Vert, trackState* predictStateVert, TVector3 measVert, TVector3 measErrVert, Double_t phiTrack){
  Bool_t magFlag;

  newTrack->setNHits(thisSystem->getNHits());

  trackSite* siteMwpc4=(trackSite*)thisSystem->At(1);
  trackState* stateMwpc4=siteMwpc4->getStatePtr(trackSite::kSmoothed);

  trackSite* siteMwpc3=(trackSite*)thisSystem->At(2);
  trackState* stateMwpc3=siteMwpc3->getStatePtr(trackSite::kSmoothed);

  trackSite* siteMwpc2=(trackSite*)thisSystem->At(3);
  trackState* stateMwpc2=siteMwpc2->getStatePtr(trackSite::kSmoothed);


  Double_t theta,phi;
  Double_t nx,ny,nz;

  //C3 information
  theta=(*stateMwpc3)(kIdxTheta,0);
  phi=(*stateMwpc3)(kIdxPhi,0);
  nx=sin(theta)*cos(phi);
  ny=sin(theta)*sin(phi);
  nz=cos(theta);
  newTrack->setMwpc3PathLocation(siteMwpc3->getPathLocation());
  newTrack->setMwpc3PathTime(siteMwpc3->getPathTime());
  newTrack->setMwpc3MX(siteMwpc3->getM()(kIdxX,0));
  newTrack->setMwpc3MY(siteMwpc3->getM()(kIdxY,0));
  newTrack->setMwpc3MZ(siteMwpc3->getM()(kIdxZ,0));
  newTrack->setMwpc3MXErr(sqrt(siteMwpc3->getV()(kIdxX,kIdxX)));
  newTrack->setMwpc3MYErr(sqrt(siteMwpc3->getV()(kIdxY,kIdxY)));
  newTrack->setMwpc3MZErr(sqrt(siteMwpc3->getV()(kIdxZ,kIdxZ)));
  newTrack->setMwpc3SX((*stateMwpc3)(kIdxX,0));
  newTrack->setMwpc3SY((*stateMwpc3)(kIdxY,0));
  newTrack->setMwpc3SZ((*stateMwpc3)(kIdxZ,0));
  newTrack->setMwpc3SNx(nx);
  newTrack->setMwpc3SNy(ny);
  newTrack->setMwpc3SNz(nz);
  newTrack->setMwpc3SP(fParticleCharge/(*stateMwpc3)(kIdxQP,0));

  //C2 information
  theta=(*stateMwpc2)(kIdxTheta,0);
  phi=(*stateMwpc2)(kIdxPhi,0);
  nx=sin(theta)*cos(phi);
  ny=sin(theta)*sin(phi);
  nz=cos(theta);
  newTrack->setMwpc2PathLocation(siteMwpc2->getPathLocation());
  newTrack->setMwpc2PathTime(siteMwpc2->getPathTime());
  newTrack->setMwpc2MX(siteMwpc2->getM()(kIdxX,0));
  newTrack->setMwpc2MY(siteMwpc2->getM()(kIdxY,0));
  newTrack->setMwpc2MZ(siteMwpc2->getM()(kIdxZ,0));
  newTrack->setMwpc2MXErr(sqrt(siteMwpc2->getV()(kIdxX,kIdxX)));
  newTrack->setMwpc2MYErr(sqrt(siteMwpc2->getV()(kIdxY,kIdxY)));
  newTrack->setMwpc2MZErr(sqrt(siteMwpc2->getV()(kIdxZ,kIdxZ)));
  newTrack->setMwpc2SX((*stateMwpc2)(kIdxX,0));
  newTrack->setMwpc2SY((*stateMwpc2)(kIdxY,0));
  newTrack->setMwpc2SZ((*stateMwpc2)(kIdxZ,0));
  newTrack->setMwpc2SNx(nx);
  newTrack->setMwpc2SNy(ny);
  newTrack->setMwpc2SNz(nz);
  newTrack->setMwpc2SP(fParticleCharge/(*stateMwpc2)(kIdxQP,0));

  //C3 propagate back to C4, then get C4 information
  trackState *stateMwpc4Prop=stateMwpc3->predictSVatFirstPlane(fMwpc3->getBackPlaneHit(thisSystem->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc3Mat);
  stateMwpc4Prop->predictSVatNextPlane(fMwpc4->getFrontPlaneHit(thisSystem->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge, magFlag, kAir);
  planeHit *hitC4=new planeHit(siteMwpc4->getPlaneHit().getPlaneId(),siteMwpc4->getPlaneHit().getPlaneDir(),siteMwpc4->getPlaneHit().getZ());
  stateMwpc4Prop->predictSVatLastPlane(hitC4, fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc4Mat);
  delete hitC4;

  theta=(*stateMwpc4Prop)(kIdxTheta,0);
  phi=(*stateMwpc4Prop)(kIdxPhi,0);
  nx=sin(theta)*cos(phi);
  ny=sin(theta)*sin(phi);
  nz=cos(theta);
  newTrack->setMwpc4PathLocation(stateMwpc4Prop->getPathLocation());
  newTrack->setMwpc4PathTime(stateMwpc4Prop->getPathTime());
  newTrack->setMwpc4MX(siteMwpc4->getM()(kIdxX,0));
  newTrack->setMwpc4MY(siteMwpc4->getM()(kIdxY,0));
  newTrack->setMwpc4MZ(siteMwpc4->getM()(kIdxZ,0));
  newTrack->setMwpc4MXErr(sqrt(siteMwpc4->getV()(kIdxX,kIdxX)));
  newTrack->setMwpc4MYErr(sqrt(siteMwpc4->getV()(kIdxY,kIdxY)));
  newTrack->setMwpc4MZErr(sqrt(siteMwpc4->getV()(kIdxZ,kIdxZ)));
  newTrack->setMwpc4SX((*stateMwpc4Prop)(kIdxX,0));
  newTrack->setMwpc4SY((*stateMwpc4Prop)(kIdxY,0));
  newTrack->setMwpc4SZ((*stateMwpc4Prop)(kIdxZ,0));
  newTrack->setMwpc4SNx(nx);
  newTrack->setMwpc4SNy(ny);
  newTrack->setMwpc4SNz(nz);
  newTrack->setMwpc4SP(fParticleCharge/(*stateMwpc4Prop)(kIdxQP,0));

  //Chi2 per NDF
  Double_t deltaChi2Mwpc3=pow((*stateMwpc3)(kIdxX,0)-siteMwpc3->getM()(kIdxX,0),2)/siteMwpc3->getV()(kIdxX,kIdxX) + pow((*stateMwpc3)(kIdxY,0)-siteMwpc3->getM()(kIdxY,0),2)/siteMwpc3->getV()(kIdxY,kIdxY);
  Double_t deltaChi2Mwpc2=pow((*stateMwpc2)(kIdxZ,0)-siteMwpc2->getM()(kIdxZ,0),2)/siteMwpc2->getV()(kIdxZ,kIdxZ) + pow((*stateMwpc2)(kIdxY,0)-siteMwpc2->getM()(kIdxY,0),2)/siteMwpc2->getV()(kIdxY,kIdxY);
  Double_t deltaChi2Mwpc4=pow((*stateMwpc4Prop)(kIdxX,0)-siteMwpc4->getM()(kIdxX,0),2)/siteMwpc4->getV()(kIdxX,kIdxX) + pow((*stateMwpc4Prop)(kIdxY,0)-siteMwpc4->getM()(kIdxY,0),2)/siteMwpc4->getV()(kIdxY,kIdxY);
  Double_t chi2PerNDF=(deltaChi2Mwpc2+deltaChi2Mwpc3+deltaChi2Mwpc4+chi2Vert)/(2*3+3-5);
  newTrack->setChi2PerNDF(chi2PerNDF);

  //ToF2 information
  trackState *predictStateTof2=stateMwpc4Prop->predictSVatFirstPlane(fMwpc4->getBackPlaneHit(thisSystem->getGapNumTof2()),fParticleId, fParticleMass, fParticleCharge, magFlag, kMwpc4Mat);
  planeHit *planeHitBUpEdge=new planeHit(kBUpEdge, kPlaneDirZ, 1000);
  predictStateTof2->predictSVatNextPlane(planeHitBUpEdge,fParticleId, fParticleMass, fParticleCharge, magFlag);
  delete planeHitBUpEdge;
  predictStateTof2->predictSVatLastPlaneStraightLine(fTof2->getMiddlePlaneHit(thisSystem->getGapNumTof2()), fParticleId, fParticleMass, fParticleCharge);

  theta=(*predictStateTof2)(kIdxTheta,0);
  phi=(*predictStateTof2)(kIdxPhi,0);
  nx=sin(theta)*cos(phi);
  ny=sin(theta)*sin(phi);
  nz=cos(theta);
  newTrack->setTof2GapNum(thisSystem->getGapNumTof2());
  newTrack->setTof2PathLocation(predictStateTof2->getPathLocation());
  newTrack->setTof2PathTime(predictStateTof2->getPathTime());
  newTrack->setTof2Part(fTof2->getPart());
  newTrack->setTof2TimeI(fTof2->getTimeI());
  newTrack->setTof2TimeO(fTof2->getTimeO());
  newTrack->setTof2TimeMean(fTof2->getTimeMean());
  newTrack->setTof2SX((*predictStateTof2)(kIdxX,0));
  newTrack->setTof2SY((*predictStateTof2)(kIdxY,0));
  newTrack->setTof2SZ((*predictStateTof2)(kIdxZ,0));
  newTrack->setTof2SNx(nx);
  newTrack->setTof2SNy(ny);
  newTrack->setTof2SNz(nz);
  newTrack->setTof2SP(fParticleCharge/(*predictStateTof2)(kIdxQP,0));

  delete predictStateTof2;
  delete stateMwpc4Prop;

  //B lower edge information
  theta=(*predictStateBLowEdge)(kIdxTheta,0);
  phi=(*predictStateBLowEdge)(kIdxPhi,0);
  nx=sin(theta)*cos(phi);
  ny=sin(theta)*sin(phi);
  nz=cos(theta);
  newTrack->setBLowEdgeSX((*predictStateBLowEdge)(kIdxX,0));
  newTrack->setBLowEdgeSY((*predictStateBLowEdge)(kIdxY,0));
  newTrack->setBLowEdgeSZ((*predictStateBLowEdge)(kIdxZ,0));
  newTrack->setBLowEdgeSNx(nx);
  newTrack->setBLowEdgeSNy(ny);
  newTrack->setBLowEdgeSNz(nz);
  newTrack->setBLowEdgeSP(fParticleCharge/(*predictStateBLowEdge)(kIdxQP,0));

  //ToF1 information
  Bool_t flagToFMatch=kFALSE;
  newTrack->setFlagToFMatch(flagToFMatch);
  newTrack->setTof1GapNum(kDefaultValue);
  newTrack->setTof1PathLocation(kDefaultValue);
  newTrack->setTof1PathTime(kDefaultValue);
  newTrack->setTof1TimeU(kDefaultValue);
  newTrack->setTof1TimeD(kDefaultValue);
  newTrack->setTof1TimeMean(kDefaultValue);
  newTrack->setTof1SX(kDefaultValue);
  newTrack->setTof1SY(kDefaultValue);
  newTrack->setTof1SZ(kDefaultValue);
  newTrack->setTof1SNx(kDefaultValue);
  newTrack->setTof1SNy(kDefaultValue);
  newTrack->setTof1SNz(kDefaultValue);
  newTrack->setTof1SP(kDefaultValue);

  //SFT information
  newTrack->setSftSX(kDefaultValue);
  newTrack->setSftSY(kDefaultValue);
  newTrack->setSftSZ(kDefaultValue);
  newTrack->setSftSNx(kDefaultValue);
  newTrack->setSftSNy(kDefaultValue);
  newTrack->setSftSNz(kDefaultValue);
  newTrack->setSftSP(kDefaultValue);

  //Vert information
  theta=(*predictStateVert)(kIdxTheta,0);
  phi=(*predictStateVert)(kIdxPhi,0);
  nx=sin(theta)*cos(phi);
  ny=sin(theta)*sin(phi);
  nz=cos(theta);

  newTrack->setVertMX(measVert.X());
  newTrack->setVertMY(measVert.Y());
  newTrack->setVertMPhi(measVert.Z());
  newTrack->setVertMXErr(measErrVert.X());
  newTrack->setVertMYErr(measErrVert.Y());
  newTrack->setVertMPhiErr(measErrVert.Z());
  newTrack->setVertSX((*predictStateVert)(kIdxX,0));
  newTrack->setVertSY((*predictStateVert)(kIdxY,0));
  newTrack->setVertSZ((*predictStateVert)(kIdxZ,0));
  newTrack->setVertSNx(nx);
  newTrack->setVertSNy(ny);
  newTrack->setVertSNz(nz);
  newTrack->setVertSP(kDefaultValue);
  newTrack->setVertSPhi(phiTrack);

  //Calculate energy loss from AC to Target, and update information at detector planes, such as ToF1, SFT and Vertex
  // Normalize momentum direction at B lower edge
  Double_t thetaACtoTarget=(*predictStateBLowEdge)(kIdxTheta,0);
  Double_t phiACtoTarget=(*predictStateBLowEdge)(kIdxPhi,0);
  Double_t nxACtoTarget=sin(thetaACtoTarget)*cos(phiACtoTarget);
  Double_t nyACtoTarget=sin(thetaACtoTarget)*sin(phiACtoTarget);
  Int_t hitFlagACToTarget;
  Double_t xLayerAcToTarget;
  Bool_t flag;
  Double_t x0=(*predictStateBLowEdge)(kIdxX,0);
  Double_t y0=(*predictStateBLowEdge)(kIdxY,0);
  trackState *stateVertProp;
  trackState *stateVertPropPionPlus;
  trackState *stateVertPropPositron;

  //AC housing top
  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fAc->getRadiusTopHousingOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=13;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    return;
  }
  else{
    planeHit *hitAcTopHousingOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp=predictStateBLowEdge->predictSVatFirstPlaneStraightLine(hitAcTopHousingOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    stateVertPropPionPlus=predictStateBLowEdge->predictSVatFirstPlaneStraightLine(hitAcTopHousingOuter,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
    stateVertPropPositron=predictStateBLowEdge->predictSVatFirstPlaneStraightLine(hitAcTopHousingOuter,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
    delete hitAcTopHousingOuter;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fAc->getRadiusTopHousingInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=12;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitAcTopHousingInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcTopHousingInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitAcTopHousingInner,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAluminium);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitAcTopHousingInner,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAluminium);
    delete hitAcTopHousingInner;
  }

  //AC
  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fAc->getRadiusOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=11;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitAcOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitAcOuter,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitAcOuter,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
    delete hitAcOuter;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fAc->getRadiusInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=10;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitAcInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcInner,fParticleId, fParticleMass, fParticleCharge, kAerogel);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitAcInner,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAerogel);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitAcInner,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAerogel);
    delete hitAcInner;
  }

  //AC housing bottom
  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fAc->getRadiusBottomHousingOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=9;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitAcBottomHousingOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcBottomHousingOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitAcBottomHousingOuter,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitAcBottomHousingOuter,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
    delete hitAcBottomHousingOuter;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fAc->getRadiusBottomHousingInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=8;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitAcBottomHousingInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcBottomHousingInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitAcBottomHousingInner,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAluminium);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitAcBottomHousingInner,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAluminium);
    delete hitAcBottomHousingInner;
  }

  //ToF1
  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fTof1->getRadiusOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=7;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitToF1Outer=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Outer,fParticleId, fParticleMass, fParticleCharge, kAir);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitToF1Outer,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitToF1Outer,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
    delete hitToF1Outer;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fTof1->getRadiusInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=6;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, (fTof1->getRadiusOuter()+fTof1->getRadiusInner())/2, xLayerAcToTarget);
    planeHit *hitToF1Middle=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Middle,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitToF1Middle,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kScintillator);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitToF1Middle,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kScintillator);
    delete hitToF1Middle;

    if( fTof1->getGapNum() >= 0 && fTof1->getGapNum() <= 11 && ( (thisSystem->getGapNumTof2() == fTof1->getGapNum() - 1) || ( thisSystem->getGapNumTof2() == fTof1->getGapNum() ) || (thisSystem->getGapNumTof2() == fTof1->getGapNum() + 1) || (thisSystem->getGapNumTof2() == 0 && fTof1->getGapNum() == 11) || (thisSystem->getGapNumTof2() == 11 && fTof1->getGapNum() == 0) ) ){

      flagToFMatch=kTRUE;
      theta=(*stateVertProp)(kIdxTheta,0);
      phi=(*stateVertProp)(kIdxPhi,0);
      nx=sin(theta)*cos(phi);
      ny=sin(theta)*sin(phi);
      nz=cos(theta);

      newTrack->setFlagToFMatch(flagToFMatch);
      newTrack->setTof1GapNum(fTof1->getGapNum()); 
      newTrack->setTof1PathLocation(stateVertProp->getPathLocation());
      newTrack->setTof1PathTime(stateVertProp->getPathTime());
      newTrack->setTof1TimeU(fTof1->getTimeU());
      newTrack->setTof1TimeD(fTof1->getTimeD());
      newTrack->setTof1TimeMean(fTof1->getTimeMean());
      newTrack->setTof1SX((*stateVertProp)(kIdxX,0));
      newTrack->setTof1SY((*stateVertProp)(kIdxY,0));
      newTrack->setTof1SZ((*stateVertProp)(kIdxZ,0));
      newTrack->setTof1SNx(nx);
      newTrack->setTof1SNy(ny);
      newTrack->setTof1SNz(nz);
      newTrack->setTof1SP(fParticleCharge/(*stateVertProp)(kIdxQP,0));

      newTrack->setTof1PathTimePionPlus(stateVertPropPionPlus->getPathTime());
      newTrack->setTof1SPPionPlus(kParticleCharge[kPionPlusId]/(*stateVertPropPionPlus)(kIdxQP,0));

      newTrack->setTof1PathTimePositron(stateVertPropPositron->getPathTime());
      newTrack->setTof1SPPositron(kParticleCharge[kPositronId]/(*stateVertPropPositron)(kIdxQP,0));
    }

    xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fTof1->getRadiusInner(), xLayerAcToTarget);
    planeHit *hitToF1Inner=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Inner,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitToF1Inner,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kScintillator);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitToF1Inner,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kScintillator);
    delete hitToF1Inner;
  }

  //SFT
  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fSft->getRadiusOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=5;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitSFTOuter=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitSFTOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitSFTOuter,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitSFTOuter,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
    delete hitSFTOuter;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fSft->getRadiusInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=4;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fSft->getRadiusMiddle(), xLayerAcToTarget);
    planeHit *hitSFTMiddle=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitSFTMiddle,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitSFTMiddle,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kScintillator);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitSFTMiddle,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kScintillator);
    delete hitSFTMiddle;

    theta=(*stateVertProp)(kIdxTheta,0);
    phi=(*stateVertProp)(kIdxPhi,0);
    nx=sin(theta)*cos(phi);
    ny=sin(theta)*sin(phi);
    nz=cos(theta);

    newTrack->setSftSX((*stateVertProp)(kIdxX,0));
    newTrack->setSftSY((*stateVertProp)(kIdxY,0));
    newTrack->setSftSZ((*stateVertProp)(kIdxZ,0));
    newTrack->setSftSNx(nx);
    newTrack->setSftSNy(ny);
    newTrack->setSftSNz(nz);
    newTrack->setSftSP(fParticleCharge/(*stateVertProp)(kIdxQP,0));

    newTrack->setSftSPPionPlus(kParticleCharge[kPionPlusId]/(*stateVertPropPionPlus)(kIdxQP,0));
    newTrack->setSftSPPositron(kParticleCharge[kPositronId]/(*stateVertPropPositron)(kIdxQP,0));

    xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fSft->getRadiusInner(), xLayerAcToTarget);
    planeHit *hitSFTInner=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitSFTInner,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitSFTInner,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kScintillator);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitSFTInner,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kScintillator);
    delete hitSFTInner;
  }

  //Target
  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fTarget->getRadiusHolderOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=3;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitTargetHolderOuter=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitTargetHolderOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitTargetHolderOuter,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitTargetHolderOuter,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
    delete hitTargetHolderOuter;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fTarget->getRadiusHolderInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=2;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    planeHit *hitTargetHolderInner=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitTargetHolderInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitTargetHolderInner,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAluminium);
    stateVertPropPositron->predictSVatNextPlaneStraightLine(hitTargetHolderInner,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAluminium);
    delete hitTargetHolderInner;
  }

  flag=xAtCylinder(x0, y0, nxACtoTarget, nyACtoTarget, fTarget->getRadius(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=1;
    newTrack->setFlagTargetToAC(hitFlagACToTarget);
    delete stateVertProp;
    delete stateVertPropPionPlus;
    delete stateVertPropPositron;
    return;
  }
  else{
    if(sqrt((*predictStateVert)(kIdxX,0)*(*predictStateVert)(kIdxX,0)+(*predictStateVert)(kIdxY,0)*(*predictStateVert)(kIdxY,0))>fTarget->getRadius()){//vertex from tracking is out of the target radius
      hitFlagACToTarget=-1;
      delete stateVertProp;
      delete stateVertPropPionPlus;
      delete stateVertPropPositron;
      return;
    }
    else{
      hitFlagACToTarget=0;
      newTrack->setFlagTargetToAC(hitFlagACToTarget);
      planeHit *hitTarget=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
      stateVertProp->predictSVatNextPlaneStraightLine(hitTarget,fParticleId, fParticleMass, fParticleCharge, kAir);
      stateVertPropPionPlus->predictSVatNextPlaneStraightLine(hitTarget,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kAir);
      stateVertPropPositron->predictSVatNextPlaneStraightLine(hitTarget,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kAir);
      delete hitTarget;

      planeHit *hitVert=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(), (*predictStateVert)(kIdxX,0));
      stateVertProp->predictSVatLastPlaneStraightLine(hitVert,fParticleId, fParticleMass, fParticleCharge, kScintillator);
      stateVertPropPionPlus->predictSVatLastPlaneStraightLine(hitVert,kPionPlusId, kParticleMass[kPionPlusId], kParticleCharge[kPionPlusId], kScintillator);
      stateVertPropPositron->predictSVatLastPlaneStraightLine(hitVert,kPositronId, kParticleMass[kPositronId], kParticleCharge[kPositronId], kScintillator);
      delete hitVert;

      newTrack->setVertPathTime(stateVertProp->getPathTime());
      newTrack->setVertSP(fParticleCharge/(*stateVertProp)(kIdxQP,0));

      newTrack->setVertPathTimePionPlus(stateVertPropPionPlus->getPathTime());
      newTrack->setVertSPPionPlus(kParticleCharge[kPionPlusId]/(*stateVertPropPionPlus)(kIdxQP,0));

      newTrack->setVertPathTimePositron(stateVertPropPositron->getPathTime());
      newTrack->setVertSPPositron(kParticleCharge[kPositronId]/(*stateVertPropPositron)(kIdxQP,0));

      delete stateVertProp;
      delete stateVertPropPionPlus;
      delete stateVertPropPositron;
      return;
    }
  }
}

//______________________________________________________________________________
Bool_t trackFinder::xAtCylinder(Double_t x0, Double_t y0, Double_t nx, Double_t ny, Double_t r, Double_t &x){
  Double_t a=nx*nx+ny*ny;
  Double_t b=2*(nx*x0+ny*y0);
  Double_t c=x0*x0+y0*y0-r*r;
  Double_t delta=b*b-4*a*c;
  Double_t t1,t2;
  Double_t x1,x2;
  if(delta<0) return kFALSE;
  else{
    t1=(-b+sqrt(delta))/2/a;
    t2=(-b-sqrt(delta))/2/a;
    x1=nx*t1+x0;
    x2=nx*t2+x0;
    if(x1>x2) x=x1;
    else x=x2;
    return kTRUE;
  }
}

//______________________________________________________________________________
Bool_t trackFinder::calcVertChi2(Int_t gapNumToF2, trackState* stateBLowEdge, Double_t &chi2){
  Double_t theta=(*stateBLowEdge)(kIdxTheta,0);
  Double_t phi=(*stateBLowEdge)(kIdxPhi,0);
  Double_t nx=sin(theta)*cos(phi);
  Double_t ny=sin(theta)*sin(phi);
  Double_t x0=(*stateBLowEdge)(kIdxX,0);
  Double_t y0=(*stateBLowEdge)(kIdxY,0);

  Double_t x1prime=fTarget->getVertX();
  Double_t y1prime=fTarget->getVertY();
  Double_t phi1prime=fTarget->getVertPhi();
  Double_t x1primeErr=fTarget->getVertXErr();
  Double_t y1primeErr=fTarget->getVertYErr();
  Double_t phi1primeErr=fTarget->getVertPhiErr();

  Double_t alpha=30*(gapNumToF2-2)/180.*TMath::Pi();

  Double_t x1=x1prime*cos(alpha)-y1prime*sin(alpha);
  Double_t y1=x1prime*sin(alpha)+y1prime*cos(alpha);
  Double_t phi1=phi1prime+alpha/TMath::Pi()*180;
  if(phi1<0) phi1=phi1+360;
  else if(phi1>=360) phi1=phi1-360;
  Double_t x1Err=sqrt(x1primeErr*x1primeErr*cos(alpha)*cos(alpha)+y1primeErr*y1primeErr*sin(alpha)*sin(alpha)+4.4*4.4);
  Double_t y1Err=sqrt(x1primeErr*x1primeErr*sin(alpha)*sin(alpha)+y1primeErr*y1primeErr*cos(alpha)*cos(alpha)+4.4*4.4);
  Double_t phi1Err=sqrt(phi1primeErr*phi1primeErr+1.5*1.5);

  Double_t deltaX=x1-x0;
  Double_t deltaY=y1-y0;

  Double_t varX=0, varY=0, varPhi=0;

  if(x1Err==0 || y1Err==0 || (nx==0 && ny ==0)) return kFALSE;
  else{
    Double_t t=(deltaX*nx/x1Err/x1Err+deltaY*ny/y1Err/y1Err)/(nx*nx/x1Err/x1Err+ny*ny/y1Err/y1Err);
    Double_t phi0=0;
    if(nx==0 && ny>0) phi0=90;
    else if(nx==0 && ny<0) phi0=270;
    else if(nx>0 && ny>=0) phi0=atan(ny/nx)/TMath::Pi()*180;
    else if(nx>0 && ny<0) phi0=atan(ny/nx)/TMath::Pi()*180+360;
    else if(nx<0) phi0=atan(ny/nx)/TMath::Pi()*180+180;

    Double_t absDeltaPhi;
    if( (phi0>=0 && phi0<90 && phi1>270 && phi1<360) || (phi1>=0 && phi1<90 && phi0>270 && phi0<360) ) absDeltaPhi=fabs(phi1+phi0-360);
    else absDeltaPhi=fabs(phi1-phi0);

    chi2=(nx*t-deltaX)*(nx*t-deltaX)/x1Err/x1Err+(ny*t-deltaY)*(ny*t-deltaY)/y1Err/y1Err+absDeltaPhi*absDeltaPhi/phi1Err/phi1Err;

    if(chi2>kChi2PerMdimCut) return kFALSE;
    else{
      Double_t xCalc=nx*t+x0;
      Double_t yCalc=ny*t+y0;
      //Double_t rCalc=sqrt(xCalc*xCalc+yCalc*yCalc);

      Int_t hitFlagACToTarget;
      Double_t xLayerAcToTarget;
      Bool_t flag;
      trackState *stateVertProp;
      Bool_t flagToFMatch=kFALSE;

      planeHit *hitVertCalc=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xCalc);

      //AC housing top
      flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusTopHousingOuter(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=13;
	stateVertProp=stateBLowEdge->predictSVatFirstPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	planeHit *hitAcTopHousingOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
	stateVertProp=stateBLowEdge->predictSVatFirstPlaneStraightLine(hitAcTopHousingOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitAcTopHousingOuter;
      }

      flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusTopHousingInner(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=12;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	goto label;
      }
      else{
	planeHit *hitAcTopHousingInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitAcTopHousingInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	delete hitAcTopHousingInner;
      }

      //AC
      flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusOuter(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=11;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	planeHit *hitAcOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitAcOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitAcOuter;
      }

      flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusInner(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=10;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAerogel);
	goto label;
      }
      else{
	planeHit *hitAcInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitAcInner,fParticleId, fParticleMass, fParticleCharge, kAerogel);
	delete hitAcInner;
      }

      //AC housing bottom
      flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusBottomHousingOuter(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=9;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	planeHit *hitAcBottomHousingOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitAcBottomHousingOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitAcBottomHousingOuter;
      }
      
      flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusBottomHousingInner(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=8;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	goto label;
      }
      else{
	planeHit *hitAcBottomHousingInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitAcBottomHousingInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	delete hitAcBottomHousingInner;
      }

      //ToF1
      flag=xAtCylinder(x0, y0, nx, ny, fTof1->getRadiusOuter(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=7;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	planeHit *hitToF1Outer=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Outer,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitToF1Outer;
      }

      flag=xAtCylinder(x0, y0, nx, ny, fTof1->getRadiusInner(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=6;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	goto label;
      }
      else{
	xAtCylinder(x0, y0, nx, ny, (fTof1->getRadiusOuter()+fTof1->getRadiusInner())/2, xLayerAcToTarget);
	planeHit *hitToF1Middle=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Middle,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	delete hitToF1Middle;

	if( fTof1->getGapNum() >= 0 && fTof1->getGapNum() <= 11 && ( (gapNumToF2 == fTof1->getGapNum() - 1) || ( gapNumToF2 == fTof1->getGapNum() ) || (gapNumToF2 == fTof1->getGapNum() + 1) || (gapNumToF2 == 0 && fTof1->getGapNum() == 11) || (gapNumToF2 == 11 && fTof1->getGapNum() == 0) ) ){
	  flagToFMatch=kTRUE;	  
	}

	xAtCylinder(x0, y0, nx, ny, fTof1->getRadiusInner(), xLayerAcToTarget);
	planeHit *hitToF1Inner=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Inner,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	delete hitToF1Inner;
      }
      
      //SFT
      flag=xAtCylinder(x0, y0, nx, ny, fSft->getRadiusOuter(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=5;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	planeHit *hitSFTOuter=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitSFTOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitSFTOuter;
      }

      flag=xAtCylinder(x0, y0, nx, ny, fSft->getRadiusInner(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=4;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	goto label;
      }
      else{
	xAtCylinder(x0, y0, nx, ny, (fSft->getRadiusOuter() + fSft->getRadiusInner())/2, xLayerAcToTarget);
	planeHit *hitSFTMiddle=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitSFTMiddle,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	delete hitSFTMiddle;

	xAtCylinder(x0, y0, nx, ny, fSft->getRadiusInner(), xLayerAcToTarget);
	planeHit *hitSFTInner=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitSFTInner,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	delete hitSFTInner;
      }

      //Target
      flag=xAtCylinder(x0, y0, nx, ny, fTarget->getRadiusHolderOuter(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=3;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	planeHit *hitTargetHolderOuter=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitTargetHolderOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitTargetHolderOuter;
      }

      flag=xAtCylinder(x0, y0, nx, ny, fTarget->getRadiusHolderInner(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=2;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	goto label;
      }
      else{
	planeHit *hitTargetHolderInner=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitTargetHolderInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	delete hitTargetHolderInner;
      }

      flag=xAtCylinder(x0, y0, nx, ny, fTarget->getRadius(), xLayerAcToTarget);
      if(flag==kFALSE){
	hitFlagACToTarget=1;
	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
	goto label;
      }
      else{
	hitFlagACToTarget=0;
	planeHit *hitTarget=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
	stateVertProp->predictSVatNextPlaneStraightLine(hitTarget,fParticleId, fParticleMass, fParticleCharge, kAir);
	delete hitTarget;

	stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kScintillator);
	goto label;
      }
    label:
      {
	stateVertProp->predictSVatLastPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
	varX=(stateVertProp->getCovMat())(kIdxX,kIdxX);
	varY=(stateVertProp->getCovMat())(kIdxY,kIdxY);
	varPhi=(stateVertProp->getCovMat())(kIdxPhi,kIdxPhi);
	delete hitVertCalc;
	delete stateVertProp;

	x1Err=sqrt(x1primeErr*x1primeErr*cos(alpha)*cos(alpha)+y1primeErr*y1primeErr*sin(alpha)*sin(alpha));
	y1Err=sqrt(x1primeErr*x1primeErr*sin(alpha)*sin(alpha)+y1primeErr*y1primeErr*cos(alpha)*cos(alpha));
	phi1Err=phi1primeErr;

	x1Err=sqrt(x1Err*x1Err+varX);
	y1Err=sqrt(y1Err*x1Err+varY);
	phi1Err=sqrt(phi1Err*phi1Err+varPhi);

	t=(deltaX*nx/x1Err/x1Err+deltaY*ny/y1Err/y1Err)/(nx*nx/x1Err/x1Err+ny*ny/y1Err/y1Err);

	chi2=(nx*t-deltaX)*(nx*t-deltaX)/x1Err/x1Err+(ny*t-deltaY)*(ny*t-deltaY)/y1Err/y1Err+absDeltaPhi*absDeltaPhi/phi1Err/phi1Err;

	if(chi2>kChi2PerMdimCut) return kFALSE;
	else return kTRUE;
      }
    }
  }
}
//______________________________________________________________________________
Double_t trackFinder::calcVertChi2(Int_t gapNumToF2, trackState* stateBLowEdge, trackState* stateVert, TVector3 &measVert, TVector3 &measErrVert, Double_t &phiTrack){
  Double_t theta=(*stateBLowEdge)(kIdxTheta,0);
  Double_t phi=(*stateBLowEdge)(kIdxPhi,0);
  Double_t nx=sin(theta)*cos(phi);
  Double_t ny=sin(theta)*sin(phi);
  Double_t nz=cos(theta);
  Double_t x0=(*stateBLowEdge)(kIdxX,0);
  Double_t y0=(*stateBLowEdge)(kIdxY,0);
  Double_t z0=(*stateBLowEdge)(kIdxZ,0);

  Double_t x1prime=fTarget->getVertX();
  Double_t y1prime=fTarget->getVertY();
  Double_t phi1prime=fTarget->getVertPhi();
  Double_t x1primeErr=fTarget->getVertXErr();
  Double_t y1primeErr=fTarget->getVertYErr();
  Double_t phi1primeErr=fTarget->getVertPhiErr();

  Double_t alpha=30*(gapNumToF2-2)/180.*TMath::Pi();

  Double_t x1=x1prime*cos(alpha)-y1prime*sin(alpha);
  Double_t y1=x1prime*sin(alpha)+y1prime*cos(alpha);
  Double_t phi1=phi1prime+alpha/TMath::Pi()*180;
  if(phi1<0) phi1=phi1+360;
  else if(phi1>=360) phi1=phi1-360;
  Double_t x1Err=sqrt(x1primeErr*x1primeErr*cos(alpha)*cos(alpha)+y1primeErr*y1primeErr*sin(alpha)*sin(alpha)+4.4*4.4);
  Double_t y1Err=sqrt(x1primeErr*x1primeErr*sin(alpha)*sin(alpha)+y1primeErr*y1primeErr*cos(alpha)*cos(alpha)+4.4*4.4);
  Double_t phi1Err=sqrt(phi1primeErr*phi1primeErr+1.5*1.5);

  Double_t deltaX=x1-x0;
  Double_t deltaY=y1-y0;

  Double_t t=(deltaX*nx/x1Err/x1Err+deltaY*ny/y1Err/y1Err)/(nx*nx/x1Err/x1Err+ny*ny/y1Err/y1Err);

  Double_t varX=0, varY=0, varPhi=0;

  Double_t phi0=0;
  if(nx==0 && ny>0) phi0=90;
  else if(nx==0 && ny<0) phi0=270;
  else if(nx>0 && ny>=0) phi0=atan(ny/nx)/TMath::Pi()*180;
  else if(nx>0 && ny<0) phi0=atan(ny/nx)/TMath::Pi()*180+360;
  else if(nx<0) phi0=atan(ny/nx)/TMath::Pi()*180+180;

  Double_t absDeltaPhi;
  if( (phi0>=0 && phi0<90 && phi1>270 && phi1<360) || (phi1>=0 && phi1<90 && phi0>270 && phi0<360) ) absDeltaPhi=fabs(phi1+phi0-360);
  else absDeltaPhi=fabs(phi1-phi0);

  Double_t xCalc=nx*t+x0;
  Double_t yCalc=ny*t+y0;
  //Double_t rCalc=sqrt(xCalc*xCalc+yCalc*yCalc);

  Int_t hitFlagACToTarget;
  Double_t xLayerAcToTarget;
  Bool_t flag;
  trackState *stateVertProp;
  Bool_t flagToFMatch=kFALSE;

  planeHit *hitVertCalc=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xCalc);

  //AC housing top
  flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusTopHousingOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=13;
    stateVertProp=stateBLowEdge->predictSVatFirstPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    planeHit *hitAcTopHousingOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp=stateBLowEdge->predictSVatFirstPlaneStraightLine(hitAcTopHousingOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitAcTopHousingOuter;
  }

  flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusTopHousingInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=12;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    goto label;
  }
  else{
    planeHit *hitAcTopHousingInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcTopHousingInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    delete hitAcTopHousingInner;
  }

  //AC
  flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=11;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    planeHit *hitAcOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitAcOuter;
  }

  flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=10;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAerogel);
    goto label;
  }
  else{
    planeHit *hitAcInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcInner,fParticleId, fParticleMass, fParticleCharge, kAerogel);
    delete hitAcInner;
  }

  //AC housing bottom
  flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusBottomHousingOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=9;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    planeHit *hitAcBottomHousingOuter=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcBottomHousingOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitAcBottomHousingOuter;
  }
      
  flag=xAtCylinder(x0, y0, nx, ny, fAc->getRadiusBottomHousingInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=8;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    goto label;
  }
  else{
    planeHit *hitAcBottomHousingInner=new planeHit(fAc->getPlaneId(),fAc->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitAcBottomHousingInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    delete hitAcBottomHousingInner;
  }

  //ToF1
  flag=xAtCylinder(x0, y0, nx, ny, fTof1->getRadiusOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=7;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    planeHit *hitToF1Outer=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Outer,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitToF1Outer;
  }

  flag=xAtCylinder(x0, y0, nx, ny, fTof1->getRadiusInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=6;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    goto label;
  }
  else{
    xAtCylinder(x0, y0, nx, ny, (fTof1->getRadiusOuter()+fTof1->getRadiusInner())/2, xLayerAcToTarget);
    planeHit *hitToF1Middle=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Middle,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    delete hitToF1Middle;

    if( fTof1->getGapNum() >= 0 && fTof1->getGapNum() <= 11 && ( (gapNumToF2 == fTof1->getGapNum() - 1) || ( gapNumToF2 == fTof1->getGapNum() ) || (gapNumToF2 == fTof1->getGapNum() + 1) || (gapNumToF2 == 0 && fTof1->getGapNum() == 11) || (gapNumToF2 == 11 && fTof1->getGapNum() == 0) ) ){
      flagToFMatch=kTRUE;	  
    }

    xAtCylinder(x0, y0, nx, ny, fTof1->getRadiusInner(), xLayerAcToTarget);
    planeHit *hitToF1Inner=new planeHit(fTof1->getPlaneId(),fTof1->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitToF1Inner,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    delete hitToF1Inner;
  }
      
  //SFT
  flag=xAtCylinder(x0, y0, nx, ny, fSft->getRadiusOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=5;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    planeHit *hitSFTOuter=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitSFTOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitSFTOuter;
  }

  flag=xAtCylinder(x0, y0, nx, ny, fSft->getRadiusInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=4;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    goto label;
  }
  else{
    xAtCylinder(x0, y0, nx, ny, (fSft->getRadiusOuter() + fSft->getRadiusInner())/2, xLayerAcToTarget);
    planeHit *hitSFTMiddle=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitSFTMiddle,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    delete hitSFTMiddle;

    xAtCylinder(x0, y0, nx, ny, fSft->getRadiusInner(), xLayerAcToTarget);
    planeHit *hitSFTInner=new planeHit(fSft->getPlaneId(),fSft->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitSFTInner,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    delete hitSFTInner;
  }

  //Target
  flag=xAtCylinder(x0, y0, nx, ny, fTarget->getRadiusHolderOuter(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=3;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    planeHit *hitTargetHolderOuter=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitTargetHolderOuter,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitTargetHolderOuter;
  }

  flag=xAtCylinder(x0, y0, nx, ny, fTarget->getRadiusHolderInner(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=2;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    goto label;
  }
  else{
    planeHit *hitTargetHolderInner=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitTargetHolderInner,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    delete hitTargetHolderInner;
  }

  flag=xAtCylinder(x0, y0, nx, ny, fTarget->getRadius(), xLayerAcToTarget);
  if(flag==kFALSE){
    hitFlagACToTarget=1;
    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAir);
    goto label;
  }
  else{
    hitFlagACToTarget=0;
    planeHit *hitTarget=new planeHit(fTarget->getPlaneId(),fTarget->getPlaneDir(),xLayerAcToTarget);
    stateVertProp->predictSVatNextPlaneStraightLine(hitTarget,fParticleId, fParticleMass, fParticleCharge, kAir);
    delete hitTarget;

    stateVertProp->predictSVatNextPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kScintillator);
    goto label;
  }
 label:
  {
    stateVertProp->predictSVatLastPlaneStraightLine(hitVertCalc,fParticleId, fParticleMass, fParticleCharge, kAluminium);
    varX=(stateVertProp->getCovMat())(kIdxX,kIdxX);
    varY=(stateVertProp->getCovMat())(kIdxY,kIdxY);
    varPhi=(stateVertProp->getCovMat())(kIdxPhi,kIdxPhi);
    delete hitVertCalc;
    delete stateVertProp;

    x1Err=sqrt(x1primeErr*x1primeErr*cos(alpha)*cos(alpha)+y1primeErr*y1primeErr*sin(alpha)*sin(alpha));
    y1Err=sqrt(x1primeErr*x1primeErr*sin(alpha)*sin(alpha)+y1primeErr*y1primeErr*cos(alpha)*cos(alpha));
    phi1Err=phi1primeErr;

    measVert.SetXYZ(x1,y1,phi1);
    measErrVert.SetXYZ(x1Err,y1Err,phi1Err);
    phiTrack=phi0;

    x1Err=sqrt(x1Err*x1Err+varX);
    y1Err=sqrt(y1Err*x1Err+varY);
    phi1Err=sqrt(phi1Err*phi1Err+varPhi);

    t=(deltaX*nx/x1Err/x1Err+deltaY*ny/y1Err/y1Err)/(nx*nx/x1Err/x1Err+ny*ny/y1Err/y1Err);

    Double_t chi2=(nx*t-deltaX)*(nx*t-deltaX)/x1Err/x1Err+(ny*t-deltaY)*(ny*t-deltaY)/y1Err/y1Err+absDeltaPhi*absDeltaPhi/phi1Err/phi1Err;

    (*stateVert)(kIdxX,0)=nx*t+x0;
    (*stateVert)(kIdxY,0)=ny*t+y0;
    (*stateVert)(kIdxZ,0)=nz*t+z0;
    (*stateVert)(kIdxTheta,0)=(*stateBLowEdge)(kIdxTheta,0);
    (*stateVert)(kIdxPhi,0)=(*stateBLowEdge)(kIdxPhi,0);
    (*stateVert)(kIdxQP,0)=(*stateBLowEdge)(kIdxQP,0);
    return chi2;
  }
}

ClassImp(trackFinder)
