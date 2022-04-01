#ifndef ROOT_TRACK_FINDER
#define ROOT_TRACK_FINDER
//ROOT
#include <vector>
#include <map>
using namespace std;

#include "TVector3.h"
#include "TClonesArray.h"
#include "TObject.h"

#include "utility.h"
#include "particle.h"
#include "planeHit.h"
#include "mwpc.h"
#include "target.h"
#include "sft.h"
#include "tof1.h"
#include "ac.h"
#include "tof2.h"
#include "trackState.h"
#include "trackSite.h"
#include "track.h"
#include "trackSystem.h"
#include "TVector3.h"

class planeHit;
class trackState;
class trackSite;
class trackSystem;

class trackFinder: public particle{
 public:
  trackFinder(particleId id, Double_t mass, Double_t charge, target *uTarget, sft *uSft, tof1* uTof1, ac* uAc, mwpc* mwpc2, mwpc* mwpc3, mwpc* mwpc4, tof2 *uTof2);
  ~trackFinder();
  void Clear( Option_t* opt="");



  Bool_t processHits(TClonesArray* theTracks);
  Bool_t findSeed();
  Bool_t trackFilter();
  void trackSmoothBack();
  Bool_t finalSelection(TClonesArray* theTracks);


  struct seed{
    Int_t fGapNumTof2;
    seedType type;
    vector<planeHit> hits; //Last two planes must be MWPC3 and MWPC4
    Bool_t isActive;
    TVector3 initPos;
    Double_t initTheta;
    Double_t initPhi;
    Double_t initMom;
    Double_t charge;
    Double_t initCov[5]; // covariance between x, y, z, theta and phi
  
    seed() {}
  seed(Int_t gapNumTof2, seedType& t, vector<planeHit> h, TVector3 pos, Double_t theta, Double_t phi, Double_t& mom, Double_t& q, Double_t cov[5]): fGapNumTof2(gapNumTof2), type(t), hits(h), isActive(kTRUE), initPos(pos), initTheta(theta), initPhi(phi), initMom(mom), charge(q){
      for(Int_t i=0;i<5;i++){
	initCov[i]=cov[i];
      }
    }

    ~seed(){hits.clear();}
    void deactive() { isActive = kFALSE; }
    Int_t getGapNumTof2() const{return fGapNumTof2;}
  };


  Double_t fastTracker(vector<planeHit> hits);
  Bool_t determineInitMom(Int_t gapNumToF2, seedType type, vector<planeHit> hits, TVector3 initPos, Double_t initTheta, Double_t initPhi, Double_t* initCov, Double_t &mom); // Determine initial momentum for each seed. Seeds are eliminated when momentum or Chi2 over limits 
  Bool_t seedTrackTry(seed* gSeed, Double_t &chi2PerNDF); // Outout Chi2/NDF after track passes through filter with a given initial momentum
  trackSite* siteInit(seed* thisSeed); //Set inital state and covariance for each seed
  void calcInitCoV(planeHit* hitMwpc3, planeHit* hitMwpc4, Double_t initCov[]); // Calculate initial covariance

  Bool_t chooseBestHit(trackState* currentState, TClonesArray* hits, Int_t &bestHitIndex); 

  void copyTrack(track* newTrack, trackSystem* thisSystem, trackState *predictStateBLowEdge, Double_t chi2Vert, trackState* predictStateVert, TVector3 measVert, TVector3 measErrVert, Double_t phiTrack);

  Bool_t calcVertChi2(Int_t gapNumToF2, trackState* stateBLowEdge, Double_t &chi2);
  Double_t calcVertChi2(Int_t gapNumToF2, trackState* stateBLowEdge, trackState* stateVert, TVector3 &measVert, TVector3 &measErrVert, Double_t &phiTrack);

  Bool_t xAtCylinder(Double_t x0, Double_t y0, Double_t nx, Double_t ny, Double_t r, Double_t &x);

  TClonesArray* getCoarseTracks() const {return fCoarseTracks;}

 private:
  target* fTarget;
  sft* fSft;
  tof1* fTof1;
  ac* fAc;
  mwpc* fMwpc2; 
  mwpc* fMwpc3; 
  mwpc* fMwpc4; 
  tof2* fTof2;

  TClonesArray* fCoarseTracks;
  int fNSeeds;
  Double_t fChi2PerMdimCut;
  map<seedType, vector<seed> > fSeedPool;
  map<Int_t, vector<planeHit> > fGoodHits;
  Int_t fNGoodTracks;

  ClassDef(trackFinder,1) 

};





#endif
