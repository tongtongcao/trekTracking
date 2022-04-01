#ifndef ROOT_TRACK_SYSTEM
#define ROOT_TRACK_SYSTEM

#include "TObjArray.h"
#include "TMath.h"

#include "utility.h"
#include "trackSite.h"
#include "trackState.h"

//class matrix;
class trackState;

class trackSystem: public TObjArray{
  friend class trackSite;
  
  public:
  trackSystem(Int_t n = 1);
  ~trackSystem();
  void Add(TObject *obj);
  void smoothBackTo(Int_t k);


  inline void      setSeedType(seedType type) { fSeedType = type; }
  inline seedType  getSeedType() const { return fSeedType; }

  void checkTrackStatus();
  inline Bool_t    getTrackStatus() const   { return fIsGood; }
  inline void      setTrackStatus(Bool_t isGood) { fIsGood = isGood; }
  
  static  trackSystem *getCurInstancePtr() { return fgCurInstancePtr; }
  void setCurInstancePtr(trackSystem *ksp) { fgCurInstancePtr = ksp; }

  inline trackSite* getCurSite() { return fCurSitePtr;}
  void increaseChi2(Double_t deltaChi2) { fChi2 += deltaChi2; }

  inline void addMissingHits() { fNMissingHits++; } 
  inline void setMissingHits(Int_t n) { fNMissingHits = n;  }
  inline Double_t  getMissingHits() const   { return fNMissingHits; }

  inline void setNHits(Int_t n) { fNHits = n;  }
  inline Int_t getNHits() const   { return fNHits; }

  inline void setChi2(Double_t chi2) { fChi2 = chi2;  }
  inline Double_t getChi2() const   { return fChi2; }
  inline void setNDF(Int_t ndf ) { fNDF = ndf;  }
  inline Double_t getNDF() const   { return fNDF; }

  inline Double_t getChi2PerNDF() const {return fChi2/fNDF;}

  inline void setGapNumTof2(Int_t gapNumTof2 ) { fGapNumTof2 = gapNumTof2;  }
  inline Double_t getGapNumTof2() const   { return fGapNumTof2; }

 inline void setTrackId(trackId id ) { fTrackId = id;  }
  inline trackId getTrackId() const   { return fTrackId; }

  virtual Int_t Compare( const TObject* obj ) const;
  virtual Bool_t IsSortable () const { return kTRUE; }
       
 private:
  static trackSystem *fgCurInstancePtr;  //! currently active instance

  seedType fSeedType;
  Int_t fNMissingHits;
  Int_t fNHits;

  Double_t fChi2;  // current total chi2
  Int_t fNDF;

  Double_t     fMomentum;

  trackSite* fCurSitePtr;  // pointer to current site

  Bool_t       fIsGood;      // for track fitting monitering 

  Int_t fGapNumTof2;

  trackId fTrackId;



  ClassDef(trackSystem,1)  // Base class for Kalman Filter

};
#endif
