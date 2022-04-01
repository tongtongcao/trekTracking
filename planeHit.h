#ifndef ROOT_PLANEHit
#define ROOT_PLANEHit

#include "TObject.h"
#include "TVector3.h"

#include "utility.h"

class planeHit:public TObject{
public:
 planeHit(){}
 planeHit(planeId id, planeDir dir, Double_t xOrZ);
 planeHit(planeId id, planeDir dir, int hitLocalXId, int hitLocalYId, Double_t x, Double_t y, Double_t z, Double_t xErr=1, Double_t yErr=1, Double_t zErr=1):fPlaneId(id), fPlaneDir(dir), fHitLocalXId(hitLocalXId), fHitLocalYId(hitLocalYId), fX(x), fY(y), fZ(z), fXErr(xErr), fYErr(yErr), fZErr(zErr){}
  planeId getPlaneId() const {return fPlaneId;}
  planeDir getPlaneDir() const {return fPlaneDir;}
  int getLocalXHitId () const {return fHitLocalXId;}
  int getLocalYHitId () const {return fHitLocalYId;}
  TVector3 getPosition () const;
  TVector3 getPositionError () const;

  void setX (Double_t x){ fX = x; }
  void setY (Double_t y){ fY = y; }
  void setZ (Double_t z){ fZ = z; }

  void setXError (Double_t xErr){ fXErr = xErr; }
  void setYError (Double_t yErr){ fYErr = yErr; }
  void setZError (Double_t zErr){ fZErr = zErr; }

  Double_t getX () const{ return fX; }
  Double_t getY () const{ return fY; }
  Double_t getZ () const{ return fZ; }

  Double_t getXError () const{ return fXErr; }
  Double_t getYError () const{ return fYErr; }
  Double_t getZError () const{ return fZErr; }


  virtual ~planeHit() {}

private:
  planeId fPlaneId; 
  planeDir fPlaneDir; 
  int fHitLocalXId; // hit # in local X direction; start from 0
  int fHitLocalYId; // hit # in local Y direction; start from 0
  Double_t fX;
  Double_t fY;
  Double_t fZ;
  Double_t fXErr;
  Double_t fYErr;
  Double_t fZErr;

  ClassDef(planeHit,1)    

};

#endif
