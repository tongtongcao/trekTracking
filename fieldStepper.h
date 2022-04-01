#ifndef ROOT_FIELD_STEPPER
#define ROOT_FIELD_STEPPER

#include "TVector3.h"

#include "utility.h"
#include "fieldMap.h"
#include "matrix.h"
#include "planeHit.h"
#include "trackSite.h"
#include "trackState.h"
class trackSite;
class trackState;

class fieldStepper:public TObject{
 public:
  fieldStepper();
  ~fieldStepper();

  static fieldStepper* getInstance() {
    if (fFieldStepper == NULL) fFieldStepper = new fieldStepper();
    return fFieldStepper;
  }

  void turnOnMS()       { fIsMSOn = kTRUE;  }
  void turnOffMS()      { fIsMSOn = kFALSE; }
  void turnOnDEDX ()    { fIsDEDXOn = kTRUE;  }
  void turnOffDEDX ()   { fIsDEDXOn = kFALSE; }
  
  inline Bool_t isMSOn()   const { return fIsMSOn;     }
  inline Bool_t isDEDXOn() const { return fIsDEDXOn;   }
  inline Double_t getEnergyLoss()  const { return energyLoss; }
  inline Double_t getStepLength() const { return stepLength; }
  inline Double_t getTrackLength() const { return trackLength; }
  inline Double_t getTrackTime() const { return trackTime; }
  inline Double_t getPathLocation() const { return pathLocation; }
  inline Double_t getPathTime() const { return pathTime; }

  void useDefaultStep();

  void transport(trackState* svFrom, planeHit *hitTo, matrix &sv, matrix &F, matrix &Q, particleId particleId, Double_t particleMass, Double_t particleCharge, Bool_t & magFlag, matType type = kAir); 
  Double_t rKPropagation(matrix &stateVec, matrix &fPropStep, Double_t stepSize, Bool_t bCalcJac, Bool_t dir, Double_t particleMass, Double_t particleCharge);
  Bool_t propagateStraightLine(matrix &stateVec, matrix &fPropChange, const planeHit *hitTo, Double_t particleMass, Double_t particleCharge);
  void propagateStraightLine(matrix &stateVec, matrix &fPropChange, Double_t step);
  Double_t estimateDistanceToNextPlane(const trackState* svFrom, const planeHit *hitTo);
  Double_t estimateDistanceToNextPlane(const matrix svFrom, const planeHit *hitTo);

  void transportStraightLine(trackState *svFrom, planeHit *hitTo, matrix &sv, matrix &F, matrix &Q, particleId particleId, Double_t particleMass, Double_t particleCharge, matType type);

  Double_t calcMultScat(matrix &Q, matrix &svFrom, Double_t length, Double_t beta, matType type = kAir);
  Double_t calcEnergyLoss(matrix &Q, matrix &svTo, Double_t length, Double_t qp, Double_t beta,particleId particleId, Double_t particleMass, Double_t particleCharge, matType type = kAir);
  Double_t calcDEDXIonPositron(Double_t qp, Double_t ZoverA, Double_t density, Double_t I);
  Double_t calcDEDXBetheBloch(Double_t beta, Double_t ZoverA, Double_t density, Double_t I, Double_t particleMass, Double_t particleCharge);
 
 protected:  
  void initDetMaterial();
  fieldMap* fFieldMap;
  
  static fieldStepper* fFieldStepper;

  Bool_t    fIsBackward;     //kTRUE: backward propagation; kFALSE: forward propagation
  Bool_t    fIsMSOn;         //! switch for multiple scattering
  Bool_t    fIsDEDXOn;       //! switch for energy loss

  Double_t  minPrecision; // minimum precision for the Runge-Kutta method
  Double_t  maxPrecision; // maximum precision for the Runge-Kutta method
  Double_t  minStepSize;  //minimum step size for the Runge-Kutta method
  Double_t  maxStepSize;  //maximum step size for the Runge-Kutta method
  Double_t  maxNumSteps;  //maximum number of steps for the Runge-Kutta method
  Double_t  stepSizeDec;  //step size reduction factor for the Runge-Kutta method
  Double_t  stepSizeInc;  //step size increasment factor for the Runge-Kutta method
  Double_t  maxDist;      //maximum distance for straight line approximation
  Double_t maxDistStraight; //maximum distance for straight line approximation in the area with no B
  Double_t  minLengthCalcQ; //when the step size of propagation is larger than this value, process noise will be calculated

  Double_t  pathLocation;  //path position of the track; mm
  Double_t pathTime; // time of the track; ns
  Double_t  stepLength;   //distance between the Runge-Kutta start and end points; mm
  Double_t stepTime; //time between the Runge-Kutta start and end points; ns
  Double_t  trackLength;  //total track length; mm
  Double_t trackTime; //total track time; ns
  Double_t  energyLoss;   //total energy loss; MeV


  Int_t     jstep;        //Runge-Kutta step number
  Double_t  initialStepSize; //initial step size for the Runge-Kutta method
  
  Double_t  fDetMatProperties[kNMatType][4];

 private:
  ClassDef(fieldStepper,1) 
  
};

#endif
