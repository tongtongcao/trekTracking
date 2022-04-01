#include <iostream>
#include <cassert>
#include <cmath>

#include "TMath.h"

#include "fieldStepper.h"
using namespace std;

fieldStepper *fieldStepper::fFieldStepper = 0;

//_________________________________________________________________
fieldStepper::fieldStepper(){
  fIsBackward = kTRUE;
  fIsMSOn = kTRUE;
  fIsDEDXOn = kTRUE;
  useDefaultStep();
  fFieldMap = fieldMap::getInstance();
  initDetMaterial();
}
//_________________________________________________________________
fieldStepper::~fieldStepper(){
  if(fFieldMap!=0) delete fFieldMap;
}
//_________________________________________________________________
void fieldStepper::useDefaultStep(){
  initialStepSize = 10;     //mm
  stepSizeDec     = 0.75;
  stepSizeInc     = 1.25;
  maxStepSize     = 100;     //mm
  minStepSize     = 1;    //mm
  minPrecision    = 1.e-2;
  maxPrecision    = 1.e-3;
  maxNumSteps     = 3e3/initialStepSize; //maximum number of steps between two hits
  maxDist         = 5.e-2;   //mm
  minLengthCalcQ  = 1.e-3;    //mm
  maxDistStraight = 1.e-4; // mm 
}

//////////////////////////////////// Propagation in B ///////////////////////////////
//_________________________________________________________________
void fieldStepper::transport(trackState *svFrom, planeHit *hitTo, matrix &sv, matrix &F, matrix &Q, particleId particleId, Double_t particleMass, Double_t particleCharge, Bool_t &magFlag, matType type){
  //-----------------------initialization begin--------------------------//
  pathLocation = svFrom->getPathLocation();
  pathTime = svFrom->getPathTime();
  trackLength    = 0.;
  trackTime = 0.;
  stepLength     = 0.;
  stepTime = 0.;
  jstep          = 0;
  energyLoss     = 0.;
  Double_t beta;
  Bool_t bCalcJac = kTRUE;
  
  if(hitTo->getPlaneDir() == kPlaneDirX){
    if(hitTo->getPosition().X() <= (*svFrom)(kIdxX, 0)) fIsBackward = kTRUE;
    else fIsBackward = kFALSE;
  }
  else if(hitTo->getPlaneDir() == kPlaneDirZ){
    if(hitTo->getPosition().Z() <= (*svFrom)(kIdxZ, 0)) fIsBackward = kTRUE;
    else fIsBackward = kFALSE;
  }

  matrix svTo(kSdim, 1);
  matrix svPreStep(kSdim, 1);
  matrix DF(kSdim, kSdim);                  // propagator matrix segment
  for (Int_t i=0; i<kSdim; i++){
    svTo(i,0) = (*svFrom)(i,0);
  }
  F.UnitMatrix(); // initialize F to unity
  Q.Zero();       // initialize Q to zero
  TVector3 posPreStep;
  TVector3 posAt; 
  posAt.SetXYZ(svTo(kIdxX,0),svTo(kIdxY,0),svTo(kIdxZ,0));
  
  Double_t step = initialStepSize; 
  Double_t nextStep = step;
  Double_t stepFac = 1.;
  Bool_t doNextStep = kTRUE;

  Double_t d = estimateDistanceToNextPlane(svFrom, hitTo);
  if(d < step){
    step = d;
    nextStep = step;
  }
  //-----------------------initialization end--------------------------//
  
  //-----------------------begin Runge-Kutta stepping------------------//
  magFlag = kTRUE;
  planeId hitToPlaneId;

  while (d >= maxDist && doNextStep && jstep < maxNumSteps){
    // Store some properties before stepping.
    posPreStep.SetXYZ(svTo(kIdxX, 0), svTo(kIdxY, 0), svTo(kIdxZ, 0));
    svPreStep = svTo;
    stepFac = rKPropagation(svTo, DF, step, bCalcJac, fIsBackward, particleMass, particleCharge); // do one step 
    posAt.SetXYZ(svTo(kIdxX, 0), svTo(kIdxY, 0), svTo(kIdxZ, 0));//update position vector after RK propagation

    Double_t sd = 0;
    if(hitTo->getPlaneDir() == kPlaneDirX){
      sd = svTo(kIdxX, 0) - hitTo->getPosition().X();
    }
    else if(hitTo->getPlaneDir() == kPlaneDirZ){
      sd = svTo(kIdxZ, 0) - hitTo->getPosition().Z();
    }
    if(((sd > 0) && fIsBackward == kFALSE) || ((sd < 0) && fIsBackward == kTRUE)){
      // Track went past target plane during propagtion. Repeating last step with reduced step size.
      step       *= stepSizeDec;
      nextStep   *= stepSizeDec;
      svTo  = svPreStep;
      posAt.SetXYZ(posPreStep.X(), posPreStep.Y(), posPreStep.Z());
      trackLength -= stepLength;
      trackTime -= stepTime;
      jstep--;
      if(fIsBackward == kFALSE) {
	pathLocation-=stepLength;
	pathTime-=stepTime;
      }
      else{
	pathLocation+=stepLength;
	pathTime+=stepTime;
      }
      continue;
    }

    hitToPlaneId=hitTo->getPlaneId();
    if(hitToPlaneId == kMwpc2 || hitToPlaneId == kMwpc3 || hitToPlaneId == kMwpc4 || hitToPlaneId == kBLowEdge || hitToPlaneId == kBUpEdge){ // Propagation between C4 and C2 must be within the magnetic area

      Double_t xPrime = svTo(kIdxX, 0);
      Double_t yPrime = svTo(kIdxY, 0);

      Double_t x = xPrime;
      Double_t y = yPrime;
      Double_t z = svTo(kIdxZ, 0);

      Double_t phi = atan(y/x); // x > 0 in the B area, so phi = [-pi/2, pi/2]
      const Double_t phi30 = 30./180. * TMath::Pi();
      Double_t deltaPhi = ((int)(phi/phi30)) * phi30;

      if(deltaPhi != 0){
	x = cos(deltaPhi) * xPrime + sin(deltaPhi) * yPrime;
	y = -sin(deltaPhi) * xPrime + cos(deltaPhi) * yPrime; 
      }

      if(x < 300 || x > 1600 || y < -100 || y > 100 || z < -400 || z > 1000){ 
	magFlag = kFALSE; 
	break;
      }
    }
  
    // After each Runge-Kutta step i the covariance matrix C_i would have to be
    // transported with the propagator matrix F_i:
    // C_1 = F_1 * C_0   * F_1^T
    // C_2 = F_2 * C_1   * F_2^T
    // ...
    // C_n = F_n * C_n-1 * F_n^T
    // This can be transformed to:
    // C_n = F_n * ... * F_2 * F_1 * C_0 * F_1^T * F_2^T * ... * F_n^T
    //     =             F         * C_0 *         F^T 
       
    F = DF * F; // update propagator

    //calculating process noice and energy loss 
    if ((isMSOn() || isDEDXOn()) && stepLength > minLengthCalcQ ){
	       
      beta = 1. / TMath::Sqrt( 1. + (particleMass*particleMass) * (svPreStep(kIdxQP, 0)/particleCharge) * (svPreStep(kIdxQP, 0)/particleCharge) ); 
      if ( isMSOn() ) calcMultScat(Q, svPreStep, stepLength, beta, type);
      if ( isDEDXOn() )   energyLoss += calcEnergyLoss(Q, svTo, stepLength, svPreStep(kIdxQP, 0), beta, particleId, particleMass, particleCharge, type);
    }
    //end calculating process noice and energy loss       

    // Decide if more steps must be done and calculate new step size
    nextStep = step * stepFac;
    d = estimateDistanceToNextPlane(svTo, hitTo); // if d < maxDist, loop ends   
    if(d > nextStep) step = nextStep;
    else step = d;

    if(hitTo->getPlaneDir() == kPlaneDirX){
      if(fIsBackward == kFALSE){
	if(posAt.X() < hitTo->getPosition().X()){ doNextStep = kTRUE; }
	else{ doNextStep = kFALSE; }
      }
      else{
	if(posAt.X() > hitTo->getPosition().X()){ doNextStep = kTRUE; }
	else { doNextStep = kFALSE;}
      }
    }
    if(hitTo->getPlaneDir() == kPlaneDirZ){
      if(fIsBackward == kFALSE){
	if(posAt.Z() < hitTo->getPosition().Z()){ doNextStep = kTRUE; }
	else{ doNextStep = kFALSE; }
      }
      else{
	if(posAt.Z() > hitTo->getPosition().Z()){ doNextStep = kTRUE; }
	else { doNextStep = kFALSE;}
      }
    }
  }

  //-----------------------end Runge-Kutta stepping--------------------//
  
  // To make sure the track position is on the target layer propagate to the target plane
  // using a straight line.
  if(magFlag == kTRUE){
    svPreStep = svTo;
    posPreStep.SetXYZ(svPreStep(kIdxX, 0), svPreStep(kIdxY, 0), svPreStep(kIdxZ, 0));

    assert(propagateStraightLine(svTo, DF, hitTo, particleMass, particleCharge));

    //calculate the energy loss and multiple scattering for this straight line propagation
    if ((isMSOn() || isDEDXOn()) && stepLength > minLengthCalcQ ){ 
      beta = 1. / TMath::Sqrt( 1. + (particleMass*particleMass) * (svPreStep(kIdxQP, 0)/particleCharge) * (svPreStep(kIdxQP, 0)/particleCharge) );
      if ( isMSOn() )     calcMultScat(Q, svPreStep, stepLength, beta);
      if ( isDEDXOn() )   energyLoss += calcEnergyLoss(Q, svTo, stepLength, svPreStep(kIdxQP, 0), beta, particleId, particleMass, particleCharge, type); 
    }

    F = DF * F; // final propagator matrix
  }
  sv = svTo;
}

Double_t fieldStepper::rKPropagation(matrix &stateVec, matrix &fPropStep, Double_t stepSize, Bool_t bCalcJac, Bool_t dir, Double_t particleMass, Double_t particleCharge){
  // One step of track tracing from track state

  const Int_t numPars = 5; // x, y, z, theta, phi
  const Double_t kappa = TMath::C() / 1.e9; // in MeV/(c * mm * T)

  // for rk stepping
  //Int_t step4;
  const Int_t rksteps = 4; // The 4 points used in Runge-Kutta stepping: start, 2 mid points and end
  const Int_t rkStart = 0;
  const Int_t rkMid1  = 1;
  const Int_t rkMid2  = 2;
  const Int_t rkEnd   = 3;

  // Constants for RK stepping.
  static Double_t a[rksteps] = { 0.0    , 0.5    , 0.5    , 1.0     };
  static Double_t b[rksteps] = { 0.0    , 0.5    , 0.5    , 1.0     };
  static Double_t c[rksteps] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };

  // Change of track parameters in each step.
  // First index: step point (0 = start; 1,2 = mid, 3 = end), second index: state parameter.
  Double_t k[rksteps][numPars];

  // for propagator matrix
  Double_t FXTheta[rksteps], FXPhi[rksteps], FYTheta[rksteps], FYPhi[rksteps], FZTheta[rksteps], FThetaPhi[rksteps], FThetaQP[rksteps], FPhiTheta[rksteps], FPhiPhi[rksteps], FPhiQP[rksteps];
  //----------------------------------------------------------------
  Double_t est     = 0.; // error estimation
  Double_t stepFac = 1.;

  TVector3 B;                  // B-field
  Double_t h = stepSize;       // step size
  if(fIsBackward == kTRUE){
    h *= -1;                 // stepping in negative z-direction
  }

  Double_t half    = h * 0.5;  // half step interval for fourth order RK
    
  Double_t qpIn   = stateVec(kIdxQP, 0);
  TVector3 posFrom = TVector3(stateVec(kIdxX, 0), stateVec(kIdxY, 0), stateVec(kIdxZ, 0));
  TVector3 posAt   = posFrom;
  Double_t pathIn = pathLocation;

  // Input track state vector and state vector during stepping.
  Double_t svIn[numPars], svStep[numPars];
  svIn[kIdxX] = stateVec(kIdxX, 0);
  svIn[kIdxY] = stateVec(kIdxY, 0);
  svIn[kIdxZ] = stateVec(kIdxZ, 0);
  svIn[kIdxTheta] = stateVec(kIdxTheta, 0);
  svIn[kIdxPhi] = stateVec(kIdxPhi, 0);

  //------------------------------------------------------------------------
  //   Runge-Kutta step
  //
  Int_t istep;
  Int_t ipar;
  
  do{
    half = h * 0.5;
    for(istep = 0; istep < rksteps; ++istep){  // k1,k2,k3,k4  (k1=start, k2,k3 = half, k4=end of interval)
      for(ipar = 0; ipar < numPars; ++ipar){  // 5 track parameters
	if(istep == 0){
	  svStep[ipar] = svIn[ipar];     // in first step copy input track parameters (x, y, z, theta, phi)
	} 
	else{
	  svStep[ipar] = svIn[ipar] + b[istep] * k[istep-1][ipar];     // do step
	}
      }
      pathLocation = pathIn + a[istep] * h; // move z along with track
      posAt.SetXYZ(svStep[kIdxX], svStep[kIdxY], svStep[kIdxZ]); // update z value for current position

            
      //get the magnatic field 
      B = fFieldMap->getBField(posAt.X(), posAt.Y(), posAt.Z());

      Double_t theta=svStep[kIdxTheta];
      Double_t phi=svStep[kIdxPhi];

      Double_t gX        = sin(theta)*cos(phi);
      Double_t gY        = sin(theta)*sin(phi);
      Double_t gZ        = cos(theta);
  
      Double_t wx = kappa * qpIn * B.X();
      Double_t wy = kappa * qpIn * B.Y();
      Double_t wz = kappa * qpIn * B.Z();

      Double_t gTheta = wx * sin(phi) - wy * cos(phi);
      Double_t gPhi = 1/tan(theta) * (wx * cos(phi) + wy * sin(phi)) - wz;

      // Change of track parameters in each step.
      k[istep][kIdxX]       = gX * h;              // dX
      k[istep][kIdxY]       = gY * h;              // dY
      k[istep][kIdxZ]       = gZ * h;              // dZ
      k[istep][kIdxTheta]   = gTheta * h; // dTheta
      k[istep][kIdxPhi]   = gPhi * h; // dPhi


      //------------------------------------------------------------------------
      // for transport matrix
      FXTheta[istep] = cos(theta)*cos(phi);
      FXPhi[istep] = -sin(theta)*sin(phi);

      FYTheta[istep] = cos(theta)*sin(phi);
      FYPhi[istep] = sin(theta)*cos(phi);

      FZTheta[istep] = -sin(theta);

      FThetaPhi[istep] = wx*cos(phi) + wy*sin(phi);
      FThetaQP[istep] = gTheta/qpIn;

      FPhiTheta[istep] = -1/sin(theta)/sin(theta) * (cos(phi)*wx +sin(phi)*wy);
      FPhiPhi[istep] = 1/tan(theta) * (-sin(phi)*wx + cos(phi)*wy);
      FPhiQP[istep] = gPhi/qpIn;
    }  // end of Runge-Kutta steps
    //------------------------------------------------------------------------

    //------------------------------------------------------------------------
    // error estimation ala Geant
    est = 0.;

    est += fabs(k[rkStart][kIdxX] + k[rkEnd][kIdxX] - k[rkMid1][kIdxX] - k[rkMid2][kIdxX]) * half;
    est += fabs(k[rkStart][kIdxY] + k[rkEnd][kIdxY] - k[rkMid1][kIdxY] - k[rkMid2][kIdxY]) * half;
    est += fabs(k[rkStart][kIdxZ] + k[rkEnd][kIdxZ] - k[rkMid1][kIdxZ] - k[rkMid2][kIdxZ]) * half;
    est += fabs(k[rkStart][kIdxTheta] + k[rkEnd][kIdxTheta] - k[rkMid1][kIdxTheta] - k[rkMid2][kIdxTheta]) * half;
    est += fabs(k[rkStart][kIdxPhi] + k[rkEnd][kIdxPhi] - k[rkMid1][kIdxPhi] - k[rkMid2][kIdxPhi]) * half;
    //------------------------------------------------------------------------
    //cout<<"h: "<<h <<" B: "<<B.Mag()<<"  "<<" est: "<<est<<endl;
    if(fabs(est) < minPrecision || fabs(h) <= minStepSize){
      // Find a step size with good precision
      jstep ++;
      break;
    } 
    else{
      // precision not good enough. make smaller step
      stepFac *= stepSizeDec;
      h *= stepSizeDec;
    }
  }while (jstep < maxNumSteps);
    
  if (fabs(est) < maxPrecision && fabs(h) < maxStepSize) {
    stepFac *= stepSizeInc;
  }
    
  //------------------------------------------------------------------------
  // set output track parameters
  for(ipar = 0; ipar < numPars; ++ ipar) {
    // yn+1        = yn         + 1/6*k1                    + 1/3*k2                  + 1/3*k3                  + 1/6*k4
    stateVec(ipar, 0) = svIn[ipar]+c[rkStart]*k[rkStart][ipar]+c[rkMid1]*k[rkMid1][ipar]+c[rkMid2]*k[rkMid2][ipar]+c[rkEnd]*k[rkEnd][ipar];
  }

  //------------------------------------------------------------------------
   
  stepLength = fabs(h);
  Double_t betaGamma=(particleCharge/qpIn+particleCharge/(stateVec(kIdxQP,0)))/particleMass/2;
  Double_t beta=1/sqrt(1/betaGamma/betaGamma+1);
  stepTime=stepLength*pow(10,6)/beta/TMath::C();
  trackLength += stepLength;  // calculate track length
  trackTime += stepTime;
  pathTime +=stepTime*h/fabs(h); 

  if(!bCalcJac) {
    return stepFac;
  }
    
  //------------------------------------------------------------------------
  //F propagate
  matrix I(kSdim,kSdim), *F[rksteps], *G[rksteps];
  I.UnitMatrix();
  for(istep = 0; istep < rksteps; ++istep){
    F[istep]=new matrix(kSdim,kSdim);
    F[istep]->Zero();
    G[istep]=new matrix(kSdim,kSdim);
    G[istep]->Zero();
  }
  for(istep = 0; istep < rksteps; ++istep){
    (*G[istep])(kIdxX,kIdxTheta) = FXTheta[istep];
    (*G[istep])(kIdxX,kIdxPhi) = FXPhi[istep];

    (*G[istep])(kIdxY,kIdxTheta) = FYTheta[istep];
    (*G[istep])(kIdxY,kIdxPhi) = FYPhi[istep];

    (*G[istep])(kIdxZ,kIdxTheta) = FZTheta[istep];

    (*G[istep])(kIdxTheta,kIdxPhi) = FThetaPhi[istep];
    (*G[istep])(kIdxTheta,kIdxQP) = FThetaQP[istep];

    (*G[istep])(kIdxPhi,kIdxTheta) = FPhiTheta[istep];
    (*G[istep])(kIdxPhi,kIdxPhi) = FPhiPhi[istep];
    (*G[istep])(kIdxPhi,kIdxQP) = FPhiQP[istep];

  }

  (*F[0]) = h * (*G[0]);
  (*F[1]) = h * (*G[1]) * (I + b[rkMid1] * (*F[0]));
  (*F[2]) = h * (*G[2]) * (I + b[rkMid1] * (*F[1]));
  (*F[3]) = h * (*G[3]) * (I + b[rkEnd] * (*F[2]));


  fPropStep.UnitMatrix();
  for(istep = 0; istep < rksteps; ++istep){
    fPropStep += c[istep] * (*F[istep]);
  }

  for(istep = 0; istep < rksteps; ++istep){
    delete F[istep];
    F[istep]=NULL;
    delete G[istep];
    G[istep]=NULL;
  }
  return stepFac;
}
//______________________________________________________________________________________________________________
void fieldStepper::propagateStraightLine(matrix &stateVec, matrix &fPropChange, Double_t step){
  Int_t h;
  if(fIsBackward == kFALSE) h = 1;
  else h = -1;

  // Update state vector.
  Double_t theta=stateVec(kIdxTheta, 0);
  Double_t phi=stateVec(kIdxPhi, 0);
  Double_t nx=sin(theta)*cos(phi);
  Double_t ny=sin(theta)*sin(phi);
  Double_t nz=cos(theta);
  stateVec(kIdxX, 0) += nx *step * h;
  stateVec(kIdxY, 0) += ny *step * h;
  stateVec(kIdxZ, 0) += nz *step * h;

  // Update propagator matrix.
  fPropChange.UnitMatrix();
  fPropChange(kIdxX, kIdxTheta) = cos(theta)*cos(phi)*step * h;
  fPropChange(kIdxX, kIdxPhi) = -sin(theta)*sin(phi)*step * h;

  fPropChange(kIdxY, kIdxTheta) = cos(theta)*sin(phi)*step * h;
  fPropChange(kIdxY, kIdxPhi) = sin(theta)*cos(phi)*step * h;

  fPropChange(kIdxZ, kIdxTheta) = -sin(theta)*step * h;

}
//______________________________________________________________________________________________________________
Bool_t fieldStepper::propagateStraightLine(matrix &stateVec, matrix &fPropChange, const planeHit *hitTo, Double_t particleMass, Double_t particleCharge){
  // From the position and direction stored in the track state vector, propagate the track
  // to a target plane using a straight line. The track state and reference layer are updated
  // and the propagator matrix is calculated. 
  //
  // Output:
  // fPropChange: Change in propagator matrix
  //
  // Input and output:
  // stateVec:    Current Track state vector.
  //
  // Input:
  // propDir:     Propagation direction.

  stepLength = estimateDistanceToNextPlane(stateVec, hitTo);
  Double_t betaGamma=particleCharge/(stateVec(kIdxQP,0))/particleMass;
  Double_t beta=1/sqrt(1/betaGamma/betaGamma+1);
  stepTime=fabs(stepLength/beta/TMath::C())*pow(10,6);

  Double_t sd = 0;
  if(hitTo->getPlaneDir() == kPlaneDirX){
    sd = stateVec(kIdxX, 0) - hitTo->getPosition().X();
  }
  else if(hitTo->getPlaneDir() == kPlaneDirZ){
    sd = stateVec(kIdxZ, 0) - hitTo->getPosition().Z();
  }

  if((sd <= 0. && fIsBackward == kFALSE) || (sd >= 0. && fIsBackward == kTRUE)) {
    propagateStraightLine(stateVec, fPropChange, stepLength);
    trackLength += stepLength;
    trackTime += stepTime;
    if(fIsBackward == kTRUE) {
      pathLocation-=stepLength;
      pathTime-=stepTime;
    }
    else{
      pathLocation+=stepLength;
      pathTime+=stepTime;
    }
    return kTRUE;
  } 
  else {
    fPropChange.UnitMatrix();
    return kFALSE;
  }
}
//_____________________________________________________________________________________________________________
Double_t fieldStepper::estimateDistanceToNextPlane(const trackState* svFrom, const planeHit *hitTo){
  Double_t dist=0;
  Double_t nx=sin((*svFrom)(kIdxTheta,0))*cos((*svFrom)(kIdxPhi,0));
  Double_t nz=cos((*svFrom)(kIdxTheta,0));
  if(hitTo->getPlaneDir() == kPlaneDirX) dist=fabs(((*svFrom)(kIdxX,0) - hitTo->getPosition().X())/nx);
  else   if(hitTo->getPlaneDir() == kPlaneDirZ) dist=fabs(((*svFrom)(kIdxZ,0) - hitTo->getPosition().Z())/nz);
  return dist;
}
//_____________________________________________________________________________________________________________
Double_t fieldStepper::estimateDistanceToNextPlane(const matrix svFrom, const planeHit *hitTo){
  Double_t dist=0;
  Double_t nx=sin(svFrom(kIdxTheta,0))*cos(svFrom(kIdxPhi,0));
  Double_t nz=cos(svFrom(kIdxTheta,0));
  if(hitTo->getPlaneDir() == kPlaneDirX) dist=fabs((svFrom(kIdxX,0) - hitTo->getPosition().X())/nx);
  else   if(hitTo->getPlaneDir() == kPlaneDirZ) dist=fabs((svFrom(kIdxZ,0) - hitTo->getPosition().Z())/nz);
  return dist;
}
//////////////////////////////////// Propagation in B ///////////////////////////////


//////////////////////////////////// Propagation without B ///////////////////////////////
//_____________________________________________________________________________________________________________
void fieldStepper::transportStraightLine(trackState *svFrom, planeHit *hitTo, matrix &sv, matrix &F, matrix &Q, particleId particleId, Double_t particleMass, Double_t particleCharge, matType type){
  //-----------------------initialization begin--------------------------//
  pathLocation = svFrom->getPathLocation();
  pathTime = svFrom->getPathTime();
  trackLength    = 0.;
  trackTime = 0.;
  stepLength     = 0.;
  stepTime = 0.;
  jstep          = 0;
  energyLoss     = 0.;
  Double_t beta;
  Double_t betaGamma;  

  if(hitTo->getPlaneDir() == kPlaneDirX){
    if(hitTo->getPosition().X() <= (*svFrom)(kIdxX, 0)) fIsBackward = kTRUE;
    else fIsBackward = kFALSE;
  }
  else if(hitTo->getPlaneDir() == kPlaneDirZ){
    if(hitTo->getPosition().Z() <= (*svFrom)(kIdxZ, 0)) fIsBackward = kTRUE;
    else fIsBackward = kFALSE;
  }

  matrix svTo(kSdim, 1);
  matrix svPreStep(kSdim, 1);
  matrix DF(kSdim, kSdim);                  // propagator matrix segment
  for (Int_t i=0; i<kSdim; i++){
    svTo(i,0) = (*svFrom)(i,0);
  }
  F.UnitMatrix(); // initialize F to unity
  Q.Zero();       // initialize Q to zero

  Double_t step;
  if(type == kAir) step = 10 * initialStepSize;
  else if( (type == kAluminium) || (type == kScintillator) ) step = 0.1 * initialStepSize;
  else step = initialStepSize;

  Double_t d = estimateDistanceToNextPlane(svFrom, hitTo);
  if(d < step) step = d;

  //-----------------------initialization end--------------------------//
  
  while (d >= maxDistStraight && jstep < maxNumSteps){
    svPreStep = svTo;
   
    stepLength = step;
    betaGamma=particleCharge/(svPreStep(kIdxQP,0))/particleMass;
    beta=1/sqrt(1/betaGamma/betaGamma+1);
    stepTime=fabs(stepLength/beta/TMath::C())*pow(10,6);
    trackLength += stepLength;
    trackTime += stepTime;
    if(fIsBackward == kTRUE) {
      pathLocation-=stepLength;
      pathTime-=stepTime;
    }
    else{
      pathLocation+=stepLength;
      pathTime+=stepTime;
    }

    propagateStraightLine(svTo, DF, step);
   
    //calculate the energy loss and multiple scattering for this straight line propagation
    if ((isMSOn() || isDEDXOn()) && stepLength > minLengthCalcQ){ 
      if ( isMSOn() )     calcMultScat(Q, svPreStep, stepLength, beta);
      if ( isDEDXOn() )   energyLoss += calcEnergyLoss(Q, svTo, stepLength, svPreStep(kIdxQP, 0), beta, particleId, particleMass, particleCharge, type); 
    }


    F = DF * F; // update propagator matrix

    // Decide if more steps must be done and calculate new step size
    d = estimateDistanceToNextPlane(svTo, hitTo); // d < maxDistStraight, loop ends    
    if(d < step) step = d; 
    jstep++; 
  }

  sv = svTo;
}
//////////////////////////////////// Propagation without B ///////////////////////////////

//////////////////////////////////// Energy loss and multi-scattering ///////////////////////////////
//_____________________________________________________________________________________________________________
Double_t fieldStepper::calcMultScat(matrix &Q, matrix &svFrom, Double_t length, Double_t beta, matType type){

  // Add multiple scattering to the process noise covariance.
  //
  // Input and output:
  // fProc:     Process noise covariance. Contribution of multiple scattering is added to this matrix.
  //
  // Input:
  // stateVec:  Track state vector at start of an RK step.
  // length:    Track length in cm.
  // radLength: Radiation length of passed material in cm.
  // beta:      v/c of particle.
  // pid:       Geant particle ID.
  
  Double_t lx0 = length * 1.e-3 / fDetMatProperties[type][kRadLength];
  Double_t cms2=13.6*13.6 /beta/beta *svFrom(kIdxQP,0)*svFrom(kIdxQP,0) * lx0 * TMath::Power((1 + 0.038 * TMath:: Log(lx0) ),2); 

  Double_t theta = svFrom(kIdxTheta, 0);

  Q(kIdxTheta, kIdxTheta) += cms2;
  Q(kIdxPhi, kIdxPhi) += cms2 / sin(theta) / sin(theta);

  return cms2;
}
//_____________________________________________________________________________________________________________
Double_t fieldStepper::calcEnergyLoss(matrix &Q, matrix &svTo, Double_t length, Double_t qp, Double_t beta, particleId particleId, Double_t particleMass, Double_t particleCharge, matType type)
{
  Double_t ZoverA = fDetMatProperties[type][kProtonOverAtomNum];
  Double_t p    = particleCharge / qp;
  Double_t elossIon = 0.;
  
  if(particleId == kPositronId) {

    // Critical energy for gases:
    // E_c = 710 MeV / (Z + 0.92)
    // From "Passage of particles through matter", Particle Data Group, 2009

    elossIon = length * calcDEDXIonPositron(qp, ZoverA, fDetMatProperties[type][kDensity],fDetMatProperties[type][kExcitEnergy]);
  } 
  else { // Energy loss for heavy particles.
    // Energy loss due to ionization.
    elossIon = length * calcDEDXBetheBloch(beta, ZoverA, fDetMatProperties[type][kDensity],fDetMatProperties[type][kExcitEnergy], particleMass, particleCharge);
    
  }
  if(fIsBackward == kTRUE) {
    elossIon *= -1.;
  }
  Double_t eloss = elossIon; // delta(E)
  Double_t p2 = p * p;
  Double_t e = TMath::Sqrt(p2 + particleMass*particleMass);

  Double_t pnew;
  if(p2 + 2.*e*eloss + eloss*eloss<0) pnew=-TMath::Sqrt(-(p2 + 2.*e*eloss + eloss*eloss)); // pnew is negative, means it is impossible; this kind of cases will be cancelled during tracking package
  else pnew = TMath::Sqrt(p2 + 2.*e*eloss + eloss*eloss);
  svTo(kIdxQP, 0) = particleCharge / pnew;

  return eloss;
}
//______________________________________________________________________________________________________________
Double_t fieldStepper::calcDEDXIonPositron(Double_t qp, Double_t ZoverA, Double_t density, Double_t I)
{
  // Calculates energy loss dE/dx in MeV/mm due to ionization for relativistic electrons/positrons.
  //
  // For formula used, see:
  // Track fitting with energy loss
  // Stampfer, Regler and Fruehwirth
  // Computer Physics Communication 79 (1994), 157-164
  //
  // qp:         particle charge divided by momentum
  // ZoverA:     atomic number / atomic mass of passed material
  // density:    density of material in g/mm^3
  // I : mean excitation energy in GeV


  Double_t K = 0.307075; // Mev * cm^2/g

  Double_t me       = 0.510998918;   // electron mass in MeV/c^2
  Double_t E        = sqrt(1/qp * 1/qp + me * me);             // Energy for relativistic electrons.
  Double_t gamma    = E / me;                          // Relativistic gamma-factor.

  // Formula is slightly different for electrons and positrons.
  // Factors for positrons:
  Double_t gammaFac = 4.;
  Double_t corr     = 2;

  Double_t dedx = 0.5 * K * density * ZoverA * (2 * TMath::Log(2*me/I) + gammaFac * TMath::Log(gamma) - corr);
  
  return (- dedx*0.1);//convert from MeV/cm to Mev/mm
}
//________________________________________________________________________________________________________________
Double_t fieldStepper::calcDEDXBetheBloch(Double_t beta, Double_t ZoverA, Double_t density, Double_t I, Double_t particleMass, Double_t particleCharge){
  // Returns the ionization energy loss for thin materials with Bethe-Bloch formula.
  // - dE/dx = K/A * Z * z^2 / beta^2 * (0.5 * ln(2*me*beta^2*gamma^2*T_max/I^2) - beta^2)
  // T_max = (2*me*beta^2*gamma^2) / (1 + 2*gamma*me/mass + me^2/mass^2)
  //
  // beta:       v/c of particle
  // mass:       mass of the particle in MeV/c^2 traversing through the material
  // ZoverA:     atomic number / atomic mass of passed material
  // density:    density of passed material in g/cm^3
  // exciteEner: mean excitation energy of passed material in MeV
  // z:          atomic number of particle

  Double_t me         = 0.510998918; // electron mass in MeV/c^2
  Double_t K          = 0.307075;               // in MeV * cm^2 / g for A = 1 g/mol

  Double_t beta2      = beta*beta;
  Double_t gamma      = 1./TMath::Sqrt(1 - beta2);
  Double_t betagamma  = beta * gamma;
  Double_t betagamma2 = betagamma * betagamma;

  if(betagamma < 0.1 || betagamma > 1000.) {
    //cout<<"Bethe-Bloch formula is only good for values between 0.1 and 1000."<<" beta= "<<beta<<" gamma="<<gamma<<endl;
    return 0.;
  }
  // Maximum kinetic energy that can be imparted on a free electron in a single collision.
  Double_t tmax = (2. * me * betagamma2) / (1 + 2*gamma*me/(particleMass) + (me*me)/(particleMass*particleMass));

  Double_t dedx = K *density * particleCharge * particleCharge * ZoverA / beta2 * (0.5 * TMath::Log((2. * me * betagamma2 * tmax)/(I*I)) - beta2);

  return (- dedx*0.1);//convert from MeV/cm to Mev/mm

}
//________________________________________________________________________________________________________________
void fieldStepper::initDetMaterial()
{
  //property array: effective proton number Z over effective atom number A, mean excitation energy (Mev), density (g/cm^3), radiation length(m)
  fDetMatProperties[kAir][kProtonOverAtomNum] = 0.49919;
  fDetMatProperties[kAir][kExcitEnergy] = 85.7e-6;  //MeV
  fDetMatProperties[kAir][kDensity] = 1.205e-3; //g/cm^3
  fDetMatProperties[kAir][kRadLength] = 303.9; //m
  
  fDetMatProperties[kMwpc2Mat][kProtonOverAtomNum] = 0.46481;
  fDetMatProperties[kMwpc2Mat][kExcitEnergy] = 153.1e-6;  //MeV
  fDetMatProperties[kMwpc2Mat][kDensity] = 9.19e-3; //g/cm^3
  fDetMatProperties[kMwpc2Mat][kRadLength] = 20.8; //m

  fDetMatProperties[kMwpc3Mat][kProtonOverAtomNum] = 0.46476;
  fDetMatProperties[kMwpc3Mat][kExcitEnergy] = 153.5e-6;  //MeV
  fDetMatProperties[kMwpc3Mat][kDensity] = 9.25e-3; //g/cm^3
  fDetMatProperties[kMwpc3Mat][kRadLength] = 20.6; //m

  fDetMatProperties[kMwpc4Mat][kProtonOverAtomNum] = 0.46476;
  fDetMatProperties[kMwpc4Mat][kExcitEnergy] = 153.5e-6;  //MeV
  fDetMatProperties[kMwpc4Mat][kDensity] = 9.25e-3; //g/cm^3
  fDetMatProperties[kMwpc4Mat][kRadLength] = 20.6; //m

  fDetMatProperties[kAerogel][kProtonOverAtomNum] = 0.50034;
  fDetMatProperties[kAerogel][kExcitEnergy] = 98.880e-6;  //MeV
  fDetMatProperties[kAerogel][kDensity] = 0.2; //g/cm^3
  fDetMatProperties[kAerogel][kRadLength] = 1.493; //m

  fDetMatProperties[kAluminium][kProtonOverAtomNum] = 13.0/26.9815385;
  fDetMatProperties[kAluminium][kExcitEnergy] = 166.0e-6;  //MeV
  fDetMatProperties[kAluminium][kDensity] = 2.699; //g/cm^3
  fDetMatProperties[kAluminium][kRadLength] = 8.896e-2; //m

  fDetMatProperties[kScintillator][kProtonOverAtomNum] = 0.50341;
  fDetMatProperties[kScintillator][kExcitEnergy] = 64.684e-6;  //MeV
  fDetMatProperties[kScintillator][kDensity] = 1.032; //g/cm^3
  fDetMatProperties[kScintillator][kRadLength] = 42.549e-2; //m
}
//////////////////////////////////// Energy loss and multi-scattering ///////////////////////////////

ClassImp(fieldStepper)


















