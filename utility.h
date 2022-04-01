#ifndef ROOT_Utility
#define ROOT_Utility

//#include "THashTable.h"
#include "TObject.h"
#include <iostream>
#include <algorithm>
#include <cassert>
using namespace std;
#define kDefaultValue -10000 //Some output information is available with satisfaction of specified conditions. If the information's is equal to kDefaultValue, the information is not available
#define kSdim 6
#define kMdim 3
#define MAXNTRACKS 1000
#define kNTotPlanes 5  // Total number of planes are used for tracking
#define kNMissingHitsLimit 2 //tracking process terminates if number of missing hits is over limit
#define kNMatType 7 //Total number of material types
enum kalIndex{ kIdxX = 0, kIdxY, kIdxZ, kIdxTheta, kIdxPhi, kIdxQP};
enum matType{kAir = 0, kMwpc2Mat,kMwpc3Mat,kMwpc4Mat,kAerogel,kAluminium,kScintillator};
enum matPropertyIdx{ kProtonOverAtomNum = 0, kExcitEnergy, kDensity, kRadLength }; 


enum planeId{kTarget=0, kSft, kTof1, kAc, kBLowEdge, kMwpc2, kMwpc3, kMwpc4, kBUpEdge, kTof2}; // plane == detector
enum planeDir{kPlaneDirX=0,kPlaneDirZ}; // detector is normal to x or z direction
enum seedType{kTargetMwpc34 = 0, kSftMwpc34, kMwpc2Mwpc34};
const Double_t kChi2PerMdimCut=1e8; // For each measurement in the filter,  DeltaChi2 limit = kChi2PerMdimCut * kMdim

enum trackId{kMST=0, kMS, kMT, kST, kM, kS, kT};

enum particleId{kMuonId=0,kPositronId,kPionPlusId};
const Double_t kParticleMass[3]={105.6583745, 0.510998918, 139.57018};
const Double_t kParticleCharge[3]={1., 1., 1.};
const Double_t kElectronMass = 0.51099891;
const Double_t kPionPlusMass = 139.57018;
const Double_t kMuonMass = 105.6583745;

const Double_t kMomLowerLimit=92.5;
const Double_t kMomUpperLimit=347.5;

#endif















