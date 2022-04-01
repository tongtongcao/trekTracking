#ifndef Field_Map
#define Field_Map

#define XSIZE 170 //cm
#define YSIZE 10  //cm
#define ZSIZE 100  //cm
#define STEP 1  //cm
#define XSTART 30  //cm

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TVector3.h"
#include "TMath.h"
#include "TObject.h"
using namespace std;

class fieldMap: public TObject{
  public: 
  ~fieldMap() {};
  static fieldMap * getInstance() {
    if (fInstance == NULL) fInstance = new fieldMap();
    return fInstance;
  }
  
  TVector3 getBField(Double_t x, Double_t y, Double_t z);  //mm
  
  fieldMap();
  static fieldMap *fInstance;
  void loadFieldMap();
  
  Double_t  bx[XSIZE+1][YSIZE+1][ZSIZE+1]; 
  Double_t  by[XSIZE+1][YSIZE+1][ZSIZE+1]; 
  Double_t  bz[XSIZE+1][YSIZE+1][ZSIZE+1]; 

   private:
   ClassDef(fieldMap,1)    
};  

#endif
