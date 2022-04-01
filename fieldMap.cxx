#include "fieldMap.h"

fieldMap * fieldMap::fInstance = NULL;

//__________________________________________________________________
fieldMap::fieldMap(){


  for (Int_t i=0; i<XSIZE+1; i++){
    for (Int_t j=0; j<YSIZE+1; j++) { 
      for (Int_t k=0; k<ZSIZE+1; k++) { 
	bx[i][j][k] = 0;
	by[i][j][k] = 0;
	bz[i][j][k] = 0;
      }
    }
  } 
  loadFieldMap();
}
//__________________________________________________________________
void fieldMap::loadFieldMap(){
  ifstream infile("B-map-I1150.txt");
  if (!infile.is_open()){
    cout<<"cannot open field map file"<<endl;
    exit(0);
  }
  Double_t x,y,z;
  for (Int_t i=0; i<XSIZE; i++){
    for (Int_t j=0; j<YSIZE; j++) { 
      for (Int_t k=0; k<ZSIZE; k++) { 
	infile>>x>>y>>z>>bx[i][j][k]>>by[i][j][k]>>bz[i][j][k]; 
	bx[i][j][k]=bx[i][j][k]/10000.; //convert from Gauss to Tesla
	by[i][j][k]=by[i][j][k]/10000.; //convert from Gauss to Tesla
	bz[i][j][k]=bz[i][j][k]/10000.; //convert from Gauss to Tesla
      }
    }
  }  
  infile.close();
}
//___________________________________________________________________

TVector3 fieldMap::getBField(Double_t x, Double_t y, Double_t z){
  //here use cm, other place use mm
  Double_t xPrime = 0.1*x;
  Double_t yPrime = 0.1*y;

  x = xPrime;
  y = yPrime;
  z = 0.1*z;

  Double_t phi = atan(y/x); // x > 0 in the B area, so phi = [-pi/2, pi/2]
  const Double_t phi30 = 30./180. * TMath::Pi();
  Double_t deltaPhi = ((int)(phi/phi30)) * phi30;

  if(deltaPhi != 0){
    x = cos(deltaPhi) * xPrime + sin(deltaPhi) * yPrime;
    y = -sin(deltaPhi) * xPrime + cos(deltaPhi) * yPrime; 
  }

  TVector3  fField;

  if(x<XSTART || x>XSTART+XSIZE || y<-YSIZE || y>YSIZE || z<-ZSIZE || z>ZSIZE){
    fField.Clear();
    return fField;
  }
  else{

    Int_t xMin,xMax,yMin,yMax,zMin,zMax;
    xMin=(Int_t)x;
    xMax=(Int_t)x+STEP;
    yMin=(Int_t)fabs(y);
    yMax=(Int_t)fabs(y)+STEP;
    zMin=(Int_t)fabs(z);
    zMax=(Int_t)fabs(z)+STEP;

    Double_t f111Bx=bx[xMin-XSTART][yMin][zMin];
    Double_t f112Bx=bx[xMin-XSTART][yMin][zMax];
    Double_t f121Bx=bx[xMin-XSTART][yMax][zMin];
    Double_t f122Bx=bx[xMin-XSTART][yMax][zMax];
    Double_t f211Bx=bx[xMax-XSTART][yMin][zMin];
    Double_t f212Bx=bx[xMax-XSTART][yMin][zMax];
    Double_t f221Bx=bx[xMax-XSTART][yMax][zMin];
    Double_t f222Bx=bx[xMax-XSTART][yMax][zMax];

    Double_t f111By=by[xMin-XSTART][yMin][zMin];
    Double_t f112By=by[xMin-XSTART][yMin][zMax];
    Double_t f121By=by[xMin-XSTART][yMax][zMin];
    Double_t f122By=by[xMin-XSTART][yMax][zMax];
    Double_t f211By=by[xMax-XSTART][yMin][zMin];
    Double_t f212By=by[xMax-XSTART][yMin][zMax];
    Double_t f221By=by[xMax-XSTART][yMax][zMin];
    Double_t f222By=by[xMax-XSTART][yMax][zMax];

    Double_t f111Bz=bz[xMin-XSTART][yMin][zMin];
    Double_t f112Bz=bz[xMin-XSTART][yMin][zMax];
    Double_t f121Bz=bz[xMin-XSTART][yMax][zMin];
    Double_t f122Bz=bz[xMin-XSTART][yMax][zMax];
    Double_t f211Bz=bz[xMax-XSTART][yMin][zMin];
    Double_t f212Bz=bz[xMax-XSTART][yMin][zMax];
    Double_t f221Bz=bz[xMax-XSTART][yMax][zMin];
    Double_t f222Bz=bz[xMax-XSTART][yMax][zMax];

    Double_t bxi=f111Bx*(xMax-x)*(yMax-fabs(y))*(zMax-fabs(z))+f112Bx*(xMax-x)*(yMax-fabs(y))*(fabs(z)-zMin)+f121Bx*(xMax-x)*(fabs(y)-yMin)*(zMax-fabs(z))+f122Bx*(xMax-x)*(fabs(y)-yMin)*(fabs(z)-zMin)+f211Bx*(x-xMin)*(yMax-fabs(y))*(zMax-fabs(z))+f212Bx*(x-xMin)*(yMax-fabs(y))*(fabs(z)-zMin)+f221Bx*(x-xMin)*(fabs(y)-yMin)*(zMax-fabs(z))+f222Bx*(x-xMin)*(fabs(y)-yMin)*(fabs(z)-zMin);
    Double_t byi=f111By*(xMax-x)*(yMax-fabs(y))*(zMax-fabs(z))+f112By*(xMax-x)*(yMax-fabs(y))*(fabs(z)-zMin)+f121By*(xMax-x)*(fabs(y)-yMin)*(zMax-fabs(z))+f122By*(xMax-x)*(fabs(y)-yMin)*(fabs(z)-zMin)+f211By*(x-xMin)*(yMax-fabs(y))*(zMax-fabs(z))+f212By*(x-xMin)*(yMax-fabs(y))*(fabs(z)-zMin)+f221By*(x-xMin)*(fabs(y)-yMin)*(zMax-fabs(z))+f222By*(x-xMin)*(fabs(y)-yMin)*(fabs(z)-zMin);
    Double_t bzi=f111Bz*(xMax-x)*(yMax-fabs(y))*(zMax-fabs(z))+f112Bz*(xMax-x)*(yMax-fabs(y))*(fabs(z)-zMin)+f121Bz*(xMax-x)*(fabs(y)-yMin)*(zMax-fabs(z))+f122Bz*(xMax-x)*(fabs(y)-yMin)*(fabs(z)-zMin)+f211Bz*(x-xMin)*(yMax-fabs(y))*(zMax-fabs(z))+f212Bz*(x-xMin)*(yMax-fabs(y))*(fabs(z)-zMin)+f221Bz*(x-xMin)*(fabs(y)-yMin)*(zMax-fabs(z))+f222Bz*(x-xMin)*(fabs(y)-yMin)*(fabs(z)-zMin);
    
    if(deltaPhi == 0){
      if(y>=0 && z>=0) fField.SetXYZ(bxi,byi,bzi);
      else if(y<0 && z>=0) fField.SetXYZ(-bxi,byi,-bzi);
      else if(y>=0 && z<0) fField.SetXYZ(bxi,byi,-bzi);
      else fField.SetXYZ(-bxi,byi,bzi);
    }
    else{
      if(y>=0 && z>=0) fField.SetXYZ(cos(deltaPhi) * bxi - sin(deltaPhi) * byi, sin(deltaPhi) * bxi + cos(deltaPhi) * byi, bzi);
      else if(y<0 && z>=0) fField.SetXYZ(cos(deltaPhi) * (-bxi) - sin(deltaPhi) * byi, sin(deltaPhi) * (-bxi) + cos(deltaPhi) * byi, -bzi);
      else if(y>=0 && z<0) fField.SetXYZ(cos(deltaPhi) * bxi - sin(deltaPhi) * byi, sin(deltaPhi) * bxi + cos(deltaPhi) * byi, -bzi);
      else fField.SetXYZ(cos(deltaPhi) * (-bxi) - sin(deltaPhi) * byi, sin(deltaPhi) * (-bxi) + cos(deltaPhi) * byi, bzi);
    }

    return fField;
  }
}

ClassImp(fieldMap)












