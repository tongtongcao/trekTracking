#include "ac.h"
#include <vector>
using namespace std;

ac::ac(planeId id, planeDir dir):plane(id,dir){
  //radiusTopHousingOuter=119;
  //radiusTopHousingInner=118.5;
  //radiusOuter=90.5;
  //radiusInner=50.5;
  //radiusBottomHousingOuter=50.5;
  //radiusBottomHousingInner=50;

  radiusTopHousingOuter=120.39;
  radiusTopHousingInner=119.884;
  radiusOuter=91.5568;
  radiusInner=51.0897;
  radiusBottomHousingOuter=51.0897;
  radiusBottomHousingInner=50.5839;
}

ac::~ac(){
}

ClassImp(ac);


