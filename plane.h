#ifndef ROOT_PLANE
#define ROOT_PLANE

#include "utility.h"

class plane: public TObject{
 public:
  plane(planeId id, planeDir dir):fPlaneId(id), fPlaneDir(dir){}
  virtual ~plane(){}


  planeId getPlaneId() const{return fPlaneId;}
  planeDir getPlaneDir() const{return fPlaneDir;}


 protected:
  planeId fPlaneId;
  planeDir fPlaneDir;

 private:
  ClassDef(plane,1)      

};

#endif
