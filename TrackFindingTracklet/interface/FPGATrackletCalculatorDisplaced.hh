//This class implementes the tracklet engine
#ifndef FPGATRACKLETCALCULATORDISPLACED_H
#define FPGATRACKLETCALCULATORDISPLACED_H

#include "IMATH_TrackletCalculator.hh"
#include "IMATH_TrackletCalculatorDisk.hh"
#include "IMATH_TrackletCalculatorOverlap.hh"

#include "FPGAProcessBase.hh"
#include "FPGATrackletProjections.hh"
#include "FPGAStubTriplets.hh"
#include "FPGAConstants.hh"

using namespace std;

class FPGATrackletCalculatorDisplaced:public FPGAProcessBase{

public:

  FPGATrackletCalculatorDisplaced(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    double dphi=two_pi/NSector;
    double dphiHG=0.0;
    if (hourglass) {
      dphiHG=0.5*(dphisectorHG-two_pi/NSector);
    }
    phimin_=iSector_*dphi-dphiHG;
    phimax_=phimin_+dphi+2*dphiHG;
    if (hourglass) {
      phimin_-=0.5*two_pi/NSector;
      phimax_-=0.5*two_pi/NSector;
    }
    if (phimin_>0.5*two_pi) phimin_-=two_pi;
    if (phimax_>0.5*two_pi) phimax_-=two_pi;
    if (phimin_>phimax_)  phimin_-=two_pi;
    if (hourglass) {
      phioffset_=phimin_;
    } else {
      phioffset_=phimin_-dphi/6.0;
    }

    maxtracklet_=63;
    if (hourglass) maxtracklet_=127;
    
   trackletproj_L1PHI1_=0;
   trackletproj_L1PHI2_=0;
   trackletproj_L1PHI3_=0;
   trackletproj_L1PHI4_=0;
   trackletproj_L1PHI5_=0;
   trackletproj_L1PHI6_=0;
   trackletproj_L1PHI7_=0;
   trackletproj_L1PHI8_=0;
   

   trackletproj_L2PHI1_=0;
   trackletproj_L2PHI2_=0;
   trackletproj_L2PHI3_=0;
   trackletproj_L2PHI4_=0;

   trackletproj_L3PHI1_=0;
   trackletproj_L3PHI2_=0;
   trackletproj_L3PHI3_=0;
   trackletproj_L3PHI4_=0;

   trackletproj_L4PHI1_=0;
   trackletproj_L4PHI2_=0;
   trackletproj_L4PHI3_=0;
   trackletproj_L4PHI4_=0;

   trackletproj_L5PHI1_=0;
   trackletproj_L5PHI2_=0;
   trackletproj_L5PHI3_=0;
   trackletproj_L5PHI4_=0;

   trackletproj_L6PHI1_=0;
   trackletproj_L6PHI2_=0;
   trackletproj_L6PHI3_=0;
   trackletproj_L6PHI4_=0;

   trackletproj_L1Plus_=0; 
   trackletproj_L1Minus_=0;
                         
   trackletproj_L2Plus_=0; 
   trackletproj_L2Minus_=0;
                         
   trackletproj_L3Plus_=0; 
   trackletproj_L3Minus_=0;
                         
   trackletproj_L4Plus_=0; 
   trackletproj_L4Minus_=0;
                         
   trackletproj_L5Plus_=0; 
   trackletproj_L5Minus_=0;
                         
   trackletproj_L6Plus_=0; 
   trackletproj_L6Minus_=0;

   trackletproj_D1PHI1_=0;
   trackletproj_D1PHI2_=0;
   trackletproj_D1PHI3_=0;
   trackletproj_D1PHI4_=0;

   trackletproj_D2PHI1_=0;
   trackletproj_D2PHI2_=0;
   trackletproj_D2PHI3_=0;
   trackletproj_D2PHI4_=0;

   trackletproj_D3PHI1_=0;
   trackletproj_D3PHI2_=0;
   trackletproj_D3PHI3_=0;
   trackletproj_D3PHI4_=0;

   trackletproj_D4PHI1_=0;
   trackletproj_D4PHI2_=0;
   trackletproj_D4PHI3_=0;
   trackletproj_D4PHI4_=0;

   trackletproj_D5PHI1_=0;
   trackletproj_D5PHI2_=0;
   trackletproj_D5PHI3_=0;
   trackletproj_D5PHI4_=0;

   trackletproj_D1Plus_=0; 
   trackletproj_D1Minus_=0;
                         
   trackletproj_D2Plus_=0; 
   trackletproj_D2Minus_=0;
                         
   trackletproj_D3Plus_=0; 
   trackletproj_D3Minus_=0;
                         
   trackletproj_D4Plus_=0; 
   trackletproj_D4Minus_=0;
                         
   trackletproj_D5Plus_=0; 
   trackletproj_D5Minus_=0;

  
   layer_=0;
   disk_=0;

   string name1 = name.substr(1);//this is to correct for "TCD" having one more letter then "TC"
   if (name1[3]=='L') layer_=name1[4]-'0';    
   if (name1[3]=='D') disk_=name1[4]-'0';    


   // set TC index
   int iTC = -1;
   int iSeed = -1;
   
   if      (name1[9]=='A') iTC =0;
   else if (name1[9]=='B') iTC =1;
   else if (name1[9]=='C') iTC =2;
   else if (name1[9]=='D') iTC =3;
   else if (name1[9]=='E') iTC =4;
   else if (name1[9]=='F') iTC =5;
   else if (name1[9]=='G') iTC =6;
   else if (name1[9]=='H') iTC =7;
   else if (name1[9]=='I') iTC =8;
   else if (name1[9]=='J') iTC =9;
   else if (name1[9]=='K') iTC =10;
   else if (name1[9]=='L') iTC =11;
   else if (name1[9]=='M') iTC =12;
   else if (name1[9]=='N') iTC =13;
   else if (name1[9]=='O') iTC =14;

   assert(iTC!=-1);
   
   if (name1.substr(3,6)=="L3L4L2") iSeed = 8;
   else if (name1.substr(3,6)=="L5L6L4") iSeed = 9;
   else if (name1.substr(3,6)=="L2L3D1") iSeed = 10;
   else if (name1.substr(3,6)=="D1D2L2") iSeed = 11;

   assert(iSeed!=-1);

   if (hourglass) {
     TCIndex_ = (iSeed<<4) + iTC;
     assert(TCIndex_>=128 && TCIndex_<191);
   } else {
     TCIndex_ = (iSeed<<3) + iTC;
     assert(TCIndex_>=0 && TCIndex_<64);
   }
     
   //if (hourglass) {
   //if (iSeed!=0)  TCIndex_+=8;
     //cout << "iTC iSeed TCIndex_ "<<iTC<<" "<<iSeed<<" "<<TCIndex_<<endl;
   //}
   
   
   assert((layer_!=0)||(disk_!=0));

   
   if (iSeed==8||iSeed==9) {
     if (layer_==3) {
       rproj_[0]=rmeanL1;
       rproj_[1]=rmeanL2;
       rproj_[2]=rmeanL5;
       rproj_[3]=rmeanL6;
       lproj_[0]=1;
       lproj_[1]=2;
       lproj_[2]=5;
       lproj_[3]=6;
     }
      
     if (layer_==5) {
       rproj_[0]=rmeanL1;
       rproj_[1]=rmeanL2;
       rproj_[2]=rmeanL3;
       rproj_[3]=rmeanL4;
       lproj_[0]=1;
       lproj_[1]=2;
       lproj_[2]=3;
       lproj_[3]=4;
     }
   }


   if (iSeed==10||iSeed==11) {
     if (layer_==2) {
       zprojoverlap_[0]=zmeanD1;
       zprojoverlap_[1]=zmeanD2;
       zprojoverlap_[2]=zmeanD3;
       zprojoverlap_[3]=zmeanD4;
     }

     if (disk_==1) {
       zprojoverlap_[0]=zmeanD2;
       zprojoverlap_[1]=zmeanD3;
       zprojoverlap_[2]=zmeanD4;
       zprojoverlap_[3]=zmeanD5;
     }
   }
      
  }

  void addOutputProjection(FPGATrackletProjections* &outputProj, FPGAMemoryBase* memory){
      outputProj=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(outputProj!=0);
  }
  
  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="trackpar"){
      FPGATrackletParameters* tmp=dynamic_cast<FPGATrackletParameters*>(memory);
      assert(tmp!=0);
      trackletpars_=tmp;
      return;
    }


    if (output=="projoutL1PHI1"||output=="projoutL1PHIA") {
      addOutputProjection(trackletproj_L1PHI1_,memory);
      return;
    }
    
    if (output=="projoutL1PHI2"||output=="projoutL1PHIB") {
      addOutputProjection(trackletproj_L1PHI2_,memory);
      return;
    }

    if (output=="projoutL1PHI3"||output=="projoutL1PHIC"){
      addOutputProjection(trackletproj_L1PHI3_,memory);
      return;
    }

    if (output=="projoutL1PHID"){
      addOutputProjection(trackletproj_L1PHI4_,memory);
      return;
    }

    if (output=="projoutL1PHIE"){
      addOutputProjection(trackletproj_L1PHI5_,memory);
      return;
    }

    if (output=="projoutL1PHIF"){
      addOutputProjection(trackletproj_L1PHI6_,memory);
      return;
    }

    if (output=="projoutL1PHIG"){
      addOutputProjection(trackletproj_L1PHI7_,memory);
      return;
    }

    if (output=="projoutL1PHIH"){
      addOutputProjection(trackletproj_L1PHI8_,memory);
      return;
    }

    if (output=="projoutL2PHI1"||output=="projoutL2PHIA"){
      addOutputProjection(trackletproj_L2PHI1_,memory);
      return;
    }

    if (output=="projoutL2PHI2"||output=="projoutL2PHIB"){
      addOutputProjection(trackletproj_L2PHI2_,memory);
      return;
    }

    if (output=="projoutL2PHI3"||output=="projoutL2PHIC"){
      addOutputProjection(trackletproj_L2PHI3_,memory);
      return;
    }

    if (output=="projoutL2PHI4"||output=="projoutL2PHID"){
      addOutputProjection(trackletproj_L2PHI4_,memory);
      return;
    }

    if (output=="projoutL3PHI1"||output=="projoutL3PHIA"){
      addOutputProjection(trackletproj_L3PHI1_,memory);
      return;
    }

    if (output=="projoutL3PHI2"||output=="projoutL3PHIB"){
      addOutputProjection(trackletproj_L3PHI2_,memory);
      return;
    }

    if (output=="projoutL3PHI3"||output=="projoutL3PHIC"){
      addOutputProjection(trackletproj_L3PHI3_,memory);
      return;
    }

    if (output=="projoutL3PHI4"||output=="projoutL3PHID"){
      addOutputProjection(trackletproj_L3PHI4_,memory);
      return;
    }

    if (output=="projoutL4PHI1"||output=="projoutL4PHIA"){
      addOutputProjection(trackletproj_L4PHI1_,memory);
      return;
    }

    if (output=="projoutL4PHI2"||output=="projoutL4PHIB"){
      addOutputProjection(trackletproj_L4PHI2_,memory);
      return;
    }

    if (output=="projoutL4PHI3"||output=="projoutL4PHIC"){
      addOutputProjection(trackletproj_L4PHI3_,memory);
      return;
    }

    if (output=="projoutL4PHI4"||output=="projoutL4PHID"){
      addOutputProjection(trackletproj_L4PHI4_,memory);
      return;
    }

    if (output=="projoutL5PHI1"||output=="projoutL5PHIA"){
      addOutputProjection(trackletproj_L5PHI1_,memory);
      return;
    }

    if (output=="projoutL5PHI2"||output=="projoutL5PHIB"){
      addOutputProjection(trackletproj_L5PHI2_,memory);
      return;
    }

    if (output=="projoutL5PHI3"||output=="projoutL5PHIC"){
      addOutputProjection(trackletproj_L5PHI3_,memory);
      return;
    }

    if (output=="projoutL5PHI4"||output=="projoutL5PHID"){
      addOutputProjection(trackletproj_L5PHI4_,memory);
      return;
    }

    if (output=="projoutL6PHI1"||output=="projoutL6PHIA"){
      addOutputProjection(trackletproj_L6PHI1_,memory);
      return;
    }

    if (output=="projoutL6PHI2"||output=="projoutL6PHIB"){
      addOutputProjection(trackletproj_L6PHI2_,memory);
      return;
    }

    if (output=="projoutL6PHI3"||output=="projoutL6PHIC"){
      addOutputProjection(trackletproj_L6PHI3_,memory);
      return;
    }

    if (output=="projoutL6PHI4"||output=="projoutL6PHID"){
      addOutputProjection(trackletproj_L6PHI4_,memory);
      return;
    }

    if (output=="projoutD1PHI1"||output=="projoutD1PHIA"){
      addOutputProjection(trackletproj_D1PHI1_,memory);
      return;
    }

    if (output=="projoutD1PHI2"||output=="projoutD1PHIB"){
      addOutputProjection(trackletproj_D1PHI2_,memory);
      return;
    }

    if (output=="projoutD1PHI3"||output=="projoutD1PHIC"){
      addOutputProjection(trackletproj_D1PHI3_,memory);
      return;
    }

    if (output=="projoutD1PHI4"||output=="projoutD1PHID"){
      addOutputProjection(trackletproj_D1PHI4_,memory);
      return;
    }

    if (output=="projoutD2PHI1"||output=="projoutD2PHIA"){
      addOutputProjection(trackletproj_D2PHI1_,memory);
      return;
    }

    if (output=="projoutD2PHI2"||output=="projoutD2PHIB"){
      addOutputProjection(trackletproj_D2PHI2_,memory);
      return;
    }

    if (output=="projoutD2PHI3"||output=="projoutD2PHIC"){
      addOutputProjection(trackletproj_D2PHI3_,memory);
      return;
    }

    if (output=="projoutD2PHI4"||output=="projoutD2PHID"){
      addOutputProjection(trackletproj_D2PHI4_,memory);
      return;
    }



    if (output=="projoutD3PHI1"||output=="projoutD3PHIA"){
      addOutputProjection(trackletproj_D3PHI1_,memory);
      return;
    }

    if (output=="projoutD3PHI2"||output=="projoutD3PHIB"){
      addOutputProjection(trackletproj_D3PHI2_,memory);
      return;
    }

    if (output=="projoutD3PHI3"||output=="projoutD3PHIC"){
      addOutputProjection(trackletproj_D3PHI3_,memory);
      return;
    }
    
    if (output=="projoutD3PHI4"||output=="projoutD3PHID"){
      addOutputProjection(trackletproj_D3PHI4_,memory);
      return;
    }


    if (output=="projoutD4PHI1"||output=="projoutD4PHIA"){
      addOutputProjection(trackletproj_D4PHI1_,memory);
      return;
    }

    if (output=="projoutD4PHI2"||output=="projoutD4PHIB"){
      addOutputProjection(trackletproj_D4PHI2_,memory);
      return;
    }

    if (output=="projoutD4PHI3"||output=="projoutD4PHIC"){
      addOutputProjection(trackletproj_D4PHI3_,memory);
      return;
    }

    if (output=="projoutD4PHI4"||output=="projoutD4PHID"){
      addOutputProjection(trackletproj_D4PHI4_,memory);
      return;
    }
    


    if (output=="projoutD5PHI1"||output=="projoutD5PHIA"){
      addOutputProjection(trackletproj_D5PHI1_,memory);
      return;
    }

    if (output=="projoutD5PHI2"||output=="projoutD5PHIB"){
      addOutputProjection(trackletproj_D5PHI2_,memory);
      return;
    }

    if (output=="projoutD5PHI3"||output=="projoutD5PHIC"){
      addOutputProjection(trackletproj_D5PHI3_,memory);
      return;
    }

    if (output=="projoutD5PHI4"||output=="projoutD5PHID"){
      addOutputProjection(trackletproj_D5PHI4_,memory);
      return;
    }


    
    if (output=="projoutL1ToMinus"){
      addOutputProjection(trackletproj_L1Minus_,memory);
      return;
    }

    if (output=="projoutL1ToPlus"){
      addOutputProjection(trackletproj_L1Plus_,memory);
      return;
    }

    if (output=="projoutL2ToMinus"){
      addOutputProjection(trackletproj_L2Minus_,memory);
      return;
    }

    if (output=="projoutL2ToPlus"){
      addOutputProjection(trackletproj_L2Plus_,memory);
      return;
    }

    if (output=="projoutL3ToMinus"){
      addOutputProjection(trackletproj_L3Minus_,memory);
      return;
    }

    if (output=="projoutL3ToPlus"){
      addOutputProjection(trackletproj_L3Plus_,memory);
      return;
    }

    if (output=="projoutL4ToMinus"){
      addOutputProjection(trackletproj_L4Minus_,memory);
      return;
    }

    if (output=="projoutL4ToPlus"){
      addOutputProjection(trackletproj_L4Plus_,memory);
      return;
    }

    if (output=="projoutL5ToMinus"){
      addOutputProjection(trackletproj_L5Minus_,memory);
      return;
    }

    if (output=="projoutL5ToPlus"){
      addOutputProjection(trackletproj_L5Plus_,memory);
      return;
    }

    if (output=="projoutL6ToMinus"){
      addOutputProjection(trackletproj_L6Minus_,memory);
      return;
    }

    if (output=="projoutL6ToPlus"){
      addOutputProjection(trackletproj_L6Plus_,memory);
      return;
    }

    if (output=="projoutL3D4ToMinus"){
      addOutputProjection(trackletproj_L3Minus_,memory);
      return;
    }

    if (output=="projoutL3D4ToPlus"){
      addOutputProjection(trackletproj_L3Plus_,memory);
      return;
    }

    if (output=="projoutL4D3ToMinus"){
      addOutputProjection(trackletproj_L4Minus_,memory);
      return;
    }

    if (output=="projoutL4D3ToPlus"){
      addOutputProjection(trackletproj_L4Plus_,memory);
      return;
    }

    if (output=="projoutL5D2ToMinus"){
      addOutputProjection(trackletproj_L5Minus_,memory);
      return;
    }

    if (output=="projoutL5D2ToPlus"){
      addOutputProjection(trackletproj_L5Plus_,memory);
      return;
    }

    if (output=="projoutL6D1ToMinus"){
      addOutputProjection(trackletproj_L6Minus_,memory);
      return;
    }

    if (output=="projoutL6D1ToPlus"){
      addOutputProjection(trackletproj_L6Plus_,memory);
      return;
    }


    if (output=="projoutD1ToPlus"){
      addOutputProjection(trackletproj_D1Plus_,memory);
      return;
    }

    if (output=="projoutD2ToPlus"){
      addOutputProjection(trackletproj_D2Plus_,memory);
      return;
    }

    if (output=="projoutD3ToPlus"){
      addOutputProjection(trackletproj_D3Plus_,memory);
      return;
    }

    if (output=="projoutD4ToPlus"){
      addOutputProjection(trackletproj_D4Plus_,memory);
      return;
    }

    if (output=="projoutD5ToPlus"){
      addOutputProjection(trackletproj_D5Plus_,memory);
      return;
    }    
    

    if (output=="projoutD1ToMinus"){
      addOutputProjection(trackletproj_D1Minus_,memory);
      return;
    }

    if (output=="projoutD2ToMinus"){
      addOutputProjection(trackletproj_D2Minus_,memory);
      return;
    }

    if (output=="projoutD3ToMinus"){
      addOutputProjection(trackletproj_D3Minus_,memory);
      return;
    }

    if (output=="projoutD4ToMinus"){
      addOutputProjection(trackletproj_D4Minus_,memory);
      return;
    }

    if (output=="projoutD5ToMinus"){
      addOutputProjection(trackletproj_D5Minus_,memory);
      return;
    }    
    

    cout << "Could not find output : "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="thirdallstubin"){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      innerallstubs_.push_back(tmp);
      return;
    }
    if (input=="firstallstubin"){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      middleallstubs_.push_back(tmp);
      return;
    }
    if (input=="secondallstubin"){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      outerallstubs_.push_back(tmp);
      return;
    }
    if (input.find("stubtriplet")==0){
      FPGAStubTriplets* tmp=dynamic_cast<FPGAStubTriplets*>(memory);
      assert(tmp!=0);
      stubtriplets_.push_back(tmp);
      return;
    }
    assert(0);
  }


  void exacttracklet(double r1, double z1, double phi1,
		     double r2, double z2, double phi2, double sigmaz,
		     double& rinv, double& phi0,
		     double& t, double& z0,
		     double phiproj[4], double zproj[4], 
		     double phider[4], double zder[4],
		     double phiprojdisk[5], double rprojdisk[5], 
		     double phiderdisk[5], double rderdisk[5]) {

    double deltaphi=phi1-phi2;

    if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
    if (deltaphi<-0.5*two_pi) deltaphi+=two_pi;
    assert(fabs(deltaphi)<0.5*two_pi);

    
    if (sigmaz<-10.0) {
      cout << "Negative sigmaz"<<endl;
    }

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    double phi1tmp=phi1-phimin_+(phimax_-phimin_)/6.0;    

    phi0=phi1tmp+asin(0.5*r1*rinv);
    
    if (phi0>0.5*two_pi) phi0-=two_pi;
    if (phi0<-0.5*two_pi) phi0+=two_pi;
    assert(fabs(phi0)<0.5*two_pi);
    
    double rhopsi1=2*asin(0.5*r1*rinv)/rinv;
	    
    double rhopsi2=2*asin(0.5*r2*rinv)/rinv;
    
    t=(z1-z2)/(rhopsi1-rhopsi2);
    
    z0=z1-t*rhopsi1;

    for (int i=0;i<4;i++) {
      exactproj(rproj_[i],rinv,phi0,t,z0,
		phiproj[i],zproj[i],phider[i],zder[i]);
    }

    for (int i=0;i<5;i++) {
      //int sign=1;
      //if (t<0) sign=-1;
      exactprojdisk(zmean[i],rinv,phi0,t,z0,
		phiprojdisk[i],rprojdisk[i],phiderdisk[i],rderdisk[i]);
    }



  }

  void exacttracklet(double r1, double z1, double phi1,
		     double r2, double z2, double phi2, double sigmaz,
		     double r3, double z3, double phi3,
		     double& rinv, double& phi0, double &d0,
		     double& t, double& z0,
		     double phiproj[4], double zproj[4], 
		     double phider[4], double zder[4],
		     double phiprojdisk[5], double rprojdisk[5], 
		     double phiderdisk[5], double rderdisk[5]) {


    //two lines perpendicular to the 1->2 and 2->3
    double x1 = r1*cos(phi1);
    double x2 = r2*cos(phi2);
    double x3 = r3*cos(phi3);
    
    double y1 = r1*sin(phi1);
    double y2 = r2*sin(phi2);
    double y3 = r3*sin(phi3);
    
    double k1 = - (x2-x1)/(y2-y1);
    double k2 = - (x3-x2)/(y3-y2);
    double b1 = 0.5*(y2+y1)-0.5*(x1+x2)*k1;
    double b2 = 0.5*(y3+y2)-0.5*(x2+x3)*k2;
    //their intersection gives the center of the circle
    double y0 = (b1*k2-b2*k1)/(k2-k1);
    double x0 = (b1-b2)/(k2-k1);
    //get the radius three ways:
    double R1 = sqrt(pow(x1-x0,2)+pow(y1-y0,2));
    double R2 = sqrt(pow(x2-x0,2)+pow(y2-y0,2));
    double R3 = sqrt(pow(x3-x0,2)+pow(y3-y0,2));
    //check if the same
    double eps1 = fabs(R1/R2-1);
    double eps2 = fabs(R3/R2-1);
    if(eps1>1e-10 || eps2>1e-10)
      cout<<"&&&&&&&&&&&& bad circle! "<<R1<<"\t"<<R2<<"\t"<<R3<<"\n";

    //results
    rinv = 1./R1;
    phi0 = atan(1.)*2 + atan2(y0,x0);
    phi0 += -phimin_+(phimax_-phimin_)/6.0; 
    d0 = R1 - sqrt(x0*x0+y0*y0);
    //sign of rinv:
    double dphi = phi3 - atan2(y0,x0);
    if(dphi> 3.1415927) dphi -= 6.283185;
    if(dphi<-3.1415927) dphi += 6.283185;
    if(dphi<0) {
      rinv = -rinv;
      d0 = -d0;
      phi0 = phi0 + 4*atan(1);
      if(phi0 >8*atan(1)) phi0 = phi0-8*atan(1);
    }

    // //comparison
    // double deltaphi=phi2-phi3;
    // if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
    // if (deltaphi<-0.5*two_pi) deltaphi+=two_pi;
    // assert(fabs(deltaphi)<0.5*two_pi);
    // double dist=sqrt(r3*r3+r2*r2-2*r2*r3*cos(deltaphi));
    // double comp_rinv=2*sin(deltaphi)/dist;
    // double comp_phi0=phi2+asin(0.5*r2*rinv);
    
    // cout<<rinv<<"\t"<<comp_rinv<<"\t: "<<phi0<<"\t"<<comp_phi0<<"\n";
      
    //now in RZ:
    //turning angle
    double beta1 = atan2(y1-y0,x1-x0)-atan2(-y0,-x0);
    double beta2 = atan2(y2-y0,x2-x0)-atan2(-y0,-x0);
    double beta3 = atan2(y3-y0,x3-x0)-atan2(-y0,-x0);

    if(beta1> 3.1415927) beta1 -= 6.283185;
    if(beta1<-3.1415927) beta1 += 6.283185;
    if(beta2> 3.1415927) beta2 -= 6.283185;
    if(beta2<-3.1415927) beta2 += 6.283185;
    if(beta3> 3.1415927) beta3 -= 6.283185;
    if(beta3<-3.1415927) beta3 += 6.283185;

    double t12 = (z2-z1)/fabs(beta2-beta1)/R1;
    double z12 = (z1*beta2-z2*beta1)/(beta2-beta1);
    double t13 = (z3-z1)/fabs(beta3-beta1)/R1;
    double z13 = (z1*beta3-z3*beta1)/(beta3-beta1);

    // cout<<"::::: "<<sigmaz<<" "<<beta1<<"\t"<<beta2<<"\t"<<beta3<<"\n";
    // cout<<"::::: "<<t12<<"\t"<<t13<<"\n";
    // cout<<"::::: "<<z12<<"\t"<<z13<<"\n";
    
    if(sigmaz<0){
      //in L456, take 13 (large lever arm)
      t = t13;
      z0 = z13;
    }
    else{
      //in L234, take 12 (pixel layers)
      t = t12;
      z0 = z12;
    }
    
    for (int i=0;i<4;i++) {
      exactproj(rproj_[i],rinv,phi0,d0,t,z0,sqrt(x0*x0+y0*y0),
		phiproj[i],zproj[i],phider[i],zder[i]);
    }    

    for (int i=0;i<5;i++) {
      int sign=1;
      if (t<0) sign=-1;
      exactprojdisk(sign*zmean[i],rinv,phi0,t,z0,
		phiprojdisk[i],rprojdisk[i],phiderdisk[i],rderdisk[i]);
    }
  }

  void exacttrackletdisk(double r1, double z1, double phi1,
			 double r2, double z2, double phi2, double sigmaz,
			 double& rinv, double& phi0, double& d0,
			 double& t, double& z0,
			 double phiprojLayer[3], double zprojLayer[3], 
			 double phiderLayer[3], double zderLayer[3],
			 double phiproj[3], double rproj[3], 
			 double phider[3], double rder[3]) {

    double deltaphi=phi1-phi2;

    if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
    if (deltaphi<-0.5*two_pi) deltaphi+=two_pi;
    assert(fabs(deltaphi)<0.5*two_pi);

    if (sigmaz<-10.0) {
      cout << "Negative sigmaz"<<endl;
    }

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    double phi1tmp=phi1-phimin_+(phimax_-phimin_)/6.0;    

    //cout << "phi1 phi2 phi1tmp : "<<phi1<<" "<<phi2<<" "<<phi1tmp<<endl;

    phi0=phi1tmp+asin(0.5*r1*rinv);
    
    if (phi0>0.5*two_pi) phi0-=two_pi;
    if (phi0<-0.5*two_pi) phi0+=two_pi;
    if (!(fabs(phi0)<0.5*two_pi)) {
      cout << "phi1tmp r1 rinv phi0 deltaphi dist: "
	   <<phi1tmp<<" "<<r1<<" "<<rinv<<" "<<phi0
	   <<" "<<deltaphi<<" "<<dist<<endl;
      exit(1);
    }
    
    double rhopsi1=2*asin(0.5*r1*rinv)/rinv;
	    
    double rhopsi2=2*asin(0.5*r2*rinv)/rinv;
    
    t=(z1-z2)/(rhopsi1-rhopsi2);
    
    z0=z1-t*rhopsi1;

    if (disk_==1) {
      if (dumppars) {
	cout << "------------------------------------------------"<<endl;
	cout << "DUMPPARS0:" 
	     <<" dz= "<<z2-z1
	     <<" rinv= "<<rinv
	     <<" phi0= "<<phi0
	     <<" t= "<<t
	     <<" z0= "<<z0
	     <<endl;
      }
    }


    for (int i=0;i<3;i++) {
      exactprojdisk(zproj_[i],rinv,phi0,t,z0,
		    phiproj[i],rproj[i],
		    phider[i],rder[i]);
    }


    for (int i=0;i<3;i++) {
      exactproj(rmean[i],rinv,phi0,t,z0,
		    phiprojLayer[i],zprojLayer[i],
		    phiderLayer[i],zderLayer[i]);
    }

    d0 = 0.0;
  }

  void exacttrackletOverlap(double r1, double z1, double phi1,
			    double r2, double z2, double phi2, double sigmaz,
			    double& rinv, double& phi0, double& d0,
			    double& t, double& z0,
			    double phiprojLayer[3], double zprojLayer[3], 
			    double phiderLayer[3], double zderLayer[3],
			    double phiproj[3], double rproj[3], 
			    double phider[3], double rder[3]) {

    double deltaphi=phi1-phi2;

    if (deltaphi>0.5*two_pi) deltaphi-=two_pi;
    if (deltaphi<-0.5*two_pi) deltaphi+=two_pi;
    assert(fabs(deltaphi)<0.5*two_pi);

    if (sigmaz<-10.0) {
      cout << "Negative sigmaz"<<endl;
    }

    double dist=sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
    
    rinv=2*sin(deltaphi)/dist;

    if (r1>r2) rinv=-rinv;

    double phi1tmp=phi1-phimin_+(phimax_-phimin_)/6.0;    

    //cout << "phi1 phi2 phi1tmp : "<<phi1<<" "<<phi2<<" "<<phi1tmp<<endl;

    phi0=phi1tmp+asin(0.5*r1*rinv);
    
    if (phi0>0.5*two_pi) phi0-=two_pi;
    if (phi0<-0.5*two_pi) phi0+=two_pi;
    if (!(fabs(phi0)<0.5*two_pi)) {
      cout << "phi1tmp r1 rinv phi0 deltaphi dist: "
	   <<phi1tmp<<" "<<r1<<" "<<rinv<<" "<<phi0
	   <<" "<<deltaphi<<" "<<dist<<endl;
      exit(1);
    }
    
    double rhopsi1=2*asin(0.5*r1*rinv)/rinv;
	    
    double rhopsi2=2*asin(0.5*r2*rinv)/rinv;
    
    t=(z1-z2)/(rhopsi1-rhopsi2);
    
    z0=z1-t*rhopsi1;


    if (disk_==1) {
      if (dumppars) {
	cout << "------------------------------------------------"<<endl;
	cout << "DUMPPARS0:" 
	     <<" dz= "<<z2-z1
	     <<" rinv= "<<rinv
	     <<" phi0= "<<phi0
	     <<" t= "<<t
	     <<" z0= "<<z0
	     <<endl;
      }
    }


    for (int i=0;i<4;i++) {
      exactprojdisk(zprojoverlap_[i],rinv,phi0,t,z0,
		    phiproj[i],rproj[i],
		    phider[i],rder[i]);
    }


    for (int i=0;i<1;i++) {
      exactproj(rmean[i],rinv,phi0,t,z0,
		    phiprojLayer[i],zprojLayer[i],
		    phiderLayer[i],zderLayer[i]);
    }

    d0 = 0.0;
  }


  void exactproj(double rproj,double rinv,double phi0,
		  double t, double z0,
		  double &phiproj, double &zproj,
		  double &phider, double &zder) {

    phiproj=phi0-asin(0.5*rproj*rinv);
    zproj=z0+(2*t/rinv)*asin(0.5*rproj*rinv);

    phider=-0.5*rinv/sqrt(1-pow(0.5*rproj*rinv,2));
    zder=t/sqrt(1-pow(0.5*rproj*rinv,2));

  }

  void exactproj(double rproj,double rinv, double phi0, double d0,
		 double t, double z0, double r0,
		 double &phiproj, double &zproj,
		 double &phider, double &zder) {

    double rho = 1/rinv;
    if(rho<0) r0 = -r0;
    phiproj=phi0-asin((rproj*rproj+r0*r0-rho*rho)/(2*rproj*r0));
    double beta = acos((rho*rho+r0*r0-rproj*rproj)/(2*r0*rho));
    zproj=z0+ t*rho*beta;

    //not exact, but close
    phider=-0.5*rinv/sqrt(1-pow(0.5*rproj*rinv,2))-d0/(rproj*rproj);
    zder=t/sqrt(1-pow(0.5*rproj*rinv,2));

  }

  void exactprojdisk(double zproj,double rinv,double phi0,
		     double t, double z0,
		     double &phiproj, double &rproj,
		     double &phider, double &rder) {

    if (t<0) zproj=-zproj;
    
    double tmp=rinv*(zproj-z0)/(2.0*t);
    rproj=(2.0/rinv)*sin(tmp);
    phiproj=phi0-tmp;


    //if (fabs(1.0/rinv)>180.0&&fabs(z0)<15.0) {
    //  cout << "phiproj phi0 tmp zproj z0 t: "<<phiproj<<" "<<phi0
    //	   <<" "<<tmp<<" "<<zproj<<" "<<z0<<" "<<t<<endl;
    //}

    if (dumpproj) {
      if (fabs(zproj+300.0)<10.0) {
	cout << "DUMPPROJDISK1: "
	       << " phi="<<phiproj
	       << " r="<<rproj
	       << endl;
	}
      }


    phider=-rinv/(2*t);
    rder=cos(tmp)/t;

    //assert(fabs(phider)<0.1);

  } 

  void execute() {

    unsigned int countall=0;
    unsigned int countsel=0;

    //cout << "FPGATrackletCalculatorDisplaced execute "<<getName()<<" "<<stubtriplets_.size()<<endl;
    
    for(unsigned int l=0;l<stubtriplets_.size();l++){
      if (trackletpars_->nTracklets()>=maxtracklet_) {
	cout << "Will break on too many tracklets in "<<getName()<<endl;
	break;
      }
      for(unsigned int i=0;i<stubtriplets_[l]->nStubTriplets();i++){

	//if(stubtriplets_.size()>0)
	//  cout << "FPGATrackletCalculatorDisplaced execute "<<getName()<<" "<<stubtriplets_[l]->getName()<<" "<<stubtriplets_[l]->nStubTriplets()<<" "<<layer_<<endl;
	
	countall++;

	L1TStub* innerStub=stubtriplets_[l]->getL1TStub1(i);
	FPGAStub* innerFPGAStub=stubtriplets_[l]->getFPGAStub1(i);

	L1TStub* middleStub=stubtriplets_[l]->getL1TStub2(i);
	FPGAStub* middleFPGAStub=stubtriplets_[l]->getFPGAStub2(i);

	L1TStub* outerStub=stubtriplets_[l]->getL1TStub3(i);
	FPGAStub* outerFPGAStub=stubtriplets_[l]->getFPGAStub3(i);

	if (debug1) {
	  cout << "FPGATrackletCalculatorDisplaced execute "<<getName()<<"["<<iSector_<<"]"<<endl;
	}
	
        if (innerFPGAStub->isBarrel()&&middleFPGAStub->isBarrel()&&outerFPGAStub->isBarrel()){
	    //barrel+barrel seeding	  
	    bool accept = barrelSeeding(innerFPGAStub,innerStub,middleFPGAStub,middleStub,outerFPGAStub,outerStub);
	    
	    if (accept) countsel++;
        }
        else if (innerFPGAStub->isDisk()&&middleFPGAStub->isDisk()&&outerFPGAStub->isDisk()){
            assert(0);
        }
        else{
	    //layer+disk seeding
	    
            if (innerFPGAStub->isBarrel() && middleFPGAStub->isDisk() && outerFPGAStub->isDisk()){ //D1D2L2
                bool accept = overlapSeeding(middleFPGAStub,middleStub,outerFPGAStub,outerStub,innerFPGAStub,innerStub);

                if (accept) countsel++;
            }
            else if (innerFPGAStub->isDisk() && middleFPGAStub->isBarrel() && outerFPGAStub->isBarrel()){ //L2L3D1
                bool accept = overlapSeeding(innerFPGAStub,innerStub,middleFPGAStub,middleStub,outerFPGAStub,outerStub);

                if (accept) countsel++;
            }
            else{
                assert(0);
            }
        }

	if (trackletpars_->nTracklets()>=maxtracklet_) {
	  cout << "Will break on number of tracklets in "<<getName()<<endl;
	  break;
	}
	
	if (countall>=MAXTC) {
	  if (debug1) cout << "Will break on MAXTC 1"<<endl;
	  break;
	}
	if (debug1) {
	  cout << "FPGATrackletCalculatorDisplaced execute done"<<endl;
	}

      }
      if (countall>=MAXTC) {
	if (debug1) cout << "Will break on MAXTC 2"<<endl;
	break;
      }
    }

    if (writeTrackletCalculatorDisplaced) {
      static ofstream out("trackletcalculatordisplaced.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }


  }


  void addDiskProj(FPGATracklet* tracklet, int disk){

    
    FPGAWord fpgar=tracklet->fpgarprojdisk(disk);

    if (fpgar.value()*krprojshiftdisk<12.0) return;
    if (fpgar.value()*krprojshiftdisk>112.0) return;

    //cout << "addDiskProj neighbor: "<<tracklet->plusNeighborDisk(disk)<<" "<<tracklet->minusNeighborDisk(disk)<<endl;

    if (tracklet->plusNeighborDisk(disk)) {
      if (getName().find("L1L2")!=std::string::npos) {
	if (abs(disk)==1) addNeighborProjectionDisk(disk,trackletproj_L6Plus_,tracklet);
	if (abs(disk)==2) addNeighborProjectionDisk(disk,trackletproj_L5Plus_,tracklet);
	if (abs(disk)==3) addNeighborProjectionDisk(disk,trackletproj_L4Plus_,tracklet);
	if (abs(disk)==4) addNeighborProjectionDisk(disk,trackletproj_L3Plus_,tracklet);
	return;
      }
      if (getName().find("L3L4")!=std::string::npos) {
	if (abs(disk)==1) addNeighborProjectionDisk(disk,trackletproj_L6Plus_,tracklet);
	if (abs(disk)==2) addNeighborProjectionDisk(disk,trackletproj_L5Plus_,tracklet);
	return;
      }
      if (abs(disk)==1) addNeighborProjectionDisk(disk,trackletproj_D1Plus_,tracklet);
      if (abs(disk)==2) addNeighborProjectionDisk(disk,trackletproj_D2Plus_,tracklet);
      if (abs(disk)==3) addNeighborProjectionDisk(disk,trackletproj_D3Plus_,tracklet);
      if (abs(disk)==4) addNeighborProjectionDisk(disk,trackletproj_D4Plus_,tracklet);
      if (abs(disk)==5) addNeighborProjectionDisk(disk,trackletproj_D5Plus_,tracklet);
      return;
    }
      
    if (tracklet->minusNeighborDisk(disk)) {
      if (getName().find("L1L2")!=std::string::npos) {
	if (abs(disk)==1) addNeighborProjectionDisk(disk,trackletproj_L6Minus_,tracklet);
	if (abs(disk)==2) addNeighborProjectionDisk(disk,trackletproj_L5Minus_,tracklet);
	if (abs(disk)==3) addNeighborProjectionDisk(disk,trackletproj_L4Minus_,tracklet);
	if (abs(disk)==4) addNeighborProjectionDisk(disk,trackletproj_L3Minus_,tracklet);
	return;
      }
      if (getName().find("L3L4")!=std::string::npos) {
	if (abs(disk)==1) addNeighborProjectionDisk(disk,trackletproj_L6Minus_,tracklet);
	if (abs(disk)==2) addNeighborProjectionDisk(disk,trackletproj_L5Minus_,tracklet);
	return;
      }

      if (abs(disk)==1) addNeighborProjectionDisk(disk,trackletproj_D1Minus_,tracklet);
      if (abs(disk)==2) addNeighborProjectionDisk(disk,trackletproj_D2Minus_,tracklet);
      if (abs(disk)==3) addNeighborProjectionDisk(disk,trackletproj_D3Minus_,tracklet);
      if (abs(disk)==4) addNeighborProjectionDisk(disk,trackletproj_D4Minus_,tracklet);
      if (abs(disk)==5) addNeighborProjectionDisk(disk,trackletproj_D5Minus_,tracklet);
      return;
    }
      
    
    FPGAWord fpgaphi=tracklet->fpgaphiprojdisk(disk);
    
    int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);

    int iphi=-1;


    if (hourglass) {

      iphi=iphivmRaw/(32/nallstubsdisks[abs(disk)-1]);
      
    } else {
    
      assert(iphivmRaw>=4);
      assert(iphivmRaw<=27);

      iphi=(iphivmRaw-4)>>3;

      assert(iphi>=0);
      assert(iphi<=2);

    }
    
    if (abs(disk)==1) {
      if (iphi==0) addProjectionDisk(disk,iphi,trackletproj_D1PHI1_,tracklet);
      if (iphi==1) addProjectionDisk(disk,iphi,trackletproj_D1PHI2_,tracklet);
      if (iphi==2) addProjectionDisk(disk,iphi,trackletproj_D1PHI3_,tracklet);
      if (iphi==3) addProjectionDisk(disk,iphi,trackletproj_D1PHI4_,tracklet);
    }
    
    if (abs(disk)==2) {
      if (iphi==0) addProjectionDisk(disk,iphi,trackletproj_D2PHI1_,tracklet);
      if (iphi==1) addProjectionDisk(disk,iphi,trackletproj_D2PHI2_,tracklet);
      if (iphi==2) addProjectionDisk(disk,iphi,trackletproj_D2PHI3_,tracklet);
      if (iphi==3) addProjectionDisk(disk,iphi,trackletproj_D2PHI4_,tracklet);
    }

    if (abs(disk)==3) {
      if (iphi==0) addProjectionDisk(disk,iphi,trackletproj_D3PHI1_,tracklet);
      if (iphi==1) addProjectionDisk(disk,iphi,trackletproj_D3PHI2_,tracklet);
      if (iphi==2) addProjectionDisk(disk,iphi,trackletproj_D3PHI3_,tracklet);
      if (iphi==3) addProjectionDisk(disk,iphi,trackletproj_D3PHI4_,tracklet);
    }

    if (abs(disk)==4) {
      if (iphi==0) addProjectionDisk(disk,iphi,trackletproj_D4PHI1_,tracklet);
      if (iphi==1) addProjectionDisk(disk,iphi,trackletproj_D4PHI2_,tracklet);
      if (iphi==2) addProjectionDisk(disk,iphi,trackletproj_D4PHI3_,tracklet);
      if (iphi==3) addProjectionDisk(disk,iphi,trackletproj_D4PHI4_,tracklet);
    }

    if (abs(disk)==5) {
      if (iphi==0) addProjectionDisk(disk,iphi,trackletproj_D5PHI1_,tracklet);
      if (iphi==1) addProjectionDisk(disk,iphi,trackletproj_D5PHI2_,tracklet);
      if (iphi==2) addProjectionDisk(disk,iphi,trackletproj_D5PHI3_,tracklet);
      if (iphi==3) addProjectionDisk(disk,iphi,trackletproj_D5PHI4_,tracklet);
    }

    
  }


  bool addLayerProj(FPGATracklet* tracklet, int layer){

    
    assert(layer>0);

    FPGAWord fpgaz=tracklet->fpgazproj(layer);
    FPGAWord fpgaphi=tracklet->fpgaphiproj(layer);


    if(fpgaphi.atExtreme()) cout<<"at extreme! "<<fpgaphi.value()<<"\n";

    assert(!fpgaphi.atExtreme());
    
    if (fpgaz.atExtreme()) return false;

    if (fabs(fpgaz.value()*kz)>zlength) return false;

    if (hourglass) {
      assert(!tracklet->plusNeighbor(layer));
      assert(!tracklet->minusNeighbor(layer));
    }
    
    if (tracklet->plusNeighbor(layer)) {
      if (layer==1) addNeighborProjection(layer,trackletproj_L1Plus_,tracklet);
      if (layer==2) addNeighborProjection(layer,trackletproj_L2Plus_,tracklet);
      if (layer==3) addNeighborProjection(layer,trackletproj_L3Plus_,tracklet);
      if (layer==4) addNeighborProjection(layer,trackletproj_L4Plus_,tracklet);
      if (layer==5) addNeighborProjection(layer,trackletproj_L5Plus_,tracklet);
      if (layer==6) addNeighborProjection(layer,trackletproj_L6Plus_,tracklet);
      return true;
    }
      
    if (tracklet->minusNeighbor(layer)) {
      if (layer==1) addNeighborProjection(layer,trackletproj_L1Minus_,tracklet);
      if (layer==2) addNeighborProjection(layer,trackletproj_L2Minus_,tracklet);
      if (layer==3) addNeighborProjection(layer,trackletproj_L3Minus_,tracklet);
      if (layer==4) addNeighborProjection(layer,trackletproj_L4Minus_,tracklet);
      if (layer==5) addNeighborProjection(layer,trackletproj_L5Minus_,tracklet);
      if (layer==6) addNeighborProjection(layer,trackletproj_L6Minus_,tracklet);
      return true;
    }
      


    int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);

    int iphi=-1;
    
    if (hourglass) {
      
      iphi=iphivmRaw/(32/nallstubslayers[layer-1]);
      
    } else {

      assert(iphivmRaw>=4);
      assert(iphivmRaw<=27);

      iphi=(iphivmRaw-4)>>3;

      if (layer==2||layer==4||layer==6) {
	iphi=(iphivmRaw>>3);
      }
      assert(iphi>=0);
      assert(iphi<=7);
      

    }
      
 
    //cout << "layer fpgaphi iphivmRaw iphi : "<<layer<<" "<<fpgaphi.value()<<" "<<iphivmRaw<<" "<<iphi<<endl;

    

    if (layer==1) {
      if (iphi==0) addProjection(layer,iphi,trackletproj_L1PHI1_,tracklet);
      if (iphi==1) addProjection(layer,iphi,trackletproj_L1PHI2_,tracklet);
      if (iphi==2) addProjection(layer,iphi,trackletproj_L1PHI3_,tracklet);
      if (iphi==3) addProjection(layer,iphi,trackletproj_L1PHI4_,tracklet);
      if (iphi==4) addProjection(layer,iphi,trackletproj_L1PHI5_,tracklet);
      if (iphi==5) addProjection(layer,iphi,trackletproj_L1PHI6_,tracklet);
      if (iphi==6) addProjection(layer,iphi,trackletproj_L1PHI7_,tracklet);
      if (iphi==7) addProjection(layer,iphi,trackletproj_L1PHI8_,tracklet);
    }
    
    if (layer==2) {
      if (iphi==0) addProjection(layer,iphi,trackletproj_L2PHI1_,tracklet);
      if (iphi==1) addProjection(layer,iphi,trackletproj_L2PHI2_,tracklet);
      if (iphi==2) addProjection(layer,iphi,trackletproj_L2PHI3_,tracklet);
      if (iphi==3) addProjection(layer,iphi,trackletproj_L2PHI4_,tracklet);
    }

    if (layer==3) {
      if (iphi==0) addProjection(layer,iphi,trackletproj_L3PHI1_,tracklet);
      if (iphi==1) addProjection(layer,iphi,trackletproj_L3PHI2_,tracklet);
      if (iphi==2) addProjection(layer,iphi,trackletproj_L3PHI3_,tracklet);
      if (iphi==3) addProjection(layer,iphi,trackletproj_L3PHI4_,tracklet);
    }

    if (layer==4) {
      if (iphi==0) addProjection(layer,iphi,trackletproj_L4PHI1_,tracklet);
      if (iphi==1) addProjection(layer,iphi,trackletproj_L4PHI2_,tracklet);
      if (iphi==2) addProjection(layer,iphi,trackletproj_L4PHI3_,tracklet);
      if (iphi==3) addProjection(layer,iphi,trackletproj_L4PHI4_,tracklet);
    }

    if (layer==5) {
      if (iphi==0) addProjection(layer,iphi,trackletproj_L5PHI1_,tracklet);
      if (iphi==1) addProjection(layer,iphi,trackletproj_L5PHI2_,tracklet);
      if (iphi==2) addProjection(layer,iphi,trackletproj_L5PHI3_,tracklet);
      if (iphi==3) addProjection(layer,iphi,trackletproj_L5PHI4_,tracklet);
    }

    if (layer==6) {
      if (iphi==0) addProjection(layer,iphi,trackletproj_L6PHI1_,tracklet);
      if (iphi==1) addProjection(layer,iphi,trackletproj_L6PHI2_,tracklet);
      if (iphi==2) addProjection(layer,iphi,trackletproj_L6PHI3_,tracklet);
      if (iphi==3) addProjection(layer,iphi,trackletproj_L6PHI4_,tracklet);
    }

    return true;

  }

  void addProjection(int layer,int iphi,FPGATrackletProjections* trackletprojs, FPGATracklet* tracklet){
    if (trackletprojs==0) {
      if (warnNoMem) {
	cout << "No projection memory exists in "<<getName()<<" for layer = "<<layer<<" iphi = "<<iphi+1<<endl;
      }
      return;
    }
    assert(trackletprojs!=0);
    trackletprojs->addProj(tracklet);
  }

  void addProjectionDisk(int disk,int iphi,FPGATrackletProjections* trackletprojs, FPGATracklet* tracklet){
    if (trackletprojs==0) {
      if (layer_==3&&abs(disk)==3) return; //L3L4 projections to D3 are not used.
      if (warnNoMem) {       
	cout << "No projection memory exists in "<<getName()<<" for disk = "<<abs(disk)<<" iphi = "<<iphi+1<<endl;
      }
      return;
    }
    assert(trackletprojs!=0);
    trackletprojs->addProj(tracklet);
  }

  void addNeighborProjection(int layer, FPGATrackletProjections* trackletprojs,FPGATracklet* tracklet){
    if (trackletprojs==0) {
      if (warnNoMem) {
	string str="";
	if (tracklet->minusNeighbor(layer)){
	  str="Minus";
	}
	if (tracklet->plusNeighbor(layer)){
	  str="Plus";
	}
	assert(str!="");
	cout << "Error no projection memory exists in "<<getName()<<" for layer = "<<layer<<" to "
	     <<str<<" neighbor"<<endl;
      }
      return;
    }
    assert(trackletprojs!=0);
    trackletprojs->addProj(tracklet);
  

  }

  
  void addNeighborProjectionDisk(int disk, FPGATrackletProjections* trackletprojs,FPGATracklet* tracklet){

    if (trackletprojs==0) {
      if (warnNoMem) {
	string str="";
	if (tracklet->minusNeighborDisk(disk)){
	  str="Minus";
	}
	if (tracklet->plusNeighborDisk(disk)){
	  str="Plus";
	}
	assert(str!="");
	cout << "Error no projection memory exists in "<<getName()<<" for disk = "<<abs(disk)<<" to "
	     <<str<<" neighbor"<<endl;
      }
      return;
    }
    assert(trackletprojs!=0);

    trackletprojs->addProj(tracklet);

    
  }


  bool barrelSeeding(FPGAStub* innerFPGAStub, L1TStub* innerStub, FPGAStub* middleFPGAStub, L1TStub* middleStub, FPGAStub* outerFPGAStub, L1TStub* outerStub){
	  
    if (debug1) {
      cout << "FPGATrackletCalculatorDisplaced "<<getName()<<" "<<layer_<<" trying stub pair in layer (inner outer): "
	   <<innerFPGAStub->layer().value()<<" "<<outerFPGAStub->layer().value()<<endl;
    }
	    
    assert(outerFPGAStub->isBarrel());
    
    //assert(layer_==innerFPGAStub->layer().value()+1);
    
    //assert(layer_==1||layer_==3||layer_==5);

    	  
    double r1=innerStub->r();
    double z1=innerStub->z();
    double phi1=innerStub->phi();
    
    double r2=middleStub->r();
    double z2=middleStub->z();
    double phi2=middleStub->phi();

    double r3=outerStub->r();
    double z3=outerStub->z();
    double phi3=outerStub->phi();
    
    
    double rinv,phi0,d0,t,z0;
    
    double phiproj[4],zproj[4],phider[4],zder[4];
    double phiprojdisk[5],rprojdisk[5],phiderdisk[5],rderdisk[5];
    
    exacttracklet(r1,z1,phi1,r2,z2,phi2,outerStub->sigmaz(),
                  r3,z3,phi3,
		  rinv,phi0,d0,t,z0,
		  phiproj,zproj,phider,zder,
		  phiprojdisk,rprojdisk,phiderdisk,rderdisk);

    if (useapprox) {
      phi1=innerFPGAStub->phiapprox(phimin_,phimax_);
      z1=innerFPGAStub->zapprox();
      r1=innerFPGAStub->rapprox();

      phi2=outerFPGAStub->phiapprox(phimin_,phimax_);
      z2=outerFPGAStub->zapprox();
      r2=outerFPGAStub->rapprox();
    }
    
    double rinvapprox,phi0approx,d0approx,tapprox,z0approx;
    double phiprojapprox[4],zprojapprox[4],phiderapprox[4],zderapprox[4];
    double phiprojdiskapprox[5],rprojdiskapprox[5];
    double phiderdiskapprox[5],rderdiskapprox[5];

    //FIXME: do the actual integer calculation
    
    phi0 -= 0.171;

    //store the approcximate results
    rinvapprox = rinv;
    phi0approx = phi0;
    d0approx   = d0;
    tapprox    = t;
    z0approx   = z0;
    
    for(int i=0; i<4; ++i){
      phiproj[i] -= 0.171;
      phiprojapprox[i] = phiproj[i];
      zprojapprox[i]   = zproj[i];
      phiderapprox[i]  = phider[i];
      zderapprox[i]    = zder[i];
    }

    for(int i=0; i<5; ++i){
      phiprojdisk[i] -= 0.171;
      phiprojdiskapprox[i] = phiprojdisk[i];
      rprojdiskapprox[i]   = rprojdisk[i];
      phiderdiskapprox[i]  = phiderdisk[i];
      rderdiskapprox[i]    = rderdisk[i];
    }

    //now binary
    double krinv = kphi1/kr*pow(2,rinv_shift),
           kphi0 = kphi1*pow(2,phi0_shift),
           kt = kz/kr*pow(2,t_shift),
           kz0 = kz*pow(2,z0_shift),
           kphiproj = kphi1*pow(2,SS_phiL_shift),
           kphider = kphi1/kr*pow(2,SS_phiderL_shift),
           kzproj = kz*pow(2,PS_zL_shift),
           kzder = kz/kr*pow(2,PS_zderL_shift),
           kphiprojdisk = kphi1*pow(2,SS_phiD_shift),
           kphiderdisk = kphi1/kr*pow(2,SS_phiderD_shift),
           krprojdisk = kr*pow(2,PS_rD_shift),
           krderdisk = kr/kz*pow(2,PS_rderD_shift);
    
    int irinv,iphi0,id0,it,iz0;
    bool validproj[4];
    int iphiproj[4],izproj[4],iphider[4],izder[4];
    bool minusNeighbor[4],plusNeighbor[4];
    bool validprojdisk[5];
    int iphiprojdisk[5],irprojdisk[5],iphiderdisk[5],irderdisk[5];
    bool minusNeighborDisk[5],plusNeighborDisk[5];
      
    //store the binary results
    irinv = rinvapprox / krinv;
    iphi0 = phi0approx / kphi0;
    id0 = d0approx / kd0;
    it    = tapprox / kt;
    iz0   = z0approx / kz0;

    for(int i=0; i<4; ++i){
      iphiproj[i] = phiprojapprox[i] / kphiproj;
      izproj[i]   = zprojapprox[i] / kzproj;

      iphider[i] = phiderapprox[i] / kphider;
      izder[i]   = zderapprox[i] / kzder;

      validproj[i] = true;
      if (izproj[i]<-(1<<(nbitszprojL123-1))) validproj[i]=false;
      if (izproj[i]>=(1<<(nbitszprojL123-1))) validproj[i]=false;
      
      minusNeighbor[i]=false;
      plusNeighbor[i]=false;
      if (!hourglass) {
	if (iphiproj[i]<(1<<nbitsphistubL456)/8) {
	  minusNeighbor[i]=true;
	  iphiproj[i]+=3*(1<<nbitsphistubL456)/4;
	}
	if (iphiproj[i]>=7*(1<<nbitsphistubL456)/8) {
	  plusNeighbor[i]=true;
	  iphiproj[i]-=3*(1<<nbitsphistubL456)/4;
	}
      }

      //this is left from the original....
      if (iphiproj[i]>=(1<<nbitsphistubL456)-1) {
	iphiproj[i]=(1<<nbitsphistubL456)-2; //-2 not to hit atExtreme
	validproj[i] = false;
      }
      
      if (rproj_[i]<60.0) {
	iphiproj[i]>>=(nbitsphistubL456-nbitsphistubL123);
	if (iphiproj[i]>=(1<<nbitsphistubL123)-1) iphiproj[i]=(1<<nbitsphistubL123)-2; //-2 not to hit atExtreme
      }
      else {
	izproj[i]>>=(nbitszprojL123-nbitszprojL456);
      }

      if (iphiproj[i]<=0) {
	iphiproj[i]=1;
	validproj[i] = false;
      }

      
    }

    if(fabs(it * kt)<1.0) {
      for(int i=0; i<5; ++i) {
	//do not bother with central tracks; the calculation there is wrong anyway.
	validprojdisk[i]=false;
      }
    } else {
      for(int i=0; i<5; ++i){
	validprojdisk[i]=true;

        iphiprojdisk[i] = phiprojdiskapprox[i] / kphiprojdisk;
        irprojdisk[i]   = rprojdiskapprox[i] / krprojdisk;

	iphiderdisk[i] = phiderdiskapprox[i] / kphiderdisk;
	irderdisk[i]   = rderdiskapprox[i] / krderdisk;
      
	minusNeighborDisk[i]=false;
	plusNeighborDisk[i]=false;
	if (!hourglass) {
	  if (iphiprojdisk[i]<(1<<nbitsphistubL123)/8) {
	    minusNeighborDisk[i]=true;
	    iphiprojdisk[i]+=3*(1<<nbitsphistubL123)/4;
	  }
	  if (iphiprojdisk[i]>=7*(1<<nbitsphistubL123)/8) {
	    plusNeighborDisk[i]=true;
	    iphiprojdisk[i]-=3*(1<<nbitsphistubL123)/4;
	  }
	}
     
	if (iphiprojdisk[i]<0) {
	  iphiprojdisk[i]=0;
	  validprojdisk[i]=false;
	}
	if (iphiprojdisk[i]>=(1<<nbitsphistubL123)) {
	  iphiprojdisk[i]=(1<<nbitsphistubL123)-1;
	  validprojdisk[i]=false;
	}
      
	//"protection" from the original
	if (iphiprojdisk[i]<0) {
	  iphiprojdisk[i]=0;
	  validprojdisk[i]=false;
	}
	if (iphiprojdisk[i]>=(1<<nbitsphistubL123)) {
	  iphiprojdisk[i]=(1<<nbitsphistubL123)-1;
	  validprojdisk[i]=false;
	}
	
	if(irprojdisk[i]< 20. / krprojdisk || irprojdisk[i] > 120. / krprojdisk ){
	  validprojdisk[i]=false;
	  irprojdisk[i] = 0;
	  iphiprojdisk[i] = 0;
	  iphiderdisk[i]  = 0;
	  irderdisk[i]    = 0;
	}
      }
    }

    bool success = true;
    if(fabs(rinvapprox)>rinvcut){
      if (debug1) 
	cout << "FPGATrackletCalculator::BarrelSeeding irinv too large: "
	     <<rinvapprox<<"("<<irinv<<")\n";
      success = false;
    }
    if (layer_==1&&fabs(z0approx)>z0cut) {
      if (debug1) cout << "Failed tracklet z0 cut "<<z0approx<<" in layer 1"<<endl;
      success = false;
    }
    if (layer_>=2&&fabs(z0approx)>1.5*z0cut) { 
      if (debug1) cout << "Failed tracklet z0 cut "<<z0approx<<" in layer "<<layer_<<endl;
      success = false;
    }
    if (fabs(d0approx)>maxd0) {
      if (debug1) cout << "Failed tracklet d0 cut "<<d0approx<<endl;
      success = false;
    }
    
    if (!success) return false;

    if (hourglass) {
      double phicrit=phi0approx-asin(0.5*rcrit*rinvapprox);
      bool keep=(phicrit>phicritminmc)&&(phicrit<phicritmaxmc);
      if (!keep) return false;
    }
    
    for(unsigned int j=0;j<5;j++){
      if (minusNeighborDisk[j]) {
	phiprojdiskapprox[j]+=dphisector;
	phiprojdisk[j]+=dphisector;
      }
      if (plusNeighborDisk[j]) {
	phiprojdiskapprox[j]-=dphisector;
	phiprojdisk[j]-=dphisector;
      }
    }
	  
    for(unsigned int j=0;j<4;j++){
      if (minusNeighbor[j]) {
	phiprojapprox[j]+=dphisector;
	phiproj[j]+=dphisector;
      }
      if (plusNeighbor[j]) {
	phiprojapprox[j]-=dphisector;
	phiproj[j]-=dphisector;
      }	    
    }
    
    
    if (writeTrackletPars) {
      static ofstream out("trackletpars.txt");
      out <<"Trackpars "<<layer_
	  <<"   "<<rinv<<" "<<rinvapprox<<" "<<rinvapprox
	  <<"   "<<phi0<<" "<<phi0approx<<" "<<phi0approx
	  <<"   "<<t<<" "<<tapprox<<" "<<tapprox
	  <<"   "<<z0<<" "<<z0approx<<" "<<z0approx
	  <<endl;
    }	        
        
    FPGATracklet* tracklet=new FPGATracklet(innerStub,middleStub,outerStub,
					    innerFPGAStub,middleFPGAStub,outerFPGAStub,
					    iSector_,
					    phioffset_,
					    rinv,phi0,d0,z0,t,
					    rinvapprox,phi0approx,d0approx,
					    z0approx,tapprox,
					    irinv,iphi0,id0,iz0,it,validproj,
					    iphiproj,izproj,iphider,izder,
					    minusNeighbor,plusNeighbor,
					    phiproj,zproj,phider,zder,
					    phiprojapprox,zprojapprox,
					    phiderapprox,zderapprox,
					    validprojdisk,
					    iphiprojdisk,irprojdisk,
					    iphiderdisk,irderdisk,
					    minusNeighborDisk,
					    plusNeighborDisk,
					    phiprojdisk,rprojdisk,
					    phiderdisk,rderdisk,
					    phiprojdiskapprox,
					    rprojdiskapprox,
					    phiderdiskapprox,
					    rderdiskapprox,
					    false);
    
    if (debug1) {
      cout << "FPGATrackletCalculatorDisplaced "<<getName()<<" Found tracklet in layer = "<<layer_<<" "
	   <<iSector_<<" phi0 = "<<phi0<<endl;
    }
        

    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);

    trackletpars_->addTracklet(tracklet);
    ofstream fout ("seedsTC.txt", ofstream::app);
    fout << tracklet->eta () << " " << tracklet->getISeed () << endl;
    fout.close ();
    
    bool addL3=false;
    bool addL4=false;
    bool addL5=false;
    bool addL6=false;
    for(unsigned int j=0;j<4;j++){
      //	    cout<<" LL to L "<<lproj_[j]<<"\n";
      bool added=false;
      if (tracklet->validProj(lproj_[j])) {
	added=addLayerProj(tracklet,lproj_[j]);
	//cout << "Add tracklet proj for layer "<<lproj_[j]<<": "<<phiproj[j]<<" "<<iphiproj[j]<<" added = "
	//     <<added<<endl;
	if (added&&lproj_[j]==3) addL3=true;
	if (added&&lproj_[j]==4) addL4=true;
	if (added&&lproj_[j]==5) addL5=true;
	if (added&&lproj_[j]==6) addL6=true;
      }
    }
    
    
    for(unsigned int j=0;j<4;j++){ //no projections to 5th disk!!
      int disk=j+1;
      if (disk==4&&addL3) continue;
      if (disk==3&&addL4) continue;
      if (disk==2&&addL5) continue;
      if (disk==1&&addL6) continue;
      if (it<0) disk=-disk;
      //	    cout<<" LL to disk "<<disk<<"\n";
      if (tracklet->validProjDisk(abs(disk))) {
	//cout << "Add tracklet "<<tracklet<<" for disk "<<disk<<endl;
	addDiskProj(tracklet,disk);
      }
    }
    
    return true;

  }
    

  bool diskSeeding(FPGAStub* innerFPGAStub,L1TStub* innerStub,FPGAStub* middleFPGAStub,L1TStub* middleStub,FPGAStub* outerFPGAStub,L1TStub* outerStub){

	    
    if (debug1) {
      cout <<  "FPGATrackletCalculatorDisplaced::execute calculate disk seeds" << endl;
    }
	      
    int sign=1;
    if (innerFPGAStub->disk().value()<0) sign=-1;
    
    disk_=innerFPGAStub->disk().value();
    assert(abs(disk_)==1||abs(disk_)==3);
    
    
    assert(innerStub->isPSmodule());
    assert(outerStub->isPSmodule());
	    
    double r1=innerStub->r();
    double z1=innerStub->z();
    double phi1=innerStub->phi();
    
    double r2=middleStub->r();
    //double z2=middleStub->z();
    //double phi2=middleStub->phi();

    double r3=outerStub->r();
    double z3=outerStub->z();
    double phi3=outerStub->phi();
	    
    
    if (r2<r1+2.0) {
      //assert(0);
      return false; //Protection... Should be handled cleaner
      //to avoid problem with floating point 
      //calculation
    }
    
    double rinv,phi0,d0,t,z0;
    
    double phiproj[3],zproj[3],phider[3],zder[3];
    double phiprojdisk[3],rprojdisk[3],phiderdisk[3],rderdisk[3];
    
    exacttrackletdisk(r1,z1,phi1,r3,z3,phi3,outerStub->sigmaz(),
		      rinv,phi0,d0,t,z0,
		      phiproj,zproj,phider,zder,
		      phiprojdisk,rprojdisk,phiderdisk,rderdisk);


    //Truncates floating point positions to integer
    //representation precision
    if (useapprox) {
      phi1=innerFPGAStub->phiapprox(phimin_,phimax_);
      z1=innerFPGAStub->zapprox();
      r1=innerFPGAStub->rapprox();
      
      //phi2=outerFPGAStub->phiapprox(phimin_,phimax_);
      //z2=outerFPGAStub->zapprox();
      r2=outerFPGAStub->rapprox();
    }
    
    double rinvapprox,phi0approx,d0approx,tapprox,z0approx;
    double phiprojapprox[3],zprojapprox[3],phiderapprox[3],zderapprox[3];
    double phiprojdiskapprox[3],rprojdiskapprox[3],
      phiderdiskapprox[3],rderdiskapprox[3];
	    
    //store the approximate results
    rinvapprox = rinv;
    phi0approx = phi0;
    d0approx = d0;
    tapprox    = t;
    z0approx   = z0;


    for(int i=0; i<3; ++i){
      phiprojapprox[i] = phiproj[i];
      zprojapprox[i]   = zproj[i];
      phiderapprox[i] = phider[i];
      zderapprox[i]   = zder[i];
    }


    for(int i=0; i<3; ++i){
      phiprojdiskapprox[i] = phiprojdisk[i];
      rprojdiskapprox[i] = rprojdisk[i];
      phiderdiskapprox[i] = phiderdisk[i];
      rderdiskapprox[i]   = rderdisk[i];
    }

    //now binary
    double krinv = kphi1/kr*pow(2,rinv_shift),
           kphi0 = kphi1*pow(2,phi0_shift),
           kt = kz/kr*pow(2,t_shift),
           kz0 = kz*pow(2,z0_shift),
           kphiproj = kphi1*pow(2,SS_phiL_shift),
           kphider = kphi1/kr*pow(2,SS_phiderL_shift),
           kzproj = kz*pow(2,PS_zL_shift),
           kzder = kz/kr*pow(2,PS_zderL_shift),
           kphiprojdisk = kphi1*pow(2,SS_phiD_shift),
           kphiderdisk = kphi1/kr*pow(2,SS_phiderD_shift),
           krprojdisk = kr*pow(2,PS_rD_shift),
           krderdisk = kr/kz*pow(2,PS_rderD_shift);
    
    int irinv,iphi0,id0,it,iz0;
    bool validproj[3];
    int iphiproj[3],izproj[3],iphider[3],izder[3];
    bool minusNeighbor[3],plusNeighbor[3];
    
    bool validprojdisk[3];
    int iphiprojdisk[3],irprojdisk[3],iphiderdisk[3],irderdisk[3];
    bool minusNeighborDisk[3],plusNeighborDisk[3];
    
    //store the binary results
    irinv = rinvapprox / krinv;
    iphi0 = phi0approx / kphi0;
    id0 = d0approx / kd0;
    it    = tapprox / kt;
    iz0   = z0approx / kz0;

    for(int i=0; i<3; ++i){
      iphiproj[i] = phiprojapprox[i] / kphiproj;
      izproj[0]   = zprojapprox[i] / kzproj;

      iphider[i] = phiderapprox[i] / kphider;
      izder[i]   = zderapprox[i] / kzder;

      validproj[i] = true;
      if (izproj[i]<-(1<<(nbitszprojL123-1))) validproj[i]=false;
      if (izproj[i]>=(1<<(nbitszprojL123-1))) validproj[i]=false;
      
      minusNeighbor[i]=false;
      plusNeighbor[i]=false;
      if (!hourglass) {
	if (iphiproj[i]<(1<<nbitsphistubL456)/8) {
	  minusNeighbor[i]=true;
	  iphiproj[i]+=3*(1<<nbitsphistubL456)/4;
	}
	if (iphiproj[i]>=7*(1<<nbitsphistubL456)/8) {
	  plusNeighbor[i]=true;
	  iphiproj[i]-=3*(1<<nbitsphistubL456)/4;
	}
      }

      //this is left from the original....
      if (iphiproj[i]>=(1<<nbitsphistubL456)-1){
	iphiproj[i]=(1<<nbitsphistubL456)-2; //-2 not to hit atExtreme
	validproj[i] = false;
      }

      if (rproj_[i]<60.0)
	iphiproj[i]>>=(nbitsphistubL456-nbitsphistubL123);
      else {
	izproj[i]>>=(nbitszprojL123-nbitszprojL456);
      }
      
      if (iphiproj[i]<=0) {
	iphiproj[i]=1;
	validproj[i] = false;
      }

    }


    for(int i=0; i<3; ++i){
      iphiprojdisk[i] = phiprojdiskapprox[i] / kphiprojdisk;
      irprojdisk[i]   = rprojdiskapprox[i] / krprojdisk;

      iphiderdisk[i] = phiderdiskapprox[i] / kphiderdisk;
      irderdisk[i]   = rderdiskapprox[i] / krderdisk;
      
      minusNeighborDisk[i]=false;
      plusNeighborDisk[i]=false;
      if (!hourglass) {
	if (iphiprojdisk[i]<(1<<nbitsphistubL123)/8) {
	  minusNeighborDisk[i]=true;
	  iphiprojdisk[i]+=3*(1<<nbitsphistubL123)/4;
	}
	if (iphiprojdisk[i]>=7*(1<<nbitsphistubL123)/8) {
	  plusNeighborDisk[i]=true;
	  iphiprojdisk[i]-=3*(1<<nbitsphistubL123)/4;
	}
      }
      
      //"protection" from the original
      if (iphiprojdisk[i]<0) iphiprojdisk[i]=0;
      if (iphiprojdisk[i]>=(1<<nbitsphistubL123)) iphiprojdisk[i]=(1<<nbitsphistubL123)-1;
      
      validprojdisk[i]=true;
      if(irprojdisk[i]<=0 || irprojdisk[i] > 120. / krprojdisk ){
	validprojdisk[i]=false;
	irprojdisk[i] = 0;
	iphiprojdisk[i] = 0;
	iphiderdisk[i]  = 0;
	irderdisk[i]    = 0;
      }
    }

    bool success = true;
    if(fabs(rinvapprox)>rinvcut){
      if (debug1) 
	cout << "FPGATrackletCalculator::DiskSeeding irinv too large: "<<rinvapprox<<endl;
      success = false;
    }
    if (fabs(z0approx)>z0cut) {
      if (debug1) cout << "Failed tracklet z0 cut "<<z0approx<<" in layer 1"<<endl;
      success = false;
    }
    if (fabs(d0approx)>maxd0) {
      if (debug1) cout << "Failed tracklet d0 cut "<<d0approx<<endl;
      success = false;
    }
   
    if (!success) return false;

    if (hourglass) {
      double phicrit=phi0approx-asin(0.5*rcrit*rinvapprox);
      bool keep=(phicrit>phicritminmc)&&(phicrit<phicritmaxmc);
      if (!keep) return false;
    }
    
    for(unsigned int j=0;j<3;j++){
      if (minusNeighborDisk[j]) {
	phiprojdiskapprox[j]+=dphisector;
	phiprojdisk[j]+=dphisector;
      }
      if (plusNeighborDisk[j]) {
	phiprojdiskapprox[j]-=dphisector;
	phiprojdisk[j]-=dphisector;
      }	    
    }
	    
    for(unsigned int j=0;j<3;j++){
      if (minusNeighbor[j]) {
	phiprojapprox[j]+=dphisector;
	phiproj[j]+=dphisector;
      }
      if (plusNeighbor[j]) {
	phiprojapprox[j]-=dphisector;
	phiproj[j]-=dphisector;
      }
    }
	    
    
    if (writeTrackletParsDisk) {
      static ofstream out("trackletparsdisk.txt");
      out <<"Trackpars         "<<disk_
	  <<"   "<<rinv<<" "<<rinvapprox<<" "<<rinvapprox
	  <<"   "<<phi0<<" "<<phi0approx<<" "<<phi0approx
	  <<"   "<<t<<" "<<tapprox<<" "<<tapprox
	  <<"   "<<z0<<" "<<z0approx<<" "<<z0approx
	  <<endl;
    }
	    
    FPGATracklet* tracklet=new FPGATracklet(innerStub,middleStub,outerStub,
					    innerFPGAStub,middleFPGAStub,outerFPGAStub,
					    iSector_,
					    phioffset_,
					    rinv,phi0,d0,z0,t,
					    rinvapprox,phi0approx,d0approx,
					    z0approx,tapprox,
					    irinv,iphi0,id0,iz0,it,
					    validproj,
					    iphiproj,izproj,iphider,izder,
					    minusNeighbor,plusNeighbor,	
					    phiproj,zproj,phider,zder,
					    phiprojapprox,zprojapprox,
					    phiderapprox,zderapprox,
					    validprojdisk,
					    iphiprojdisk,irprojdisk,
					    iphiderdisk,irderdisk,
					    minusNeighborDisk,
					    plusNeighborDisk,
					    phiprojdisk,rprojdisk,
					    phiderdisk,rderdisk,
					    phiprojdiskapprox,
					    rprojdiskapprox,
					    phiderdiskapprox,
					    rderdiskapprox,
					    true);
    
    if (debug1) {
      cout << "Found tracklet in disk = "<<disk_<<" "<<tracklet
	   <<" "<<iSector_<<endl;
    }
        
    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);

    trackletpars_->addTracklet(tracklet);
    ofstream fout ("seedsTC.txt", ofstream::app);
    fout << tracklet->eta () << " " << tracklet->getISeed () << endl;
    fout.close ();
    
    if (tracklet->validProj(1)) {
      addLayerProj(tracklet,1);
    }
    
    if (tracklet->validProj(2)) {
      addLayerProj(tracklet,2);
    }
    
    for(unsigned int j=0;j<3;j++){
      if (tracklet->validProjDisk(sign*dproj_[j])) {
	addDiskProj(tracklet,sign*dproj_[j]);
      }
    }

    return true;
    
  }
  

  bool overlapSeeding(FPGAStub* innerFPGAStub, L1TStub* innerStub, FPGAStub* middleFPGAStub, L1TStub* middleStub, FPGAStub* outerFPGAStub, L1TStub* outerStub){
    
    //Deal with overlap stubs here
    assert(outerFPGAStub->isBarrel());
    
    assert(innerFPGAStub->isDisk());
    
    disk_=innerFPGAStub->disk().value();
    
    if (debug1) {
      cout << "trying to make overlap tracklet disk_ = "<<disk_<<" "<<getName()<<endl;
    }
    
    //int sign=1;
    //if (disk_<0) sign=-1;
    
    double r1=innerStub->r();
    double z1=innerStub->z();
    double phi1=innerStub->phi();
    
    double r2=middleStub->r();
    double z2=middleStub->z();
    //double phi2=middleStub->phi();

    double r3=outerStub->r();
    double z3=outerStub->z();
    double phi3=outerStub->phi();
    
    //Protection... Should be handled cleaner
    //to avoid problem with floating point 
    //calculation and with overflows
    //in the integer calculation
    if (fabs(z2)>120&&(r1<r3+1.5 || r2<r3+1.5)) {
      cout << "in overlap tracklet: radii wrong"<<endl;
      return false;
    }
    if (fabs(z2)<120&&(r1<r2+1.5)) {
      cout << "in overlap tracklet: radii wrong"<<endl;
      return false;
    }
    

    double rinv,phi0,d0,t,z0;
	    
    double phiproj[3],zproj[3],phider[3],zder[3];
    double phiprojdisk[4],rprojdisk[4],phiderdisk[4],rderdisk[4];
    
    exacttrackletOverlap(r1,z1,phi1,r3,z3,phi3,outerStub->sigmaz(),
			 rinv,phi0,d0,t,z0,
			 phiproj,zproj,phider,zder,
			 phiprojdisk,rprojdisk,phiderdisk,rderdisk);
    
    
    //Truncates floating point positions to integer
    //representation precision
    if (useapprox) {
      phi1=innerFPGAStub->phiapprox(phimin_,phimax_);
      z1=innerFPGAStub->zapprox();
      r1=innerFPGAStub->rapprox();
	      
      //phi2=outerFPGAStub->phiapprox(phimin_,phimax_);
      z2=outerFPGAStub->zapprox();
      r2=outerFPGAStub->rapprox();
    }

    double rinvapprox,phi0approx,d0approx,tapprox,z0approx;
    double phiprojapprox[3],zprojapprox[3],phiderapprox[3],zderapprox[3];
    double phiprojdiskapprox[4],rprojdiskapprox[4],
      phiderdiskapprox[4],rderdiskapprox[4];

    phi0 -= 0.171;
    
    //store the approximate results
    rinvapprox = rinv;
    phi0approx = phi0;
    d0approx = d0;
    tapprox    = t;
    z0approx   = z0;

    for(int i=0; i<3; ++i){
      phiproj[i] -= 0.171;
      phiprojapprox[i] = phiproj[i];
      zprojapprox[i]   = zproj[i];
      phiderapprox[i] = phider[i];
      zderapprox[i]   = zder[i];
    }

    for(int i=0; i<4; ++i){
      phiprojdisk[i] -= 0.171;
      phiprojdiskapprox[i] = phiprojdisk[i];
      rprojdiskapprox[i] = rprojdisk[i];
      phiderdiskapprox[i] = phiderdisk[i];
      rderdiskapprox[i]   = rderdisk[i];
    }

    //now binary
    double krinv = kphi1/kr*pow(2,rinv_shift),
           kphi0 = kphi1*pow(2,phi0_shift),
           kt = kz/kr*pow(2,t_shift),
           kz0 = kz*pow(2,z0_shift),
           kphiproj = kphi1*pow(2,SS_phiL_shift),
           kphider = kphi1/kr*pow(2,SS_phiderL_shift),
           kzproj = kz*pow(2,PS_zL_shift),
           kzder = kz/kr*pow(2,PS_zderL_shift),
           kphiprojdisk = kphi1*pow(2,SS_phiD_shift),
           kphiderdisk = kphi1/kr*pow(2,SS_phiderD_shift),
           krprojdisk = kr*pow(2,PS_rD_shift),
           krderdisk = kr/kz*pow(2,PS_rderD_shift);

    int irinv,iphi0,id0,it,iz0;
    bool validproj[3];
    int iphiproj[3],izproj[3],iphider[3],izder[3];
    bool minusNeighbor[3],plusNeighbor[3];
    
    bool validprojdisk[4];
    int iphiprojdisk[4],irprojdisk[4],iphiderdisk[4],irderdisk[4];
    bool minusNeighborDisk[4],plusNeighborDisk[4];
    
    //store the binary results
    irinv = rinvapprox / krinv;
    iphi0 = phi0approx / kphi0;
    id0 = d0approx / kd0;
    it    = tapprox / kt;
    iz0   = z0approx / kz0;

    //"protection" from the original, reinterpreted
    if (iz0>= 1<<(nbitsz0-1)) iz0=(1<<(nbitsz0-1))-1;
    if (iz0<=-(1<<(nbitsz0-1))) iz0=1-(1<<(nbitsz0-1))-1; 
    if (irinv>= (1<<(nbitsrinv-1))) irinv=(1<<(nbitsrinv-1))-1;
    if (irinv<=-(1<<(nbitsrinv-1))) irinv=1-(1<<(nbitsrinv-1))-1; 


    for(int i=0; i<3; ++i){
      iphiproj[i] = phiprojapprox[i] / kphiproj;
      izproj[i]   = zprojapprox[i] / kzproj;

      iphider[i] = phiderapprox[i] / kphider;
      izder[i]   = zderapprox[i] / kzder;

      validproj[i] = true;
      if (izproj[i]<-(1<<(nbitszprojL123-1))) validproj[i]=false;
      if (izproj[i]>=(1<<(nbitszprojL123-1))) validproj[i]=false;
      
      minusNeighbor[i]=false;
      plusNeighbor[i]=false;
      if (!hourglass) {
	if (iphiproj[i]<(1<<nbitsphistubL456)/8) {
	  minusNeighbor[i]=true;
	  iphiproj[i]+=3*(1<<nbitsphistubL456)/4;
	}
	if (iphiproj[i]>=7*(1<<nbitsphistubL456)/8) {
	  plusNeighbor[i]=true;
	  iphiproj[i]-=3*(1<<nbitsphistubL456)/4;
	}
      }

      //this is left from the original....
      if (iphiproj[i]>=(1<<nbitsphistubL456)) {
	iphiproj[i]=(1<<nbitsphistubL456)-2; //-2 not to hit atExtreme
	validproj[i] = false;
      }

      if (rproj_[i]<60.0)
	iphiproj[i]>>=(nbitsphistubL456-nbitsphistubL123);
      else {
	izproj[i]>>=(nbitszprojL123-nbitszprojL456);
      }

      if (iphiproj[i]<=0) {
	iphiproj[i]=1;
	validproj[i] = false;
      }

      
    }


    for(int i=0; i<4; ++i){
      iphiprojdisk[i] = phiprojdiskapprox[i] / kphiprojdisk;
      irprojdisk[i]   = rprojdiskapprox[i] / krprojdisk;

      iphiderdisk[i] = phiderdiskapprox[i] / kphiderdisk;
      irderdisk[i]   = rderdiskapprox[i] / krderdisk;

      minusNeighborDisk[i]=false;
      plusNeighborDisk[i]=false;
      if (!hourglass) {
	if (iphiprojdisk[i]<(1<<nbitsphistubL123)/8) {
	  minusNeighborDisk[i]=true;
	  iphiprojdisk[i]+=3*(1<<nbitsphistubL123)/4;
	}
	if (iphiprojdisk[i]>=7*(1<<nbitsphistubL123)/8) {
	  plusNeighborDisk[i]=true;
	  iphiprojdisk[i]-=3*(1<<nbitsphistubL123)/4;
	}
      }
	
      //"protection" from the original
      if (iphiprojdisk[i]<0) iphiprojdisk[i]=0;
      if (iphiprojdisk[i]>=(1<<nbitsphistubL123)) iphiprojdisk[i]=(1<<nbitsphistubL123)-1;
      
      validprojdisk[i]=true;
      if(irprojdisk[i]<=0 || irprojdisk[i] > 120. / krprojdisk ){
	validprojdisk[i]=false;
	irprojdisk[i] = 0;
	iphiprojdisk[i] = 0;
	iphiderdisk[i]  = 0;
	irderdisk[i]    = 0;
      }
    }
    
    bool success = true;
    if(fabs(tapprox)>10)
      success = false;
     if(fabs(rinvapprox)>rinvcut){
      if (debug1) 
	cout << "FPGATrackletCalculator::OverlapSeeding irinv too large: "<<rinvapprox<<endl;
      success = false;
    }
    if (fabs(z0approx)>z0cut) {
      if (debug1) cout << "Failed tracklet z0 cut "<<z0approx<<" in layer 1"<<endl;
      success = false;
    }
    if (fabs(d0approx)>maxd0) {
      if (debug1) cout << "Failed tracklet d0 cut "<<d0approx<<endl;
      success = false;
    }

    if (!success) {
      return false;
    }

    if (hourglass) {
      double phicrit=phi0approx-asin(0.5*rcrit*rinvapprox);
      bool keep=(phicrit>phicritminmc)&&(phicrit<phicritmaxmc);
      if (!keep) return false;
    }
    

    for(unsigned int j=0;j<3;j++){
      if (minusNeighborDisk[j]) {
	phiprojdiskapprox[j]+=dphisector;
	phiprojdisk[j]+=dphisector;
      }
      if (plusNeighborDisk[j]) {
	phiprojdiskapprox[j]-=dphisector;
	phiprojdisk[j]-=dphisector;
      }	    
    }
    
    for(unsigned int j=0;j<3;j++){
      if (minusNeighbor[j]) {
	phiprojapprox[j]+=dphisector;
	phiproj[j]+=dphisector;
      }
      if (plusNeighbor[j]) {
	phiprojapprox[j]-=dphisector;
	phiproj[j]-=dphisector;
      }
    }
    
    
    if (writeTrackletParsOverlap) {
      static ofstream out("trackletparsoverlap.txt");
      out <<"Trackpars "<<disk_
	  <<"   "<<rinv<<" "<<irinv<<" "<<rinvapprox
	  <<"   "<<phi0<<" "<<iphi0<<" "<<phi0approx
	  <<"   "<<t<<" "<<it<<" "<<tapprox
	  <<"   "<<z0<<" "<<iz0<<" "<<z0approx
	  <<endl;
    }
	      
    FPGATracklet* tracklet=new FPGATracklet(innerStub,middleStub,outerStub,
					    innerFPGAStub,middleFPGAStub,outerFPGAStub,
					    iSector_,
					    phioffset_,
					    rinv,phi0,d0,z0,t,
					    rinvapprox,phi0approx,d0approx,
					    z0approx,tapprox,
					    irinv,iphi0,id0,iz0,it,
					    validproj,
					    iphiproj,izproj,iphider,izder,
					    minusNeighbor,plusNeighbor,
					    
					    
					    phiproj,zproj,phider,zder,
					    phiprojapprox,zprojapprox,
					    phiderapprox,zderapprox,
					    validprojdisk,
					    iphiprojdisk,irprojdisk,
					    iphiderdisk,irderdisk,
					    minusNeighborDisk,
					    plusNeighborDisk,
					    phiprojdisk,rprojdisk,
					    phiderdisk,rderdisk,
					    phiprojdiskapprox,
					    rprojdiskapprox,
					    phiderdiskapprox,
					    rderdiskapprox,
					    false,true);
    
    if (debug1) {
      cout << "Found tracklet in overlap = "<<layer_<<" "<<disk_
	   <<" "<<tracklet<<" "<<iSector_<<endl;
    }
    
        
    tracklet->setTrackletIndex(trackletpars_->nTracklets());
    tracklet->setTCIndex(TCIndex_);
    
    trackletpars_->addTracklet(tracklet);
    ofstream fout ("seedsTC.txt", ofstream::app);
    fout << tracklet->eta () << " " << tracklet->getISeed () << endl;
    fout.close ();
    
    int layer=outerFPGAStub->layer().value()+1;
    
    if (layer==2) {
      if (tracklet->validProj(1)) {
	addLayerProj(tracklet,1);
      }
    }
    
    
    for(unsigned int disk=2;disk<6;disk++){
      if (layer==2 && disk==5 ) continue;
      if (tracklet->validProjDisk(disk)) {
	addDiskProj(tracklet,disk);
      }
    }

    return true;
    
  }
  
  
  int round_int( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
  }
 
    
private:

  int TCIndex_;
  int layer_;
  int disk_;
  double phimin_;
  double phimax_;
  double phioffset_;
  double rproj_[4];
  int lproj_[4];
  double zproj_[3];
  int dproj_[3];

  unsigned int maxtracklet_; //maximum numbor of tracklets that be stored
  
  double zprojoverlap_[4];

  vector<FPGAAllStubs*> innerallstubs_;
  vector<FPGAAllStubs*> middleallstubs_;
  vector<FPGAAllStubs*> outerallstubs_;
  vector<FPGAStubTriplets*> stubtriplets_;

  FPGATrackletParameters* trackletpars_;

  FPGATrackletProjections* trackletproj_L1PHI1_;
  FPGATrackletProjections* trackletproj_L1PHI2_;
  FPGATrackletProjections* trackletproj_L1PHI3_;
  FPGATrackletProjections* trackletproj_L1PHI4_;
  FPGATrackletProjections* trackletproj_L1PHI5_;
  FPGATrackletProjections* trackletproj_L1PHI6_;
  FPGATrackletProjections* trackletproj_L1PHI7_;
  FPGATrackletProjections* trackletproj_L1PHI8_;

  FPGATrackletProjections* trackletproj_L2PHI1_;
  FPGATrackletProjections* trackletproj_L2PHI2_;
  FPGATrackletProjections* trackletproj_L2PHI3_;
  FPGATrackletProjections* trackletproj_L2PHI4_;

  FPGATrackletProjections* trackletproj_L3PHI1_;
  FPGATrackletProjections* trackletproj_L3PHI2_;
  FPGATrackletProjections* trackletproj_L3PHI3_;
  FPGATrackletProjections* trackletproj_L3PHI4_;

  FPGATrackletProjections* trackletproj_L4PHI1_;
  FPGATrackletProjections* trackletproj_L4PHI2_;
  FPGATrackletProjections* trackletproj_L4PHI3_;
  FPGATrackletProjections* trackletproj_L4PHI4_;

  FPGATrackletProjections* trackletproj_L5PHI1_;
  FPGATrackletProjections* trackletproj_L5PHI2_;
  FPGATrackletProjections* trackletproj_L5PHI3_;
  FPGATrackletProjections* trackletproj_L5PHI4_;

  FPGATrackletProjections* trackletproj_L6PHI1_;
  FPGATrackletProjections* trackletproj_L6PHI2_;
  FPGATrackletProjections* trackletproj_L6PHI3_;
  FPGATrackletProjections* trackletproj_L6PHI4_;

  FPGATrackletProjections* trackletproj_D1PHI1_;
  FPGATrackletProjections* trackletproj_D1PHI2_;
  FPGATrackletProjections* trackletproj_D1PHI3_;
  FPGATrackletProjections* trackletproj_D1PHI4_;

  FPGATrackletProjections* trackletproj_D2PHI1_;
  FPGATrackletProjections* trackletproj_D2PHI2_;
  FPGATrackletProjections* trackletproj_D2PHI3_;
  FPGATrackletProjections* trackletproj_D2PHI4_;

  FPGATrackletProjections* trackletproj_D3PHI1_;
  FPGATrackletProjections* trackletproj_D3PHI2_;
  FPGATrackletProjections* trackletproj_D3PHI3_;
  FPGATrackletProjections* trackletproj_D3PHI4_;

  FPGATrackletProjections* trackletproj_D4PHI1_;
  FPGATrackletProjections* trackletproj_D4PHI2_;
  FPGATrackletProjections* trackletproj_D4PHI3_;
  FPGATrackletProjections* trackletproj_D4PHI4_;

  FPGATrackletProjections* trackletproj_D5PHI1_;
  FPGATrackletProjections* trackletproj_D5PHI2_;
  FPGATrackletProjections* trackletproj_D5PHI3_;
  FPGATrackletProjections* trackletproj_D5PHI4_;


  
  FPGATrackletProjections* trackletproj_L1Plus_; 
  FPGATrackletProjections* trackletproj_L1Minus_;
			                         
  FPGATrackletProjections* trackletproj_L2Plus_; 
  FPGATrackletProjections* trackletproj_L2Minus_;
			                         
  FPGATrackletProjections* trackletproj_L3Plus_; 
  FPGATrackletProjections* trackletproj_L3Minus_;
			                         
  FPGATrackletProjections* trackletproj_L4Plus_; 
  FPGATrackletProjections* trackletproj_L4Minus_;
			                         
  FPGATrackletProjections* trackletproj_L5Plus_; 
  FPGATrackletProjections* trackletproj_L5Minus_;
			                         
  FPGATrackletProjections* trackletproj_L6Plus_; 
  FPGATrackletProjections* trackletproj_L6Minus_;


  FPGATrackletProjections* trackletproj_D1Plus_; 
  FPGATrackletProjections* trackletproj_D1Minus_;
			                         
  FPGATrackletProjections* trackletproj_D2Plus_; 
  FPGATrackletProjections* trackletproj_D2Minus_;
			                         
  FPGATrackletProjections* trackletproj_D3Plus_; 
  FPGATrackletProjections* trackletproj_D3Minus_;
			                         
  FPGATrackletProjections* trackletproj_D4Plus_; 
  FPGATrackletProjections* trackletproj_D4Minus_;
			                         
  FPGATrackletProjections* trackletproj_D5Plus_; 
  FPGATrackletProjections* trackletproj_D5Minus_;
};

#endif
