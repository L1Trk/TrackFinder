#ifndef FPGAVARIANCE_H
#define FPGAVARIANCE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <assert.h>
#include <math.h>


using namespace std;

class FPGAVariance{

public:

  FPGAVariance(SLHCEvent& ev) {
    process(ev);
  }

  void process(SLHCEvent& ev) {

    cout << "Process variance:"<<ev.nsimtracks()<<endl;
    assert(ev.nsimtracks()==1);

    L1SimTrack simtrk=ev.simtrack(0);

    double eta=simtrk.eta();
    double phi0=simtrk.phi();
    double pt=simtrk.pt();
    double z0=simtrk.vz();

    double t=sinh(eta);
    double rinv=-0.01*0.3*3.8/pt;

    double layerresidphi[6];
    double layerresidz[6];
    double diskresidphi[5];
    double diskresidr[5];

    for(unsigned int i=0;i<6;i++) {
      layerresidphi[i]=-9999.9;
      layerresidz[i]=-9999.9;
      if (i<5) {
	diskresidphi[i]=-9999.9;
	diskresidr[i]=-9999.9;
      }
    }

    for(int i=0;i<ev.nstubs();i++){

      L1TStub stub=ev.stub(i);

      int disk=-1;
      int layer=stub.layer()+1;
      if (layer>999) {
	disk=stub.module();
	layer=-1;
      }

      cout << "layer disk : "<<layer<<" "<<disk<<endl;
      
      assert(disk>0||layer>0);
      assert(!((disk>0)&&(layer>0)));

      if (disk>0) {
	double zproj=stub.z();
	double tmp=rinv*(zproj-z0)/(2.0*t);
	//double rproj=(2.0/rinv)*sin(tmp);
	double phiproj=phi0-tmp;
	double dphi=phiproj-stub.phi();
	if (dphi>0.5*two_pi) dphi-=two_pi;
	if (dphi<-0.5*two_pi) dphi+=two_pi;
	diskresidphi[disk-1]=dphi*stub.r();
	//diskresidr[disk-1]=rproj-stub.r();
      }

      if (layer>0) {
	double rproj=stub.r();
	double phiproj=phi0-asin(0.5*rproj*rinv);
	//double zproj=z0+(2*t/rinv)*asin(0.5*rproj*rinv);
	double dphi=phiproj-stub.phi();
	if (dphi>0.5*two_pi) dphi-=two_pi;
	if (dphi<-0.5*two_pi) dphi+=two_pi;
	layerresidphi[layer-1]=dphi*stub.r();
	//layerresidz[layer-1]=zproj-stub.z();
      }

    }

    static ofstream out("variance.txt");

    out << pt;

    for(unsigned int i=0;i<6;i++) {
      out <<" "<<layerresidphi[i];
    }
    
    for(unsigned int i=0;i<5;i++) {
      out <<" "<<diskresidphi[i];
    }
    out << endl;
      
      
    
    
  }
  
};



#endif



