#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TArrayD.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"


class TelescopeTracks {
    
    double p0_, p1_, p2_, p3_, p4_, p5_;
    
    double xOnPlane_, yOnPlane_, zOnPlane_; 

    int iLayer;
    
    std::vector<SiPixelCluster> clusterList_;
    
    std::vector<TVector3> globalPoints;
    std::vector<TVector3> globalPoints_err;
    std::vector<TVector3> pseudolocalPoints;
    std::vector<TVector3> pseudolocalPoints_err;
    std::vector<int>      modDetId_;
    std::vector<double>   param_ass_a_;
    std::vector<double>   param_ass_b_;
    std::vector<double>   param_ass_c_;
    std::vector<double>   param_ass_d_;
    std::vector<double>   param_ass_x0_;
    std::vector<double>   param_ass_y0_;
    std::vector<double>   param_ass_z0_;
    double normChi2_;
    double chi2_;
    
    TMatrixD *rotationX, *rotationY;
    
  public:
    //Default Constructor 
    TelescopeTracks() 
    { 
        p0_ = -1000;
        p1_ = -1000;
        p2_ = -1000;
        p3_ = -1000;
	p4_ = -1000;
	p5_ = -1000;
	xOnPlane_= -1000;
	yOnPlane_= -1000;
	zOnPlane_= -1000;
	
	
/*
	rotationX = new TMatrixD(3, 3);
	rotationY = new TMatrixD(3, 3);
	
	TArrayD dataRx(9);
	TArrayD dataRy(9);
	
//	double costheta = cos(30*3.14159265/180.);
//	double sintheta = sin(30*3.14159265/180.);
	double costheta = cos(-30*3.14159265/180.);
	double sintheta = sin(-30*3.14159265/180.);
	
	
	double cosphi = cos(20*3.14159265/180.);
	double sinphi = sin(20*3.14159265/180.);
	
	dataRx[0] = 1; dataRx[1] = 0;        dataRx[2] = 0;
	dataRx[3] = 0; dataRx[4] = costheta; dataRx[5] = -1*sintheta;
	dataRx[6] = 0; dataRx[7] = sintheta; dataRx[8] = costheta;
	
	rotationX->SetMatrixArray(dataRx.GetArray());
	
	dataRy[0] = cosphi; dataRy[1] = 0; dataRy[2] = sinphi;
	dataRy[3] = 0;      dataRy[4] = 1; dataRy[5] = 0;
	dataRy[6] = -1.*sinphi; dataRy[7] = 0; dataRy[8] = cosphi;
	
	rotationY->SetMatrixArray(dataRy.GetArray());
	
*/	
	
	
	
    } 
    

   int GetLayerNumber(int theDetID)  {
       int iLayer = -1;
       if ((theDetID == 3173) || (theDetID == 3164))   {
	   iLayer = 0;
	}
	if ((theDetID == 3239) || (theDetID == 3023))   {
	   iLayer = 1;
	}
	if ((theDetID == 3265) || (theDetID == 3226))   {
	   iLayer = 2;
	}
	if ((theDetID == 3204) || (theDetID == 3192))   {
	   iLayer = 3;
	}
	if ((theDetID == 3090) || (theDetID == 3124))   {
	   iLayer = 4;
	}
	if ((theDetID == 3082) || (theDetID == 3175))   {
	   iLayer = 5;
	}
	if ((theDetID == 3009) || (theDetID == 3057))   {
	   iLayer = 6;
	}
	if ((theDetID == 3027) || (theDetID == 3074))   {
	   iLayer = 7;
	}
	return iLayer;
   }


    TelescopeTracks(double p0, double p1, double p2, double p3, double p4, double p5) 
    { 
        p0_ = p0;
        p1_ = p1;
        p2_ = p2;
        p3_ = p3;
	p4_ = p4;
	p5_ = p5;
	xOnPlane_= -1000;
	yOnPlane_= -1000;
	zOnPlane_= -1000;
	 
    } 
    void setTrackParam(int, double);
    
    double getParameter(int );
//    void propagateToPlane(std::vector<double> );
    
    double getXonPlane(){return xOnPlane_; };
    double getYonPlane(){return yOnPlane_; };
    double getZonPlane(){return zOnPlane_; };
    
    void addCluster(SiPixelCluster );
    void cleanClusterList();
    std::vector<SiPixelCluster> getclusterList();
    
    void addGlobalPoint(TVector3 );
    void cleanGlobalPoints();
    std::vector<TVector3> getGlobalPoints();
    
    void addGlobalPointErr(TVector3 );
    void cleanGlobalPointsErr();
    std::vector<TVector3> getGlobalPointsErr();
    
    void addPseudoLocalPoint(TVector3 );
    void cleanPseudoLocalPoints();
    std::vector<TVector3> getPseudoLocalPoints();
    
    void addPseudoLocalPointErr(TVector3 );
    void cleanPseudoLocalPointsErr();
    std::vector<TVector3> getPseudoLocalPointsErr();

    void addAssPlaneA(double aval){param_ass_a_.push_back(aval);};
    void addAssPlaneB(double bval){param_ass_b_.push_back(bval);};
    void addAssPlaneC(double cval){param_ass_c_.push_back(cval);};
    void addAssPlaneD(double dval){param_ass_d_.push_back(dval);};
    void addAssPlaneX0(double x0){param_ass_x0_.push_back(x0);};
    void addAssPlaneY0(double y0){param_ass_y0_.push_back(y0);};
    void addAssPlaneZ0(double z0){param_ass_z0_.push_back(z0);};
    void cleanAssPlaneA(){param_ass_a_.clear();};
    void cleanAssPlaneB(){param_ass_b_.clear();};
    void cleanAssPlaneC(){param_ass_c_.clear();};
    void cleanAssPlaneD(){param_ass_d_.clear();};
    void cleanAssPlaneX0(){param_ass_x0_.clear();};
    void cleanAssPlaneY0(){param_ass_y0_.clear();};
    void cleanAssPlaneZ0(){param_ass_z0_.clear();};
    std::vector<double> getAssPlaneA(){return param_ass_a_;};
    std::vector<double> getAssPlaneB(){return param_ass_b_;};
    std::vector<double> getAssPlaneC(){return param_ass_c_;};
    std::vector<double> getAssPlaneD(){return param_ass_d_;};
    std::vector<double> getAssPlaneX0(){return param_ass_x0_;};
    std::vector<double> getAssPlaneY0(){return param_ass_y0_;};
    std::vector<double> getAssPlaneZ0(){return param_ass_z0_;};
    
    void addModDetId(int thedetid){modDetId_.push_back(thedetid);};
    void cleanModDetId(){modDetId_.clear();};
    std::vector<int>  getModDetId(){return modDetId_;};
    
    double getNormChi2(){return normChi2_;};
    double getChi2(){return chi2_;};
    
    void line(double , double *, double &, double &, double &);
    void intersection(double *, double *, double &, double &, double &);
    
    void fitTrack();
    
    
} ;
