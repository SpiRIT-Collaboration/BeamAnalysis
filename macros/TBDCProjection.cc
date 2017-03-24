#include "TBDCProjection.h"
#include <iostream>
#include <cmath>
#include <fstream>

TBDCProjection::TBDCProjection(){TLoadField();}
TBDCProjection::TBDCProjection(double bdc1z, double bdc2z, double TGTz){
  m_bdc1z=bdc1z;
  m_bdc2z=bdc2z;
  m_end_z=TGTz;
  TLoadField();
  //std::cout << "TBDCProjection" << std::endl;
}
TBDCProjection::~TBDCProjection(){}
void TBDCProjection::TLoadField(){
    std::ifstream Bfield;
    Bfield.open("../ReducedBMap.txt");
    int ii=0;
    double tx,ty,tz;
    double tbx,tby,tbz;
    while(ii<=300){
      Bfield>>tx>>ty>>tz>>tbx>>tby>>tbz;
      m_Bxx[ii] = tbx;
      m_Byy[ii] = tby;
      m_Bzz[ii] = tbz;
      ii++;
      if (!Bfield.good()) break;
    }
    Bfield.close();
}

void TBDCProjection::ProjectParticle(double sx, double sy, double sz, double sa, double sb, double charge, double energy, double endz, double mass){
    m_mass=mass;
    m_energy=energy;
    //m_mass=131.8904*931.494;//calibration check
    //m_energy=m_mass+37345;//calibration check
    setKE();
    double Elost=0;
    double E132[9]={0.998375283,0.9962917511,0.9998661133,0.9921531869,0.9554079952,0.9522248905,0.9997032993,0.9995844958,0.9509204276};//percent of KE left after total material, Sn132
    double E124[9]={0.9982873828,0.9960907838,0.9998870631,0.9917546733,0.9521382609,0.9484465177,0.9997162405,0.9995269333,0.9478749251};//percent of KE left after total material, Sn124
    double E112[9]={0.9981455152,0.9957267604,0.999844508,0.9909489596,0.9486189579,0.9441154088,0.9996495409,0.9995091852,0.9418098913};//percent of KE left after total material, Sn112
    double E108[9]={0.9980763682,0.9955349973,0.9998709344,0.9905769975,0.9452371645,0.9403412028,0.9996334848,0.9994866906,0.9388114453};//percent of KE left after total material, Sn108
    for(int i =0;i<9;i++){
        if(m_beam>130 && m_beam < 135) Eloss[i]=(m_kinetic_energy-Elost)*(1-E132[i]);
        if(m_beam>120 && m_beam < 130) Eloss[i]=(m_kinetic_energy-Elost)*(1-E124[i]);
        if(m_beam>110 && m_beam < 120) Eloss[i]=(m_kinetic_energy-Elost)*(1-E112[i]);
        if(m_beam>100 && m_beam < 110) Eloss[i]=(m_kinetic_energy-Elost)*(1-E108[i]);
        Elost=Elost+Eloss[i];
    }
    setMomentum();
    setBeta();
    m_x=sx;
    m_y=sy;
    m_z=sz;
    m_a=sa;
    m_b=sb;
    m_end_z=endz;
    m_charge=charge;
    m_brho=3.3356*m_momentum/m_charge/1000.;//given charge in multiples of proton charge, momentum in MeV/c, provides Brho in T*m
    //TLoadField();
    m_By=0;
    if(m_brho>0){
    while(m_z<m_end_z && std::abs(m_x)<430.){
      m_By=getByField(m_x,m_y,m_z);
      dz = 1.;//step size in mm
      if(m_z+dz>m_end_z) dz=m_end_z-m_z;
      MagStep();
    }
    }

    m_py=m_momentum*std::sin(m_b/1000.);
    m_px=std::sqrt(m_momentum*m_momentum-m_py*m_py)*std::sin(m_a/1000.);
    m_pz=std::sqrt(m_momentum*m_momentum-m_py*m_py-m_px*m_px);

}
void TBDCProjection::MagStep(){
    double da;
    m_brho=3.3356*m_momentum/m_charge/1000.;
    da=0.;//no change in angle unless |B|>0
    dx=-dz*std::tan(m_a/1000.);//linear projection used unless |B|>0
    dy=dz*std::tan(m_b/1000.);
    if(std::abs(m_By)>0.){
        double mrho=m_brho/m_By*1000.;
        da=-std::asin(dz/mrho)*1000.;//gives da in mrad
        dx=mrho*(std::cos((m_a+da)/1000.)-std::cos(m_a/1000.));//gives dx in mm
    }
    dr=std::sqrt(dz*dz+dy*dy+dx*dx);
    m_x=m_x+dx;
    m_y=m_y+dy;
    m_a=m_a+da;
    m_z=m_z+dz;
    updateEnergy2();
    setMomentum();
    setBeta();

}
void TBDCProjection::updateEnergy2(){//use a table of energy loss through material (from LISE++) to apply energy loss
    double thickness[9]={0.05,141.75,0.004,274.184,0.418,0.418,8.316,0.012,1447.};//thickness of material
    double bound1[9]={-1009.05,-1009.,-867.25,-867.246,-593.062,-592.644,-592.226,-583.91,-583.898};//start of each material
    double bound2[9]={-1009.,-867.25,-867.246,-593.062,-592.644,-592.226,-583.91,-583.898,863.102};//end of each material
    for(int i=0; i<9; i++){//for each material, determine how much of the material the step includes, and apply appropriate correction
        //if(m_z>bound2[i]) break;
        double eff_thick=0.;
        eff_thick=std::min(bound2[i],m_z)-std::max(bound1[i],m_z-dz);
        if(eff_thick>0){
            if((eff_thick/thickness[i]*Eloss[i])>0) m_kinetic_energy=m_kinetic_energy-(eff_thick*Eloss[i]/thickness[i]*dr/dz);
            m_energy=m_kinetic_energy+m_mass;
            setBeta();
            setMomentum();
        }

    }
}
void TBDCProjection::setMomentum(){
  m_momentum=std::min(std::sqrt(m_energy*m_energy-m_mass*m_mass),m_mass*m_beta/std::sqrt(std::abs(1-m_beta*m_beta)));//I assume that momentum is always lost
  //m_momentum=m_mass*m_beta/std::sqrt(std::abs(1-m_beta*m_beta));
}
void TBDCProjection::setBeta(){
  m_beta=std::sqrt(1-m_mass*m_mass/m_energy/m_energy);
}
double TBDCProjection::getByField(double x, double y, double z){
    m_By=0.;
    if(std::sqrt(x*x+z*z)<3000){
    int t_index = (int)(std::sqrt(z*z+x*x)/10.+0.5);
    m_By=m_Byy[t_index];
    }
    return m_By;
}
double TBDCProjection::getX(){ return m_x; }
double TBDCProjection::getY(){ return m_y; }
double TBDCProjection::getZ(){ return m_z; }
double TBDCProjection::getA(){ return m_a; }
double TBDCProjection::getB(){ return m_b; }
double TBDCProjection::getPX(){ return m_px; }
double TBDCProjection::getPY(){ return m_py; }
double TBDCProjection::getPZ(){ return m_pz; }
void TBDCProjection::setE(){ m_energy=m_mass/std::sqrt(1-m_beta*m_beta); }
double TBDCProjection::getE(){ return m_energy; }
void TBDCProjection::setKE(){ m_kinetic_energy=m_energy-m_mass; }
double TBDCProjection::getKE(){ return m_kinetic_energy; }
double TBDCProjection::getMeVu(){ return m_kinetic_energy/(m_mass/931.494); }
double TBDCProjection::getP(){ return m_momentum; }
double TBDCProjection::getBeta(){ return m_beta; }
void TBDCProjection::setBeam(int runNo){
    if(runNo>=2174 && runNo<=2509) m_beam=108;
    if(runNo>=2520 && runNo<=2653) m_beam=112;
    if(runNo>=3044 && runNo<=3184) m_beam=124;
    if(runNo>=2819 && runNo<=3039) m_beam=132;
}
int TBDCProjection::getBeam(){
    return m_beam;
}
