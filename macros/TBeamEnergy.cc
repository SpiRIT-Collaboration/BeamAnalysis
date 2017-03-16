#include "TBeamEnergy.h"
#include <iostream>
#include <cmath>

TBeamEnergy::TBeamEnergy(){m_corr=false;}
TBeamEnergy::TBeamEnergy(double z, double aoq, double beta78){
    m_corr=false;
    setBeta78(beta78);
    setZ(z);
    setAoq(aoq);
    setMass();
    setEnergy78();
    //std::cout << "TBeamEnergy" << std::endl;
}
TBeamEnergy::~TBeamEnergy(){}

void TBeamEnergy::setBeta78(double beta78){
    m_beta78=beta78;
}
double TBeamEnergy::getBeta78(){
    return m_beta78;
}
void TBeamEnergy::setCorrection(){//use calibrate energy loss
    double dE132=0.9374555242;
    double dE124=0.9384493375;
    double dE112=0.9328597049;
    double dE108=0.9308606678;
    //setAoq(2.64);
    //m_mass=160440.1;
    //setEnergy78();
    //double kinetic_energy=39903;//m_energy78-m_mass;
    double kinetic_energy=m_energy78-m_mass;
    if(m_beam>130 && m_beam<140) kinetic_energy=kinetic_energy*dE132;
    if(m_beam>120 && m_beam<130) kinetic_energy=kinetic_energy*dE124;
    if(m_beam>110 && m_beam<120) kinetic_energy=kinetic_energy*dE112;
    if(m_beam>100 && m_beam<110) kinetic_energy=kinetic_energy*dE108;
    m_energy=kinetic_energy+m_mass;
    m_beta=std::sqrt(1-(m_mass/(m_energy))*(m_mass/(m_energy)));
    //m_momentum=m_mass*m_beta/std::sqrt(std::abs(1-m_beta*m_beta));
    m_momentum=std::sqrt(m_energy*m_energy-m_mass*m_mass);
    //m_energy=m_mass*(1/std::sqrt(1-m_beta*m_beta));
    m_corr=true;
    //std::cout << kinetic_energy << std::endl;
}

double TBeamEnergy::getBeta(){
    if(not m_corr) setCorrection();
    return m_beta;
}
double TBeamEnergy::getBrho(){
  if(not m_corr) setCorrection();
  m_brho=3.3356*m_momentum/m_z/1000.;
  return m_brho;
}
void TBeamEnergy::setZ(double z){
    m_z=z;
}
double TBeamEnergy::getZ(){
    return m_z;
}
void TBeamEnergy::setAoq(double aoq){
    m_aoq=aoq;
}
double TBeamEnergy::getAoq(){
    return m_aoq;
}
double TBeamEnergy::setMass(){
    m_mass = /*(double)((int)*/((m_aoq*m_z)*931.494*0.9993774);
}
double TBeamEnergy::getMass(){
    return m_mass;
}
double TBeamEnergy::setEnergy78(){
    m_energy78=m_mass*(1/std::sqrt(1-m_beta78*m_beta78));
}
double TBeamEnergy::getEnergy78(){
    return m_energy78;
}
double TBeamEnergy::getCorrectedEnergy(){
    if(not m_corr) setCorrection();
    return m_energy;
}
double TBeamEnergy::getCorrectedMomentum(){
    if(not m_corr) setCorrection();
    return m_momentum;
}
void TBeamEnergy::setBeam(int runNo){
    if(runNo>=2174 && runNo<=2509) m_beam=108;
    if(runNo>=2520 && runNo<=2653) m_beam=112;
    if(runNo>=3044 && runNo<=3184) m_beam=124;
    if(runNo>=2819 && runNo<=3039) m_beam=132;
}
int TBeamEnergy::getBeam(){
    return m_beam;
}
