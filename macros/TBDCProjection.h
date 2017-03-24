#ifndef TBDCProjection_h
#define TBDCProjection_h

class TBDCProjection{
public:
    TBDCProjection();
    TBDCProjection(double bdc1z, double bdc2z, double TGTz);
    virtual ~TBDCProjection();
    void ProjectParticle(double sx, double sy, double sz, double sa, double sb, double charge, double energy, double endz, double mass);
    double getX();
    double getY();
    double getZ();
    double getA();
    double getB();
    double getPX();
    double getPY();
    double getPZ();
    double getE();
    double getKE();
    double getMeVu();
    double getP();
    double getBeta();
    void setBeam(int runNo);
    int getBeam();
private:
    int m_beam;
    double m_beta;
    double m_brho;
    double m_charge;
    double m_aoq;
    double m_mass;
    double m_energy;
    double m_kinetic_energy;
    double m_momentum;
    double m_bdc1x;
    double m_bdc1y;
    double m_bdc1z;
    double m_bdc2x;
    double m_bdc2y;
    double m_bdc2z;
    double m_x;
    double m_y;
    double m_z;
    double m_px;
    double m_py;
    double m_pz;
    double m_a;
    double m_b;
    double m_Bxx[300];
    double m_Byy[300];
    double m_Bzz[300];
    double Eloss[9];
    double m_Bx;
    double m_By;
    double m_Bz;
    double m_end_z;
    double dx;
    double dy;
    double dz;
    double dr;
    void MagStep();
    void TLoadField();
    void updateEnergy2();
    void setMomentum();
    void setBeta();
    void setE();
    void setKE();
    double getByField(double x, double y, double z);
};
#endif
