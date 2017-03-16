#ifndef TBeamEnergy_h
#define TBeamEnergy_h

class TBeamEnergy{
public:
    TBeamEnergy();
    TBeamEnergy(double beta78, double z, double aoq);
    virtual ~TBeamEnergy();
    void setBeta78(double beta78);
    void setZ(double z);
    void setAoq(double aoq);
    void setCorrection();
    double getBeta78();
    double getZ();
    double getAoq();
    double setMass();
    double getMass();//returns estimated particle mass
    void setBeam(int runNo);
    int getBeam();
    double setEnergy78();
    double getEnergy78();//returns initial energy
    double getBeta();
    double getBrho();
    double getCorrectedEnergy();//returns energy after BDC2
    double getCorrectedMomentum();

private:
    int m_beam;
    double m_beta;
    double m_brho;
    double m_beta78;
    double m_z;
    double m_aoq;
    double m_mass;
    double m_energy78;
    double m_energy;
    double m_momentum;
    bool m_corr;
};
#endif
