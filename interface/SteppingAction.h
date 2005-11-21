#ifndef SimG4Core_SteppingAction_H
#define SimG4Core_SteppingAction_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4UserSteppingAction.hh"
#include "boost/signal.hpp"

class SteppingAction: public G4UserSteppingAction
{
public:
    SteppingAction(const edm::ParameterSet & ps);
    ~SteppingAction();
    void UserSteppingAction(const G4Step * aStep);

    boost::signal< void(const G4Step*)> m_g4StepSignal;
private:
    void catchLowEnergyInVacuumHere(const G4Step * aStep);
    void catchLowEnergyInVacuumNext(const G4Step * aStep);
private:
    bool   killBeamPipe;
    double theCriticalEnergyForVacuum;
    double theCriticalDensity;
    int    verbose;
};

#endif
