#include "RimatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"
#include "ElectricFieldBodyForce.h"
#include "ElectricFieldBodyForceExplicit.h"
#include "BodyForceDummy.h"
#include "ElectricFieldBodyForceExplicitRamp.h"
#include "ApparentDynamicViscosityAux.h"
// #include "INSMomentumBaseLES.h"
#include "INSMomentumLaplaceFormLES.h"
#include "INSMomentumLaplaceFormRZLES.h"
#include "FluidTurbulent.h"
#include "ApparentDynamicViscosityWALEAux.h"
#include "ElemVolumeAux.h"
#include "ApparentDynamicViscosityMLAux.h"
#include "ApparentDynamicViscosityNikuradseAux.h"
// #include "INSMomentumBaseWALE.h"
#include "INSMomentumLaplaceFormWALE.h"
#include "INSMomentumLaplaceFormRZWALE.h"
#include "LengthScaleAux.h"
#include "INSkTransport.h"
#include "INSepsilonTransport.h"
#include "INSMomentumLaplaceFormTurbulent.h"
#include "INSMomentumLaplaceFormRZTurbulent.h"
#include "ApparentDynamicViscosityWALEAverageAux.h"
#include "SmoothFunction.h"
#include "TestAux.h"
#include "ApparentDynamicViscosityWALENonmeshAverageAux.h"
#include "InverseWallDistance.h"
#include "OnePointBoundedInverseDistanceDirichletBC.h"
#include "TwoPointsMinInverseDistanceDirichletBC.h"
#include "WallDistanceAux.h"
#include "InverseWallDistanceRZ.h"
#include "InverseWallDistanceGrad.h"
#include "SpalartAllmaras.h"
#include "SpalartAllmarasAux.h"

template<>
InputParameters validParams<RimatApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

RimatApp::RimatApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  RimatApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  RimatApp::associateSyntax(_syntax, _action_factory);
}

RimatApp::~RimatApp()
{
}

// External entry point for dynamic application loading
extern "C" void RimatApp__registerApps() { RimatApp::registerApps(); }
void
RimatApp::registerApps()
{
  registerApp(RimatApp);
}

// External entry point for dynamic object registration
extern "C" void RimatApp__registerObjects(Factory & factory) { RimatApp::registerObjects(factory); }
void
RimatApp::registerObjects(Factory & factory)
{
  registerKernel(ElectricFieldBodyForce);
  registerKernel(ElectricFieldBodyForceExplicit);
  registerKernel(BodyForceDummy);
  registerKernel(ElectricFieldBodyForceExplicitRamp);
  registerAux(ApparentDynamicViscosityAux);
  // registerKernel(INSMomentumBaseLES);
  registerKernel(INSMomentumLaplaceFormLES);
  registerKernel(INSMomentumLaplaceFormRZLES);
  registerMaterial(FluidTurbulent);
  registerAux(ApparentDynamicViscosityWALEAux);
  registerAux(ElemVolumeAux);
  registerAux(ApparentDynamicViscosityMLAux);
  registerAux(ApparentDynamicViscosityNikuradseAux);
  // registerKernel(INSMomentumBaseWALE);
  registerKernel(INSMomentumLaplaceFormWALE);
  registerKernel(INSMomentumLaplaceFormRZWALE);
  registerAux(LengthScaleAux);
  registerKernel(INSkTransport);
  registerKernel(INSepsilonTransport);
  registerKernel(INSMomentumLaplaceFormTurbulent);
  registerKernel(INSMomentumLaplaceFormRZTurbulent);
  registerAux(ApparentDynamicViscosityWALEAverageAux);
  registerUserObject(SmoothFunction);
  registerAux(TestAux);
  registerAux(ApparentDynamicViscosityWALENonmeshAverageAux);
  registerKernel(InverseWallDistance);
  registerBoundaryCondition(OnePointBoundedInverseDistanceDirichletBC);
  registerBoundaryCondition(TwoPointsMinInverseDistanceDirichletBC);
  registerAux(WallDistanceAux);
  registerKernel(InverseWallDistanceRZ);
  registerKernel(InverseWallDistanceGrad);
  registerKernel(SpalartAllmaras);
  registerAux(SpalartAllmarasAux);

}

// External entry point for dynamic syntax association
extern "C" void RimatApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { RimatApp::associateSyntax(syntax, action_factory); }
void
RimatApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
