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
// #include "FluidLES.h"
#include "ApparentDynamicViscosityWALEAux.h"
#include "ElemVolumeAux.h"

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
  // registerMaterial(FluidLES);
  registerAux(ApparentDynamicViscosityWALEAux);
  registerAux(ElemVolumeAux);
}

// External entry point for dynamic syntax association
extern "C" void RimatApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { RimatApp::associateSyntax(syntax, action_factory); }
void
RimatApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
