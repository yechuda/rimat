#ifndef RIMATAPP_H
#define RIMATAPP_H

#include "MooseApp.h"

class RimatApp;

template<>
InputParameters validParams<RimatApp>();

class RimatApp : public MooseApp
{
public:
  RimatApp(InputParameters parameters);
  virtual ~RimatApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* RIMATAPP_H */
