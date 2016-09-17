[GlobalParams]
  gravity = '0 0 0'
  rho = 1.105
  mu = 0.00001989
  integrate_p_by_parts = false
[]

[Mesh]
  file = pin_to_plate.e
[]

[Variables]

  [./vel_x]
    order = SECOND
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./vel_y]
    order = SECOND
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./p]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0
    [../]
  [../]
[]

[AuxVariables]
  [./potential]
    order = FIRST
    family = LAGRANGE
  [../]

  [./space_charge_density]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]

  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    p = p
  [../]

  [./x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  [../]

  [./y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  [../]

  [./x_momentum_space]
    type = INSMomentum
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
  [../]

  [./x_body_force]
    type = ElectricFieldBodyForce
    variable = vel_x
    potential = potential
    space_charge_density = space_charge_density
    component = 0
  [../]

  [./y_momentum_space]
    type = INSMomentum
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
  [../]

  [./y_body_force]
    type = ElectricFieldBodyForce
    variable = vel_y
    potential = potential
    space_charge_density = space_charge_density
    component = 1
  [../]
[]

[AuxKernels]
  [./potential_aux]
    type = SolutionAux
    solution = corona_discharge_solution
    variable = potential
    from_variable = voltage
  [../]

  [./space_charge_density_aux]
    type = SolutionAux
    solution = corona_discharge_solution
    variable = space_charge_density
    from_variable = density
  [../]
[]

[UserObjects]
  [./corona_discharge_solution]
    type = SolutionUserObject
    mesh = pin_to_plate_metamoq_out.e
    system_variables = 'voltage density'
  [../]
[]

[BCs]
  active = 'x_no_slip y_no_slip outlet_p'

  [./x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'ground anode'
    value = 0.0
  [../]

  [./y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'ground anode'
    value = 0.0
  [../]

  [./outlet_p]
    type = DirichletBC
    variable = p
    boundary = 'far back'
    value = 0.0
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  solve_type = 'PJFNK'
  type = Transient
  end_time = 1
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dtmin = 1e-12
  trans_ss_check = true
  ss_check_tol = 5e-5

  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-4
    growth_factor = 1.2
    optimal_iterations = 15
  [../]
[]

[Outputs]
  exodus = true
[]

[Debug]
  show_var_residual_norms = true
  show_var_residual = 'vel_x vel_y p'
[]
