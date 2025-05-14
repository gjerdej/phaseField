// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
customAttributeLoader::loadVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "c");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "c, psi, grad(psi), phie");
  set_dependencies_gradient_term_RHS(0, "grad(mu), psi, grad(psi)");

  // Variable 1
  set_variable_name(1, "mu");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "c, grad(c), psi, grad(psi)");
  set_dependencies_gradient_term_RHS(1, "grad(c), psi, grad(psi)");

  // Variable 2
  set_variable_name(2, "psi");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, AUXILIARY);

  set_dependencies_value_term_RHS(2, "psi");
  set_dependencies_gradient_term_RHS(2, "");

  // Variable 3
  set_variable_name(3, "phie");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(3, "grad(psi), phie, c");
  set_dependencies_gradient_term_RHS(3, "psi, grad(phie)");
  set_dependencies_value_term_LHS(3, "");
  set_dependencies_gradient_term_LHS(3, "psi, grad(change(phie))");

  // // Variable 4
  // set_variable_name(4, "phie");
  // set_variable_type(4, SCALAR);
  // set_variable_equation_type(4, TIME_INDEPENDENT);

  // set_dependencies_value_term_RHS(4, "psi, grad(psi)");
  // set_dependencies_gradient_term_RHS(4, "psi, grad(phie)");
  // set_dependencies_value_term_LHS(4, "");
  // set_dependencies_gradient_term_LHS(4, "psi, grad(change(phie))");
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---
  scalarvalueType c = variable_list.get_scalar_value(0);

  scalargradType  mux = variable_list.get_scalar_gradient(1);

  scalarvalueType psi = variable_list.get_scalar_value(2);
  scalargradType  psix = variable_list.get_scalar_gradient(2);

  scalarvalueType phie = variable_list.get_scalar_value(3);

  scalarvalueType j0 = 24.0;
  scalarvalueType RT = 2.48e3;
  scalarvalueType F = 9.65e4;
  scalarvalueType hc = c * c * (3.0 - 2.0 * c);
  scalarvalueType Bnc = 2.0 * j0 / (2.0 * RT) * phie * hc;

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType eq_c  = c + constV(userInputs.dtValue) * (McV * psix * mux + std::sqrt(psix * psix) * Bnc) / psi;
  scalargradType  eqx_c = constV(-McV * userInputs.dtValue) * mux;

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_c);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_c);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  scalarvalueType c  = variable_list.get_scalar_value(0);
  scalargradType  cx = variable_list.get_scalar_gradient(0);

  scalarvalueType psi  = variable_list.get_scalar_value(2);
  scalargradType  psix = variable_list.get_scalar_gradient(2);

  // scalargradType  philix = variable_list.get_scalar_gradient(3);

  scalarvalueType phie  = variable_list.get_scalar_value(3);
  scalargradType  phiex = variable_list.get_scalar_gradient(3);

  scalarvalueType Bnmu = 0.0;
  scalarvalueType eps = 30.0;
  scalarvalueType F = 9.65e4;
  scalarvalueType sigma = 2.0e-1;
  scalarvalueType j0 = 24.0;
  scalarvalueType RT = 2.48e3;
  scalarvalueType hc  = c * c * (3.0 - 2.0 * c);

  // --- Setting the expressions for the terms in the governing equations ---

  // The derivative of the local free energy
  scalarvalueType fcV = 30.0 * constV(0.5) * (c * (1.0 - c) * (1.0 - c) - c * c * (1.0 - c));

  // The terms for the governing equations
  scalarvalueType eq_mu  = fcV - constV(KcV) * (psix * cx + Bnmu) / psi;
  scalargradType  eqx_mu = constV(KcV) * cx;

  // scalarvalueType eq_phili  = -psix * phiex;
  // scalargradType  eqx_phili = -psi * philix;

  scalarvalueType eq_phie  = 1.0 / sigma * std::sqrt(psix * psix) * j0 * F / (2.0 * RT) * phie * hc;
  scalargradType  eqx_phie = -(1.0 - psi) * phiex;

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(1, eq_mu);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_mu);

  variable_list.set_scalar_value_term_RHS(2, psi);

  // variable_list.set_scalar_value_term_RHS(3, eq_phili);
  // variable_list.set_scalar_gradient_term_RHS(3, eqx_phili);

  variable_list.set_scalar_value_term_RHS(3, eq_phie);
  variable_list.set_scalar_gradient_term_RHS(3, eqx_phie);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{

  // if (this -> currentFieldIndex == 3){
  //   scalarvalueType psi  = variable_list.get_scalar_value(2);
  //   scalargradType Dphilix = variable_list.get_change_in_scalar_gradient(3);

  //   scalargradType eqx_phili = psi * Dphilix;

  //   variable_list.set_scalar_gradient_term_LHS(3,eqx_phili);
  // }
  
  if (this -> currentFieldIndex == 3){
    scalarvalueType psi  = variable_list.get_scalar_value(2);
    scalargradType Dphiex = variable_list.get_change_in_scalar_gradient(3);

    scalargradType eqx_phie = (1.0 - psi) * Dphiex;
  
    variable_list.set_scalar_gradient_term_LHS(3,eqx_phie); 
  }
}
