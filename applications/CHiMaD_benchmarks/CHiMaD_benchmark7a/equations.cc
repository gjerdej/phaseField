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
variableAttributeLoader::loadVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n");
  set_dependencies_gradient_term_RHS(0, "grad(n)");
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

  // The order parameter and its derivatives
  scalarvalueType n  = variable_list.get_scalar_value(0);
  scalargradType  nx = variable_list.get_scalar_gradient(0);

  // --- Setting the expressions for the terms in the governing equations ---

  scalarvalueType fnV   = (4.0 * n * (n - 1.0) * (n - 0.5));
  scalarvalueType eq_n  = (n - constV(userInputs.dtValue) * fnV);
  scalargradType  eqx_n = (-constV(userInputs.dtValue * kappa) * nx);

  scalarvalueType source_term;

  scalarvalueType alpha = 0.25 + A1 * this->currentTime * std::sin(B1 * q_point_loc(0)) +
                          A2 * std::sin(B2 * q_point_loc(0) + C2 * this->currentTime);
  scalarvalueType alpha_t =
    A1 * std::sin(B1 * q_point_loc(0)) +
    A2 * C2 * std::cos(B2 * q_point_loc(0) + C2 * this->currentTime);
  scalarvalueType alpha_y =
    A1 * B1 * this->currentTime * std::cos(B1 * q_point_loc(0)) +
    A2 * B2 * std::cos(B2 * q_point_loc(0) + C2 * this->currentTime);
  scalarvalueType alpha_yy =
    -A1 * B1 * B1 * this->currentTime * std::sin(B1 * q_point_loc(0)) -
    A2 * B2 * B2 * std::sin(B2 * q_point_loc(0) + C2 * this->currentTime);

  for (unsigned i = 0; i < n.size(); i++)
    {
      source_term[i] =
        (-2.0 * std::sqrt(kappa) *
           std::tanh((q_point_loc(1)[i] - alpha[i]) / std::sqrt(2.0 * kappa)) *
           (alpha_y[i] * alpha_y[i]) +
         std::sqrt(2.0) * (alpha_t[i] - kappa * alpha_yy[i])) /
        (4.0 * std::sqrt(kappa)) /
        Utilities::fixed_power<2>(
          std::cosh((q_point_loc(1)[i] - alpha[i]) / std::sqrt(2.0 * kappa)));
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(0, eq_n + userInputs.dtValue * source_term);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_n);
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
{}

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
{}