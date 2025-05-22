// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index
  
  double int_width = 4*(userInputs.domain_size[0] / ((double) userInputs.subdivisions[0]) /
    std::pow(2.0, userInputs.refine_factor));

  if (index == 0){
    
        // int_width = 3*(userInputs.domain_size[0] / ((double) userInputs.subdivisions[0]) /
        //   std::pow(2.0, userInputs.refine_factor));
        // dealii::Utilities::fixed_power<2>()
    scalar_IC = 1 - 0.5*(1.0+std::tanh(-(userInputs.domain_size[1]/2.0 - p[1])/(int_width/2))) * 0.5 * (1.0+std::tanh((dealii::Utilities::fixed_power<2>(p[0] - userInputs.domain_size[0]/2.0)
      + dealii::Utilities::fixed_power<2>(p[1] - userInputs.domain_size[1]/2.0) - dealii::Utilities::fixed_power<2>(userInputs.domain_size[0]/6.0))/(30*int_width/2.0)));
      // * 0.5 * (1.0+std::tanh((dealii::Utilities::fixed_power<2>(p[0] - 2*userInputs.domain_size[0]/3.0) + dealii::Utilities::fixed_power<2>(p[1] 
      //   - userInputs.domain_size[1]/2.0) - dealii::Utilities::fixed_power<2>(userInputs.domain_size[0]/8.0))/(30*int_width/2.0)));
      // * 0.5 * (1.0+std::tanh((dealii::Utilities::fixed_power<2>(p[0] - 3*userInputs.domain_size[0]/4.0) + dealii::Utilities::fixed_power<2>(p[1] 
      //   - userInputs.domain_size[1]/2.0) - dealii::Utilities::fixed_power<2>(userInputs.domain_size[0]/10.0))/(10*int_width/2.0)));
    
    // double r = std::sqrt((p[0]-userInputs.domain_size[0]/2.0)*(p[0]-userInputs.domain_size[0]/2.0)+(p[1]-userInputs.domain_size[1]/2.0)*(p[1]-userInputs.domain_size[1]/2.0));
    // scalar_IC = 0.5*(1.0+std::tanh((r-std::sqrt(userInputs.domain_size[0]/4.0))/(0.5*int_width)));
  }

  if (index == 1){
    // 4/int_width**2*((np.tanh(2/int_width*(y-ysize/2)))**3 - np.tanh(2/int_width*(y-ysize/2)))
    // scalarvalueType c = 0.5*(1.0+std::tanh((userInputs.domain_size[1]/2.0 - p[1])/(int_width/2)));
    scalar_IC = 0.0; //constV(0.25) * (constV(2.0) * std::pow(2.0, 0.5*(1.0+std::tanh((userInputs.domain_size[1]/2.0 - p[1])/(int_width/2)))) - constV(4.0) * std::pow(3.0, 0.5*(1.0+std::tanh((userInputs.domain_size[1]/2.0 - p[1])/(int_width/2))))) - constV(1.0e-3) * 4/std::pow(2.0, int_width) * (std::pow(3.0, std::tanh((userInputs.domain_size[1]/2.0 - p[1])/(int_width/2))) - std::tanh((userInputs.domain_size[1]/2.0 - p[1])/(int_width/2)));
  }

  if (index == 2){
    // scalar_IC = 0.5*(1.0+std::tanh(-(p[1]-userInputs.domain_size[1]/2.0)/(int_width/2))) + 1.0e-6;
    scalar_IC = 1 - 0.5*(1.0+std::tanh(-(userInputs.domain_size[1]/2.0 - p[1])/(int_width/2))) * 0.5 * (1.0+std::tanh((dealii::Utilities::fixed_power<2>(p[0] - userInputs.domain_size[0]/2.0)
      + dealii::Utilities::fixed_power<2>(p[1] - userInputs.domain_size[1]/2.0) - dealii::Utilities::fixed_power<2>(userInputs.domain_size[0]/6.0))/(30*int_width/2.0)));
  }

  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(
  [[maybe_unused]] const Point<dim>  &p,
  [[maybe_unused]] const unsigned int index,
  [[maybe_unused]] const unsigned int direction,
  [[maybe_unused]] const double       time,
  [[maybe_unused]] double            &scalar_BC,
  [[maybe_unused]] Vector<double>    &vector_BC)
{
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
  // --------------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the boundary condition for each variable
  // according to its variable index. This function can be left blank if there
  // are no non-uniform Dirichlet boundary conditions. For BCs that change in
  // time, you can access the current time through the variable "time". The
  // boundary index can be accessed via the variable "direction", which starts
  // at zero and uses the same order as the BC specification in parameters.in
  // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).

  // -------------------------------------------------------------------------
}
