// All of the methods for the 'variableContainer' class
#include "../../include/variableContainer.h"

template <int dim, int degree, typename T>
variableContainer<dim, degree, T>::variableContainer(
  const dealii::MatrixFree<dim, double> &data,
  std::vector<variable_info>             _varInfoList,
  std::vector<variable_info>             _varChangeInfoList)
{
  varInfoList       = _varInfoList;
  varChangeInfoList = _varChangeInfoList;

  num_var = varInfoList.size();

  for (unsigned int i = 0; i < num_var; i++)
    {
      if (varInfoList[i].var_needed)
        {
          if (varInfoList[i].is_scalar)
            {
              dealii::FEEvaluation<dim, degree, degree + 1, 1, double> var(data, i);
              scalar_vars.push_back(var);
            }
          else
            {
              dealii::FEEvaluation<dim, degree, degree + 1, dim, double> var(data, i);
              vector_vars.push_back(var);
            }
        }

      if (varChangeInfoList[i].var_needed)
        {
          if (varChangeInfoList[i].is_scalar)
            {
              dealii::FEEvaluation<dim, degree, degree + 1, 1, double> var(data, i);
              scalar_change_in_vars.push_back(var);
            }
          else
            {
              dealii::FEEvaluation<dim, degree, degree + 1, dim, double> var(data, i);
              vector_change_in_vars.push_back(var);
            }
        }
    }
}

template <int dim, int degree, typename T>
variableContainer<dim, degree, T>::variableContainer(
  const dealii::MatrixFree<dim, double> &data,
  std::vector<variable_info>             _varInfoList)
{
  varInfoList = _varInfoList;

  num_var = varInfoList.size();

  for (unsigned int i = 0; i < num_var; i++)
    {
      if (varInfoList[i].var_needed)
        {
          if (varInfoList[i].is_scalar)
            {
              dealii::FEEvaluation<dim, degree, degree + 1, 1, double> var(data, i);
              scalar_vars.push_back(var);
            }
          else
            {
              dealii::FEEvaluation<dim, degree, degree + 1, dim, double> var(data, i);
              vector_vars.push_back(var);
            }
        }
    }
}

// Variant of the constructor where it reads from a fixed index of "data", used
// for post-processing
template <int dim, int degree, typename T>
variableContainer<dim, degree, T>::variableContainer(
  const dealii::MatrixFree<dim, double> &data,
  std::vector<variable_info>             _varInfoList,
  unsigned int                           fixed_index)
{
  varInfoList = _varInfoList;

  num_var = varInfoList.size();

  for (unsigned int i = 0; i < num_var; i++)
    {
      if (varInfoList[i].var_needed)
        {
          if (varInfoList[i].is_scalar)
            {
              dealii::FEEvaluation<dim, degree, degree + 1, 1, double> var(data,
                                                                           fixed_index);
              scalar_vars.push_back(var);
            }
          else
            {
              dealii::FEEvaluation<dim, degree, degree + 1, dim, double> var(data,
                                                                             fixed_index);
              vector_vars.push_back(var);
            }
        }
    }
}

template <int dim, int degree, typename T>
unsigned int
variableContainer<dim, degree, T>::get_num_q_points()
{
  if (scalar_vars.size() > 0)
    {
      return scalar_vars[0].n_q_points;
    }
  else if (vector_vars.size() > 0)
    {
      return vector_vars[0].n_q_points;
    }
  else if (scalar_change_in_vars.size() > 0)
    {
      return scalar_change_in_vars[0].n_q_points;
    }
  else
    {
      return vector_change_in_vars[0].n_q_points;
    }
}

template <int dim, int degree, typename T>
dealii::Point<dim, T>
variableContainer<dim, degree, T>::get_q_point_location()
{
  if (scalar_vars.size() > 0)
    {
      return scalar_vars[0].quadrature_point(q_point);
    }
  else if (vector_vars.size() > 0)
    {
      return vector_vars[0].quadrature_point(q_point);
    }
  else if (scalar_change_in_vars.size() > 0)
    {
      return scalar_change_in_vars[0].quadrature_point(q_point);
    }
  else
    {
      return vector_change_in_vars[0].quadrature_point(q_point);
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit_and_eval(const std::vector<vectorType *> &src,
                                                   unsigned int                     cell)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      if (varInfoList[i].var_needed)
        {
          if (varInfoList[i].is_scalar)
            {
              scalar_vars[varInfoList[i].variable_index].reinit(cell);
              scalar_vars[varInfoList[i].variable_index].read_dof_values(*src[i]);
              scalar_vars[varInfoList[i].variable_index].evaluate(
                varInfoList[i].evaluation_flags);
            }
          else
            {
              vector_vars[varInfoList[i].variable_index].reinit(cell);
              vector_vars[varInfoList[i].variable_index].read_dof_values(*src[i]);
              vector_vars[varInfoList[i].variable_index].evaluate(
                varInfoList[i].evaluation_flags);
            }
        }
    }
}

/**
 * This is specialized for the LHS where there will be only one change in the
 * solution needed. The RHS method takes the src as a vector of vectorTypes.
 */
template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit_and_eval_change_in_solution(
  const vectorType &src,
  unsigned int      cell,
  unsigned int      var_being_solved)
{
  if (varChangeInfoList[var_being_solved].is_scalar)
    {
      scalar_change_in_vars[0].reinit(cell);
      scalar_change_in_vars[0].read_dof_values(src);
      scalar_change_in_vars[0].evaluate(
        varChangeInfoList[var_being_solved].evaluation_flags);
    }
  else
    {
      vector_change_in_vars[0].reinit(cell);
      vector_change_in_vars[0].read_dof_values(src);
      vector_change_in_vars[0].evaluate(
        varChangeInfoList[var_being_solved].evaluation_flags);
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit_and_eval_LHS(
  const vectorType               &src,
  const std::vector<vectorType *> solutionSet,
  unsigned int                    cell,
  unsigned int                    var_being_solved)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      if (varInfoList[i].var_needed)
        {
          if (varInfoList[i].is_scalar)
            {
              scalar_vars[varInfoList[i].variable_index].reinit(cell);
              if (i == var_being_solved)
                {
                  scalar_vars[varInfoList[i].variable_index].read_dof_values(src);
                }
              else
                {
                  scalar_vars[varInfoList[i].variable_index].read_dof_values(
                    *solutionSet[i]);
                }
              scalar_vars[varInfoList[i].variable_index].evaluate(
                varInfoList[i].evaluation_flags);
            }
          else
            {
              vector_vars[varInfoList[i].variable_index].reinit(cell);
              if (i == var_being_solved)
                {
                  vector_vars[varInfoList[i].variable_index].read_dof_values(src);
                }
              else
                {
                  vector_vars[varInfoList[i].variable_index].read_dof_values(
                    *solutionSet[i]);
                }
              vector_vars[varInfoList[i].variable_index].evaluate(
                varInfoList[i].evaluation_flags);
            }
        }
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::reinit(unsigned int cell)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      if (varInfoList[i].var_needed)
        {
          if (varInfoList[i].is_scalar)
            {
              scalar_vars[varInfoList[i].variable_index].reinit(cell);
            }
          else
            {
              vector_vars[varInfoList[i].variable_index].reinit(cell);
            }
        }
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::integrate_and_distribute(
  std::vector<vectorType *> &dst)
{
  for (unsigned int i = 0; i < num_var; i++)
    {
      if (!(varInfoList[i].residual_flags & dealii::EvaluationFlags::nothing))
        {
          if (varInfoList[i].is_scalar)
            {
              scalar_vars[varInfoList[i].variable_index].integrate(
                varInfoList[i].residual_flags);
              scalar_vars[varInfoList[i].variable_index].distribute_local_to_global(
                *dst[i]);
            }
          else
            {
              vector_vars[varInfoList[i].variable_index].integrate(
                varInfoList[i].residual_flags);
              vector_vars[varInfoList[i].variable_index].distribute_local_to_global(
                *dst[i]);
            }
        }
    }
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::integrate_and_distribute_change_in_solution_LHS(
  vectorType        &dst,
  const unsigned int var_being_solved)
{
  // integrate
  if (varChangeInfoList[var_being_solved].is_scalar)
    {
      scalar_change_in_vars[0].integrate(
        varChangeInfoList[var_being_solved].residual_flags);
      scalar_change_in_vars[0].distribute_local_to_global(dst);
    }
  else
    {
      vector_change_in_vars[0].integrate(
        varChangeInfoList[var_being_solved].residual_flags);
      vector_change_in_vars[0].distribute_local_to_global(dst);
    }
}

// Need to add index checking to these functions so that an error is thrown if
// the index wasn't set
template <int dim, int degree, typename T>
T
variableContainer<dim, degree, T>::get_scalar_value(
  unsigned int global_variable_index) const
{
  if (varInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::values)
    {
      return scalar_vars[varInfoList[global_variable_index].variable_index].get_value(
        q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a variable value that "
                   "was not marked as needed in 'equations.cc'. The attempted "
                   "access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_scalar_gradient(
  unsigned int global_variable_index) const
{
  if (varInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients)
    {
      return scalar_vars[varInfoList[global_variable_index].variable_index].get_gradient(
        q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a variable gradient "
                   "that was not marked as needed in 'equations.cc'. The "
                   "attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_scalar_hessian(
  unsigned int global_variable_index) const
{
  if (varInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians)
    {
      return scalar_vars[varInfoList[global_variable_index].variable_index].get_hessian(
        q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a variable hessian "
                   "that was not marked as needed in 'equations.cc'. The "
                   "attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_vector_value(
  unsigned int global_variable_index) const
{
  if (varInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::values)
    {
      return vector_vars[varInfoList[global_variable_index].variable_index].get_value(
        q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a variable value that "
                   "was not marked as needed in 'equations.cc'. The attempted "
                   "access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_vector_gradient(
  unsigned int global_variable_index) const
{
  if (varInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients)
    {
      return vector_vars[varInfoList[global_variable_index].variable_index].get_gradient(
        q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a variable gradient "
                   "that was not marked as needed in 'equations.cc'. The "
                   "attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<3, dim, T>
variableContainer<dim, degree, T>::get_vector_hessian(
  unsigned int global_variable_index) const
{
  if (varInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians)
    {
      return vector_vars[varInfoList[global_variable_index].variable_index].get_hessian(
        q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a variable hessian "
                   "that was not marked as needed in 'equations.cc'. The "
                   "attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

// Need to add index checking to these functions so that an error is thrown if
// the index wasn't set
template <int dim, int degree, typename T>
T
variableContainer<dim, degree, T>::get_change_in_scalar_value(
  unsigned int global_variable_index) const
{
  if (varChangeInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::values)
    {
      return scalar_change_in_vars[varChangeInfoList[global_variable_index]
                                     .variable_index]
        .get_value(q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a change in variable "
                   "value that was not marked as needed in 'equations.cc'. The "
                   "attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_change_in_scalar_gradient(
  unsigned int global_variable_index) const
{
  if (varChangeInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients)
    {
      return scalar_change_in_vars[varChangeInfoList[global_variable_index]
                                     .variable_index]
        .get_gradient(q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a change in variable "
                   "gradient that was not marked as needed in 'equations.cc'. "
                   "The attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_change_in_scalar_hessian(
  unsigned int global_variable_index) const
{
  if (varChangeInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians)
    {
      return scalar_change_in_vars[varChangeInfoList[global_variable_index]
                                     .variable_index]
        .get_hessian(q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a change in variable "
                   "hessian that was not marked as needed in 'equations.cc'. "
                   "The attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<1, dim, T>
variableContainer<dim, degree, T>::get_change_in_vector_value(
  unsigned int global_variable_index) const
{
  if (varChangeInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::values)
    {
      return vector_change_in_vars[varChangeInfoList[global_variable_index]
                                     .variable_index]
        .get_value(q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a change in variable "
                   "value that was not marked as needed in 'equations.cc'. The "
                   "attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<2, dim, T>
variableContainer<dim, degree, T>::get_change_in_vector_gradient(
  unsigned int global_variable_index) const
{
  if (varChangeInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::gradients)
    {
      return vector_change_in_vars[varChangeInfoList[global_variable_index]
                                     .variable_index]
        .get_gradient(q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a change in variable "
                   "gradient that was not marked as needed in 'equations.cc'. "
                   "The attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

template <int dim, int degree, typename T>
dealii::Tensor<3, dim, T>
variableContainer<dim, degree, T>::get_change_in_vector_hessian(
  unsigned int global_variable_index) const
{
  if (varChangeInfoList[global_variable_index].evaluation_flags &
      dealii::EvaluationFlags::hessians)
    {
      return vector_change_in_vars[varChangeInfoList[global_variable_index]
                                     .variable_index]
        .get_hessian(q_point);
    }
  else
    {
      std::cerr << "PRISMS-PF Error: Attempted access of a change in variable "
                   "hessian that was not marked as needed in 'equations.cc'. "
                   "The attempted access was for variable with index "
                << global_variable_index << " ." << std::endl;
      abort();
    }
}

// The methods to set the residual terms
template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_value_term_RHS(
  unsigned int global_variable_index,
  T            val)
{
  scalar_vars[varInfoList[global_variable_index].variable_index].submit_value(val,
                                                                              q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_gradient_term_RHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> grad)
{
  scalar_vars[varInfoList[global_variable_index].variable_index].submit_gradient(grad,
                                                                                 q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_value_term_RHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> val)
{
  vector_vars[varInfoList[global_variable_index].variable_index].submit_value(val,
                                                                              q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_gradient_term_RHS(
  unsigned int              global_variable_index,
  dealii::Tensor<2, dim, T> grad)
{
  vector_vars[varInfoList[global_variable_index].variable_index].submit_gradient(grad,
                                                                                 q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_value_term_LHS(
  unsigned int global_variable_index,
  T            val)
{
  scalar_change_in_vars[varChangeInfoList[global_variable_index].variable_index]
    .submit_value(val, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_scalar_gradient_term_LHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> grad)
{
  scalar_change_in_vars[varChangeInfoList[global_variable_index].variable_index]
    .submit_gradient(grad, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_value_term_LHS(
  unsigned int              global_variable_index,
  dealii::Tensor<1, dim, T> val)
{
  vector_change_in_vars[varChangeInfoList[global_variable_index].variable_index]
    .submit_value(val, q_point);
}

template <int dim, int degree, typename T>
void
variableContainer<dim, degree, T>::set_vector_gradient_term_LHS(
  unsigned int              global_variable_index,
  dealii::Tensor<2, dim, T> grad)
{
  vector_change_in_vars[varChangeInfoList[global_variable_index].variable_index]
    .submit_gradient(grad, q_point);
}

template class variableContainer<2, 1, dealii::VectorizedArray<double>>;
template class variableContainer<2, 2, dealii::VectorizedArray<double>>;
template class variableContainer<2, 3, dealii::VectorizedArray<double>>;
template class variableContainer<2, 4, dealii::VectorizedArray<double>>;
template class variableContainer<3, 1, dealii::VectorizedArray<double>>;
template class variableContainer<3, 2, dealii::VectorizedArray<double>>;
template class variableContainer<3, 3, dealii::VectorizedArray<double>>;
template class variableContainer<3, 4, dealii::VectorizedArray<double>>;
