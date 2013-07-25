/* =========================================================================
   Copyright (c) 2010-2013, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/*
*
*   Benchmark:   Vector operations (vector.cpp and vector.cu are identical, the latter being required for compilation using CUDA nvcc)
*
*/


#define VIENNACL_DEBUG_ALL
#define VIENNACL_DEBUG_BUILD

#ifndef NDEBUG
 #define NDEBUG
#endif

#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/scheduler/execute.hpp"
#include "viennacl/generator/generate.hpp"

#include <iostream>
#include <vector>
#include "benchmark-utils.hpp"

using std::cout;
using std::cin;
using std::endl;


#define BENCHMARK_VECTOR_SIZE   2
#define BENCHMARK_RUNS          1000

void execute_generator(viennacl::generator::code_generator const & generator){
  Timer t;
  t.start();
  char* program_name = (char*)malloc(256*sizeof(char));
  generator.make_program_name(program_name);
  t.start();
  if(!viennacl::ocl::current_context().has_program(program_name)){
    std::string source_code = generator.make_program_string();
#ifdef VIENNACL_DEBUG_BUILD
    std::cout << "Building " << program_name << "..." << std::endl;
    std::cout << source_code << std::endl;
#endif
    viennacl::ocl::current_context().add_program(source_code, program_name);
  }
  viennacl::ocl::program & p = viennacl::ocl::current_context().get_program(program_name);
  t.start();
  std::list<viennacl::ocl::kernel*> kernels;
  generator.configure_program(p, kernels);
  double overhead = t.get();
  std::cout << "overhead : " << overhead << std::endl;
  for(std::list<viennacl::ocl::kernel*>::iterator it = kernels.begin() ; it != kernels.end() ; ++it){
    viennacl::ocl::enqueue(**it);
  }
  viennacl::backend::finish();
}

template<typename ScalarType>
int run_benchmark()
{

  Timer timer;
  double exec_time;

  ScalarType std_result = 0;

  std::vector<ScalarType> std_vec1(BENCHMARK_VECTOR_SIZE);
  std::vector<ScalarType> std_vec2(BENCHMARK_VECTOR_SIZE);
  std::vector<ScalarType> std_vec3(BENCHMARK_VECTOR_SIZE);
  viennacl::scalar<ScalarType> vcl_scal1(0);
  viennacl::scalar<ScalarType> vcl_scal2(0);
  viennacl::vector<ScalarType> vcl_vec1(BENCHMARK_VECTOR_SIZE);
  viennacl::vector<ScalarType> vcl_vec2(BENCHMARK_VECTOR_SIZE);
  viennacl::vector<ScalarType> vcl_vec3(BENCHMARK_VECTOR_SIZE);
  viennacl::matrix<ScalarType> vcl_mat1(BENCHMARK_VECTOR_SIZE, BENCHMARK_VECTOR_SIZE);
  viennacl::matrix<ScalarType> vcl_mat2(BENCHMARK_VECTOR_SIZE, BENCHMARK_VECTOR_SIZE);
  viennacl::matrix<ScalarType> vcl_mat3(BENCHMARK_VECTOR_SIZE, BENCHMARK_VECTOR_SIZE);



  ///////////// Vector operations /////////////////

  std_vec1[0] = 1.0;
  std_vec2[0] = 1.0;
  for (int i=1; i<BENCHMARK_VECTOR_SIZE; ++i)
  {
    std_vec1[i] = std_vec1[i-1] * ScalarType(1.000001);
    std_vec2[i] = std_vec1[i-1] * ScalarType(0.999999);
  }

  viennacl::copy(std_vec1, vcl_vec1);
  viennacl::fast_copy(std_vec1, vcl_vec1);
  viennacl::copy(std_vec2, vcl_vec2);

  viennacl::backend::finish();

  {
    viennacl::generator::code_generator generator;
    generator.add_statement(viennacl::scheduler::statement(vcl_vec1, viennacl::op_assign(), viennacl::linalg::prod(vcl_mat1, vcl_vec2)));
//    generator.add_statement(viennacl::scheduler::statement(vcl_vec2, viennacl::op_assign(), viennacl::linalg::prod(vcl_mat1, vcl_vec3)));

//    generator.add_statement(viennacl::scheduler::statement(vcl_vec1, viennacl::op_assign(), vcl_vec1 + vcl_vec2));
//    generator.add_statement(viennacl::scheduler::statement(vcl_vec1, viennacl::op_assign(), vcl_vec1 + vcl_vec3));

    execute_generator(generator);
    execute_generator(generator);
  }

  return 0;
}

int main()
{
  std::cout << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "               Device Info" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;

#ifdef VIENNACL_WITH_OPENCL
  std::cout << viennacl::ocl::current_device().info() << std::endl;
#endif

  std::cout << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "## Benchmark :: Vector" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << "   -------------------------------" << std::endl;
  std::cout << "   # benchmarking single-precision" << std::endl;
  std::cout << "   -------------------------------" << std::endl;
  run_benchmark<float>();
  return 0;
}

