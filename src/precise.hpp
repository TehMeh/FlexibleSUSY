#ifndef PRECISE_H
#define PRECISE_H


#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/eigen.hpp>
//#include <complex>
#include <Eigen/Dense>

#include <string> 

#include <boost/type_traits.hpp>
#include <boost/functional.hpp>
#include <boost/function.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/exception/to_string.hpp>
#include <boost/compute/container/valarray.hpp>

#include "boost/filesystem.hpp"
#include <boost/filesystem/fstream.hpp>

namespace boost {
   namespace multiprecision{
      typedef number<mpc_complex_backend<30>> precise_complex_type;

      typedef number<mpfr_float_backend<30>> precise_real_type;   
   }
} 

using namespace boost::multiprecision;

// User defined numeric literal; correctly converts numbers to precise types;! 
precise_real_type operator ""_p(const char* a);

namespace Eigen{
	typedef Array<precise_real_type, Dynamic, 1> ArrayXdp;
	typedef Matrix< precise_real_type, Dynamic, 1 > VectorXdp;
	typedef Matrix< precise_complex_type, Dynamic, 1 > VectorXcdp;
	typedef Matrix< precise_real_type, Dynamic, Dynamic > MatrixXdp;
	typedef Matrix< precise_complex_type, Dynamic, Dynamic > MatrixXcdp;
}
#endif



