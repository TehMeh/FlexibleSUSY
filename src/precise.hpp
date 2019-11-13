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

/*typedef boost::multiprecision::number<boost::multiprecision::mpc_complex_backend<25>> precise_complex_type;

typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<25>> precise_real_type;*/

/*typedef std::complex<long double> precise_complex_type;

typedef long double precise_real_type;

using namespace std;*/


// User defined numeric literal; correctly converts numbers to precise types;! 
precise_real_type operator ""_p(const char* a);
/*precise_real_type operator ""_p(const char* a){
	precise_real_type res=0;
	std::string number=a;

	for(size_t i=number.size()-1; i>=0; i--) {
		if(number.substr(i,1).compare("0")!=0 ){
			number=number.substr(0,i+1);
			break;
		}
	}
	std::size_t found=number.find(".");

	std::size_t found2=number.find("e");
	precise_real_type fac=1;

	if(found2!=std::string::npos){
		fac=precise_real_type(std::stoi(number.substr(found2+1,number.size()-found2-1)));
		fac=pow(precise_real_type(10),fac);
		number=number.substr(0,found2);
	}

	if(found!=std::string::npos){
		for(int i=number.substr(0,found).size()-1; i>=0; i--){
			res+=precise_real_type(std::stoi(number.substr(i,1)))*pow(precise_real_type(10), number.substr(0,found).size()-1-i);
		}
		for(int i=0; i<=number.substr(found+1).size()-1;i++){
			res+=precise_real_type(std::stoi(number.substr(found+1+i,1)))*pow(precise_real_type(10),-i-1);
		}
	}
	else{
		for(int i=number.size()-1; i>=0; i--){
			res+=precise_real_type(std::stoi(number.substr(i,1)))*pow(precise_real_type(10), number.size()-1-i);
		}
	}
	return res;
}*/

namespace Eigen{
	typedef Array<precise_real_type, Dynamic, 1> ArrayXdp;
	typedef Matrix< precise_real_type, Dynamic, 1 > VectorXdp;
	typedef Matrix< precise_complex_type, Dynamic, 1 > VectorXcdp;
	typedef Matrix< precise_real_type, Dynamic, Dynamic > MatrixXdp;
	typedef Matrix< precise_complex_type, Dynamic, Dynamic > MatrixXcdp;
}
#endif



