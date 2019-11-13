// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "precise.hpp"

#ifndef MATHLINK_UTILS_H
#define MATHLINK_UTILS_H

#include <mathlink.h>

#include <complex>
#include <string>
#include <vector>
#include <Eigen/Core>

/********************* put types *********************/

inline void MLPut(MLINK link, const std::string& s)
{
   MLPutSymbol(link, s.c_str());
}

inline void MLPut(MLINK link, int c)
{
   MLPutInteger(link, c);
}

inline void MLPut(MLINK link, precise_real_type c)
{
   MLPutReal(link, (double)(c));
}

inline void MLPut(MLINK link, precise_complex_type c)
{
   if (imag(c) == 0.) {
      MLPutReal(link, (double)(real(c)));
   } else {
      MLPutFunction(link, "Complex", 2);
      MLPutReal(link, (double)(real(c)));
      MLPutReal(link, (double)(imag(c)));
   }
}

template <int M>
void MLPut(MLINK link, const Eigen::Array<precise_real_type,M,1>& a)
{
   precise_real_type v[M];
   for (int i = 0; i < M; i++)
      v[i] = a(i);
   MLPutRealList(link, (const double*)(v), M);
}

template <int M>
void MLPut(MLINK link, const Eigen::Matrix<precise_real_type,M,1>& m)
{
   const Eigen::Array<precise_real_type,M,1> a(m.array());
   MLPut(link, a);
}

template <int M, int N>
void MLPut(MLINK link, const Eigen::Matrix<precise_real_type,M,N>& m)
{
   precise_real_type mat[M][N];
   for (int i = 0; i < M; i++)
      for (int k = 0; k < N; k++)
         mat[i][k] = m(i, k);

   long dims[] = { M, N };
   MLPutDoubleArray(link, (double*)mat, dims, NULL, 2);
}

template <int M>
void MLPut(MLINK link, const Eigen::Array<precise_complex_type,M,1>& a)
{
   MLPutFunction(link, "List", M);
   for (int i = 0; i < M; i++)
      MLPut(link, a(i));
}

template <int M>
void MLPut(MLINK link, const Eigen::Matrix<precise_complex_type,M,1>& m)
{
   const Eigen::Array<precise_complex_type,M,1> a(m.array());
   MLPut(link, a);
}

template <int M, int N>
void MLPut(MLINK link, const Eigen::Matrix<precise_complex_type,M,N>& m)
{
   MLPutFunction(link, "List", M);
   for (int i = 0; i < M; i++) {
      MLPutFunction(link, "List", N);
      for (int k = 0; k < N; k++)
         MLPut(link, m(i,k));
   }
}

/********************* put single heads *********************/

inline void MLPutHeads(MLINK link, const std::vector<std::string>& heads)
{
   for (const auto& h: heads)
      MLPutFunction(link, h.c_str(), 1);
}

/********************* put rules to types *********************/

inline void MLPutRule(MLINK link, const std::string& name, const std::vector<std::string>& heads = {})
{
   MLPutFunction(link, "Rule", 2);
   MLPutHeads(link, heads);
   MLPutUTF8Symbol(link, reinterpret_cast<const unsigned char*>(name.c_str()), name.size());
}

inline void MLPutRule(MLINK link, int number, const std::vector<std::string>& heads = {})
{
   MLPutFunction(link, "Rule", 2);
   MLPutHeads(link, heads);
   MLPutInteger(link, number);
}

inline void MLPutRule(MLINK link, long number, const std::vector<std::string>& heads = {})
{
   MLPutFunction(link, "Rule", 2);
   MLPutHeads(link, heads);
   MLPutLongInteger(link, number);
}

template <class T1, class T2>
void MLPutRuleTo(MLINK link, T1 t, const T2& name, const std::vector<std::string>& heads = {})
{
   MLPutRule(link, name, heads);
   MLPut(link, t);
}

/********************* get types *********************/

inline void MLGet(MLINK link, int *c)
{
   MLGetInteger(link, c);
}

inline void MLGet(MLINK link, long *c)
{
   MLGetLongInteger(link, c);
}

inline void MLGet(MLINK link, short *c)
{
   MLGetShortInteger(link, c);
}

#endif // MATHLINK_UTILS_H
