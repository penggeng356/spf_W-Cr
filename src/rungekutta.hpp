/* ----------------------------------------------------------------------
    SPF - Stochastic Phase Field
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the README file in the top-level SPF directory.
---------------------------------------------------------------------- */
// File: rungekutta.hpp

#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

namespace SPF_NS
{
   double marcusRK4_increment( double (*integrand)(const double& theta, 
                                          const double& bb,
                                          const double& poissonJump), 
                     const double& dtheta,
                     const double& b0,
                     const double& theta0,
                     const double& jump);

   double marcusRK4_increment( double (*integrand)(const double& theta, 
                                          const double& LL,
                                          const double& bb,
                                          const double& poissonJump), 
                     const double& dtheta,
                     const double& L0, // state of driving Poisson process
                     const double& b0,
                     const double& theta0,
                     const double& jump);

   double marcusRK4_increment( double (*integrand)(
                                          const double& tau, 
                                          const double& theta, 
                                          const double& LL,
                                          const double& bb,
                                          const double& poissonJump), 
                     const double& tau,
                     const double& dtheta,
                     const double& L0, // state of driving Poisson process
                     const double& b0,
                     const double& theta0,
                     const double& jump);

   double marcusRK4_increment(
                     double (*bIntegrand)( const double& b0,
                                           const double& LL,
                                           const double& tau,
                                           const double& theta,
                                           const double& Yk
                                           ),
                     double (*integrand)(
                                          const double& bIntegrand,
                                          const double& tau, 
                                          const double& theta, 
                                          const double& LL,
                                          const double& poissonJump), 
                     const double& tau,
                     const double& dtheta,
                     const double& L0, // state of driving Poisson process
                     const double& b0,
                     const double& theta0,
                     const double& jump);


   //double marcusRK4_increment( double (*bIntegrand)(
   //                              double (*zTauTheta)(
   //                                       const double& theta,
   //                                       ),
   //                                       const double& theta, 
   //                                       const double& LL,
   //                                       const double& bb,
   //                                       const double& poissonJump), 
   //                  const double& dtheta,
   //                  const double& L0, // state of driving Poisson process
   //                  const double& b0,
   //                  const double& theta0,
   //                  const double& jump);


   int marcusIntegral( double (*marcusIntegrand)(const double& theta,
                                          const double& bb,
                                          const double& dP), 
                     const double& rk_dtheta,    // increment of unit interval
                     const double& b0,    // path.back()
                     const double& theta0,    // t_n, t_{n+1}=t_n + dt
                     const double& dP, // Y_{k} size of jump before marcus
                     double& jumpDestination);
}

#endif
