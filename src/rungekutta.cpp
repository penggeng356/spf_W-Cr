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
// File: rungekutta.cpp

#ifndef RUNGEKUTTA_CPP
#define RUNGEKUTTA_CPP

#include "rungekutta.hpp"

#include <iostream> // cout endl;
//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

#ifndef ONESIXTH 
#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666667
#endif

int SPF_NS::marcusIntegral(
                      double (*marcusIntegrand)(const double& tt,
                                          const double& yy,
                                          const double& dP), 
                     const double& rk_dt,    // increment of unit interval
                     const double& y0,    // path.back()
                     const double& t0,    // t_n, t_{n+1}=t_n + dt
                     const double& dP,// Y_{k} size of noise jump
                     double& jumpDestination)
{
   jumpDestination = y0;
   if (( rk_dt >= 1.0) || (rk_dt <= 0))
   {
      std::cout << "Error: runge-kutta increment rk_dt >=1 or <=0" 
               << std::endl;
      return EXIT_FAILURE;
   }
   for( double tt=0.0; tt < 1.0; tt += rk_dt ) 
   {
      jumpDestination += marcusRK4_increment( marcusIntegrand,
                                    rk_dt,
                                    jumpDestination, //accumulates over[0,1]
                                    tt, // varies from 0 to 1
                                    dP);
   }
   return EXIT_SUCCESS;
}

double SPF_NS::marcusRK4_increment( // returns y0 + contribution of integrand over dt // \Phi_{Y_k}
                     double (*integrand)(const double& theta, 
                                          const double& bb,
                                          const double& dP),// Y_k
                     const double& rk_dtheta, // constant
                     const double& b0, // will be reassigned the return value of this function, accumulates value of the integral as t0 varies from 0 to 1
                     const double& theta0, // will vary from 0 to 1
                     const double& dP) // will be constant while integrating theta0 over [0,1]
{
   // Returns a contribution of the integral of integrand() over [theta0,theta0+rk_dtheta] added to y0.
   double k1, k2, k3, k4; 
   k1 = dP * rk_dtheta * integrand( theta0, b0, dP); // == y0*jump
   k2 = dP * rk_dtheta * integrand( theta0 + 0.5*rk_dtheta, b0 + 0.5*k1, dP);
   k3 = dP * rk_dtheta * integrand( theta0 + 0.5*rk_dtheta, b0 + 0.5*k2, dP);
   k4 = dP * rk_dtheta * integrand( theta0 + rk_dtheta, b0 + k3, dP);

   //std::cout << "t0: " << t0 << ", y0:" << y0  << ", rk_dt: " << rk_dt << ", jump: " << jump << ", k's: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << std::endl;

   return ONESIXTH *(k1 + 2.0*k2 + 2.0*k3 + k4 );
}

double SPF_NS::marcusRK4_increment( // returns y0 + contribution of integrand over dt // \Phi_{Y_k}
                     double (*integrand)(const double& theta,
                                          const double& LL,
                                          const double& bb,
                                          const double& dP),// Y_k
                     const double& rk_dtheta, // constant
                     const double& LL, // state of driving Poisson process
                     const double& b0, // will be reassigned the return value of this function, accumulates value of the integral as t0 varies from 0 to 1
                     const double& theta0, // will vary from 0 to 1
                     const double& dP) // will be constant while integrating theta0 over [0,1]
{
   // Returns a contribution of the integral of integrand() over [theta0,theta0+rk_dtheta] added to y0.
   double k1, k2, k3, k4; 
   k1 = dP * rk_dtheta * integrand( theta0, LL, b0, dP); // == y0*jump
   k2 = dP * rk_dtheta * integrand( theta0 + 0.5*rk_dtheta, LL, b0 + 0.5*k1, dP);
   k3 = dP * rk_dtheta * integrand( theta0 + 0.5*rk_dtheta, LL, b0 + 0.5*k2, dP);
   k4 = dP * rk_dtheta * integrand( theta0 + rk_dtheta, LL, b0 + k3, dP);

   //std::cout << "t0: " << t0 << ", y0:" << y0  << ", rk_dt: " << rk_dt << ", jump: " << jump << ", k's: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << std::endl;

   return ONESIXTH *(k1 + 2.0*k2 + 2.0*k3 + k4 );
}

//double SPF_NS::marcus_increment(
//      )
//{
//   integrand( tau, theta0, LL, b0, dP); // == y0*jump
//}

double SPF_NS::marcusRK4_increment( // returns y0 + contribution of integrand over dt // \Phi_{Y_k}
                     double (*integrand)(
                                          const double& tau,
                                          const double& theta,
                                          const double& LL,
                                          const double& bb,
                                          const double& dP),// Y_k
                     const double& tau,
                     const double& rk_dtheta, // constant
                     const double& LL, // state of driving Poisson process
                     const double& b0, // will be reassigned the return value of this function, accumulates value of the integral as t0 varies from 0 to 1
                     const double& theta0, // will vary from 0 to 1
                     const double& dP) // will be constant while integrating theta0 over [0,1]
{
   // Returns a contribution of the integral of integrand() over [theta0,theta0+rk_dtheta] added to y0.
   double k1, k2, k3, k4; 
   k1 = dP * rk_dtheta * integrand( tau, theta0, LL, b0, dP); // == y0*jump
   k2 = dP * rk_dtheta * integrand( tau, theta0 + 0.5*rk_dtheta, LL, b0 + 0.5*k1, dP);
   k3 = dP * rk_dtheta * integrand( tau, theta0 + 0.5*rk_dtheta, LL, b0 + 0.5*k2, dP);
   k4 = dP * rk_dtheta * integrand( tau, theta0 + rk_dtheta, LL, b0 + k3, dP);

   //std::cout << "t0: " << t0 << ", y0:" << y0  << ", rk_dt: " << rk_dt << ", jump: " << jump << ", k's: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << std::endl;

   return ONESIXTH *(k1 + 2.0*k2 + 2.0*k3 + k4 );
}

double SPF_NS::marcusRK4_increment( // returns y0 + contribution of integrand over dt // \Phi_{Y_k}
                     double (*bIntegrand)( const double& b0,
                                           const double& LL,
                                           const double& tau,
                                           const double& theta,
                                           const double& Yk
                                           ),
                     double (*integrand)( // marcusIntegrand3
                                          //const double& xTauOfTheta, // same as b0
                                          const double& bIntegrand,
                                          const double& tau,
                                          const double& theta,
                                          const double& LL,
                                          const double& poissonJump),// Y_k
                     const double& tau,
                     const double& rk_dtheta, // constant
                     const double& L0, // state of driving Poisson process
                     const double& b0, // will be reassigned the return value of this function, accumulates value of the integral as t0 varies from 0 to 1
                     const double& theta0, // will vary from 0 to 1
                     const double& jump) // will be constant while integrating theta0 over [0,1]
{
   // Returns a contribution of the integral of integrand() over [theta0,theta0+rk_dtheta] added to y0.
   double k1, k2, k3, k4; 
   k1 = jump * rk_dtheta * integrand(
         bIntegrand( b0, L0, tau, theta0, jump),
         tau, theta0, L0, jump); // == y0*jump
   k2 = jump * rk_dtheta * integrand(
         bIntegrand( b0 + 0.5*k1, L0, tau, theta0, jump),
         tau, theta0 + 0.5*rk_dtheta, L0, jump);
   k3 = jump * rk_dtheta * integrand(
         bIntegrand( b0 + 0.5*k2, L0, tau, theta0, jump),
         tau, theta0 + 0.5*rk_dtheta, L0, jump);
   k4 = jump * rk_dtheta * integrand(
         bIntegrand( b0 + k3, L0, tau, theta0, jump),
         tau, theta0 + rk_dtheta, L0, jump);

   //std::cout << "t0: " << t0 << ", y0:" << y0  << ", rk_dt: " << rk_dt << ", jump: " << jump << ", k's: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << std::endl;

   return ONESIXTH *(k1 + 2.0*k2 + 2.0*k3 + k4 );
}


#endif
