#ifndef __LAMBERT_HPP__
#define __LAMBERT_HPP__

// Solver for Lambert's Problem
//
// Based on "A Unified Form of Lambert's Theorem"
// E.R. Lancaster and R.C. Blanchard (GSFC)
// NASA Technical Note D-5368, September 1969 

// For the vec_math namespace
#include "conic.hpp"

class Lambert { public:

    double eta;

    Lambert() {
        eta = 1e-10;
    }

    // Minimum eccentricity transfer time
    double T_minE(double q) const {
        double beta_m = 2.*asin(q);
        return M_PI - beta_m + sin(beta_m);
    }

    // Parabolic transfer time
    double T_para(double theta,double q) const {
        if (theta>0) return (4./3.)*(1. - q*q*q);
        else return (4./3.)*(1. + q*q*q);
    }

    // Calculate x from the inverse semimajor axis
    double calc_x(double ia,double s,bool neg_alpha) const {
        // Parabolic case
        if (fabs(ia)<eta) return 1;
        double E = -s*ia*0.5;
        // Elliptical case
        if (E<0) {
            double alpha = 2*asin(sqrt(-E));
            if (neg_alpha) alpha = 2.*M_PI - alpha;
            return cos(alpha/2.);
        // Hyperbolic case
        } else {
            double gamma = 2*asinh(sqrt(E));
            return cosh(gamma/2.);
        }
    }

    // Calculated the inverse of semimajor axis as a function of x
    double calc_ia(double x,double s) const {
        if (fabs(x-1.)<eta) return 0;
        double E = x*x - 1.;
        return -2.*E/s;
    }
     
    // Special functions for parabola
    double sigma(double u) const {
        double rootu = sqrt(u);
        double hnum = asin(rootu) - rootu*sqrt(1.-u);
        return 2.*hnum/(rootu*rootu*rootu);
    }
    double sigmap(double u) const {
        double root1u = sqrt(1.-u), rootu = sqrt(u);
        return (2.*(u/root1u - root1u + 1./root1u)/(rootu*rootu*rootu) + 
                6.*(rootu*root1u - asin(rootu))/(u*u));
    }

    // T as a function of x and q
    double Txq(double x,double q) const {
        double K = q*q, E = x*x - 1.;
        double rho = fabs(E);
        if (rho<eta)
            return sigma(-E) - q*K*sigma(-K*E);
        double y = sqrt(rho);
        double z = sqrt(1+K*E);
        double f = y*(z-q*x);
        double g = x*z - q*E, d;
        if (E<0) d = atan2(f,g);
        else d = log(f+g);
        return 2.*(x - q*z - d/y)/E;
    }

    // first derivative of T(x,q) with respect to x
    double dTdx(double x,double q) const { 
        double K = q*q, E = x*x - 1.;
        if (fabs(1.-x) < eta)
            return 2.*x*(q*K*K*sigmap(-K*E) - sigmap(-E));
        double z = sqrt(1. + K*E);
        double T = Txq(x,q);
        return (4. - 4.*q*K*x/z - 3.*x*T)/E;
    }

    // x as a function of T and q
    double xTq(double T,double q,double theta) const {
        double Tp=T_para(theta,q), x, d=1;
        // Parabola
        if (fabs(T-Tp)<eta) return 1;
        // Hyperbola
        else if (T<Tp) x = 1.5;
        // Ellipses
        else {
            double Tm = T_minE(q);
            // Minimum energy ellipse
            if (fabs(T-Tm)<eta) return 0;
            // Long ellipse
            else if (T>Tm) x = -0.5;
            // Short ellipse
            else x = 0.5;
        }
        // Newton's method to solve for x
        for(int ic=0;ic<100;ic++) {
            d = (Txq(x,q)-T)/dTdx(x,q);
            x -= d;
            if (fabs(d) < eta) break;
        }
        return x;
    }

    std::vector<double> transfer(double mu,const std::vector<double> &r1v,const std::vector<double> &r2v,double dt) const {
        using namespace vec_math;

        // Normalize the vectors
        double r1 = norm(r1v), r2 = norm(r2v);
        std::vector<double> r1n = div(r1v,r1), r2n = div(r2v,r2);
        
        // Define the geometry of the problem
        double theta = acos(dot(r1n,r2n));
        if (theta<0) theta += 2.*M_PI;
        double c = sqrt(r1*r1 + r2*r2 - 2.*r1*r2*cos(theta));
        double s = (r1 + r2 + c)/2.;
        double q = sqrt(1.-(c/s));
        
        // Convert to normalized time units
        double T = sqrt(8.*mu/s)*dt/s;
        
        // Solve for the inverse of semimajor axis
        double x = xTq(T,q,theta);
        double ia = calc_ia(x,s);
        double r1dot,v1,v1t;
        
        // Special case for a parabolic transfer
        if(fabs(ia)<eta) {
            r1dot = 0;
            v1 = mu*2/r1;
            v1t = sqrt(v1);
        
        // General case
        } else {
            double K = q*q, E = x*x - 1.;
            double z = sqrt(1. + K*E);
            r1dot = sqrt(2.*mu*s)*(q*z*(s-r1) - x*(s-r2))/(c*r1);
            double p = 2.*r1 - r1*r1*ia - (r1*r1*r1dot*r1dot)/mu;
            v1 = mu*((2/r1)-ia);
            v1t = sqrt(mu*p)/r1;
        }
        
        // Calculate the required velocity at r1
        double r1c = r1dot - v1t/tan(theta);
        double r2c = v1t/sin(theta);
        return add( mult(r1c,r1n), mult(r2c,r2n) );
    }

};

#endif
