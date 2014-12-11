import OrbFitLib as ofl
import SimSpice as spice
import numpy as np

# Demo script 1 for OrbFitLib
#
# Compute the best orbit and uncertainites for Neptune using five synthetic data points from HST

if __name__=='__main__':

    ############################################
    # Parameters

    # The MPC-format file that we are trying to fit
    fn = 'neptune3.mpc'

    # The NAIF SPICE meta kernel
    #   Needs at least an spk file for Earth and a leapsecond file
    mk = 'basic.tm'

    # Maximum allowed RMS error of initial solution (in milliarcseconds)
    maxrms = 100.

    # Number of walkers
    nwalkers = 256

    # Number of initial emcee iterations to burn-in the distribution
    niter1 = 1000

    # Number of emcee iterations in the final chain
    niter2 = 1000

    # Flag to turn on/off plotting the result
    plot = True

    ############################################
    # Make the cloud

    # Create the MPC file object
    m = ofl.MPC_File(fn)
                
    # Add SPICE information to MPC object
    m.add_spice(mk)

    # Use emcee to create the distributions
    #   (This is where the magic happens)
    sampler = m.make_cloud(maxrms,nwalkers,niter1,niter2,verbose=True)

    # Pull out the final flattened chain
    samples = sampler.chain.reshape((-1,6))

    ############################################
    # Display results

    # Convert the sample state vectors to orbital elements
    seles = []
    for samp in samples:
        # Get the orbit and convert it from equitorial to ecliptic
        orb = ofl.eq2ec_orbit(m.to_orbit(samp))
        seles.append([orb.rp/(m.AU*(1-orb.e)),orb.e,orb.i*180./np.pi])

    # Use SPICE to find what the right answer should be for the state vector
    spice.furnsh(mk)
    state_true,lt = spice.spkgeo(8,m.et0,"J2000",0)
    spice.kclear()

    # Pull out the orbital elements of the true orbit
    orb_true = ofl.eq2ec_orbit(m.to_orbit(state_true))
    a_true = orb_true.rp/(m.AU*(1-orb_true.e))
    e_true,i_true = (orb_true.e,orb_true.i*180./np.pi)

    # Use a few percentiles to get the uncertainty
    a_err, e_err, i_err = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
            zip(*np.percentile(seles, [16, 50, 84], axis=0)))

    # Display those uncertainties
    print('\n                    True           Fit           Uncertainty')
    print('Semimajor axis:  {:10.7f}    {:10.7f}  +{:010.7f}  -{:010.7f}'.format(a_true,*a_err))
    print('Eccentricity:    {:10.7f}    {:10.7f}  +{:010.7f}  -{:010.7f}'.format(e_true,*e_err))
    print('Inclination (d): {:10.7f}    {:10.7f}  +{:010.7f}  -{:010.7f}'.format(i_true,*i_err))

    # Plot the Variance-Covariance of the resulting distribution
    if plot:
        import matplotlib.pyplot as plt
        import triangle
        fig = triangle.corner(seles, labels=["a (AU)","Eccentricity","I (deg)"],truths=[a_true,e_true,i_true])
        plt.show()
