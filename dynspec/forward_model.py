import numpy as np
from scipy.stats import exponnorm
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
import matplotlib.pyplot as plt

import argparse


def create_time_axis(t0, dt, nbin, t0_is_ctr_of_first_bin=False):
    if t0_is_ctr_of_first_bin:
        return np.arange(nbin)*dt + t0
    else:
        return np.arange(nbin)*dt + t0 + dt/2

def create_freq_axis(fctr, bw, nchan):
    flo = fctr - bw/2
    df = bw/nchan
    return np.arange(nchan)*df + flo + df/2

def main():

    parser = argparse.ArgumentParser(description="Forward modelling for highly scattered and dispersed LPT pulses")

    parser.add_argument('coord', help="Sky coordinate of LPT in 'HH:MM:SS.S ±DD:MM:SS.S' format")
    parser.add_argument('telescope', help="Telescope name (must be a listed site in Astropy)")
    parser.add_argument('--dynspec_plot', help="The output filename for the dynamic spectrum plot")

    parser_frequency_group = parser.add_argument_group('Frequency arguments')
    parser_frequency_group.add_argument('--fctr_MHz', type=float, default=200.0, help="Centre frequency in MHz (default=200)")
    parser_frequency_group.add_argument('--bw_MHz', type=float, default=50.0, help="Bandwidth in MHz (default=50)")
    parser_frequency_group.add_argument('--nchan', type=int, help="Number of channels. If not provided, set to round(bw_MHz).")

    parser_time_group = parser.add_argument_group('Time arguments')
    parser_time_group.add_argument('--dt', type=float, default=1.0, help="Sampling time in seconds (default=1)")
    parser_time_group.add_argument('--t0', type=float, help="Start time in MJD. If not provided, set to PEPOCH - 10*width as determined from the ephemeris")
    parser_time_group.add_argument('--nbin', type=int, help="Number of time bins. If not provided, set to round(20*width/dt)")
    parser_time_group.add_argument('--t0_is_ctr_of_first_bin', action='store_true', help="If set, t0 is interpreted as the centre of the first time bin. Otherwise it is interpreted as the left edge of the first time bin")

    parser_model_group = parser.add_argument_group('LPT model parameters')
    parser_model_group.add_argument('--PEPOCH', type=float, default=60000.0, help="Reference epoch (MJD) for pulse numbers (default=60000.0)")
    parser_model_group.add_argument('--period', type=float, default=3600.0, help="Pulse period in seconds (default=3600.0)")
    parser_model_group.add_argument('--DM', type=float, default=100.0, help="Dispersion measure (default=100.0)")
    parser_model_group.add_argument('--tau_sc_1GHz', type=float, default=1e-3, help="Scattering timescale at 1 GHz in seconds (default=1e-3)")
    parser_model_group.add_argument('--sc_idx', type=float, default=-4.0, help="Scattering index (default=-4)")
    parser_model_group.add_argument('--width', type=float, default=100.0, help="Pulse width in seconds (default=100)")

    args = parser.parse_args()

    # Sort out argparse defaults that need to be calculated from other values
    fctr = args.fctr_MHz * u.MHz
    bw = args.bw_MHz * u.MHz
    nchan = args.nchan or int(np.round((bw/u.MHz).decompose()))
    dt = args.dt * u.s
    PEPOCH = Time(args.PEPOCH, scale='utc', format='mjd')
    period = args.period * u.s
    width = args.width * u.s
    tau_sc_1GHz = args.tau_sc_1GHz * u.s
    DM = args.DM * u.pc / u.cm**3

    if args.t0 is not None:
        t0 = Time(args.t0, scale='utc', format='mjd')
    else:
        t0 = PEPOCH - 10*width

    nbin = args.nbin or int(np.round((20*width/dt).decompose()))

    bary_t = create_time_axis(t0, dt, nbin, t0_is_ctr_of_first_bin=args.t0_is_ctr_of_first_bin)
    f = create_freq_axis(fctr, bw, nchan)

    #print(f"{nbin = }, {width = }, {dt = }, {t0 = }, {t[-1] = }")

    # Correct for barycentring
    coord = SkyCoord(args.coord, unit=(u.hourangle, u.deg), frame='fk4')
    telescope = EarthLocation.of_site(args.telescope)
    topo_t = bary_t + bary_t.light_travel_time(coord, ephemeris='jpl', location=telescope)

    # Convert to pulse phase
    pulse, phase = np.divmod(((topo_t - PEPOCH)/period).decompose() + 0.5, 1.0)
    phase -= 0.5

    PHASE, F = np.meshgrid(phase, f)

    # And reintroduce time units
    PHASE_time = PHASE * period

    # Correct for dispersion (scattering is dealt with later)
    D = 4148.808 * u.MHz**2 * u.cm**3 * u.s / u.pc
    PHASE_time -= D * DM / F**2

    # Make an exponentially scattered pulse. Use sigma = width/2
    TAU_SC = tau_sc_1GHz * (F / u.GHz).decompose()**args.sc_idx
    modelled_pulse = exponnorm.pdf(
        PHASE_time,
        TAU_SC.to(PHASE_time.unit),
        loc=0.0,
        scale=width.to(PHASE_time.unit)/2
    ) * np.sqrt(2*np.pi)*width.to(PHASE_time.unit).value/2 # <-- normalisation so that peak = 1

    # Plot it
    if args.dynspec_plot is not None:
        plt.pcolormesh((topo_t - topo_t[0]).to('s').value, f.to('MHz').value, modelled_pulse, shading='auto', vmax=1.0, vmin=0.0)
        plt.colorbar()
        plt.xlabel(f"Time (s) since MJD {topo_t[0].mjd:.7f} (= GPS {topo_t[0].gps:.0f})")
        plt.ylabel(f"Frequency (MHz)")
        plt.savefig(args.dynspec_plot)

if __name__ == '__main__':
    main()
