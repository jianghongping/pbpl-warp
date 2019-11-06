#!/usr/bin/env python
import sys, math, os, glob
import re as regex
import argparse
import toml
import numpy as np
import h5py
import asteval
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
from pbpl.units import *
from .core import setup_plot

def get_parser():
    parser = argparse.ArgumentParser(
        description='Make field animation from HDF5 dump')
    parser.add_argument(
        'conf', metavar='CONF',
        help='Configuration file (e.g., field-anim.toml)')
    return parser

def get_args():
    parser = get_parser()
    args = parser.parse_args()
    args.conf = toml.load(args.config_filename)
    return args

def plot_frame(output, args, aeval):
    fig = plot.figure(figsize=(244.0/72, 140.0/72))
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(
        aeval(args.xaxis_value), aeval(args.yaxis_value),
        marker='o', ls='', markersize=3, markeredgewidth=0,
        color='#0083b8', alpha=0.5)


    ax.set_xlabel(args.xaxis_title, labelpad=0.0)
    ax.set_ylabel(args.yaxis_title, labelpad=0.0)

    ax.set_xlim(args.xaxis_min, args.xaxis_max)
    ax.set_ylim(args.yaxis_min, args.yaxis_max)

    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    output.savefig(fig, transparent=True)

def get_sim_steps(path):
    A = glob.glob(path + '_[0-9]*.h5')
    m = regex.compile('.*_([0-9]+).h5$')
    result = [int(y.groups()[0]) for y in [m.match(x) for x in A] ]
    result.sort()
    return result

def main():
    args = get_args()
    sys.exit()
    
    available_steps = get_sim_steps(args.input_path)
    if args.input_i1 == -1:
        args.input_i1 = available_steps[-1]
    steps = range(args.input_i0, args.input_i1, args.input_istep)

    setup_plot()

    # create safe interpreter for evaluation of scale expressions
    aeval = asteval.Interpreter()
    import pbpl.units
    for x in pbpl.units.__all__:
        aeval.symtable[x] = pbpl.units.__dict__[x]


    output = PdfPages(args.output_filename)
    for i in steps:
        f = h5py.File('{}_{:04}.h5'.format(args.input_path, i), 'r')
        particles = f['particles'].value.T
        N = args.input_num_particles
        if N == -1:
            N = particles.shape[1]
        particles = particles[:, 0:N]

        # reference particle
        m0 = f['mass'].value * (GeV/c_light**2)
        p0 = f['pz'].value * (GeV/c_light)
        # beta0 = p0 / np.sqrt((m0 * c_light)**2 + p0**2)
        # gamma0 = 1 / np.sqrt(1 - beta0**2)
        gamma0 = np.sqrt(1 + (p0/(m0 * c_light))**2)
        beta0 = np.sqrt(1 - (1/gamma0**2))

        x = particles[0] * meter
        px = particles[1] * p0
        y = particles[2] * meter
        py = particles[3] * p0
        cdt = particles[4] * meter
        deltap = particles[5] * p0
        zeta = -beta0 * cdt
        p = p0 + deltap
        pz = np.sqrt(p**2 - px**2 - py**2)
        dpz = pz - p0
        gamma = np.sqrt(1 + (p/(m0 * c_light))**2)
        beta = np.sqrt(1 - (1/gamma**2))
        vz = (p0 + deltap) / (gamma * m0)
        v0 = beta0 * c_light
        print(px/p0)
        print(py/p0)
        print(dpz/p0)
#        sys.exit()

        aeval.symtable['m0'] = m0
        aeval.symtable['p0'] = p0
        aeval.symtable['gamma0'] = gamma0
        aeval.symtable['beta0'] = beta0

        aeval.symtable['x'] = x
        aeval.symtable['px'] = px
        aeval.symtable['y'] = y
        aeval.symtable['py'] = py
        aeval.symtable['zeta'] = zeta
        aeval.symtable['deltap'] = deltap
        aeval.symtable['p'] = p
        aeval.symtable['pz'] = pz
        aeval.symtable['dpz'] = dpz

        aeval.symtable['vz'] = vz
        aeval.symtable['v0'] = v0

        plot_frame(output, args, aeval)
        f.close()
    output.close()

if __name__ == '__main__':
    sys.exit(main())
