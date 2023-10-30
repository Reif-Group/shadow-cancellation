#!/usr/bin/env python
"""
A template file to simulate a specific system of Ordinary Differential Equations (ODEs).

... or an autogenerated script from the *crnsimulator* Python package.

Note: If this file is executable, it is *autogenerated*.
    This means it contains a system of hardcoded ODEs together with some
    default parameters. While it may be tempting to tweak some functions,
    beware that this file may be overwritten by the next execution of the
    `crnsimulator` executable.  Edits should be done at the source file:
    "crnsimulator.odelib_template.py" or you can provide an alternative
    template file. Use the option --output to avoid overwriting this file.

Usage: 
    python odesystem.py --help
"""

import logging
logger = logging.getLogger(__name__)

import argparse
import numpy as np
from scipy.integrate import odeint

class ODETemplateError(Exception):
    pass

rates = {
    
}

def odesystem(p0, t0, r):
    C, O, HelperCC, ProduceCC, ReactCC, U_10, U_22, U_33, U_37, U_9 = p0
    if not r : r = rates


    dCdt = -0.0020865*C*ReactCC - 0.00104769*C*U_22 + 0.00383326*HelperCC*U_22 + 0.00104769*ProduceCC*U_9 + 0.00175818*U_37*U_9
    dOdt = 0.045*O
    dHelperCCdt = -0.0038442886*HelperCC*U_22 + 0.000330858*U_37*U_9
    dProduceCCdt = 0.00104769*C*U_22 - 0.00104769*ProduceCC*U_9
    dReactCCdt = -0.0020865*C*ReactCC
    dU_10dt = 0.0020865*C*ReactCC
    dU_22dt = -0.00104769*C*U_22 - 0.0038442886*HelperCC*U_22 + 0.00104769*ProduceCC*U_9 + 0.000330858*U_37*U_9
    dU_33dt = 0.00383326*HelperCC*U_22 + 0.00175818*U_37*U_9
    dU_37dt = 1.10286e-5*HelperCC*U_22 - 0.002089038*U_37*U_9
    dU_9dt = 0.0020865*C*ReactCC + 0.00104769*C*U_22 + 1.10286e-5*HelperCC*U_22 - 0.00104769*ProduceCC*U_9 - 0.002089038*U_37*U_9
    return np.array([dCdt, dOdt, dHelperCCdt, dProduceCCdt, dReactCCdt, dU_10dt, dU_22dt, dU_33dt, dU_37dt, dU_9dt])

def jacobian(p0, t0, r):
    C, O, HelperCC, ProduceCC, ReactCC, U_10, U_22, U_33, U_37, U_9 = p0
    if not r : r = rates


    J = [[[] for i in range(len(p0))] for j in range(len(p0))]
    J[0][0] = -0.0020865*ReactCC - 0.00104769*U_22
    J[0][1] = 0
    J[0][2] = 0.00383326*U_22
    J[0][3] = 0.00104769*U_9
    J[0][4] = -0.0020865*C
    J[0][5] = 0
    J[0][6] = -0.00104769*C + 0.00383326*HelperCC
    J[0][7] = 0
    J[0][8] = 0.00175818*U_9
    J[0][9] = 0.00104769*ProduceCC + 0.00175818*U_37
    J[1][0] = 0
    J[1][1] = 0.0450000000000000
    J[1][2] = 0
    J[1][3] = 0
    J[1][4] = 0
    J[1][5] = 0
    J[1][6] = 0
    J[1][7] = 0
    J[1][8] = 0
    J[1][9] = 0
    J[2][0] = 0
    J[2][1] = 0
    J[2][2] = -0.0038442886*U_22
    J[2][3] = 0
    J[2][4] = 0
    J[2][5] = 0
    J[2][6] = -0.0038442886*HelperCC
    J[2][7] = 0
    J[2][8] = 0.000330858*U_9
    J[2][9] = 0.000330858*U_37
    J[3][0] = 0.00104769*U_22
    J[3][1] = 0
    J[3][2] = 0
    J[3][3] = -0.00104769*U_9
    J[3][4] = 0
    J[3][5] = 0
    J[3][6] = 0.00104769*C
    J[3][7] = 0
    J[3][8] = 0
    J[3][9] = -0.00104769*ProduceCC
    J[4][0] = -0.0020865*ReactCC
    J[4][1] = 0
    J[4][2] = 0
    J[4][3] = 0
    J[4][4] = -0.0020865*C
    J[4][5] = 0
    J[4][6] = 0
    J[4][7] = 0
    J[4][8] = 0
    J[4][9] = 0
    J[5][0] = 0.0020865*ReactCC
    J[5][1] = 0
    J[5][2] = 0
    J[5][3] = 0
    J[5][4] = 0.0020865*C
    J[5][5] = 0
    J[5][6] = 0
    J[5][7] = 0
    J[5][8] = 0
    J[5][9] = 0
    J[6][0] = -0.00104769*U_22
    J[6][1] = 0
    J[6][2] = -0.0038442886*U_22
    J[6][3] = 0.00104769*U_9
    J[6][4] = 0
    J[6][5] = 0
    J[6][6] = -0.00104769*C - 0.0038442886*HelperCC
    J[6][7] = 0
    J[6][8] = 0.000330858*U_9
    J[6][9] = 0.00104769*ProduceCC + 0.000330858*U_37
    J[7][0] = 0
    J[7][1] = 0
    J[7][2] = 0.00383326*U_22
    J[7][3] = 0
    J[7][4] = 0
    J[7][5] = 0
    J[7][6] = 0.00383326*HelperCC
    J[7][7] = 0
    J[7][8] = 0.00175818*U_9
    J[7][9] = 0.00175818*U_37
    J[8][0] = 0
    J[8][1] = 0
    J[8][2] = 1.10286e-5*U_22
    J[8][3] = 0
    J[8][4] = 0
    J[8][5] = 0
    J[8][6] = 1.10286e-5*HelperCC
    J[8][7] = 0
    J[8][8] = -0.002089038*U_9
    J[8][9] = -0.002089038*U_37
    J[9][0] = 0.0020865*ReactCC + 0.00104769*U_22
    J[9][1] = 0
    J[9][2] = 1.10286e-5*U_22
    J[9][3] = -0.00104769*U_9
    J[9][4] = 0.0020865*C
    J[9][5] = 0
    J[9][6] = 0.00104769*C + 1.10286e-5*HelperCC
    J[9][7] = 0
    J[9][8] = -0.002089038*U_9
    J[9][9] = -0.00104769*ProduceCC - 0.002089038*U_37
    return J


def add_integrator_args(parser):
    """ODE integration aruments."""
    solver = parser.add_argument_group('odeint parameters')
    plotter = parser.add_argument_group('plotting parameters')

    # required: simulation time and output settings
    solver.add_argument("--t0", type=float, default=0, metavar='<flt>',
            help="First time point of the time-course.")
    solver.add_argument("--t8", type=float, default=100, metavar='<flt>',
            help="End point of simulation time.")
    plotter.add_argument("--t-lin", type=int, default=500, metavar='<int>',
            help="Returns --t-lin evenly spaced numbers on a linear scale from --t0 to --t8.")
    plotter.add_argument("--t-log", type=int, default=None, metavar='<int>',
            help="Returns --t-log evenly spaced numbers on a logarithmic scale from --t0 to --t8.")

    # required: initial concentration vector
    solver.add_argument("--p0", nargs='+', metavar='<int/str>=<flt>',
            help="""Vector of initial species concentrations. 
            E.g. \"--p0 1=0.5 3=0.7\" stands for 1st species at a concentration of 0.5 
            and 3rd species at a concentration of 0.7. You may chose to address species
            directly by name, e.g.: --p0 C=0.5.""")
    # advanced: scipy.integrate.odeint parameters
    solver.add_argument("-a", "--atol", type=float, default=None, metavar='<flt>',
            help="Specify absolute tolerance for the solver.")
    solver.add_argument("-r", "--rtol", type=float, default=None, metavar='<flt>',
            help="Specify relative tolerance for the solver.")
    solver.add_argument("--mxstep", type=int, default=0, metavar='<int>',
            help="Maximum number of steps allowed for each integration point in t.")

    # optional: choose output formats
    plotter.add_argument("--list-labels", action='store_true',
            help="Print all species and exit.")
    plotter.add_argument("--labels", nargs='+', default=[], metavar='<str>+',
            help="""Specify the (order of) species which should appear in the pyplot legend, 
            as well as the order of species for nxy output format.""")
    plotter.add_argument("--labels-strict", action='store_true',
            help="""When using pyplot, only plot tracjectories corresponding to labels,
            when using nxy, only print the trajectories corresponding to labels.""")
 
    plotter.add_argument("--nxy", action='store_true',
            help="Print time course to STDOUT in nxy format.")
    plotter.add_argument("--header", action='store_true',
            help="Print header for trajectories.")

    plotter.add_argument("--pyplot", default='', metavar='<str>',
            help="Specify a filename to plot the ODE simulation.")
    plotter.add_argument("--pyplot-xlim", nargs=2, type=float, default=None, metavar='<flt>',
            help="Specify the limits of the x-axis.")
    plotter.add_argument("--pyplot-ylim", nargs=2, type=float, default=None, metavar='<flt>',
            help="Specify the limits of the y-axis.")
    plotter.add_argument("--pyplot-labels", nargs='+', default=[], metavar='<str>+',
            help=argparse.SUPPRESS)
    return

def flint(inp):
    return int(float(inp)) if float(inp) == int(float(inp)) else float(inp)

def set_logger(verbose, logfile):
    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(logfile) if logfile else logging.StreamHandler()
    if verbose == 0:
        handler.setLevel(logging.WARNING)
    elif verbose == 1:
        handler.setLevel(logging.INFO)
    elif verbose == 2:
        handler.setLevel(logging.DEBUG)
    elif verbose >= 3:
        handler.setLevel(logging.NOTSET)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def integrate(args, setlogger = False):
    """Main interface to solve the ODE-system.

    Args:
      args (:obj:`argparse.ArgumentParser()`): An argparse object containing all of
        the arguments of :obj:`crnsimulator.add_integrator_args()`.

    Prints:
      - plot files
      - time-course

    Returns:
      Nothing
    """
    if setlogger:
        set_logger(args.verbose, args.logfile)

    if args.pyplot_labels:
        logger.warning('Deprecated argument: --pyplot_labels.')

    svars = ["C", "O", "HelperCC", "ProduceCC", "ReactCC", "U_10", "U_22", "U_33", "U_37", "U_9"]

    p0 = [0] * len(svars)
    p0[0] = 1.0
    p0[1] = 1.0
    p0[2] = 150.0
    p0[3] = 200.0
    p0[4] = 100.0
    
    const = None
    #<&>CONSTANT_SPECIES_INFO<&>#
    if args.p0:
        for term in args.p0:
            p, o = term.split('=')
            try:
                pi = svars.index(p)
            except ValueError as e:
                pi = int(p) - 1
            finally:
                p0[pi] = flint(o)
    else:
        msg = 'Specify a vector of initial concentrations: ' + \
                'e.g. --p0 1=0.1 2=0.005 3=1e-6 (see --help)'
        if sum(p0) == 0:
            logger.warning(msg)
            args.list_labels = True
        else:
            logger.info(msg)

    if args.list_labels:
        print('List of variables and initial concentrations:')
        for e, v in enumerate(svars, 1):
            if args.labels_strict and e > len(args.labels):
                break
            print(f'{e} {v} {p0[e-1]} {"constant" if const and const[e-1] else ""}')
        raise SystemExit('Initial concentrations can be overwritten by --p0 argument')

    if not args.nxy and not args.pyplot:
        logger.warning('Use --pyplot and/or --nxy to plot your results.')

    if not args.t8:
        raise ODETemplateError('Specify a valid end-time for the simulation: --t8 <flt>')

    if args.t_log:
        if args.t0 == 0:
            raise ODETemplateError('--t0 cannot be 0 when using log-scale!')
        time = np.logspace(np.log10(args.t0), np.log10(args.t8), num=args.t_log)
    elif args.t_lin:
        time = np.linspace(args.t0, args.t8, num=args.t_lin)
    else:
        raise ODETemplateError('Please specify either --t-lin or --t-log. (see --help)')

    # It would be nice if it is possible to read alternative rates from a file instead.
    # None triggers the default-rates that are hard-coded in the (this) library file.
    rates = None

    logger.info(f'Initial concentrations: {list(zip(svars, p0))}')
    # TODO: logging should report more info on parameters.

    ny = odeint(odesystem,
        np.array(p0), time, (rates, ), Dfun = jacobian,
        atol=args.atol, rtol=args.rtol, mxstep=args.mxstep).T

    # Output
    if args.nxy and args.labels_strict:
        end = len(args.labels)
        if args.header:
            print(' '.join(['{:15s}'.format(x) for x in ['time'] + svars[:end]]))
        for i in zip(time, *ny[:end]):
            print(' '.join(map("{:.9e}".format, i)))
    elif args.nxy:
        if args.header:
            print(' '.join(['{:15s}'.format(x) for x in ['time'] + svars]))
        for i in zip(time, *ny):
            print(' '.join(map("{:.9e}".format, i)))

    if args.pyplot:
        from crnsimulator.plotting import ode_plotter
        plotfile = ode_plotter(args.pyplot, time, ny, svars,
                               log=True if args.t_log else False,
                               labels=set(args.labels),
                               xlim = args.pyplot_xlim,
                               ylim = args.pyplot_ylim,
                               labels_strict = args.labels_strict)
        logger.info(f"Plotting successfull. Wrote plot to file: {plotfile}")

    return zip(time, *ny)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--verbose", action='count', default = 0,
        help = "Print logging output. (-vv increases verbosity.)")
    parser.add_argument('--logfile', default = '', action = 'store', metavar = '<str>',
        help = """Redirect logging information to a file.""")
    add_integrator_args(parser)
    args = parser.parse_args()
    integrate(args, setlogger = True)

