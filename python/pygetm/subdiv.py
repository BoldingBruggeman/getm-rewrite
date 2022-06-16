import argparse
import logging
import cProfile
import pstats

import mpi4py
rank = mpi4py.MPI.COMM_WORLD.rank

import pygetm
import pygetm.legacy
import pygetm.parallel

def main():
    logging.basicConfig(level=logging.INFO if rank == 0 else logging.ERROR)
    logger = logging.getLogger()

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd', required=True)
    
    optimize_parser = subparsers.add_parser('optimize', help='compute the optimal subdomain division')
    optimize_parser.add_argument('path', help='path to topo file')
    optimize_parser.add_argument('ncpus', type=int, help='number of cores (active subdomains)')
    optimize_parser.add_argument('--max_protrude', type=float, default=0.5, help='maximum fraction of the subdomain that can protrude from the global domain (and thus be empty)')
    optimize_parser.add_argument('--pickle', help='path to save subdomain division to')
    optimize_parser.set_defaults(func=optimize)

    show_parser = subparsers.add_parser('show', help='describe existing subdomain division')
    show_parser.add_argument('path', help='path to load subdomain division from')
    show_parser.add_argument('--topo', help='path to topo file')
    show_parser.set_defaults(func=load)

    for p in (optimize_parser, show_parser):
        p.add_argument('--plot', action='store_true', help='plot subdomain decomposition')
        p.add_argument('--savefig', help='path to save figure to')
        p.add_argument('--profile', action='store_true')
        p.add_argument('--legacy', action='store_true')

    args = parser.parse_args()
    if args.profile and rank == 0:
        p = cProfile.Profile()
        tiling, background = p.runcall(args.func, args, logger)
        p.print_stats(sort=pstats.SortKey.TIME)
    else:
        tiling, background = args.func(args, logger)

    tiling.report(logging.getLogger())
    logger.info('Subdomain rank map:\n%s' % (tiling.map,))

    if args.plot and rank == 0:
        from matplotlib import pyplot
        ny, nx = tiling.map.shape
        fig, ax = pyplot.subplots(figsize=(0.75 * nx, 0.75 * ny))
        tiling.plot(ax=ax, background=background)
        if args.savefig:
            fig.savefig(args.savefig, dpi=300)
        else:
            pyplot.show()

def load_topo(path: bool, legacy: bool, logger: logging.Logger):
    logger.info('Reading topo from %s...' % path)
    if legacy:
        domain = pygetm.legacy.domain_from_topo(path, nlev=1, logger=logger, glob=True, z0=0.)
    else:
        domain = pygetm.domain.load(path, 1, logger=logger, glob=True)
    domain.initialize(1)
    return domain.T.mask.values.copy(), domain.T.H.ma.copy()

def optimize(args, logger):
    mask, background = load_topo(args.path, args.legacy, logger)
    tiling = pygetm.parallel.Tiling.autodetect(mask, logger=logger, ncpus=args.ncpus, max_protrude=args.max_protrude)

    if args.pickle and rank == 0:
        logger.info('Saving subdomain decomposition to %s...' % args.pickle)
        tiling.dump(args.pickle)

    return tiling, background

def load(args, logger):
    background = None
    if args.topo:
        _, background = load_topo(args.topo, args.legacy, logger)
    return pygetm.parallel.Tiling.load(args.path), background

if __name__ == '__main__':
    main()
