import argparse
import sys
import logging

from . import legacy
from . import parallel
import pygetm

def main():
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd', required=True)
    
    optimize_parser = subparsers.add_parser('optimize', help='compute the optimal subdomain division')
    optimize_parser.add_argument('--legacy', action='store_true')
    optimize_parser.add_argument('path', help='path to topo file')
    optimize_parser.add_argument('ncpus', type=int, help='number of cores (active subdomains)')
    optimize_parser.add_argument('--pickle', help='path to save subdomain division to')
    optimize_parser.set_defaults(func=optimize)

    show_parser = subparsers.add_parser('show', help='describe existing subdomain division')
    show_parser.add_argument('path', help='path to load subdomain division from')
    show_parser.set_defaults(func=load)

    for p in (optimize_parser, show_parser):
        p.add_argument('--plot', action='store_true', help='plot subdomain decomposition')

    args = parser.parse_args()
    tiling = args.func(args, logger)

    tiling.report(logging.getLogger())
    logger.info('Subdomain rank map:\n%s' % (tiling.map,))

    if args.plot:
        from matplotlib import pyplot
        fig, ax = pyplot.subplots()
        tiling.plot(ax=ax)
        pyplot.show()

def optimize(args, logger):
    if not args.legacy:
        print('Currently only legacy topo files are supported. Supply --legacy')
        sys.exit(2)

    logger.info('Reading topo from %s...' % args.path)
    domain = legacy.domain_from_topo(args.path, nlev=1, logger=logger)
    domain.initialize(1)

    tiling = parallel.Tiling.autodetect(domain.T.mask, logger=logger, ncpus=args.ncpus)

    if args.pickle:
        logger.info('Saving subdomain decomposition to %s...' % args.pickle)
        tiling.dump(args.pickle)

    return tiling

def load(args, logger):
    return parallel.Tiling.load(args.path)

if __name__ == '__main__':
    main()