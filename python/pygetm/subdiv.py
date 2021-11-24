import argparse
import sys
import logging
import pickle

from . import legacy
from . import parallel

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--legacy', action='store_true')
    parser.add_argument('path', help='path to topo file')
    parser.add_argument('ncpus', type=int, help='number of cores (active subdomains)')
    parser.add_argument('--pickle', help='path to save subdomain division to')
    parser.add_argument('--bdyinfo', help='path to bdyinfo.dat')
    args = parser.parse_args()

    if not args.legacy:
        print('Currently only legacy topo files are supported. Supply --legacy')
        sys.exit(2)

    logging.basicConfig(level=logging.INFO)

    logger = logging.getLogger()
    logger.info('Reading topo from %s...' % args.path)
    domain = legacy.domain_from_topo(args.path, nlev=1, logger=logger)
    if args.bdyinfo:
        legacy.load_bdyinfo(domain, args.bdyinfo)
    domain.initialize(1)

    tiling = parallel.Tiling.autodetect(domain.T.mask, logger=logger, ncpus=args.ncpus)

    if args.pickle:
        logger.info('Saving subdomain decomposition to %s...' % args.pickle)
        tiling.dump(args.pickle)

if __name__ == '__main__':
    main()
