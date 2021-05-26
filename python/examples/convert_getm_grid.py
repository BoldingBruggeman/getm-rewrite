import argparse
from pygetm import domain, input

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('topo_file')
    parser.add_argument('out_file')
    parser.add_argument('-g', '--grids', action='store_true', help='Save metrics per separate grid (T,U,V,X) instead of for the supergrid')
    parser.add_argument('-f', '--full', action='store_true', help='Save secondary supergrid metrics. These can be inferred from the primary metrics, but saving them can help debugging.')
    args = parser.parse_args()

    dom = domain.Domain.from_getm_topo(args.topo_file, nlev=1)
    #H = input.request_from_netcdf(r"C:\Users\jornb\OneDrive\Code\igotm\data\GEBCO\GEBCO_2020.nc", 'elevation')
    #dom.set_bathymetry(H, scale_factor=-1, minimum_depth=0.)
    dom.initialize(runtype=1)
    if args.grids:
        dom.save_grids(args.out_file)
    else:
        dom.save(args.out_file, full=args.full)