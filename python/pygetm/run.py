import datetime
import os.path
import argparse
import sys
from typing import Mapping

import pygetm
import pygetm.config
import pygetm.domain
import pygetm.input
import pygetm.airsea
import pygetm.legacy
import pygetm.input.tpxo

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('configuration', help='Path to configuration fiel in yaml format')
    parser.add_argument('-f', '--force', action='store_true', help='Continue even if unknown configuration settings are encountered')
    parser.add_argument('-r', '--report', type=int, help='Reporting interval', default=100)
    args = parser.parse_args()

    root_dir = '.'

    config = pygetm.config.configure(args.configuration)

    start = config['time/start']
    stop = config['time/stop']
    timestep = config['time/timestep']

    topo_path = config['domain/legacy_topo']
    nlev = config['domain/nlev']
    z0_const = config.get('bottom/z0_const', default=0.001)

    domain = pygetm.legacy.domain_from_topo(os.path.join(root_dir, topo_path), nlev=nlev, z0_const=z0_const)

    bdyinfo = config.get('open_boundaries/legacy_info')
    if bdyinfo is not None:
        print('Loading legacy boundary info from %s...' % bdyinfo)
        pygetm.legacy.load_bdyinfo(domain, os.path.join(root_dir, bdyinfo))

    runtype = config.get('runtype', 1)
    uv_adv_scheme = config.get('momentum/advection/scheme', 1)
    sim = pygetm.Simulation(domain, runtype=runtype, advection_scheme=uv_adv_scheme)

    inputs = []

    # Tidal boundary forcing from TPXO
    tpxo_dir = config.get('open_boundaries/elevation/tpxo_dir')
    if tpxo_dir is not None:
        print('Loading TPXO tidal constituents forcing from %s...' % tpxo_dir)
        tpxo_dir = os.path.join(root_dir, tpxo_dir)
        bdy_lon = domain.T.lon[domain.bdy_j, domain.bdy_i]
        bdy_lat = domain.T.lat[domain.bdy_j, domain.bdy_i]
        print('- elevations...')
        tidal_h = pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=tpxo_dir)
        print('- Eastward velocities...')
        tidal_u = pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='u', root=tpxo_dir)
        print('- Northward velocities...')
        tidal_v = pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='v', root=tpxo_dir)
        inputs += [tidal_h.getm, tidal_u.getm, tidal_v.getm]

    # Meteorology from ERA
    print('Loading meteorology...')
    meteo = config['meteorology']
    #tcc = meteo.get_input('tcc', domain, verbose=True)
    t2m = meteo.get_input('t2m', domain, verbose=True)
    d2m = meteo.get_input('d2m', domain, verbose=True)
    sp = meteo.get_input('sp', domain, verbose=True)
    u10 = meteo.get_input('u10', domain, verbose=True)
    v10 = meteo.get_input('v10', domain, verbose=True)
    inputs += [t2m, d2m, sp, u10, v10]

    output = config.get('output')
    if output:
        assert isinstance(output, Mapping), 'output must be a mapping'
        for path, file_info in output.items():
            interval = file_info.get('time_step', 1)
            path += '.nc'
            print('Output file %s will contain:' % path)
            f = sim.output_manager.add_netcdf_file(path, interval=interval)
            for source in file_info['variables']:
                assert isinstance(source, str), 'Currently all entries under output/FILE/variables must be variable names, not %s' % (source,)
                print('- %s' % source)
                f.request(source)

    unused = config.check()
    if unused:
        print('%s: the following setting(s) in %s are not recognized:' % ('WARNING' if args.force else 'ERROR', args.configuration))
        for path in unused:
            print('- %s' % path)
        if not args.force:
            print('If you want to ignore these settings, provide the argument -f/--force.')
            sys.exit(2)
        print('The simulation will continue because you specified -f/--force, but these settings will not be used.')

    inputs = t2m, d2m, sp, u10, v10
    timedelta = datetime.timedelta(seconds=timestep)
    date = start
    istep = 0
    while date < stop:
        date += timedelta

        for variable in inputs:
            variable.update(date)

        # Update tidal elevations and velocities at the boundary
        if tpxo_dir is not None:
            sim.zbdy[:] = tidal_h
            sim.bdyu[:] = tidal_u
            sim.bdyv[:] = tidal_v

        es, ea, qs, qa, rhoa = pygetm.airsea.humidity(domain.T, 3, d2m - 273.15, sp, t2m - 273.15, t2m - 273.15)
        tausx, tausy, qe, qh = pygetm.airsea.airsea_fluxes(domain.T, 1, t2m, t2m, u10, v10, rhoa, qs, qa)

        sim.update_sealevel_boundaries(timestep)
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
        sim.U.update_halos()
        sim.V.update_halos()
        sim.update_sealevel(timestep, sim.U, sim.V)
        sim.update_depth()

        istep += 1
        if args.report != 0 and istep % args.report == 0:
            print(date)

        sim.output_manager.save()

    sim.output_manager.close()