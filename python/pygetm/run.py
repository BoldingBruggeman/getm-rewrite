import datetime
import os.path
import argparse
import sys
from typing import Mapping
import logging

import numpy

import pygetm
import pygetm.config
import pygetm.domain
import pygetm.input
import pygetm.airsea
import pygetm.legacy
import pygetm.input.tpxo
import pygetm.mixing

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('configuration', help='Path to configuration fiel in yaml format')
    parser.add_argument('-f', '--force', action='store_true', help='Continue even if unknown configuration settings are encountered')
    parser.add_argument('-r', '--report', type=int, help='Reporting interval', default=100)
    parser.add_argument('-l', '--log', type=str, help='Log level', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'), default='INFO')
    args = parser.parse_args()

    logging.basicConfig(level=args.log.upper())
    logger = logging.getLogger()

    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        pygetm.input.debug_nc_reads(logger=logger)

    root_dir = '.'

    config = pygetm.config.configure(args.configuration)

    start = config['time/start']
    stop = config['time/stop']
    timestep = config['time/timestep']
    split_factor = config.get('time/split_factor', 1)

    topo_path = config['domain/legacy_topo']
    nlev = config['domain/nlev']
    z0_const = config.get('bottom/z0_const', default=0.001)

    domain = pygetm.legacy.domain_from_topo(os.path.join(root_dir, topo_path), nlev=nlev, z0_const=z0_const)

    bdyinfo = config.get('open_boundaries/legacy_info')
    if bdyinfo is not None:
        logger.info('Loading legacy boundary info from %s...' % bdyinfo)
        pygetm.legacy.load_bdyinfo(domain, os.path.join(root_dir, bdyinfo))
        logging.info('%i open boundaries found.' % sum(map(len, domain.open_boundaries.values())))

    runtype = config.get('runtype', 1)
    uv_adv_scheme = config.get('momentum/advection/scheme', 1)
    fabm_config = config.get('fabm/configuration', 'fabm.yaml')
    if not config.get('fabm/use', False):
        fabm_config = None
    sim = pygetm.Simulation(domain, runtype=runtype, advection_scheme=uv_adv_scheme, fabm=fabm_config)

    inputs = []

    # Tidal boundary forcing from TPXO
    tpxo_dir = config.get('open_boundaries/elevation/tpxo_dir')
    if tpxo_dir is not None:
        logging.info('Loading TPXO tidal constituents at open boundaries from %s...' % tpxo_dir)
        tpxo_dir = os.path.join(root_dir, tpxo_dir)
        bdy_lon = domain.T.lon[domain.bdy_j, domain.bdy_i]
        bdy_lat = domain.T.lat[domain.bdy_j, domain.bdy_i]
        logger.info('- elevations...')
        tidal_h = pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=tpxo_dir)
        logger.info('- Eastward velocities...')
        tidal_u = pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='u', root=tpxo_dir)
        logger.info('- Northward velocities...')
        tidal_v = pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='v', root=tpxo_dir)
        inputs += [tidal_h.getm, tidal_u.getm, tidal_v.getm]

    # Meteorology from ERA
    logger.info('Loading meteorology...')
    meteo = config['meteorology']
    #tcc = meteo.get_input('tcc', domain, logger=logger)
    t2m = meteo.get_input('t2m', domain, logger=logger)
    d2m = meteo.get_input('d2m', domain, logger=logger)
    sp = meteo.get_input('sp', domain, logger=logger)
    u10 = meteo.get_input('u10', domain, logger=logger)
    v10 = meteo.get_input('v10', domain, logger=logger)
    inputs += [t2m, d2m, sp, u10, v10]

    if runtype == 4:
        temp = domain.T.array(fill=5., is_3d=True, name='temp', units='degrees_Celsius', long_name='conservative temperature')
        salt = domain.T.array(fill=35., is_3d=True, name='salt', units='-', long_name='absolute salinity')
        sst = temp.isel(z=-1)
        temp_source = domain.T.array(fill=0., is_3d=True, units='W m-2')
        salt_source = domain.T.array(fill=0., is_3d=True)
        shf = temp_source.isel(z=-1, name='shf', long_name='surface heat flux')

        SS = domain.W.array(fill=0., is_3d=True, name='SS', units='s-2', long_name='shear frequency squared')
        NN = domain.W.array(fill=0., is_3d=True, name='NN', units='s-2', long_name='buoyancy frequency squared')
        u_taus = domain.T.array(fill=0., name='u_taus', units='m s-1', long_name='surface shear velocity')
        u_taub = domain.T.array(fill=0., name='u_taub', units='m s-1', long_name='bottom shear velocity')
        z0s = domain.T.array(fill=0.1, name='z0s', units='m', long_name='hydrodynamic surface roughness')
        z0b = domain.T.array(fill=0.1, name='z0b', units='m', long_name='hydrodynamic bottom roughness')

        vertical_diffusion = pygetm.VerticalDiffusion(domain.T, cnpar=1.)
        turbulence = pygetm.mixing.Turbulence(domain)
        hf_scale_factor = split_factor * timestep / (pygetm.rho0 * pygetm.cp)

        # ensure ho and hn are up to date and identical
        domain.do_vertical()

    output = config.get('output')
    if output:
        assert isinstance(output, Mapping), 'output must be a mapping'
        for path, file_info in output.items():
            interval = file_info.get('time_step', 1)
            path += '.nc'
            logger.info('Output file %s will contain:' % path)
            f = sim.output_manager.add_netcdf_file(path, interval=interval)
            for source in file_info['variables']:
                assert isinstance(source, str), 'Currently all entries under output/FILE/variables must be variable names, not %s' % (source,)
                logger.info('- %s' % source)
                f.request(source)

    unused = config.check()
    if unused:
        logger.log(logging.WARNING if args.force else logging.ERROR, 'The following setting(s) in %s are not recognized:' % (args.configuration,))
        for path in unused:
            logger.log(logging.WARNING if args.force else logging.ERROR, '- %s' % path)
        if not args.force:
            logger.error('If you want to ignore these settings, provide the argument -f/--force.')
            sys.exit(2)
        logger.warn('The simulation will continue because you specified -f/--force, but these settings will not be used.')

    timedelta = datetime.timedelta(seconds=timestep)
    date = start
    istep = 0

    for variable in inputs:
        variable.update(date)
    sim.output_manager.save()

    logger.info('Starting simulation at %s' % date)
    itime = 0
    while date < stop:
        date += timedelta

        for variable in inputs:
            variable.update(date)

        # Update tidal elevations and velocities at the boundary
        if tpxo_dir is not None:
            sim.zbdy[:] = tidal_h
            sim.bdyu[:] = tidal_u
            sim.bdyv[:] = tidal_v

        if runtype != 4:
            sst = t2m - 273.15
        es, ea, qs, qa, rhoa = pygetm.airsea.humidity(domain.T, 3, d2m - 273.15, sp, sst, t2m - 273.15)
        tausx, tausy, qe, qh = pygetm.airsea.airsea_fluxes(domain.T, 1, sst + 273.15, t2m, u10, v10, rhoa, qs, qa)

        sim.update_sealevel_boundaries(timestep)
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
        sim.U.update_halos()
        sim.V.update_halos()
        sim.update_sealevel(timestep, sim.U, sim.V)
        sim.update_depth()

        istep += 1
        if args.report != 0 and istep % args.report == 0:
            logger.info(date)

        if runtype == 4 and itime % split_factor == 0:
            shf.all_values[...] = qe.all_values + qh.all_values
            vertical_diffusion(turbulence.nuh, split_factor * timestep, temp, ea4=temp_source * hf_scale_factor)
            vertical_diffusion(turbulence.nuh, split_factor * timestep, salt, ea4=salt_source)

            pygetm.mixing.get_buoyancy_frequency(salt, temp, out=NN)
            u_taus.all_values[...] = (tausx.all_values**2 + tausy.all_values**2)**0.25 / numpy.sqrt(pygetm.rho0)
            turbulence(split_factor * timestep, u_taus, u_taub, z0s, z0b, NN, SS)

        sim.output_manager.save()

        itime += 1

    sim.output_manager.close()
    logger.info('Simulation complete')