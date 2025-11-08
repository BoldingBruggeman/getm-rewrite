import argparse
import numpy as np
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt

import pygetm

N = 5000
Ri, Ro = 5000.0, 10000.0
Nr, Nt = 15, 45
z0 = 0.1
lat = 53.5
H = 10.0
taux = 0.01
tauy = 0.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p = argparse.ArgumentParser(
        description="Various simulation configurations in an annulus domain setup"
    )
    p.add_argument("--N", type=int, help="Number of timesteps", default=Nr)
    p.add_argument("--Ri", type=float, help="Inner radius", default=Ri)
    p.add_argument("--Ro", type=float, help="Outer radius", default=Ro)
    p.add_argument("--Nr", type=int, help="Number of radial points", default=Nr)
    p.add_argument("--Nt", type=int, help="Number of azimutal points", default=Nt)
    p.add_argument("--z0", type=float, help="Friction factor", default=z0)
    p.add_argument("--lat", type=float, help="Latitude of setup", default=lat)
    p.add_argument(
        "--taux", type=float, help="Surface stress - x-direction", default=taux
    )
    p.add_argument(
        "--tauy", type=float, help="Surface stress - y-direction", default=tauy
    )
    p.add_argument("--H", type=float, help="Undisturbed water depth", default=H)
    p.add_argument("--plot", action="store_true", help="Show some plots")
    p.add_argument(
        "--no_cyclic",
        action="store_true",
        help="Don't use cyclic boundaries",
    )
    p.add_argument(
        "--radial_wind",
        action="store_true",
        help="Run the radial wind setup",
    )
    return p.parse_args()


def cartesian_annulus_grid(r_in, r_out, Nx, Ny, origin=(0.0, 0.0)):
    """Return (X, Y, mask) where mask=True for points inside the ring."""
    ox, oy = origin
    # Build a square domain that fully contains the outer circle
    L = max(abs(r_out), abs(r_in)) * 1.01  # a tiny margin
    x = np.linspace(-L, L, Nx)
    y = np.linspace(-L, L, Ny)
    X, Y = np.meshgrid(x, y, indexing="ij")
    # Distance from centre
    R = np.sqrt((X - ox) ** 2 + (Y - oy) ** 2)
    mask = (R >= r_in) & (R <= r_out)
    return X, Y, mask


def polar_annulus_grid(r_in, r_out, Nr, Nt, origin=(0.0, 0.0)):
    """Return (R, Θ, X, Y) on a regular r‑θ lattice."""
    ox, oy = origin
    r = np.linspace(r_in, r_out, Nr)  # radial nodes
    theta = np.linspace(0.0, 2 * np.pi, Nt, endpoint=True)  # azimuthal nodes
    R, Θ = np.meshgrid(r, theta, indexing="ij")  # shape (Nr, Nt)

    # Convert to Cartesian if you need (x, y) coordinates
    X = ox + R * np.cos(Θ)
    Y = oy + R * np.sin(Θ)
    return R, Θ, X, Y


def create_domain(Ri, Ro, Nr, Nt, cyclic=True, plot=False):
    R, theta, xx, yx = polar_annulus_grid(Ri, Ro, Nr, Nt, origin=(0.0, 0.0))
    # xx, yx, mask = cartesian_annulus_grid(Ri, Ro, Nr, Nt, origin=(0.0, 0.0))
    topo = np.zeros((Nr - 1, Nt - 1))
    topo[:, :] = H
    topo[0, :] = -10.0
    topo[Nr - 2, :] = -10.0
    if not cyclic:
        topo[:, 0] = -10.0
        topo[:, Nt - 2] = -10.0

    topo = np.ma.masked_where(topo == -10, topo)
    domain = pygetm.domain.create_cartesian(
        xx[:, :],
        yx[:, :],
        interfaces=True,
        z0=z0,
        lat=lat,
        H=topo,
        periodic_x=True,
    )
    domain.to_xarray().to_netcdf("bathymetry.nc")

    if plot:
        plt.figure(figsize=(5, 5))
        plt.scatter(xx.ravel(), yx.ravel(), s=4, c="steelblue", alpha=0.6)
        plt.gca().set_aspect("equal")
        plt.title("Polar annulus grid (Nr={}, Nt={})".format(Nr, Nt))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True, linestyle=":")

        domain.plot(show_mesh=True)
        plt.savefig("domain_mesh.png")
        domain.plot(show_mask=True)
        plt.savefig("domain_mask.png")

    return domain


def sim_analytical(domain, N, taux, tauy, plot=False):
    sim = pygetm.Simulation(
        domain,
        runtype=pygetm.RunType.BAROTROPIC_2D,
        airsea=pygetm.airsea.Fluxes(taux=taux, tauy=tauy),
        # vertical_coordinates=pygetm.vertical_coordinates.Sigma(args.nz),
    )

    dt = np.floor(0.8 * domain.cfl_check())
    sim.logger.info(f"Used timestep: {dt}")

    sim.airsea.taux.set(taux)
    sim.airsea.tauy.set(tauy)

    _path = Path("analytical.nc")
    output = sim.output_manager.add_netcdf_file(
        str(_path),
        interval_units=pygetm.TimeUnit.TIMESTEPS,
        interval=1,
        sync_interval=None,
        default_dtype=np.float32,
        save_initial=True,
    )
    output.request(
        "Ht",
        "xx",
        "yx",
        "rotationx",
        "zt",
        "U",
        "V",
        "u1",
        "v1",
        "tausx",
        "tausy",
    )

    if plot:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_aspect("equal")
        pc = ax.pcolormesh(sim.X.x, sim.X.y, sim.T.z, vmin=-0.0001, vmax=0.0001)
        cb = fig.colorbar(pc)
        cb.set_label("elevation (m)")
        u = sim.momentum.u1.interp(sim.T)
        v = sim.momentum.v1.interp(sim.T)
        Q = ax.quiver(
            sim.T.x[::1, ::1],
            sim.T.y[::1, ::1],
            u[::1, ::1],
            v[::1, ::1],
            # scale=0.005,
        )
        title = ax.set_title("time: 0 s")

    sim.start(
        datetime(2000, 1, 1),
        dt,
        report=N // 100,
        report_totals=N // 10,
    )
    for n in range(N):
        sim.advance()

    if plot:
        if True:
            u = sim.momentum.u1.interp(sim.T)
            v = sim.momentum.v1.interp(sim.T)
        else:
            u = sim.airsea.taux.interp(sim.T)
            v = sim.airsea.tauy.interp(sim.T)
        Q.set_UVC(u[::1, ::1], v[::1, ::1])
        title.set_text(f"time: {(n + 1) * dt} s")
        pc.set_array(sim.T.z[...].ravel())
        plt.show()
    return sim


if __name__ == "__main__":
    args = parse_args()
    domain = create_domain(
        args.Ri, args.Ro, args.Nr, args.Nt, cyclic=not args.no_cyclic, plot=args.plot
    )
    if args.radial_wind:
        print("radial wind")
    else:
        sim_analytical(domain, args.N, args.taux, args.tauy, plot=args.plot)
