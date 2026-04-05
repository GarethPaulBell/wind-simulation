"""
harbour_flow.py
===============
2-D external-flow validation case using XLB (Lattice Boltzmann Method).

Purpose
-------
This is the **first validation step** for the wind-simulation pipeline.
It is *not* a realistic Auckland Harbour model.  The goal is to confirm that
the solver, GPU/CPU execution, boundary conditions, obstacle masking, and
output visualisation all work correctly on a coarse grid.

Domain: 512 × 256 cells (left-to-right flow)
Obstacle: one circular bluff body near x = 0.35·Nx, y = 0.50·Ny, radius ≈ 20

Design note – geometry separation
----------------------------------
The obstacle mask is produced by ``make_circular_obstacle_mask()`` and kept
entirely separate from the solver setup.  To swap in a rasterised coastline
or topography mask later, replace that function call with one that loads an
array from a GeoTIFF or NumPy file, then pass the resulting boolean mask
array to ``_indices_from_mask()``.

References
----------
XLB library: https://github.com/Autodesk/XLB
D2Q9 velocity set, BGK collision, ZouHe velocity inlet,
ExtrapolationOutflow outlet, FullwayBounceback walls/obstacle.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless / server use
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# XLB imports
# ---------------------------------------------------------------------------
import xlb
from xlb.compute_backend import ComputeBackend
from xlb.grid import grid_factory
from xlb.operator.boundary_condition import (
    ExtrapolationOutflowBC,
    FullwayBounceBackBC,
    ZouHeBC,
)
from xlb.operator.macroscopic import Macroscopic
from xlb.operator.stepper import IncompressibleNavierStokesStepper
from xlb.precision_policy import PrecisionPolicy

# ---------------------------------------------------------------------------
# Public configuration dataclass
# ---------------------------------------------------------------------------

@dataclass
class SimConfig:
    """All tunable parameters in one place.  Edit here to change the run."""

    # Grid
    nx: int = 512          # domain width  (cells)
    ny: int = 256          # domain height (cells)

    # Obstacle (circular bluff body – the geometry that can later be replaced)
    obs_cx_frac: float = 0.35   # centre x as fraction of Nx
    obs_cy_frac: float = 0.50   # centre y as fraction of Ny
    obs_radius: int = 20        # radius in lattice cells

    # Physics
    inlet_vel: float = 0.04     # lattice-unit inlet speed (stable range: 0.03–0.08)
    omega: float = 1.7          # BGK relaxation  → ν = (1/ω - 0.5)/3

    # Simulation
    num_steps: int = 5000       # total timesteps
    print_every: int = 500      # diagnostic printout interval

    # Output
    output_dir: str = "output"  # directory for images and arrays
    save_arrays: bool = True    # also save numpy field arrays

    # Backend  (WARP runs on GPU if CUDA is available, otherwise falls back to CPU)
    backend: ComputeBackend = ComputeBackend.WARP
    precision: PrecisionPolicy = PrecisionPolicy.FP32FP32


# ---------------------------------------------------------------------------
# Helper: index extraction from boolean mask
# ---------------------------------------------------------------------------

def _indices_from_mask(mask: np.ndarray) -> list[list[int]]:
    """Return [x_list, y_list] index lists for all True cells in *mask*.

    Parameters
    ----------
    mask : 2-D boolean array of shape (ny, nx)  (row = y, col = x)

    Returns
    -------
    list of two lists  [x_indices, y_indices]  – the format expected by
    XLB boundary-condition constructors.
    """
    ys, xs = np.where(mask)
    return [xs.tolist(), ys.tolist()]


# ---------------------------------------------------------------------------
# Interface 1: make_domain
# ---------------------------------------------------------------------------

def make_domain(nx: int, ny: int) -> dict[str, Any]:
    """Return boundary-index lists for a rectangular 2-D domain.

    The four face indices are computed so that no grid node belongs to more
    than one face (corners are assigned to top/bottom walls; inlet/outlet
    span only the interior rows).

    Parameters
    ----------
    nx, ny : grid dimensions

    Returns
    -------
    dict with keys ``left``, ``right``, ``bottom``, ``top`` – each a
    ``[x_list, y_list]`` pair ready for XLB BC constructors.
    """
    # Interior row range (excludes corner nodes shared with top/bottom walls)
    interior_rows = list(range(1, ny - 1))

    domain = {
        # Inlet (left face, interior rows only – corners go to walls)
        "left": [[0] * len(interior_rows), interior_rows],
        # Outlet (right face, interior rows only)
        "right": [[nx - 1] * len(interior_rows), interior_rows],
        # No-slip bottom wall (full width including corners)
        "bottom": [list(range(nx)), [0] * nx],
        # Free-stream / no-slip top wall (full width including corners)
        "top": [list(range(nx)), [ny - 1] * nx],
    }
    return domain


# ---------------------------------------------------------------------------
# Interface 2: make_circular_obstacle_mask
# ---------------------------------------------------------------------------

def make_circular_obstacle_mask(
    nx: int,
    ny: int,
    cx: int,
    cy: int,
    r: int,
) -> np.ndarray:
    """Return a boolean 2-D mask with True inside the circular obstacle.

    Parameters
    ----------
    nx, ny : domain dimensions
    cx, cy : obstacle centre in lattice cells
    r      : obstacle radius in lattice cells

    Returns
    -------
    mask : shape (ny, nx), dtype bool – True where the obstacle is solid.

    Replacement point
    -----------------
    To use a rasterised coastline or topography mask instead, replace this
    function with one that loads your geometry array and returns it as a
    (ny, nx) boolean array.  The rest of the pipeline is unchanged.
    """
    Y, X = np.mgrid[0:ny, 0:nx]
    mask = (X - cx) ** 2 + (Y - cy) ** 2 <= r ** 2
    return mask


# ---------------------------------------------------------------------------
# Interface 3: run_simulation
# ---------------------------------------------------------------------------

def run_simulation(cfg: SimConfig) -> dict[str, np.ndarray]:
    """Run the LBM simulation and return macroscopic field arrays.

    Parameters
    ----------
    cfg : SimConfig instance

    Returns
    -------
    dict with keys:
        ``rho``       – density field, shape (ny, nx)
        ``ux``        – x-velocity field, shape (ny, nx)
        ``uy``        – y-velocity field, shape (ny, nx)
        ``vel_mag``   – velocity magnitude, shape (ny, nx)
        ``vorticity`` – 2-D vorticity (ω_z = ∂u_y/∂x − ∂u_x/∂y), shape (ny, nx)
    """
    nx, ny = cfg.nx, cfg.ny

    # ------------------------------------------------------------------
    # XLB initialisation
    # ------------------------------------------------------------------
    velocity_set = xlb.velocity_set.D2Q9(
        precision_policy=cfg.precision,
        compute_backend=cfg.backend,
    )
    xlb.init(
        velocity_set=velocity_set,
        default_backend=cfg.backend,
        default_precision_policy=cfg.precision,
    )

    # ------------------------------------------------------------------
    # Domain boundary indices
    # ------------------------------------------------------------------
    domain = make_domain(nx, ny)

    # ------------------------------------------------------------------
    # Obstacle mask → boundary indices
    # ------------------------------------------------------------------
    obs_cx = int(cfg.obs_cx_frac * nx)
    obs_cy = int(cfg.obs_cy_frac * ny)
    obs_mask = make_circular_obstacle_mask(nx, ny, obs_cx, obs_cy, cfg.obs_radius)

    # Exclude domain-boundary rows/columns from the obstacle mask so that no
    # obstacle cell overlaps with the inlet, outlet, or wall BC indices.
    obs_mask[:, 0]    = False   # left face (inlet)
    obs_mask[:, -1]   = False   # right face (outlet)
    obs_mask[0, :]    = False   # bottom wall
    obs_mask[-1, :]   = False   # top wall
    obs_idx = _indices_from_mask(obs_mask)
    print(
        f"Obstacle: centre=({obs_cx}, {obs_cy}), radius={cfg.obs_radius}, "
        f"solid cells={len(obs_idx[0])}"
    )

    # ------------------------------------------------------------------
    # Boundary conditions
    # ------------------------------------------------------------------
    # Inlet: uniform x-velocity via ZouHe velocity BC
    bc_inlet = ZouHeBC(
        bc_type="velocity",
        prescribed_value=(cfg.inlet_vel, 0.0),
        indices=domain["left"],
    )
    # Outlet: extrapolation / open outflow (avoids strong reflections)
    bc_outlet = ExtrapolationOutflowBC(indices=domain["right"])
    # Bottom: no-slip wall
    bc_bottom = FullwayBounceBackBC(indices=domain["bottom"])
    # Top: no-slip wall (use FullwayBounceBack for simplicity; swap to
    #      a slip/free-stream BC if XLB supports one in a future version)
    bc_top = FullwayBounceBackBC(indices=domain["top"])
    # Obstacle: no-slip bounce-back on circular bluff body
    bc_obs = FullwayBounceBackBC(indices=obs_idx)

    # Order matters for overlap check: walls first, then inlet/outlet, then obstacle
    all_bcs = [bc_bottom, bc_top, bc_inlet, bc_outlet, bc_obs]

    # ------------------------------------------------------------------
    # Grid and stepper
    # ------------------------------------------------------------------
    grid = grid_factory((nx, ny), compute_backend=cfg.backend)
    stepper = IncompressibleNavierStokesStepper(
        grid,
        boundary_conditions=all_bcs,
        collision_type="BGK",
    )
    f_0, f_1, bc_mask, missing_mask = stepper.prepare_fields()

    # Macroscopic operator (reused each printout)
    macro = Macroscopic(
        velocity_set=velocity_set,
        precision_policy=cfg.precision,
        compute_backend=cfg.backend,
    )

    # ------------------------------------------------------------------
    # Simulation loop
    # ------------------------------------------------------------------
    print(
        f"\nRunning {cfg.num_steps} steps on "
        f"{cfg.backend.name} backend (ω={cfg.omega}, u_inlet={cfg.inlet_vel})"
    )
    # Kinematic viscosity: ν = (1/ω − 0.5) / 3  (D2Q9 lattice units)
    # Reynolds number based on obstacle diameter as length scale: Re = U · 2r / ν
    nu = (1.0 / cfg.omega - 0.5) / 3.0
    re = cfg.inlet_vel * 2 * cfg.obs_radius / nu
    print(f"Grid: {nx}×{ny}  | ν={nu:.4f}  | Re ≈ {re:.0f}\n")

    for step in range(cfg.num_steps):
        f_0, f_1 = stepper(f_0, f_1, bc_mask, missing_mask, cfg.omega, step)
        f_0, f_1 = f_1, f_0  # swap double buffers

        if (step + 1) % cfg.print_every == 0:
            _print_diagnostics(f_0, grid, velocity_set, cfg, step + 1, macro)

    # ------------------------------------------------------------------
    # Extract final macroscopic fields
    # ------------------------------------------------------------------
    rho_field = grid.create_field(cardinality=1, fill_value=0.0)
    u_field = grid.create_field(cardinality=velocity_set.d, fill_value=0.0)
    rho_field, u_field = macro(f_0, rho_field, u_field)

    # Convert from Warp arrays (cardinality, nx, ny, 1) → numpy (ny, nx)
    rho_np = np.array(rho_field.numpy()[0, :, :, 0]).T   # (ny, nx)
    ux_np  = np.array(u_field.numpy()[0, :, :, 0]).T     # (ny, nx)
    uy_np  = np.array(u_field.numpy()[1, :, :, 0]).T     # (ny, nx)

    vel_mag = np.sqrt(ux_np ** 2 + uy_np ** 2)

    # 2-D vorticity: ω_z = ∂u_y/∂x − ∂u_x/∂y  (central differences)
    vorticity = _compute_vorticity_2d(ux_np, uy_np)

    # Zero out obstacle interior for cleaner images
    # obs_mask is (ny, nx) — matches the numpy field layout
    rho_np[obs_mask]     = np.nan
    ux_np[obs_mask]      = np.nan
    uy_np[obs_mask]      = np.nan
    vel_mag[obs_mask]    = np.nan
    vorticity[obs_mask]  = np.nan

    return {
        "rho": rho_np,
        "ux": ux_np,
        "uy": uy_np,
        "vel_mag": vel_mag,
        "vorticity": vorticity,
    }


# ---------------------------------------------------------------------------
# Diagnostics helper
# ---------------------------------------------------------------------------

def _print_diagnostics(
    f_0,
    grid,
    velocity_set,
    cfg: SimConfig,
    step: int,
    macro: Macroscopic,
) -> None:
    """Print min/max density and velocity magnitude for monitoring."""
    rho_tmp = grid.create_field(cardinality=1, fill_value=0.0)
    u_tmp   = grid.create_field(cardinality=velocity_set.d, fill_value=0.0)
    rho_tmp, u_tmp = macro(f_0, rho_tmp, u_tmp)

    rho_np = rho_tmp.numpy()
    u_np   = u_tmp.numpy()
    ux = u_np[0, :, :, 0]
    uy = u_np[1, :, :, 0]
    vmag = np.sqrt(ux ** 2 + uy ** 2)

    print(
        f"  step {step:6d}/{cfg.num_steps}  |  "
        f"ρ [{rho_np.min():.4f}, {rho_np.max():.4f}]  |  "
        f"|u| [{vmag.min():.4f}, {vmag.max():.4f}]"
    )


# ---------------------------------------------------------------------------
# Vorticity helper (2-D, pure NumPy)
# ---------------------------------------------------------------------------

def _compute_vorticity_2d(ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    """Compute 2-D vorticity ω_z = ∂u_y/∂x − ∂u_x/∂y using central diffs.

    Parameters
    ----------
    ux, uy : velocity components, shape (ny, nx)

    Returns
    -------
    vorticity : shape (ny, nx)
    """
    # np.gradient returns derivatives w.r.t. the given axis
    # axis=1 → x-direction, axis=0 → y-direction
    duy_dx = np.gradient(uy, axis=1)
    dux_dy = np.gradient(ux, axis=0)
    return duy_dx - dux_dy


# ---------------------------------------------------------------------------
# Interface 4: plot_velocity
# ---------------------------------------------------------------------------

def plot_velocity(results: dict[str, np.ndarray], out_path: str) -> None:
    """Save a colour-map image of the velocity magnitude field.

    Parameters
    ----------
    results  : dict returned by ``run_simulation()``
    out_path : full file path for the output image (e.g. ``output/velocity.png``)
    """
    vel_mag = results["vel_mag"]
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        vel_mag,
        origin="lower",
        cmap="viridis",
        interpolation="bilinear",
    )
    plt.colorbar(im, ax=ax, label="Velocity magnitude (lattice units)")
    ax.set_title("Velocity magnitude – 2-D harbour flow validation")
    ax.set_xlabel("x (cells)")
    ax.set_ylabel("y (cells)")
    _save_fig(fig, out_path)


# ---------------------------------------------------------------------------
# Interface 5: plot_vorticity
# ---------------------------------------------------------------------------

def plot_vorticity(results: dict[str, np.ndarray], out_path: str) -> None:
    """Save a diverging colour-map image of the 2-D vorticity field.

    Parameters
    ----------
    results  : dict returned by ``run_simulation()``
    out_path : full file path for the output image (e.g. ``output/vorticity.png``)
    """
    vort = results["vorticity"]
    # symmetric colour scale
    vmax = np.nanpercentile(np.abs(vort), 99)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        vort,
        origin="lower",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="bilinear",
    )
    plt.colorbar(im, ax=ax, label="Vorticity ω_z (lattice units)")
    ax.set_title("Vorticity ω_z – 2-D harbour flow validation")
    ax.set_xlabel("x (cells)")
    ax.set_ylabel("y (cells)")
    _save_fig(fig, out_path)


# ---------------------------------------------------------------------------
# Internal: figure save helper
# ---------------------------------------------------------------------------

def _save_fig(fig: plt.Figure, out_path: str) -> None:
    parent = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(parent, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")
