#!/usr/bin/env python3
"""
harbour_flow_sim.py
===================
Runnable entry point for the 2-D harbour-flow LBM validation case.

Quick start
-----------
    python harbour_flow_sim.py

With custom output directory::

    python harbour_flow_sim.py --output-dir /tmp/my_results

With fewer steps (useful for a quick smoke-test)::

    python harbour_flow_sim.py --steps 200

See ``wind_simulation/harbour_flow.py`` for all configurable parameters.

What this script does
---------------------
1. Sets up a 512×256 lattice with left-to-right flow.
2. Inserts a single circular bluff body (the placeholder geometry).
3. Runs the BGK-LBM simulation for ``num_steps`` timesteps.
4. Saves velocity-magnitude and vorticity images plus optional NumPy arrays.

Next steps (not implemented here)
----------------------------------
- Replace ``make_circular_obstacle_mask()`` with a function that loads a
  rasterised coastline or topography mask from a GeoTIFF or NumPy file.
- Tune ``inlet_vel``, ``omega``, and ``num_steps`` for the target Reynolds
  number once real geometry is available.
- Extend to 3-D (D3Q19 or D3Q27) and atmospheric boundary-layer inflow.
"""

import argparse
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="2-D harbour-flow LBM validation (XLB / D2Q9)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="output",
        help="Directory for output images and arrays (default: output/)",
    )
    parser.add_argument(
        "--steps", "-n",
        type=int,
        default=None,
        help="Override number of timesteps (default: from SimConfig)",
    )
    parser.add_argument(
        "--nx",
        type=int,
        default=None,
        help="Override grid width (default: 512)",
    )
    parser.add_argument(
        "--ny",
        type=int,
        default=None,
        help="Override grid height (default: 256)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = _parse_args()

    # Import here so the module is only loaded after potential --help exit
    from wind_simulation.harbour_flow import (
        SimConfig,
        plot_velocity,
        plot_vorticity,
        run_simulation,
    )

    # Build config (CLI overrides take priority)
    cfg = SimConfig(output_dir=args.output_dir)
    if args.steps is not None:
        cfg.num_steps = args.steps
    if args.nx is not None:
        cfg.nx = args.nx
    if args.ny is not None:
        cfg.ny = args.ny

    print("=" * 60)
    print("  2-D Harbour-Flow LBM Validation  (XLB / D2Q9 / BGK)")
    print("=" * 60)
    print(f"  Grid        : {cfg.nx} × {cfg.ny}")
    print(f"  Steps       : {cfg.num_steps}")
    print(f"  Inlet vel   : {cfg.inlet_vel}")
    print(f"  Omega       : {cfg.omega}")
    print(f"  Output dir  : {cfg.output_dir}")
    print("=" * 60)

    # ------------------------------------------------------------------
    # Run simulation
    # ------------------------------------------------------------------
    results = run_simulation(cfg)

    # ------------------------------------------------------------------
    # Save output images
    # ------------------------------------------------------------------
    os.makedirs(cfg.output_dir, exist_ok=True)

    print("\nSaving outputs …")
    plot_velocity(results, os.path.join(cfg.output_dir, "velocity_magnitude.png"))
    plot_vorticity(results, os.path.join(cfg.output_dir, "vorticity.png"))

    # ------------------------------------------------------------------
    # Optionally save NumPy arrays
    # ------------------------------------------------------------------
    if cfg.save_arrays:
        for key, arr in results.items():
            path = os.path.join(cfg.output_dir, f"{key}.npy")
            np.save(path, arr)
            print(f"  Saved → {path}")

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    vel_mag = results["vel_mag"]
    vort    = results["vorticity"]
    print("\n--- Final field summary ---")
    print(f"  |u| max : {np.nanmax(vel_mag):.5f}")
    print(f"  |u| min : {np.nanmin(vel_mag):.5f}")
    print(f"  ω_z max : {np.nanmax(np.abs(vort)):.5f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
