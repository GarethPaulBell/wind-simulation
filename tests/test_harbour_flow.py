"""
tests/test_harbour_flow.py
==========================
Unit tests for wind_simulation.harbour_flow.

These tests use a tiny grid to keep runtime short and avoid GPU requirements.
They verify the public API contracts (shapes, types, value ranges) rather
than numerical accuracy.
"""

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Helpers (import target module lazily so we can test even if XLB unavailable)
# ---------------------------------------------------------------------------

def _import_module():
    """Import harbour_flow, skipping if XLB is not installed."""
    pytest.importorskip("xlb", reason="xlb is required for harbour_flow tests")
    from wind_simulation import harbour_flow
    return harbour_flow


# ---------------------------------------------------------------------------
# make_circular_obstacle_mask
# ---------------------------------------------------------------------------

class TestMakeCircularObstacleMask:
    def test_returns_bool_array(self):
        m = _import_module()
        mask = m.make_circular_obstacle_mask(64, 32, 32, 16, 5)
        assert mask.dtype == bool

    def test_shape(self):
        m = _import_module()
        mask = m.make_circular_obstacle_mask(64, 32, 32, 16, 5)
        assert mask.shape == (32, 64)  # (ny, nx)

    def test_centre_is_solid(self):
        m = _import_module()
        mask = m.make_circular_obstacle_mask(64, 32, 32, 16, 5)
        assert mask[16, 32]  # centre cell must be True

    def test_far_corner_is_fluid(self):
        m = _import_module()
        mask = m.make_circular_obstacle_mask(64, 32, 32, 16, 5)
        assert not mask[0, 0]   # top-left corner is far from centre

    def test_radius_boundary(self):
        m = _import_module()
        # A cell exactly at distance r should be solid (≤ r²)
        mask = m.make_circular_obstacle_mask(100, 100, 50, 50, 10)
        assert mask[50, 60]   # (50-50)² + (60-50)² = 100 = r² → solid
        assert not mask[50, 61]  # (61-50)² = 121 > 100 → fluid


# ---------------------------------------------------------------------------
# make_domain
# ---------------------------------------------------------------------------

class TestMakeDomain:
    def setup_method(self):
        m = _import_module()
        self.domain = m.make_domain(32, 16)

    def test_keys(self):
        assert set(self.domain.keys()) == {"left", "right", "bottom", "top"}

    def test_left_x_all_zero(self):
        xs = self.domain["left"][0]
        assert all(x == 0 for x in xs)

    def test_right_x_all_nx_minus_one(self):
        xs = self.domain["right"][0]
        assert all(x == 31 for x in xs)

    def test_bottom_y_all_zero(self):
        ys = self.domain["bottom"][1]
        assert all(y == 0 for y in ys)

    def test_top_y_all_ny_minus_one(self):
        ys = self.domain["top"][1]
        assert all(y == 15 for y in ys)

    def test_no_overlap_between_faces(self):
        """No (x, y) pair should appear in more than one face."""
        all_pairs: set = set()
        for face_idx in self.domain.values():
            pairs = set(zip(face_idx[0], face_idx[1]))
            # Check this face has no overlap with previously seen pairs
            assert pairs.isdisjoint(all_pairs), "Duplicate indices found between faces"
            all_pairs.update(pairs)

    def test_inlet_excludes_corners(self):
        """Inlet (left face) must not contain y=0 or y=ny-1 (those belong to walls)."""
        ys = self.domain["left"][1]
        assert 0 not in ys
        assert 15 not in ys


# ---------------------------------------------------------------------------
# run_simulation  (tiny grid, few steps)
# ---------------------------------------------------------------------------

class TestRunSimulation:
    """Integration-level test on a 32×16 grid for 20 steps."""

    @pytest.fixture(scope="class")
    def results(self):
        m = _import_module()
        # Use a small obstacle radius appropriate for this tiny test grid
        cfg = m.SimConfig(
            nx=32, ny=16, num_steps=20, print_every=10,
            save_arrays=False, obs_radius=3,
        )
        return m.run_simulation(cfg)

    def test_result_keys(self, results):
        assert set(results.keys()) == {"rho", "ux", "uy", "vel_mag", "vorticity"}

    def test_shapes(self, results):
        ny, nx = 16, 32
        for key, arr in results.items():
            assert arr.shape == (ny, nx), f"{key} has wrong shape {arr.shape}"

    def test_vel_mag_non_negative(self, results):
        vel_mag = results["vel_mag"]
        assert np.nanmin(vel_mag) >= 0.0

    def test_rho_near_one(self, results):
        rho = results["rho"]
        # Density should stay close to 1 for incompressible LBM
        assert np.nanmin(rho) > 0.8
        assert np.nanmax(rho) < 1.5

    def test_vel_mag_stable(self, results):
        """Velocity magnitude must not blow up (< 0.5 in lattice units)."""
        assert np.nanmax(results["vel_mag"]) < 0.5


# ---------------------------------------------------------------------------
# plot_velocity / plot_vorticity  (file creation only)
# ---------------------------------------------------------------------------

class TestPlotFunctions:
    def _dummy_results(self):
        ny, nx = 16, 32
        ux = np.random.rand(ny, nx) * 0.05
        uy = np.random.rand(ny, nx) * 0.01
        return {
            "rho": np.ones((ny, nx)),
            "ux": ux,
            "uy": uy,
            "vel_mag": np.sqrt(ux**2 + uy**2),
            "vorticity": np.random.randn(ny, nx) * 0.01,
        }

    def test_plot_velocity_creates_file(self, tmp_path):
        m = _import_module()
        results = self._dummy_results()
        out = str(tmp_path / "vel.png")
        m.plot_velocity(results, out)
        import os
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

    def test_plot_vorticity_creates_file(self, tmp_path):
        m = _import_module()
        results = self._dummy_results()
        out = str(tmp_path / "vort.png")
        m.plot_vorticity(results, out)
        import os
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0
