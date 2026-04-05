# wind-simulation

Create visualizations of wind flows around typical yacht racing areas using Lattice Boltzmann Method (LBM) computational fluid dynamics.

> **Status**: Early Alpha (v0.0.1) — proof-of-concept stage. Core workflows are demonstrated in Jupyter notebooks; the installable Python package is a minimal scaffold that will grow as the project matures.

## Purpose

Wind simulation is a tool for modelling how wind flows over real terrain in yacht racing areas. Starting from publicly available LiDAR Digital Elevation Model (DEM) data, it builds a 3-D voxel representation of the terrain and then runs a Lattice Boltzmann fluid simulation to visualise wind patterns. The initial focus area is Auckland, New Zealand, using 1 m-resolution LiDAR data from the New Zealand Land Information System (LINZ).

The same pipeline can be applied to arbitrary 3-D objects (e.g. boat hulls, harbour structures) provided as STL CAD files, making it useful for general aerodynamic prototyping as well.

## Current State

| Area | Status |
|---|---|
| Terrain loading (LiDAR GeoTIFF) | Working (notebook `00_condition_data`) |
| Car / object aerodynamics example | Working (notebook `01_car_example`) |
| Boundary condition definitions | Documented (notebook `02_define_boundary_conditions`) |
| 2-D harbour-flow LBM validation (D2Q9 / BGK) | Working (`wind_simulation/harbour_flow.py`) |
| Python package (`wind_simulation`) | Scaffold + harbour-flow module (v0.0.1) |
| Test suite | Working (`tests/test_harbour_flow.py`, pytest) |
| Documentation site | Auto-deployed from notebooks via Quarto |
| PyPI release | Not yet published |

The project is developed using [nbdev](https://nbdev.fast.ai/) — all working code lives in the `nbs/` notebooks and is exported to the `wind_simulation` package automatically. A standalone 2-D validation module (`harbour_flow.py`) lives directly in the package and is covered by a pytest test suite.

## Technical Approach

- **Computational Fluid Dynamics**: [Lattice Boltzmann Method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) via the [XLB](https://github.com/Autodesk/XLB) library.
  - 2-D validation: D2Q9 velocity set, BGK collision, ZouHe velocity inlet, ExtrapolationOutflow outlet, FullwayBounceback walls and obstacles.
  - 3-D terrain simulation: D3Q19 velocity set / `KBCSim` (planned).
- **Terrain data**: GeoTIFF rasters loaded and merged with [rasterio](https://rasterio.readthedocs.io/) / GDAL; 1 m Auckland North LiDAR DEM (2016-2018).
- **3-D object support**: STL files voxelised with [trimesh](https://trimesh.org/) for solid boundary conditions.
- **Simulation domain**: Configurable box with named faces — North, South, East, West, Top, Bottom — where wind inlet/outlet and wall conditions are assigned.
- **Visualisation**: matplotlib + rasterio plotting for terrain and flow-field outputs.

## Install

```sh
pip install wind_simulation
```

> **Note**: The PyPI package currently contains only the package scaffold. Full functionality requires running the notebooks directly (see [How to use](#how-to-use)).

### Development install (recommended)

```sh
git clone https://github.com/GarethPaulBell/wind-simulation.git
cd wind-simulation
pip install -e ".[dev]"
```

Key dependencies (install manually if needed):

```sh
pip install rasterio numpy matplotlib trimesh xlb
```

You will also need GDAL installed in your environment for rasterio to handle GeoTIFF files.

## How to use

The main workflows are exposed as numbered notebooks in the `nbs/` folder:

### 1. Load and condition terrain data (`00_condition_data.ipynb`)

Reads one or more LiDAR GeoTIFF tiles, merges them into a single raster, and plots the resulting elevation model.

```python
import rasterio
from rasterio.merge import merge

# open tiles
src1 = rasterio.open("data/raw/.../tile_a.tif")
src2 = rasterio.open("data/raw/.../tile_b.tif")

mosaic, transform = merge([src1, src2])
```

### 2. Run an aerodynamic simulation (`01_car_example.ipynb`)

Demonstrates the full LBM pipeline using a car STL model as the solid obstacle.

```python
from wind_simulation.core import Car

sim = Car(shape=(128, 64, 64), velocity_set="D3Q19")
sim.voxelize_stl("path/to/car.stl")
sim.run(num_steps=1000)
sim.plot_streamlines()
```

### 3. Define boundary conditions (`02_define_boundary_conditions.ipynb`)

Shows how to configure inlet (wind source), outlet, and solid wall boundaries for a custom simulation domain.

### 4. Run the 2-D harbour-flow validation

A self-contained 2-D LBM validation case is provided as both a Python module and a runnable script. It simulates flow around a circular bluff body (circular obstacle) on a 512 × 256 lattice and saves velocity-magnitude and vorticity images.

**Quick start (default 5000 steps, 512 × 256 grid):**

```sh
python harbour_flow_sim.py
```

**Custom options:**

```sh
# Fewer steps (quick smoke-test)
python harbour_flow_sim.py --steps 200

# Custom grid size and output directory
python harbour_flow_sim.py --steps 1000 --nx 256 --ny 128 --output-dir /tmp/results
```

Output files are written to the `output/` directory (or the path given with `--output-dir`):

```
output/
├── velocity_magnitude.png
├── vorticity.png
├── rho.npy
├── ux.npy
├── uy.npy
├── vel_mag.npy
└── vorticity.npy
```

The module can also be used directly from Python:

```python
from wind_simulation.harbour_flow import SimConfig, run_simulation, plot_velocity, plot_vorticity

cfg = SimConfig(nx=512, ny=256, num_steps=5000, output_dir="output")
results = run_simulation(cfg)          # returns dict of numpy arrays
plot_velocity(results, "output/velocity_magnitude.png")
plot_vorticity(results, "output/vorticity.png")
```

The circular obstacle geometry is intentionally kept separate in `make_circular_obstacle_mask()`. To use a real coastline or topography mask, replace that function with one that loads a boolean array from a GeoTIFF or NumPy file — the rest of the pipeline is unchanged.

## Data

LiDAR DEM data for Auckland North (1 m resolution, 2016-2018) is available from [LINZ Data Service](https://data.linz.govt.nz/). Place downloaded GeoTIFF tiles under:

```
data/raw/lds-auckland-north-lidar-1m-dem-2016-2018-GTiff/
```

## Project Structure

```
wind-simulation/
├── nbs/                          # Jupyter notebooks (source of truth)
│   ├── index.ipynb               # Project index / documentation home
│   ├── 00_condition_data.ipynb
│   ├── 01_car_example.ipynb
│   └── 02_define_boundary_conditions.ipynb
├── wind_simulation/              # Python package
│   ├── __init__.py
│   ├── core.py                   # Auto-generated from notebooks (scaffold)
│   └── harbour_flow.py           # 2-D LBM validation module (D2Q9 / BGK)
├── tests/
│   └── test_harbour_flow.py      # pytest tests for harbour_flow module
├── harbour_flow_sim.py           # CLI entry point for 2-D validation run
├── data/raw/                     # Raw LiDAR GeoTIFF tiles (not committed)
├── settings.ini                  # nbdev project configuration
└── setup.py
```

## Testing

A pytest test suite covers the `harbour_flow` module (geometry helpers, domain builder, simulation integration, and plot functions):

```sh
pip install pytest
python -m pytest tests/test_harbour_flow.py
```

The tests use a tiny grid to keep runtime short and avoid GPU requirements. `xlb` must be installed (see [Key dependencies](#development-install-recommended)); tests are skipped automatically if it is not available.

## Documentation

Full documentation is published automatically from the notebooks to GitHub Pages:
https://GarethPaulBell.github.io/wind-simulation

## Contributing

This project uses [nbdev](https://nbdev.fast.ai/). All source changes to the notebook-driven code should be made in the `nbs/` notebooks; running `nbdev_export` updates the Python package. The `harbour_flow` module and its tests live directly in `wind_simulation/` and `tests/` respectively and are edited as regular Python files.

```sh
pip install nbdev
nbdev_install_hooks   # set up git hooks
nbdev_export          # export notebooks to package
nbdev_test            # run notebook tests

# Run the standalone pytest suite
python -m pytest tests/
```

## License

[Apache 2.0](LICENSE) © 2024 Gareth Paul Bell
