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
| Python package (`wind_simulation`) | Scaffold only (v0.0.1) |
| Documentation site | Auto-deployed from notebooks via Quarto |
| PyPI release | Not yet published |

The project is developed using [nbdev](https://nbdev.fast.ai/) — all working code lives in the `nbs/` notebooks and is exported to the `wind_simulation` package automatically.

## Technical Approach

- **Computational Fluid Dynamics**: [Lattice Boltzmann Method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) via the [XLB](https://github.com/Autodesk/XLB) library (D3Q19 velocity set / `KBCSim`).
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

## Data

LiDAR DEM data for Auckland North (1 m resolution, 2016-2018) is available from [LINZ Data Service](https://data.linz.govt.nz/). Place downloaded GeoTIFF tiles under:

```
data/raw/lds-auckland-north-lidar-1m-dem-2016-2018-GTiff/
```

## Project Structure

```
wind-simulation/
├── nbs/                    # Jupyter notebooks (source of truth)
│   ├── index.ipynb         # Project index / documentation home
│   ├── 00_condition_data.ipynb
│   ├── 01_car_example.ipynb
│   └── 02_define_boundary_conditions.ipynb
├── wind_simulation/        # Auto-generated Python package
│   ├── __init__.py
│   └── core.py
├── data/raw/               # Raw LiDAR GeoTIFF tiles (not committed)
├── settings.ini            # nbdev project configuration
└── setup.py
```

## Documentation

Full documentation is published automatically from the notebooks to GitHub Pages:
https://GarethPaulBell.github.io/wind-simulation

## Contributing

This project uses [nbdev](https://nbdev.fast.ai/). All source changes should be made in the `nbs/` notebooks; running `nbdev_export` updates the Python package.

```sh
pip install nbdev
nbdev_install_hooks   # set up git hooks
nbdev_export          # export notebooks to package
nbdev_test            # run notebook tests
```

## License

[Apache 2.0](LICENSE) © 2024 Gareth Paul Bell
