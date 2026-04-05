# Reproducing Optimization on X299-SSD

This documents the exact local workflow used so far for the KrakenOS merit-function optimization work on branch `feat/merit-function-optimization`.

## Repositories

Required local checkouts:

- `~/Projects/Kraken-Optical-Simulator`
- `~/Projects/pagmo2`
- `~/Projects/pygmo2`

## 1. Enter the Kraken development environment

From the Kraken repo:

```bash
cd ~/Projects/Kraken-Optical-Simulator
direnv allow
devenv shell
```

The current `devenv.nix` provides:

- Python 3.13
- `cmake`
- `ninja`
- `gcc`
- `pkg-config`
- `boost`
- `tbb`
- Python packages needed by the local pygmo build:
  - `cloudpickle`
  - `pybind11`

It also adds this runtime path automatically when present:

```bash
/home/thinky/Projects/pagmo2/_install/lib64
```

## 2. Build local pagmo2

Build and install pagmo into a local prefix:

```bash
nix shell \
  nixpkgs#cmake \
  nixpkgs#ninja \
  nixpkgs#gcc \
  nixpkgs#pkg-config \
  nixpkgs#boost \
  nixpkgs#boost.dev \
  nixpkgs#tbb \
  nixpkgs#tbb.dev \
  -c bash -lc '
    export CMAKE_PREFIX_PATH=/nix/store/msvfnfpsp60jlmnyddcbakky4qxjd2qg-onetbb-2022.3.0-dev:/nix/store/adv73gqmglkg48w5ql70s778ba6wpbs4-boost-1.89.0-dev${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}
    export PKG_CONFIG_PATH=/nix/store/msvfnfpsp60jlmnyddcbakky4qxjd2qg-onetbb-2022.3.0-dev/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
    cd ~/Projects/pagmo2
    rm -rf build_local _install
    cmake -S . -B build_local -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PWD/_install
    cmake --build build_local -j2
    cmake --install build_local
  '
```

Expected output install root:

```bash
~/Projects/pagmo2/_install
```

## 3. Build local pygmo2 against that pagmo install

Build and install pygmo into the Kraken venv site-packages:

```bash
nix shell \
  nixpkgs#cmake \
  nixpkgs#ninja \
  nixpkgs#gcc \
  nixpkgs#pkg-config \
  nixpkgs#boost \
  nixpkgs#boost.dev \
  nixpkgs#tbb \
  nixpkgs#tbb.dev \
  nixpkgs#python313Packages.pybind11 \
  -c bash -lc '
    cd ~/Projects/Kraken-Optical-Simulator
    VENV=$PWD/.devenv/state/venv
    $VENV/bin/python -m pip install cloudpickle

    PYBIND11_DIR=/nix/store/qsrcf2wc72sn2yvb6zx72mqz0kz390ff-python3.13-pybind11-3.0.1/lib/python3.13/site-packages/pybind11/share/cmake/pybind11
    PAGMO_PREFIX=~/Projects/pagmo2/_install
    export CMAKE_PREFIX_PATH=$PAGMO_PREFIX:/nix/store/msvfnfpsp60jlmnyddcbakky4qxjd2qg-onetbb-2022.3.0-dev:/nix/store/adv73gqmglkg48w5ql70s778ba6wpbs4-boost-1.89.0-dev:$PYBIND11_DIR${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}
    export PKG_CONFIG_PATH=/nix/store/msvfnfpsp60jlmnyddcbakky4qxjd2qg-onetbb-2022.3.0-dev/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}

    cd ~/Projects/pygmo2
    rm -rf build_local
    cmake -S . -B build_local -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DPython3_EXECUTABLE=$VENV/bin/python \
      -DPYGMO_INSTALL_PATH=$VENV/lib/python3.13/site-packages \
      -DCMAKE_INSTALL_PREFIX=$PWD/_install
    cmake --build build_local -j2
    cmake --install build_local
  '
```

## 4. Verify pygmo is importable from the Kraken environment

```bash
cd ~/Projects/Kraken-Optical-Simulator
devenv shell -- bash -lc 'python - <<\"PY\"
import pygmo as pg
print(pg.__version__)
PY'
```

Expected:

```text
2.19.7
```

## 5. Run the Kraken optimization example

```bash
cd ~/Projects/Kraken-Optical-Simulator
devenv shell -- bash -lc 'python KrakenOS/Optimization/examples/optimize_doublet_spot.py'
```

This performs:

- initial merit evaluation
- SciPy differential evolution run
- real `pygmo.de` run
- JSON result export to:
  - `doublet_merit_result.json`

## 6. Launch the Kraken layout editor with optimization controls

```bash
cd ~/Projects/Kraken-Optical-Simulator
devenv shell -- bash -lc 'python -m KrakenOS.UI.layout_editor'
```

Current optimization UI behavior:

- mark optimization variables by toggling `Rc*` and `T*`
- or select rows and use `Mark Selected`
- choose merit function:
  - `Spot RMS`
  - `Wavefront RMS`
  - `Spot + Wavefront`
- click `Start Optimization`

## Notes

- The current optimizer backend wired into the UI is `pygmo.de`.
- The current bounds are heuristic and derived from the current prescription values.
- The current UI supports optimization of:
  - `Rc`
  - `Thickness`
- `Thin Lens` focal-length optimization is not wired yet.

## Known unrelated local files

These were intentionally left untouched by the optimization work:

- `KrakenOS/Examples/mi_objeto.pkl`
- `PR_GUIDE.md`
