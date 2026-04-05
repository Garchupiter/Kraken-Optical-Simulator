# Merit Function Optimization Plan

## Goal
Build a KrakenOS optimization workflow comparable to Zemax merit-function optimization.
The initial target is multi-variable optimization of optical prescriptions such as curvature and thickness against spot size, wavefront quality, and related paraxial or pupil objectives.

## Branching
- Base branch: `feat/layout-editor-ui`
- Optimization branch: `feat/merit-function-optimization`
- Keep optimization core separate from UI work initially.

## Local dependencies
- KrakenOS repo: `~/Projects/Kraken-Optical-Simulator`
- pagmo2 C++ repo: `~/Projects/pagmo2`
- pygmo2 Python implementation repo: `~/Projects/pygmo2`

## Immediate architectural decision
Do not begin with direct optimizer coupling.
First define a Kraken merit-evaluation layer with a stable API.
Then plug optimizers into that interface.

This keeps the design backend-agnostic and allows:
- a pure Python prototype,
- later `pygmo2` integration,
- optional future `pagmo2` C++ integration if needed.

## Target optimization model

### 1. Variables
Represent design variables explicitly.

Proposed object:
- `OpticalVariable`
  - `surface_index`
  - `parameter`
  - `lower_bound`
  - `upper_bound`
  - `scale`
  - `enabled`
  - optional `metadata`

Parameters to support first:
- `Rc`
- `Thickness`

Parameters to support later:
- `k`
- `AspherData[i]`
- `Diameter`
- tilts / decenters where meaningful

### 2. Merit operands
Represent each merit term independently.

Proposed object:
- `MeritOperand`
  - `kind`
  - `weight`
  - `target`
  - `surface_index`
  - `wavelength`
  - `field`
  - operand-specific settings

Initial operand set:
- spot RMS at image plane
- image centroid offset
- effective focal length target
- thickness penalty
- invalid trace penalty

Second-phase operand set:
- wavefront RMS
- wavefront P-V
- Seidel coefficient targets
- entrance / exit pupil targets
- chief-ray constraints

### 3. Merit function
Proposed object:
- `MeritFunction`
  - list of operands
  - penalties
  - aggregation strategy

Default aggregation:
- weighted sum of scalar operand residuals
- hard penalties for invalid systems or failed traces

## Core evaluator
Build a reusable evaluator that accepts a variable vector and returns a merit value.

Proposed API:
- `evaluate_merit(base_surfaces, setup, variables, operands, x) -> MeritResult`

Responsibilities:
- deep-copy the base prescription
- apply variable vector to copied surfaces
- rebuild Kraken system
- run required traces / analyses
- compute each operand
- aggregate total merit
- return diagnostics

Proposed result object:
- `MeritResult`
  - `total`
  - `operands`
  - `valid`
  - `message`
  - `prescription`

## Suggested package layout
Create a new package:
- `KrakenOS/Optimization/`

Suggested modules:
- `variables.py`
- `operands.py`
- `merit.py`
- `evaluator.py`
- `problems.py`
- `adapters/pygmo2_adapter.py`
- `examples/optimize_doublet_spot.py`

## First implementation milestone
Implement a complete end-to-end optimization without UI.

### Example problem
Use a doublet lens prescription.

Variables:
- surface 1 `Rc`
- surface 2 `Rc`
- surface 3 `Rc`
- one or more internal `Thickness`

Objective:
- minimize spot RMS at the image plane

Optional secondary term:
- hold effective focal length near target

Output:
- best variable vector
- before/after prescription
- before/after merit breakdown
- before/after spot plot

## Analysis pipeline choices

### Spot merit
Use Kraken trace results and image-plane intersections.
Candidate path:
- trace a controlled ray set
- use `raykeeper.pick(-1)`
- compute RMS via Kraken utilities or direct calculation

### Wavefront merit
Use Kraken pupil and phase tools.
Candidate path:
- `PupilCalc(...)`
- `Phase(...)`
- compute RMS or P-V from returned phase map

### Seidel merit
Use:
- `Seidel(Pupil)`

## Optimizer integration plan

### Phase 1: backend-agnostic
Implement evaluator and test it with a simple Python-side optimization loop.
This is only for validation.

### Phase 2: Python optimizer backend
Use the local `~/Projects/pygmo2` implementation as the first serious optimization backend.
Wrap the merit evaluator as a pygmo2-compatible problem interface.

Candidate algorithms for first trials:
- DE
- SADE
- CMA-ES if available and bounds are well-scaled

### Phase 3: advanced workflows
Add support for:
- multi-objective optimization
- constrained optimization
- restart strategies
- population logging

## Constraints and penalties
Support explicit constraints early.

Examples:
- minimum center thickness
- non-negative air gaps
- edge-thickness lower bound
- diameter / clear-aperture consistency
- maximum ray failure ratio

If a candidate is invalid:
- return a large penalty merit
- preserve a diagnostic message

## Diagnostics and logging
Every optimization run should save:
- optimizer settings
- initial prescription
- final prescription
- best merit value
- merit operand breakdown
- iteration / generation history

Preferred formats:
- JSON for machine-readable logs
- Markdown summary for human review

## UI integration phase
Do not start here.
Once the optimization core is stable, add a UI layer that can:
- choose variables
- choose operands and weights
- select optimizer backend
- run / stop optimization
- display live best merit and current prescription

This should build on the existing layout editor branch, not replace it.

## Risks
- Kraken analyses may require careful setup for stop surface, aperture type, wavelength, and field.
- Some merit terms are expensive; wavefront-based optimization will cost more than spot-only merit.
- Poor variable scaling will harm optimizer performance.
- Directly integrating `pagmo2` C++ first would slow progress; Python-side integration is lower risk.

## Recommended next tasks
1. Create `KrakenOS/Optimization/` with the core dataclasses and evaluator scaffold.
2. Implement spot-RMS merit on a doublet example.
3. Add a `pygmo2` adapter using the local `~/Projects/pygmo2` checkout.
4. Add wavefront merit once spot optimization is stable.
5. Only then connect optimization controls into the Kraken UI.

## Current implementation status
- `KrakenOS/Optimization/` scaffold exists.
- A doublet merit example runs end-to-end with `scipy.optimize.differential_evolution`.
- The first wavefront merit term is implemented via `PupilCalc()` and `Phase()`.
- A pygmo-compatible UDP wrapper exists in:
  - `KrakenOS/Optimization/adapters/pygmo2_adapter.py`

## Current blocker for local pygmo2
The local Python repo at `~/Projects/pygmo2` is not importable yet in the Kraken environment.
At the moment:
- `pygmo` is not installed in the Kraken venv,
- the local `~/Projects/pygmo2` checkout does not contain a built `pygmo.core` extension,
- building it will require a working local `pagmo` installation plus the required CMake toolchain.

The example script therefore does two things:
- always runs a real optimization with SciPy,
- attempts a `pygmo` run only if `pygmo` becomes importable later.

## Local pygmo bootstrap notes
The intended local dependency chain is:
- `~/Projects/pagmo2` -> build/install C++ pagmo locally
- `~/Projects/pygmo2` -> build/install Python bindings against that pagmo install
- Kraken optimization example -> import `pygmo` and evolve the UDP

Environment support already added to `devenv.nix`:
- `cmake`
- `ninja`
- `gcc`
- `boost`
- `tbb`
- pip packages:
  - `cloudpickle`
  - `pybind11`
