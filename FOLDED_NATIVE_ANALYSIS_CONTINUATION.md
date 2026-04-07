# Folded Native Analysis Continuation

## Current state

There are now three distinct concepts in the editor:

- `2D`
  - readable folded preview
  - custom folded geometry and custom folded rays
  - not suitable as the source of truth for optical analysis
- `Native`
  - readable folded layout scaffold
  - Kraken-derived folded optical surfaces over that scaffold
  - native-hit-sequence-gated ray display
  - useful for debugging folded/off-axis geometry
- Kraken standard axial analysis backend
  - used by the existing `Spot`, `RMS`, `Wavefront`, `Pupil`, `MTF`, `FC/Dist`, `PSF`, `Seidel` paths
  - works for normal axial/sequential layouts
  - does not yet work as a proper folded native analysis backend

That is why the UI currently reports `Analysis unavailable for folded preview`.

## What is missing

The folded display is ahead of the folded analysis backend.

The existing analysis code assumes one of these:

- normal axial sequential plotting
- local pupil/image-plane analysis driven by the current axial system assumptions

A folded/off-axis system needs analysis data from Kraken's native traced system, not from the custom folded preview.

## Correct next implementation order

### 1. Folded-native spot first

Add a folded-native `Spot` / `RMS` path that:

- builds the real Kraken system
- traces a dedicated analysis bundle against that real system
- reads the true image-plane hit cloud from native trace results
- plots the spot on a normal analysis axes, not in the folded layout axes

This is the first milestone because it is the simplest reliable imaging metric for folded layouts.

### 2. Folded-native pupil

Once spot works, wire:

- entrance pupil
- exit pupil
- aperture-driven sampling

against the real folded system.

### 3. Folded-native wavefront

After pupil sampling is stable, add:

- wavefront RMS
- phase-based calculations

Only do this after the pupil path is confirmed stable for the folded system.

### 4. Folded-native MTF

Do not start here.

MTF should be built only after:

- folded-native spot is stable
- folded-native pupil is stable
- folded-native wavefront or PSF support is stable

Otherwise the result will be easy to plot and hard to trust.

## Recommended architecture

Do not route folded analysis through the custom folded display.

Instead:

- display backend:
  - `2D`
  - `Native`
- analysis backend:
  - always use the real Kraken system for folded layouts

In other words:

- display may be folded and human-readable
- analysis must be native and system-driven

## Concrete code direction

### A. Introduce a folded-analysis capability check

Add a helper such as:

- `_is_folded_native_layout()`

and use it to branch analysis behavior.

### B. Split analysis data generation from plotting

The current code mixes:

- system building
- ray tracing
- plotting

For folded work, factor each analysis mode into:

1. native data collection
2. mode-specific plotting

### C. Start with image-plane hit extraction

A reusable helper should return:

- image-plane hit coordinates
- wavelength
- surface index reached
- ray success/failure count

That helper becomes the base for:

- `Spot`
- `RMS`
- later `PSF`
- later geometric `MTF`

## Sanity checks for the next chat

The next chat should validate folded-native analysis in this order:

1. `Spot` on the default folded mirror + singlet layout
2. `RMS` numeric value matches the plotted spot cloud
3. changing `Mirror 2 Thickness` moves the folded-native spot in a physically sensible way
4. changing `Field type / value` updates folded-native spot accordingly

## Non-goals for the next chat

Do not try to solve all of these at once:

- full folded-native `Wavefront`
- full folded-native `MTF`
- all operand support for folded mode
- complete non-sequential stray-light analysis

The next chat should be disciplined and land folded-native `Spot` / `RMS` first.
