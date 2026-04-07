# KrakenOS vs Optiland Gap Closure

## Current position

KrakenOS is now much stronger in these areas:

- editable UI for optical prescriptions
- optimization variable/operand registry
- native and folded display work
- off-axis/folded visualization
- optimization integration through `pygmo`
- headless snapshot workflow for debugging

Optiland is still ahead in overall product maturity.

## Main remaining gaps

### 1. Folded/off-axis analysis

Optiland's value is not just display. It is the tighter coupling between:

- model
- analysis
- optimization
- GUI

Kraken now has a readable folded native display, but folded analysis is still incomplete.

This is the most important remaining technical gap for Kraken's off-axis story.

### 2. Problem definition UX

Kraken has:

- variable registry
- operand registry

But the UI is still more manual than Optiland.

Needed next:

- per-operand typed widgets only for relevant parameters
- presets for common optimization tasks
- saved optimization problem definitions
- cleaner constraint vs target separation

### 3. Analysis breadth

Optiland remains ahead in integrated analysis breadth.

Kraken should prioritize:

1. folded-native spot / RMS
2. folded-native pupil / wavefront
3. folded-native MTF
4. better field curvature / distortion
5. report export

### 4. Solver backend resilience

Kraken now works with local `pygmo`, but it is still fragile compared with a platform where the optimizer is a first-class dependency.

Needed next:

- detect optimizer availability clearly at startup
- show solver backend in UI
- allow fallback optimizer backend when `pygmo` is not available

### 5. UI workflow polish

Needed next:

- better scroll/overflow behavior in configuration panels
- clearer folded/native distinction in the plot legend and status
- per-mode help text
- operand presets for common workflows

## Recommended roadmap

### Short term

- finish folded-native `Spot` / `RMS`
- keep improving native folded surfaces only when directly visible issues remain
- add optimizer backend status in UI
- keep strengthening the markdown/manual workflow

### Medium term

- folded-native `Wavefront`
- folded-native `MTF`
- problem save/load
- optimization presets
- report export

### Long term

- better native off-axis geometry extraction so the folded scaffold becomes less necessary
- more complete non-sequential workflow
- broader operand catalog with real engineering constraints

## What not to do next

Do not spend the next cycle on cosmetic surface tweaks if the folded-native analysis backend is still missing.

That would close the wrong gap.

The real gap to Optiland now is:

- integrated folded analysis and optimization

not another small improvement in line styling.
