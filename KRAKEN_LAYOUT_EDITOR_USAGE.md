# Kraken Layout Editor Manual

## Launch

```bash
cd ~/Projects/Kraken-Optical-Simulator
devenv shell -- bash -lc 'python -m KrakenOS.UI.layout_editor'
```

Headless native snapshot:

```bash
cd ~/Projects/Kraken-Optical-Simulator
./.devenv/state/venv/bin/python -m KrakenOS.UI.render_layout_snapshot --mode native --output ~/Pictures/kraken_layout_headless.jpg
```

## Main UI layout

### Left side

- `Display`
  - object mode
  - wavelength
  - orientation
  - ray fan count
  - pupil factor
  - analysis stop surface
  - aperture type/value
- `Field`
  - field type
  - field value
  - field count
- editable prescription table
- plot area and analysis toolbar

### Right side

- `Information`
- `Optimization`
- `Progress`
- `Debug`

## Display panel

### Object mode

- `Infinity`
  - source is a collimated field definition
- `Finite`
  - source is a finite-distance object definition

`Field angle = 0` does not make `Finite` and `Infinity` equivalent. They remain different source models.

### Orientation

- `Vertical`
  - ordinary axial-style display
- `Horizontal`
  - folded/off-axis display orientation used for mirror-fold systems

### Analysis stop surface

- `Auto`
  - use the editor's default analysis stop
- or choose an explicit row by index and name

### Aperture

- `STOP`
- `EPD`

These settings affect analysis modes that depend on pupil construction.

## Field panel

### Field type

Available definitions:

- `Angle`
- `Object Height`
- `Paraxial Image Height`
- `Real Image Height`

Only one field definition is active at a time.

### Recommended usage

- for `Infinity`
  - use `Angle`
- for `Finite`
  - use `Object Height`

The status bar shows:

- preferred field note
- any warning such as field exceeding object radius
- converted field summary

## Prescription table

### Editing

- click a cell to select it
- double-click or type into editable numeric cells
- right-click `Surface` and `Glass` for popup choices

### Selection

- click: single selection
- `Ctrl` + click: toggle row selection
- `Shift` + click: contiguous row range
- arrow keys move the active cell

### Toolbar actions

- `Add surface`
- `Delete`
- `Duplicate`
- `Flip`
- `▲`
- `▼`
- `Common Optical Layout`
- `Examples`

### Common Optical Layout insertion

A common layout now inserts after the last selected row.

If nothing is selected, it inserts before the final `Image` row.

## Plot modes

Toolbar buttons:

- `Open 3D`
- `2D`
- `Native`
- `Spot`
- `PSF`
- `RMS`
- `FC/Dist`
- `Pupil`
- `Seidel`
- `Wavefront`
- `MTF`

### `2D`

Use this for the stable folded preview.

### `Native`

Use this for the folded-native debug/display path.

It shows:

- native-derived optical surfaces
- native-gated rays
- native surface diagnostics in `Information` and `Debug`

### Current folded limitation

For folded preview layouts, not every analysis mode is available yet.

That is expected. Folded-native analysis is still being built out.

## Optimization workflow

### Mark variables

Right-click a supported numeric cell.

Currently useful variable types include:

- `Radius`
- `Thickness`

Thickness optimization now supports:

- `Object`
- `Mirror`
- `Standard`
- `Thin Lens`
- `Grating`

`Image` remains excluded.

### Bounds

Right-click the same cell and choose:

- `Set bounds...`
- `Clear bounds`

### Operand setup

In the `Optimization` panel:

1. select one or more merit operands
2. configure per-operand fields such as:
   - `Weight`
   - `Target`
   - `Wvl`
   - `Field`
   - `Surf`
   - `Aper`
   - `AVal`
   - `Freq`
   - `Mode`
3. click `Start Optimization`

### Example: optimize a mirror distance

1. load `Double Mirror Fold`
2. right-click `Mirror 2` `Thickness`
3. select it for optimization
4. choose operand `Spot RMS`
5. click `Start Optimization`

### Example: optimize a singlet radius

1. load `Single Lens`
2. right-click front or back `Rc`
3. select it for optimization
4. set bounds
5. choose operand `Spot RMS` or `Wavefront RMS`
6. click `Start Optimization`

## Native folded workflow example

### Default folded system

The default startup layout is a folded mirror + singlet system.

Use it like this:

1. launch the editor
2. keep `Orientation = Horizontal`
3. use `Native`
4. inspect:
   - lens body
   - mirror overlays
   - rays
   - `Information`
   - `Debug`

### Insert a doublet after the singlet

1. click the row where insertion should happen
   - usually after the singlet back surface
2. choose `Doublet Lens` from `Common Optical Layout`
3. confirm the inserted rows appear immediately after the selected row
4. use `▲` / `▼` only for fine adjustment after insertion

## Files and persistence

### Save

- `File -> Save`
- `File -> Save As`

Saved Python layout files preserve optimization marks and bounds.

### Open

- `File -> Open`

## Diagnostics

### `Debug`

Use this for:

- native hit sequence information
- native overlay metrics
- fallback/error messages

### `Progress`

Use this for:

- optimization progress
- long analysis generation steps

### Headless snapshots

Use the headless renderer for reproducible image inspection when debugging visual issues.

## Known limitations

- folded-native display is ahead of folded-native analysis
- some complex example imports still need manual cleanup after insertion
- native folded view still uses a readable scaffold for placement; it is not a raw Kraken projection view
