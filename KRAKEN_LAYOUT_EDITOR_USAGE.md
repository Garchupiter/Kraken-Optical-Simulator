# Kraken Layout Editor Usage

## Launch

```bash
cd ~/Projects/Kraken-Optical-Simulator
devenv shell -- bash -lc 'python -m KrakenOS.UI.layout_editor'
```

## Main areas

- Top toolbar:
  - common optical layouts
  - examples
  - surface editing actions
  - analysis mode buttons
  - optimization start/stop
- Left panel `Display`:
  - object mode and ray controls
  - analysis surface and aperture controls
  - merit-function selection
  - toggle for `PP / EP / XP`
- Center table:
  - editable prescription
  - Excel-like cell focus border and light grid overlay
- Right panel `Information`:
  - imaging, principal-plane, pupil, and spot data
- Bottom panels:
  - `Debug`
  - `Progress`

## Layout and example loading

- Use `Common Optical Layout` to load a starter layout from `KrakenOS/common_optical_layouts`.
- Use `Examples` to load a supported file from `KrakenOS/Examples`.
- The selected item stays visible in the dropdown after loading.

## Table editing

- Double-click a normal cell to edit it.
- Right-click `Surface` or `Glass` cells to choose from popup menus.
- `Add surface` inserts before the final `Image` row.
- `Delete`, `Duplicate`, `Flip`, `Move Up`, and `Move Down` work on the current selection.
- `Flip` also reverses lens-block signs/names/media handling for standard surfaces.

## Source and ray controls

### Object mode

- `Finite`
  - object height comes from the first row `Object` diameter
  - object distance comes from the first row `Object` thickness
  - `Field count` = number of object points
  - `Ray fan count` = number of rays emitted by each object point
  - `Max field angle [deg]` = angular spread of each finite-field fan

- `Infinity`
  - `Field count` = number of field angles
  - `Ray fan count` = pupil sampling used for the preview bundle
  - `Max field angle [deg]` = maximum off-axis field angle

### Notes

- The first row `Object` thickness is the source-to-next-surface distance.
- The first row `Object` diameter sets the field extent in finite mode.
- Preview rays are for visualization; spot/RMS analysis uses a dedicated analysis ray set.

## Analysis modes

Toolbar buttons:

- `Layout`
- `Spot`
- `RMS`
- `Pupil`
- `Seidel`
- `Wavefront`

Current behavior:

- `Spot` and `RMS` use a dedicated `PupilCalc + TraceLoop` analysis bundle.
- `Pupil`, `Seidel`, and `Wavefront` use the selected analysis surface plus current aperture settings.
- `Show PP / EP / XP` overlays:
  - front principal plane
  - back principal plane
  - entrance pupil
  - exit pupil

## Information panel

Current sections:

- `Imaging`
  - EFFL
  - magnification
- `Principal Planes`
  - front principal plane
  - back principal plane
- `Pupils`
  - entrance pupil radius / diameter / z
  - exit pupil radius / diameter / z
  - Airy radius
- `Spot`
  - spot RMS
  - centroid X/Y

## Optimization workflow

### Mark variables

- Right-click an `Rc` or `Thickness` cell.
- Choose:
  - `Select to optimize`
  - `Unselect from optimize`
  - `Set bounds...`
  - `Clear bounds`
- Marked numeric cells show a trailing `*`.
- Marked rows are highlighted.

### Run optimization

- Choose a merit function in the left panel:
  - `Spot RMS`
  - `Wavefront RMS`
  - `Spot + Wavefront`
- Click `Start Optimization`.
- Click `Stop` to request cancellation after the current generation.

### Progress and debug

- `Debug` shows captured Kraken/runtime output.
- `Progress` shows:
  - generation progress
  - best-merit updates
  - final operand breakdown
  - spinner and percentage in the header

## Save and open

- `File -> Open` loads a layout `.py` file.
- `File -> Save` / `Save As` writes a Kraken-compatible Python layout script.
- Saved files preserve optimization metadata for marked `Rc` / `Thickness` cells.

## Current limitations

- The table is still `ttk.Treeview`, so the Excel behavior is approximate rather than native.
- Grid lines are drawn as overlays and depend on visible rows.
- The preview ray bundle is separate from the analysis ray bundle by design.
- Some complex examples from `KrakenOS/Examples` remain best-effort when loaded into the editor.
