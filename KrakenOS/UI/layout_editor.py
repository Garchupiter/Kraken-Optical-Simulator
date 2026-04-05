"""Simple KrakenOS layout editor.

This is an initial editor scaffold that mirrors the RayTracing workflow:
- file-backed starter layouts
- editable surface table
- embedded axial sketch with a small traced ray fan
"""

from __future__ import annotations

import importlib.util
import io
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass, asdict
import os
from pathlib import Path
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

import KrakenOS as Kos
from KrakenOS.Display import Plot2DRays, Plot2DSurf
from KrakenOS.Optimization import (
    MeritEvaluator,
    MeritFunction,
    OpticalVariable,
    SpotRMSOperand,
    WavefrontRMSOperand,
)
from KrakenOS.Optimization.adapters.pygmo2_adapter import Pygmo2MeritProblem


LAYOUTS_DIR = Path(__file__).resolve().parent.parent / "common_optical_layouts"
EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "Examples"
FIELDS = (
    "label",
    "surface",
    "name",
    "rc",
    "thickness",
    "diameter",
    "glass",
)
COLUMN_LABELS = {
    "label": "#",
    "surface": "Surface",
    "name": "Name",
    "rc": "Rc [mm]",
    "thickness": "Thickness [mm]",
    "diameter": "Diameter [mm]",
    "glass": "Glass",
}
NUMERIC_FIELDS = {"rc", "thickness", "diameter"}
SURFACE_TYPES = ("Object", "Standard", "Thin Lens", "Grating", "Image")
MERIT_MODES = ("Spot RMS", "Wavefront RMS", "Spot + Wavefront")


class _CapturedExample(Exception):
    def __init__(self, surfaces):
        super().__init__("Captured example system")
        self.surfaces = surfaces


@dataclass
class SurfaceRow:
    label: str = "0"
    surface: str = "Standard"
    name: str = "Surface"
    optimize_rc: bool = False
    optimize_rc_bounds: tuple[float, float] | None = None
    rc: float = 0.0
    optimize_thickness: bool = False
    optimize_thickness_bounds: tuple[float, float] | None = None
    thickness: float = 0.0
    diameter: float = 25.0
    glass: str = "AIR"


def _load_python_data(path: Path) -> dict:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Could not load layout file: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    title = getattr(module, "TITLE", path.stem.replace("_", " ").title())
    surfaces = getattr(module, "SURFACES", None)
    if not isinstance(surfaces, list) or not surfaces:
        raise ValueError(f"{path.name} does not define a non-empty SURFACES list.")
    return {"title": title, "surfaces": surfaces}


def _load_python_title(path: Path) -> str:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Could not load layout file: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, "TITLE", path.stem.replace("_", " ").title())


def _coerce_opt_flag(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip() in {"*", "1", "true", "True", "yes", "on"}


def _coerce_bounds(value) -> tuple[float, float] | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return (float(value[0]), float(value[1]))
    return None


class KrakenLayoutEditor(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("KrakenOS Layout Editor")
        self.geometry("1400x850")
        self.minsize(1100, 720)

        self.current_layout_file: Path | None = None
        self.layout_files: dict[str, Path] = {}
        self.example_files: dict[str, Path] = {}
        self.rows: list[SurfaceRow] = []
        self.editor: tk.Widget | None = None
        self.popup_menu: tk.Menu | None = None
        self.current_menu_row_id: str | None = None
        self.current_menu_field: str | None = None
        self.analysis_mode = "none"
        self.last_system = None
        self.last_rays = None
        self.optimization_running = False
        self.optimization_cancel_requested = False
        self.optimization_context: dict | None = None

        self._build_menu()
        self._build_ui()
        self.load_layouts()
        self.load_examples()
        if self.layout_names:
            self.load_layout_by_name(self.layout_names[0])

    def _build_menu(self) -> None:
        menubar = tk.Menu(self)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open", command=self.open_layout)
        file_menu.add_command(label="Save", command=self.save_layout)
        file_menu.add_command(label="Save As", command=self.save_layout_as)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.destroy)
        menubar.add_cascade(label="File", menu=file_menu)
        self.config(menu=menubar)

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)

        toolbar = ttk.Frame(self, padding=8)
        toolbar.grid(row=0, column=0, sticky="ew")
        toolbar.columnconfigure(0, weight=0)
        toolbar.columnconfigure(1, weight=0)
        toolbar.columnconfigure(2, weight=0)
        toolbar.columnconfigure(3, weight=1)

        self.layout_var = tk.StringVar(value="Common Optical Layout")
        self.layout_menu = ttk.Combobox(
            toolbar,
            textvariable=self.layout_var,
            state="readonly",
            width=30,
        )
        self.layout_menu.grid(row=0, column=0, padx=(0, 8), pady=4, sticky="w")
        self.layout_menu.bind("<<ComboboxSelected>>", self._on_layout_selected)

        self.example_var = tk.StringVar(value="Examples")
        self.example_menu = ttk.Combobox(
            toolbar,
            textvariable=self.example_var,
            state="readonly",
            width=34,
        )
        self.example_menu.grid(row=0, column=1, padx=(0, 8), pady=4, sticky="w")
        self.example_menu.bind("<<ComboboxSelected>>", self._on_example_selected)

        ttk.Button(toolbar, text="Add surface", command=self.add_surface).grid(
            row=0, column=2, padx=4, pady=4
        )
        ttk.Button(toolbar, text="Delete", command=self.delete_selected).grid(
            row=0, column=3, padx=4, pady=4, sticky="w"
        )
        ttk.Button(toolbar, text="Duplicate", command=self.duplicate_selected).grid(
            row=0, column=4, padx=4, pady=4, sticky="w"
        )
        ttk.Button(toolbar, text="Flip", command=self.flip_selected).grid(
            row=0, column=5, padx=4, pady=4, sticky="w"
        )
        ttk.Button(toolbar, text="▲", width=3, command=self.move_up).grid(
            row=0, column=6, padx=(12, 2), pady=4
        )
        ttk.Button(toolbar, text="▼", width=3, command=self.move_down).grid(
            row=0, column=7, padx=2, pady=4
        )
        ttk.Button(toolbar, text="Refresh Plot", command=self.refresh_plot).grid(
            row=0, column=8, padx=(12, 0), pady=4
        )
        ttk.Button(toolbar, text="Layout", command=lambda: self.set_analysis_mode("none")).grid(
            row=0, column=9, padx=(12, 2), pady=4
        )
        ttk.Button(toolbar, text="Spot", command=lambda: self.set_analysis_mode("spot")).grid(
            row=0, column=10, padx=2, pady=4
        )
        ttk.Button(toolbar, text="RMS", command=lambda: self.set_analysis_mode("rms")).grid(
            row=0, column=11, padx=2, pady=4
        )
        ttk.Button(toolbar, text="Pupil", command=lambda: self.set_analysis_mode("pupil")).grid(
            row=0, column=12, padx=2, pady=4
        )
        ttk.Button(toolbar, text="Seidel", command=lambda: self.set_analysis_mode("seidel")).grid(
            row=0, column=13, padx=2, pady=4
        )
        ttk.Button(toolbar, text="Wavefront", command=lambda: self.set_analysis_mode("wavefront")).grid(
            row=0, column=14, padx=2, pady=4
        )
        ttk.Button(toolbar, text="Start Optimization", command=self.start_optimization).grid(
            row=0, column=15, padx=(12, 0), pady=4
        )
        ttk.Button(toolbar, text="Stop", command=self.stop_optimization).grid(
            row=0, column=16, padx=(4, 0), pady=4
        )

        main = ttk.Panedwindow(self, orient=tk.VERTICAL)
        main.grid(row=1, column=0, sticky="nsew")

        top = ttk.Panedwindow(main, orient=tk.HORIZONTAL)
        main.add(top, weight=3)

        controls = ttk.LabelFrame(top, text="Display", padding=8)
        controls.columnconfigure(0, weight=1)
        top.add(controls, weight=1)

        table_frame = ttk.Frame(top, padding=8)
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)
        top.add(table_frame, weight=3)

        results = ttk.LabelFrame(top, text="Information", padding=8)
        results.columnconfigure(0, weight=1)
        results.rowconfigure(0, weight=1)
        top.add(results, weight=2)

        plot_frame = ttk.Frame(main, padding=8)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)
        main.add(plot_frame, weight=2)

        bottom = ttk.Panedwindow(main, orient=tk.HORIZONTAL)
        main.add(bottom, weight=1)

        debug_frame = ttk.LabelFrame(bottom, text="Debug", padding=8)
        debug_frame.columnconfigure(0, weight=1)
        debug_frame.rowconfigure(0, weight=1)
        bottom.add(debug_frame, weight=1)

        progress_frame = ttk.LabelFrame(bottom, text="Progress", padding=8)
        progress_frame.columnconfigure(0, weight=1)
        progress_frame.rowconfigure(0, weight=1)
        bottom.add(progress_frame, weight=1)

        self.table = ttk.Treeview(table_frame, columns=FIELDS, show="headings", selectmode="extended")
        for field in FIELDS:
            self.table.heading(field, text=COLUMN_LABELS[field])
            width = (
                55 if field == "label"
                else 140 if field == "surface"
                else 160 if field == "name"
                else 110
            )
            self.table.column(field, width=width, stretch=True, anchor="center")
        self.table.grid(row=0, column=0, sticky="nsew")
        self.table.bind("<Double-1>", self.begin_edit)
        self.table.bind("<Button-3>", self.show_context_menu)
        self.table.tag_configure("optimize", background="#fff4bf")

        yscroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.table.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        self.table.configure(yscrollcommand=yscroll.set)

        self._build_controls_panel(controls)
        self._build_results_panel(results)

        self.figure = Figure(figsize=(7, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        self.debug_text = tk.Text(debug_frame, wrap="word", height=8)
        self.debug_text.grid(row=0, column=0, sticky="nsew")
        debug_scroll = ttk.Scrollbar(debug_frame, orient="vertical", command=self.debug_text.yview)
        debug_scroll.grid(row=0, column=1, sticky="ns")
        self.debug_text.configure(yscrollcommand=debug_scroll.set)

        self.progress_text = tk.Text(progress_frame, wrap="word", height=8)
        self.progress_text.grid(row=0, column=0, sticky="nsew")
        progress_scroll = ttk.Scrollbar(progress_frame, orient="vertical", command=self.progress_text.yview)
        progress_scroll.grid(row=0, column=1, sticky="ns")
        self.progress_text.configure(yscrollcommand=progress_scroll.set)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(self, textvariable=self.status_var, anchor="w").grid(
            row=2, column=0, sticky="ew", padx=8, pady=(0, 8)
        )

    def _build_controls_panel(self, parent) -> None:
        ttk.Label(parent, text="Wavelength [um]").grid(row=0, column=0, sticky="w", pady=(0, 2))
        self.wavelength_var = tk.StringVar(value="0.55")
        ttk.Entry(parent, textvariable=self.wavelength_var, width=12).grid(row=1, column=0, sticky="ew", pady=(0, 8))

        ttk.Label(parent, text="Ray count").grid(row=2, column=0, sticky="w", pady=(0, 2))
        self.ray_count_var = tk.StringVar(value="5")
        ttk.Entry(parent, textvariable=self.ray_count_var, width=12).grid(row=3, column=0, sticky="ew", pady=(0, 8))

        ttk.Label(parent, text="Ray height factor").grid(row=4, column=0, sticky="w", pady=(0, 2))
        self.ray_height_factor_var = tk.StringVar(value="0.8")
        ttk.Entry(parent, textvariable=self.ray_height_factor_var, width=12).grid(row=5, column=0, sticky="ew", pady=(0, 8))

        ttk.Label(parent, text="Analysis stop surface").grid(row=6, column=0, sticky="w", pady=(0, 2))
        self.analysis_surface_var = tk.StringVar(value="Auto")
        self.analysis_surface_menu = ttk.Combobox(
            parent,
            textvariable=self.analysis_surface_var,
            state="readonly",
            values=["Auto"],
        )
        self.analysis_surface_menu.grid(row=7, column=0, sticky="ew", pady=(0, 8))
        self.analysis_surface_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())

        ttk.Label(parent, text="Aperture type").grid(row=8, column=0, sticky="w", pady=(0, 2))
        self.aperture_type_var = tk.StringVar(value="STOP")
        self.aperture_type_menu = ttk.Combobox(
            parent,
            textvariable=self.aperture_type_var,
            state="readonly",
            values=["STOP", "EPD"],
        )
        self.aperture_type_menu.grid(row=9, column=0, sticky="ew", pady=(0, 8))
        self.aperture_type_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())

        ttk.Label(parent, text="Aperture value").grid(row=10, column=0, sticky="w", pady=(0, 2))
        self.aperture_value_var = tk.StringVar(value="1.0")
        ttk.Entry(parent, textvariable=self.aperture_value_var, width=12).grid(row=11, column=0, sticky="ew", pady=(0, 8))

        ttk.Label(parent, text="Merit function").grid(row=12, column=0, sticky="w", pady=(0, 2))
        self.merit_mode_var = tk.StringVar(value=MERIT_MODES[0])
        self.merit_mode_menu = ttk.Combobox(
            parent,
            textvariable=self.merit_mode_var,
            state="readonly",
            values=MERIT_MODES,
        )
        self.merit_mode_menu.grid(row=13, column=0, sticky="ew", pady=(0, 8))

        ttk.Label(
            parent,
            text="Right-click Rc/Thickness cells to mark optimization variables.",
            wraplength=180,
            justify="left",
        ).grid(row=14, column=0, sticky="ew", pady=(4, 8))
        ttk.Button(parent, text="Clear Marks", command=self.clear_optimization_marks).grid(
            row=15, column=0, sticky="ew", pady=(0, 4)
        )
        ttk.Button(parent, text="Refresh", command=self.refresh_plot).grid(row=16, column=0, sticky="ew", pady=(8, 0))

    def _build_results_panel(self, parent) -> None:
        self.results_table = ttk.Treeview(parent, columns=("property", "value"), show="headings", selectmode="none")
        self.results_table.heading("property", text="Property")
        self.results_table.heading("value", text="Value")
        self.results_table.column("property", width=180, anchor="w")
        self.results_table.column("value", width=140, anchor="w")
        self.results_table.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(parent, orient="vertical", command=self.results_table.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.results_table.configure(yscrollcommand=scroll.set)

    def load_layouts(self) -> None:
        self.layout_files = {}
        for path in sorted(LAYOUTS_DIR.glob("*.py")):
            if path.name.startswith("_") or path.name == "__init__.py":
                continue
            try:
                title = _load_python_title(path)
            except Exception:
                continue
            self.layout_files[title] = path
        self.layout_names = sorted(self.layout_files)
        self.layout_menu["values"] = ["Common Optical Layout", *self.layout_names]
        self.layout_var.set("Common Optical Layout")

    def load_examples(self) -> None:
        self.example_files = {}
        for path in sorted(EXAMPLES_DIR.glob("*.py")):
            self.example_files[path.stem] = path
        self.example_names = sorted(self.example_files)
        self.example_menu["values"] = ["Examples", *self.example_names]
        self.example_var.set("Examples")

    def set_analysis_mode(self, mode: str) -> None:
        self.analysis_mode = mode
        self.refresh_plot()

    def load_layout_by_name(self, name: str) -> None:
        path = self.layout_files.get(name)
        if path is None:
            return
        self.current_layout_file = path
        try:
            info = _load_python_data(path)
            self.rows = [
                SurfaceRow(
                    surface=str(item.get("surface", self._infer_surface_type(item))),
                    name=str(item.get("name", "Surface")),
                    optimize_rc=_coerce_opt_flag(item.get("optimize_rc", item.get("opt_rc", ""))),
                    optimize_rc_bounds=_coerce_bounds(item.get("optimize_rc_bounds")),
                    rc=float(item.get("rc", 0.0)),
                    optimize_thickness=_coerce_opt_flag(item.get("optimize_thickness", item.get("opt_thickness", ""))),
                    optimize_thickness_bounds=_coerce_bounds(item.get("optimize_thickness_bounds")),
                    thickness=float(item.get("thickness", 0.0)),
                    diameter=float(item.get("diameter", 25.0)),
                    glass=str(item.get("glass", "AIR")),
                )
                for item in info["surfaces"]
            ]
        except Exception:
            surfaces = self._extract_surfaces_from_example(path)
            self.rows = [self._row_from_surface(surface, index, len(surfaces)) for index, surface in enumerate(surfaces)]
        self._normalize_special_rows()
        self._sync_table()
        self.refresh_plot()
        self.layout_var.set("Common Optical Layout")
        self.status_var.set(f"Loaded {name}")

    def load_example_by_name(self, name: str) -> None:
        path = self.example_files.get(name)
        if path is None:
            return
        try:
            surfaces = self._extract_surfaces_from_example(path)
        except Exception as exc:
            self.example_var.set("Examples")
            self.status_var.set(f"Failed to load example {name}: {exc}")
            return
        self.current_layout_file = None
        self.rows = [self._row_from_surface(surface, index, len(surfaces)) for index, surface in enumerate(surfaces)]
        self._normalize_special_rows()
        self._apply_example_display_defaults(path)
        self._sync_table()
        self.refresh_plot()
        self.example_var.set("Examples")
        self.status_var.set(f"Loaded example {name}")

    def _sync_table(self) -> None:
        self.table.delete(*self.table.get_children())
        for index, row in enumerate(self.rows):
            row.label = str(index)
            values = [
                row.label,
                row.surface,
                row.name,
                self._format_numeric_cell("rc", row),
                self._format_numeric_cell("thickness", row),
                f"{row.diameter:g}",
                row.glass,
            ]
            tags = ("optimize",) if self._row_has_optimization(row) else ()
            self.table.insert("", "end", values=values, tags=tags)
        self._refresh_analysis_surface_choices()

    def _refresh_analysis_surface_choices(self) -> None:
        options = ["Auto"]
        for index, row in enumerate(self.rows):
            options.append(f"{index}: {row.name}")
        current = self.analysis_surface_var.get()
        self.analysis_surface_menu["values"] = options
        if current not in options:
            self.analysis_surface_var.set("Auto")

    @staticmethod
    def _parse_numeric_display(value: str) -> float:
        return float(value.replace("*", "").strip())

    @staticmethod
    def _format_numeric_cell(field: str, row: SurfaceRow) -> str:
        value = row.rc if field == "rc" else row.thickness
        mark = row.optimize_rc if field == "rc" else row.optimize_thickness
        text = f"{value:g}"
        if mark:
            text += " *"
        return text

    def _read_rows_from_table(self) -> None:
        rows: list[SurfaceRow] = []
        for item in self.table.get_children():
            values = self.table.item(item, "values")
            rows.append(
                SurfaceRow(
                    label=str(values[0]),
                    surface=str(values[1]),
                    name=str(values[2]),
                    optimize_rc=self.rows[len(rows)].optimize_rc if len(rows) < len(self.rows) else False,
                    optimize_rc_bounds=self.rows[len(rows)].optimize_rc_bounds if len(rows) < len(self.rows) else None,
                    rc=self._parse_numeric_display(str(values[3])),
                    optimize_thickness=self.rows[len(rows)].optimize_thickness if len(rows) < len(self.rows) else False,
                    optimize_thickness_bounds=self.rows[len(rows)].optimize_thickness_bounds if len(rows) < len(self.rows) else None,
                    thickness=self._parse_numeric_display(str(values[4])),
                    diameter=float(values[5]),
                    glass=str(values[6]),
                )
            )
        self.rows = rows

    def add_surface(self) -> None:
        insert_at = len(self.rows)
        if self.rows and self.rows[-1].surface == "Image":
            insert_at -= 1
        self.rows.insert(insert_at, SurfaceRow())
        self._sync_table()
        self.refresh_plot()

    def delete_selected(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        indices = sorted(self.table.index(item) for item in selected)
        for index in reversed(indices):
            del self.rows[index]
        self._sync_table()
        self.refresh_plot()

    def duplicate_selected(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        indices = sorted(self.table.index(item) for item in selected)
        insert_at = indices[-1] + 1
        duplicates = [SurfaceRow(**asdict(self.rows[index])) for index in indices]
        for offset, row in enumerate(duplicates):
            self.rows.insert(insert_at + offset, row)
        self._normalize_special_rows()
        self._sync_table()
        new_items = self.table.get_children()[insert_at:insert_at + len(duplicates)]
        self.table.selection_set(new_items)
        self.refresh_plot()

    def flip_selected(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        indices = sorted(self.table.index(item) for item in selected)
        if len(indices) < 2:
            return
        selected_rows = [SurfaceRow(**asdict(self.rows[index])) for index in indices]
        selected_thicknesses = [row.thickness for row in selected_rows]
        selected_glasses = [row.glass for row in selected_rows]
        flipped_rows = list(reversed(selected_rows))

        for row in flipped_rows:
            if row.surface == "Standard" and row.rc != 0.0:
                row.rc = -row.rc
            row.name = self._flipped_name(row.name)

        if len(flipped_rows) >= 2:
            remapped_thicknesses = list(reversed(selected_thicknesses[:-1])) + [selected_thicknesses[-1]]
            remapped_glasses = list(reversed(selected_glasses[:-1])) + [selected_glasses[-1]]
        else:
            remapped_thicknesses = selected_thicknesses
            remapped_glasses = selected_glasses

        for row, thickness, glass in zip(flipped_rows, remapped_thicknesses, remapped_glasses):
            row.thickness = thickness
            row.glass = glass

        for index, row in zip(indices, flipped_rows):
            self.rows[index] = row
        self._normalize_special_rows()
        self._sync_table()
        items = self.table.get_children()
        self.table.selection_set([items[index] for index in indices])
        self.refresh_plot()

    def move_up(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        index = min(self.table.index(item) for item in selected)
        if index == 0:
            return
        self.rows[index - 1], self.rows[index] = self.rows[index], self.rows[index - 1]
        self._sync_table()
        self.table.selection_set(self.table.get_children()[index - 1])
        self.refresh_plot()

    def move_down(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        index = max(self.table.index(item) for item in selected)
        if index >= len(self.rows) - 1:
            return
        self.rows[index + 1], self.rows[index] = self.rows[index], self.rows[index + 1]
        self._sync_table()
        self.table.selection_set(self.table.get_children()[index + 1])
        self.refresh_plot()

    def begin_edit(self, event: tk.Event) -> None:
        row_id = self.table.identify_row(event.y)
        column_id = self.table.identify_column(event.x)
        if not row_id or not column_id:
            return
        column_index = int(column_id.replace("#", "")) - 1
        field = FIELDS[column_index]
        if field == "label":
            return
        x, y, width, height = self.table.bbox(row_id, column_id)
        current_value = self.table.set(row_id, field)
        if field in {"rc", "thickness"}:
            current_value = current_value.replace("*", "").strip()

        if self.editor is not None:
            self.editor.destroy()

        if field == "surface":
            self._show_choice_menu(row_id, field, SURFACE_TYPES, event.x_root, event.y_root)
            return
        elif field == "glass":
            self._show_choice_menu(
                row_id,
                field,
                ("AIR", "BK7", "F2", "MIRROR"),
                event.x_root,
                event.y_root,
            )
            return
        else:
            editor = ttk.Entry(self.table)
            editor.insert(0, current_value)
            editor.bind("<Return>", lambda e: self._finish_edit(row_id, field))

        editor.place(x=x, y=y, width=width, height=height)
        editor.focus_set()
        editor.bind("<FocusOut>", lambda e: self._finish_edit(row_id, field, quiet=True))
        self.editor = editor

    def show_context_menu(self, event: tk.Event) -> None:
        row_id = self.table.identify_row(event.y)
        column_id = self.table.identify_column(event.x)
        if not row_id or not column_id:
            return
        column_index = int(column_id.replace("#", "")) - 1
        field = FIELDS[column_index]
        if field not in {"rc", "thickness"}:
            return
        row_index = self.table.index(row_id)
        row = self.rows[row_index]
        if row.surface in {"Object", "Image"}:
            return
        if self.popup_menu is not None:
            self.popup_menu.destroy()
        self.current_menu_row_id = row_id
        self.current_menu_field = field
        marked = row.optimize_rc if field == "rc" else row.optimize_thickness
        bounds = row.optimize_rc_bounds if field == "rc" else row.optimize_thickness_bounds
        menu = tk.Menu(self, tearoff=0)
        menu.add_command(
            label="Unselect from optimize" if marked else "Select to optimize",
            command=self.toggle_current_optimization_cell,
        )
        menu.add_separator()
        menu.add_command(label="Set bounds...", command=self.edit_current_bounds)
        menu.add_command(label="Clear bounds", command=self.clear_current_bounds, state=("normal" if bounds else "disabled"))
        self.popup_menu = menu
        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()

    def _finish_edit(self, row_id: str, field: str, quiet: bool = False) -> None:
        if self.editor is None:
            return
        value = self.editor.get().strip()
        self.editor.destroy()
        self.editor = None
        if not value:
            return
        if field in NUMERIC_FIELDS:
            try:
                float(value)
            except ValueError:
                if not quiet:
                    messagebox.showerror("Invalid value", f"{COLUMN_LABELS[field]} expects a number.")
                return
        self.table.set(row_id, field, value)
        self._read_rows_from_table()
        self._normalize_special_rows()
        self._sync_table()
        self.refresh_plot()

    def _cancel_edit(self) -> None:
        if self.editor is None:
            return
        self.editor.destroy()
        self.editor = None

    def _show_choice_menu(
        self,
        row_id: str,
        field: str,
        values: tuple[str, ...],
        x_root: int,
        y_root: int,
    ) -> None:
        if self.popup_menu is not None:
            self.popup_menu.destroy()
        menu = tk.Menu(self, tearoff=0)
        for value in values:
            menu.add_command(
                label=value,
                command=lambda selected=value: self._apply_choice(row_id, field, selected),
            )
        self.popup_menu = menu
        try:
            menu.tk_popup(x_root, y_root)
        finally:
            menu.grab_release()

    def _apply_choice(self, row_id: str, field: str, value: str) -> None:
        self.table.set(row_id, field, value)
        self._read_rows_from_table()
        self._normalize_special_rows()
        self._sync_table()
        self.refresh_plot()
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None

    def _on_layout_selected(self, _event: tk.Event) -> None:
        selected = self.layout_var.get().strip()
        if selected == "Common Optical Layout":
            return
        self.load_layout_by_name(selected)

    def _on_example_selected(self, _event: tk.Event) -> None:
        selected = self.example_var.get().strip()
        if selected == "Examples":
            return
        self.load_example_by_name(selected)

    @staticmethod
    def _row_has_optimization(row: SurfaceRow) -> bool:
        return row.optimize_rc or row.optimize_thickness

    def toggle_current_optimization_cell(self) -> None:
        if self.current_menu_row_id is None or self.current_menu_field is None:
            return
        index = self.table.index(self.current_menu_row_id)
        row = self.rows[index]
        if self.current_menu_field == "rc":
            row.optimize_rc = not row.optimize_rc
        elif self.current_menu_field == "thickness":
            row.optimize_thickness = not row.optimize_thickness
        self._sync_table()
        self.refresh_plot()
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None
        self.current_menu_row_id = None
        self.current_menu_field = None

    def edit_current_bounds(self) -> None:
        if self.current_menu_row_id is None or self.current_menu_field is None:
            return
        index = self.table.index(self.current_menu_row_id)
        row = self.rows[index]
        current = row.optimize_rc_bounds if self.current_menu_field == "rc" else row.optimize_thickness_bounds
        current_value = row.rc if self.current_menu_field == "rc" else row.thickness
        default_bounds = current or self._optimization_bounds(self.current_menu_field.capitalize() if self.current_menu_field == "rc" else "Thickness", current_value)

        dialog = tk.Toplevel(self)
        dialog.title(f"Bounds for {row.name} {self.current_menu_field}")
        dialog.transient(self)
        dialog.grab_set()
        dialog.resizable(False, False)

        ttk.Label(dialog, text="Lower").grid(row=0, column=0, padx=12, pady=(12, 4), sticky="w")
        lower_var = tk.StringVar(value=f"{default_bounds[0]:g}")
        ttk.Entry(dialog, textvariable=lower_var, width=16).grid(row=1, column=0, padx=12, pady=(0, 8))

        ttk.Label(dialog, text="Upper").grid(row=2, column=0, padx=12, pady=(0, 4), sticky="w")
        upper_var = tk.StringVar(value=f"{default_bounds[1]:g}")
        ttk.Entry(dialog, textvariable=upper_var, width=16).grid(row=3, column=0, padx=12, pady=(0, 12))

        def accept():
            try:
                lower = float(lower_var.get())
                upper = float(upper_var.get())
            except ValueError:
                self.append_debug("Invalid optimization bounds entry.")
                return
            if lower >= upper:
                self.append_debug("Optimization bounds rejected: lower must be less than upper.")
                return
            if self.current_menu_field == "rc":
                row.optimize_rc_bounds = (lower, upper)
            else:
                row.optimize_thickness_bounds = (lower, upper)
            self.append_progress(
                f"Bounds set for row {index} {self.current_menu_field}: [{lower:g}, {upper:g}]"
            )
            dialog.destroy()

        buttons = ttk.Frame(dialog)
        buttons.grid(row=4, column=0, padx=12, pady=(0, 12), sticky="w")
        ttk.Button(buttons, text="Save", command=accept).pack(side="left")
        ttk.Button(buttons, text="Cancel", command=dialog.destroy).pack(side="left", padx=(8, 0))

        self.wait_window(dialog)
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None
        self.current_menu_row_id = None
        self.current_menu_field = None

    def clear_current_bounds(self) -> None:
        if self.current_menu_row_id is None or self.current_menu_field is None:
            return
        index = self.table.index(self.current_menu_row_id)
        row = self.rows[index]
        if self.current_menu_field == "rc":
            row.optimize_rc_bounds = None
        else:
            row.optimize_thickness_bounds = None
        self.append_progress(f"Bounds cleared for row {index} {self.current_menu_field}.")
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None
        self.current_menu_row_id = None
        self.current_menu_field = None

    def clear_optimization_marks(self) -> None:
        for row in self.rows:
            row.optimize_rc = False
            row.optimize_thickness = False
        self._sync_table()

    def build_system(self):
        surfaces = []
        wavelength = self._current_wavelength()
        for row in self.rows:
            surface = Kos.surf()
            display_name = row.name if row.surface in {"Object", "Image"} else ""
            surface.Name = display_name.replace(" ", "\n")
            surface.Rc = row.rc
            surface.Thickness = row.thickness
            surface.Diameter = row.diameter
            surface.Glass = row.glass
            surface.Nm_Pos = self._name_offset(row)
            if row.surface == "Thin Lens":
                surface.Thin_Lens = row.rc if row.rc != 0 else 100.0
                surface.Rc = 0.0
            elif row.surface == "Grating":
                surface.Diff_Ord = 1.0
                surface.Grating_D = 1.0
            surfaces.append(surface)
        return Kos.system(surfaces, Kos.Setup(), build=1)

    def _current_wavelength(self) -> float:
        try:
            return float(self.wavelength_var.get())
        except ValueError:
            return 0.55

    def _current_aperture_type(self) -> str:
        value = self.aperture_type_var.get().strip().upper()
        if value in {"STOP", "EPD"}:
            return value
        return "STOP"

    def _current_aperture_value(self) -> float:
        try:
            value = float(self.aperture_value_var.get())
        except ValueError:
            return 1.0
        if value == 0.0:
            return 1.0
        return value

    @staticmethod
    def _name_offset(row: SurfaceRow) -> tuple[float, float]:
        if row.surface not in {"Object", "Image"}:
            return (0.0, 0.0)
        name = row.name.lower()
        base_y = max(row.diameter * 0.08, 0.0)
        if "front" in name:
            return (-max(row.diameter * 0.35, 8.0), base_y)
        if "back" in name:
            return (max(row.diameter * 0.15, 2.0), base_y)
        return (0.0, base_y)

    def refresh_plot(self) -> None:
        if not self.rows:
            self.ax.clear()
            self.ax.set_title("Axial Layout")
            self.canvas.draw_idle()
            return

        max_radius = 1.0
        for row in self.rows:
            radius = max(row.diameter / 2.0, 0.5)
            max_radius = max(max_radius, radius)

        self.figure.clear()
        if self.analysis_mode == "none":
            self.ax = self.figure.add_subplot(111)
            analysis_ax = None
        else:
            gs = self.figure.add_gridspec(1, 2, width_ratios=[3, 2], wspace=0.28)
            self.ax = self.figure.add_subplot(gs[0])
            analysis_ax = self.figure.add_subplot(gs[1])

        try:
            wavelength = self._current_wavelength()
            capture = io.StringIO()
            with redirect_stdout(capture), redirect_stderr(capture):
                system = self.build_system()
                if getattr(system.Pr3D, "ExistSolid", 0) == 0:
                    original_build = system.BUILD
                    system.BUILD = 1
                    system.build()
                    system.BUILD = original_build
                rays = Kos.raykeeper(system)
                for y0 in self._sample_ray_heights(max_radius):
                    system.Trace([0.0, y0, 0.0], [0.0, 0.0, 1.0], wavelength)
                    rays.push()
            self.append_debug(capture.getvalue())
            self.last_system = system
            self.last_rays = rays
            Plot2DSurf(system, 0, self.ax)
            surf_line_count = len(self.ax.lines)
            Plot2DRays(rays, 0, 0, self.ax, 0)
            self._style_embedded_plot(surf_line_count)
            self._plot_analysis(analysis_ax, system, rays, wavelength)
            self._update_results(system, rays, wavelength)
            self.status_var.set("Plot refreshed")
        except Exception as exc:
            self.last_system = None
            self.last_rays = None
            self._plot_fallback_preview(max_radius)
            if analysis_ax is not None:
                analysis_ax.clear()
                analysis_ax.text(0.5, 0.5, "Analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
            self._set_results([("Status", "Unavailable"), ("Error", str(exc))])
            self.status_var.set(f"Plot refreshed with fallback preview: {exc}")
            self.append_debug(f"Plot refresh error: {exc}")

        self.ax.grid(True, alpha=0.2)
        self.ax.set_xlabel("Z [mm]")
        self.ax.set_ylabel("Y [mm]")
        self.ax.set_title("KrakenOS Layout")
        self.figure.subplots_adjust(left=0.07, right=0.98, bottom=0.11, top=0.92, wspace=0.28)
        self.canvas.draw_idle()

    def _plot_analysis(self, analysis_ax, system, rays, wavelength: float) -> None:
        if analysis_ax is None:
            return
        analysis_ax.clear()
        try:
            X, Y, Z, L, M, N = rays.pick(-1)
        except Exception:
            X = Y = Z = L = M = N = np.asarray([])

        if X.size == 0:
            analysis_ax.text(0.5, 0.5, "No ray data", ha="center", va="center")
            analysis_ax.set_axis_off()
            return

        if self.analysis_mode == "spot":
            analysis_ax.scatter(X, Y, s=18, c="#c0392b", alpha=0.8)
            analysis_ax.axhline(0.0, color="#2c3e50", linewidth=0.6, alpha=0.5)
            analysis_ax.axvline(0.0, color="#2c3e50", linewidth=0.6, alpha=0.5)
            analysis_ax.set_title("Spot Diagram")
            analysis_ax.set_xlabel("X [mm]")
            analysis_ax.set_ylabel("Y [mm]")
            analysis_ax.set_aspect("equal", adjustable="box")
            analysis_ax.grid(True, alpha=0.2)
            return

        if self.analysis_mode == "rms":
            rms, cenX, cenY = Kos.RMS(X, Y, Z, L, M, N)
            radii = np.sqrt((X - cenX) ** 2 + (Y - cenY) ** 2)
            bins = min(max(5, int(np.sqrt(max(len(radii), 1)))), 20)
            analysis_ax.hist(radii, bins=bins, color="#4f81bd", edgecolor="white")
            analysis_ax.set_title(f"Spot Radius Histogram  |  RMS = {float(rms):.4g} mm")
            analysis_ax.set_xlabel("Radius [mm]")
            analysis_ax.set_ylabel("Count")
            analysis_ax.grid(True, axis="y", alpha=0.2)
            return

        if self.analysis_mode == "pupil":
            try:
                pupil = Kos.PupilCalc(
                    system,
                    self._analysis_surface_index(),
                    wavelength,
                    self._current_aperture_type(),
                    self._current_aperture_value(),
                )
                labels = [
                    ("Input Radius", float(pupil.RadPupInp)),
                    ("Input Z", float(pupil.PosPupInp[2])),
                    ("Output Radius", float(pupil.RadPupOut)),
                    ("Output Z", float(pupil.PosPupOut[2])),
                    ("Airy Radius", float(pupil.FocusAiryRadius)),
                ]
                y_pos = np.arange(len(labels))
                values = [item[1] for item in labels]
                analysis_ax.barh(y_pos, values, color="#16a085")
                analysis_ax.set_yticks(y_pos, [item[0] for item in labels])
                analysis_ax.set_title("Pupil Summary")
                analysis_ax.set_xlabel("mm")
                analysis_ax.grid(True, axis="x", alpha=0.2)
            except Exception:
                analysis_ax.text(0.5, 0.5, "Pupil analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
            return

        if self.analysis_mode == "seidel":
            try:
                pupil = Kos.PupilCalc(
                    system,
                    self._analysis_surface_index(),
                    wavelength,
                    self._current_aperture_type(),
                    self._current_aperture_value(),
                )
                seidel = Kos.Seidel(pupil)
                values = np.asarray(seidel.SCW_TOTAL, dtype=float)
                labels = seidel.SCW_NM
                analysis_ax.bar(np.arange(len(values)), values, color="#8e44ad")
                analysis_ax.set_xticks(np.arange(len(values)), labels, rotation=25, ha="right")
                analysis_ax.set_title("Seidel Coefficients in Waves")
                analysis_ax.set_ylabel("Waves")
                analysis_ax.grid(True, axis="y", alpha=0.2)
            except Exception:
                analysis_ax.text(0.5, 0.5, "Seidel analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
            return

        if self.analysis_mode == "wavefront":
            try:
                pupil = Kos.PupilCalc(
                    system,
                    self._analysis_surface_index(),
                    wavelength,
                    self._current_aperture_type(),
                    self._current_aperture_value(),
                )
                pupil.Samp = max(6, min(16, int(np.sqrt(max(1, self._current_ray_count())) * 3)))
                px, py, phase, p2v = Kos.Phase(pupil)
                scatter = analysis_ax.scatter(py, px, c=phase, cmap="RdBu_r", s=20)
                analysis_ax.set_title(f"Wavefront  |  P2V = {float(p2v):.4g}")
                analysis_ax.set_xlabel("X pupil")
                analysis_ax.set_ylabel("Y pupil")
                analysis_ax.set_aspect("equal", adjustable="box")
                analysis_ax.grid(True, alpha=0.2)
                self.figure.colorbar(scatter, ax=analysis_ax, fraction=0.046, pad=0.04, label="Waves")
            except Exception:
                analysis_ax.text(0.5, 0.5, "Wavefront analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
            return

    def _analysis_surface_index(self) -> int:
        selected = self.analysis_surface_var.get().strip()
        if selected and selected != "Auto":
            try:
                return int(selected.split(":", 1)[0])
            except ValueError:
                pass
        if len(self.rows) <= 2:
            return max(0, len(self.rows) - 1)
        candidate_indices = [i for i, row in enumerate(self.rows[1:-1], start=1)]
        if not candidate_indices:
            return 1
        return min(candidate_indices, key=lambda i: max(self.rows[i].diameter, 1e-9))

    def _style_embedded_plot(self, surf_line_count: int) -> None:
        neon_green = "#39FF14"
        for index, line in enumerate(self.ax.lines):
            if index < surf_line_count:
                line.set_linewidth(max(line.get_linewidth(), 1.25))
            else:
                line.set_linewidth(max(line.get_linewidth(), 1.8))
                line.set_color(neon_green)
                line.set_alpha(0.95)

    def _apply_example_display_defaults(self, path: Path) -> None:
        code = path.read_text(encoding="utf-8", errors="ignore")

        wavelength_match = re.search(r"\bW\s*=\s*([0-9]*\.?[0-9]+)", code)
        if wavelength_match:
            self.wavelength_var.set(wavelength_match.group(1))

        aperture_type_match = re.search(r"\b(?:AperType|ApType)\s*=\s*['\"](STOP|EPD)['\"]", code)
        if aperture_type_match:
            self.aperture_type_var.set(aperture_type_match.group(1))

        aperture_value_match = re.search(r"\b(?:AperVal|ApVal)\s*=\s*([0-9]*\.?[0-9]+)", code)
        if aperture_value_match:
            self.aperture_value_var.set(aperture_value_match.group(1))

        surf_match = re.search(r"\b(?:Surf|sup)\s*=\s*([0-9]+)", code)
        if surf_match:
            surf_index = surf_match.group(1)
            label = None
            for option in self.analysis_surface_menu["values"]:
                if option.startswith(f"{surf_index}:"):
                    label = option
                    break
            if label is not None:
                self.analysis_surface_var.set(label)
            else:
                self.analysis_surface_var.set("Auto")
        else:
            self.analysis_surface_var.set("Auto")

    def _set_results(self, items) -> None:
        self.results_table.delete(*self.results_table.get_children())
        for key, value in items:
            self.results_table.insert("", "end", values=(key, value))

    def append_debug(self, message: str) -> None:
        if not message:
            return
        self.debug_text.insert("end", message.rstrip() + "\n")
        self.debug_text.see("end")
        self.update_idletasks()

    def append_progress(self, message: str) -> None:
        if not message:
            return
        self.progress_text.insert("end", message.rstrip() + "\n")
        self.progress_text.see("end")
        self.update_idletasks()

    def _update_results(self, system, rays, wavelength: float) -> None:
        items = []
        items.append(("Surface count", str(len(self.rows))))
        items.append(("Optimized vars", str(len(self._build_optimization_variables()))))
        items.append(("Wavelength [um]", f"{wavelength:.4g}"))
        items.append(("Analysis surface", str(self._analysis_surface_index())))
        items.append(("Aperture type", self._current_aperture_type()))
        items.append(("Aperture value", f"{self._current_aperture_value():.4g}"))

        total_length = sum(max(float(row.thickness), 0.0) for row in self.rows)
        items.append(("Total length [mm]", f"{total_length:.4g}"))

        try:
            _, _, _, _, _, _, _, effl, ppa, ppp, _, _, _ = system.Parax(wavelength)
            items.append(("EFFL [mm]", f"{float(effl):.4g}"))
            items.append(("PPA [mm]", f"{float(ppa):.4g}"))
            items.append(("PPP [mm]", f"{float(ppp):.4g}"))
        except Exception:
            items.append(("Paraxial data", "Unavailable"))

        try:
            X, Y, Z, L, M, N = rays.pick(-1)
            if X.size:
                rms, cenX, cenY = Kos.RMS(X, Y, Z, L, M, N)
                items.append(("Spot RMS [mm]", f"{float(rms):.4g}"))
                items.append(("Spot centroid X [mm]", f"{float(cenX):.4g}"))
                items.append(("Spot centroid Y [mm]", f"{float(cenY):.4g}"))
            else:
                items.append(("Spot RMS [mm]", "No rays"))
        except Exception:
            items.append(("Spot RMS [mm]", "Unavailable"))

        try:
            pupil = Kos.PupilCalc(
                system,
                self._analysis_surface_index(),
                wavelength,
                self._current_aperture_type(),
                self._current_aperture_value(),
            )
            items.append(("EP radius [mm]", f"{float(pupil.RadPupInp):.4g}"))
            items.append(("EP z [mm]", f"{float(pupil.PosPupInp[2]):.4g}"))
            items.append(("XP radius [mm]", f"{float(pupil.RadPupOut):.4g}"))
            items.append(("XP z [mm]", f"{float(pupil.PosPupOut[2]):.4g}"))
            items.append(("Airy radius [mm]", f"{float(pupil.FocusAiryRadius):.4g}"))
        except Exception:
            items.append(("Pupil data", "Unavailable"))

        self._set_results(items)

    @staticmethod
    def _optimization_bounds(parameter: str, value: float) -> tuple[float, float]:
        if parameter == "Rc":
            if abs(value) < 1e-6:
                return (-100.0, 100.0)
            scale = max(abs(value) * 0.5, 5.0)
            return (value - scale, value + scale)
        if parameter == "Thickness":
            if value <= 0.0:
                return (0.01, 10.0)
            lower = max(0.01, value * 0.5)
            upper = max(lower + 0.5, value * 1.5)
            return (lower, upper)
        raise ValueError(f"Unsupported optimization parameter: {parameter}")

    def _build_optimization_variables(self) -> list[OpticalVariable]:
        variables: list[OpticalVariable] = []
        for index, row in enumerate(self.rows):
            if row.surface in {"Object", "Image"}:
                continue
            if row.optimize_rc and row.surface == "Standard":
                lower, upper = row.optimize_rc_bounds or self._optimization_bounds("Rc", row.rc)
                variables.append(
                    OpticalVariable(index, "Rc", lower, upper, name=f"{row.name} Rc")
                )
            if row.optimize_thickness:
                lower, upper = row.optimize_thickness_bounds or self._optimization_bounds("Thickness", row.thickness)
                variables.append(
                    OpticalVariable(index, "Thickness", lower, upper, name=f"{row.name} Thickness")
                )
        return variables

    def _build_merit_function(self, merit_mode: str) -> MeritFunction:
        operands = []
        if merit_mode in {"Spot RMS", "Spot + Wavefront"}:
            operands.append(
                SpotRMSOperand(
                    name="Spot RMS",
                    weight=1.0,
                    target=0.0,
                    surface_index=-1,
                    wavelength=self._current_wavelength(),
                    ray_count=max(5, self._current_ray_count()),
                    ray_height_factor=self._current_ray_height_factor(),
                )
            )
        if merit_mode in {"Wavefront RMS", "Spot + Wavefront"}:
            operands.append(
                WavefrontRMSOperand(
                    name="Wavefront RMS",
                    weight=1e-2 if merit_mode == "Spot + Wavefront" else 1.0,
                    target=0.0,
                    surface_index=self._analysis_surface_index(),
                    wavelength=self._current_wavelength(),
                    aperture_type=self._current_aperture_type(),
                    aperture_value=self._current_aperture_value(),
                    sample_size=9,
                )
            )
        return MeritFunction(operands=operands)

    def start_optimization(self) -> None:
        if self.optimization_running:
            return
        self._read_rows_from_table()
        variables = self._build_optimization_variables()
        if not variables:
            self.append_progress("Optimization skipped: no Rc/Thickness cells marked.")
            return

        try:
            system = self.build_system()
        except Exception as exc:
            self.append_progress(f"Optimization aborted: system build failed: {exc}")
            return

        merit_mode = self.merit_mode_var.get().strip()
        merit = self._build_merit_function(merit_mode)
        if not merit.operands:
            self.append_progress("Optimization aborted: no merit operands selected.")
            return

        x0 = []
        for variable in variables:
            row = self.rows[variable.surface_index]
            x0.append(row.rc if variable.parameter == "Rc" else row.thickness)

        evaluator = MeritEvaluator(system.SDT, setup=system.SETUP, merit_function=merit)
        initial = evaluator.evaluate(variables, x0)
        self.status_var.set(f"Optimization running: initial merit = {initial.total:.6g}")
        self.append_progress(f"Optimization start | merit mode: {merit_mode}")
        self.append_progress(f"Variables: {', '.join(v.normalized_name() for v in variables)}")
        self.append_progress(f"Initial merit: {initial.total:.6g}")
        self.update_idletasks()

        try:
            import ctypes
            pagmo_lib = Path(os.path.expanduser("~/Projects/pagmo2/_install/lib64/libpagmo.so"))
            if pagmo_lib.exists():
                try:
                    ctypes.CDLL(str(pagmo_lib), mode=ctypes.RTLD_GLOBAL)
                except OSError:
                    pass
            import pygmo as pg  # type: ignore
        except Exception as exc:
            self.append_progress(f"Optimization aborted: pygmo unavailable: {exc}")
            return

        udp = Pygmo2MeritProblem(evaluator=evaluator, variables=variables)
        problem = pg.problem(udp)
        population = pg.population(problem, size=12, seed=42)
        population.push_back(x0)

        self.optimization_running = True
        self.optimization_cancel_requested = False
        self.optimization_context = {
            "pg": pg,
            "variables": variables,
            "evaluator": evaluator,
            "initial": initial,
            "population": population,
            "generations_total": 12,
            "generation_done": 0,
            "verbosity_every": 3,
        }
        self.after(0, self._optimization_step)

    def stop_optimization(self) -> None:
        if not self.optimization_running:
            self.append_progress("Stop ignored: no optimization is running.")
            return
        self.optimization_cancel_requested = True
        self.append_progress("Stop requested. Optimization will stop after the current generation.")

    def _optimization_step(self) -> None:
        if not self.optimization_running or self.optimization_context is None:
            return

        ctx = self.optimization_context
        if self.optimization_cancel_requested:
            self._finish_optimization(cancelled=True)
            return

        pg = ctx["pg"]
        algorithm = pg.algorithm(pg.de(gen=1, seed=42 + ctx["generation_done"]))
        algorithm.set_verbosity(1)
        capture = io.StringIO()
        try:
            with redirect_stdout(capture), redirect_stderr(capture):
                ctx["population"] = algorithm.evolve(ctx["population"])
        except Exception as exc:
            self.append_debug(capture.getvalue())
            self.append_progress(f"Optimization failed at generation {ctx['generation_done'] + 1}: {exc}")
            self._finish_optimization(cancelled=True)
            return

        self.append_debug(capture.getvalue())
        ctx["generation_done"] += 1
        logs = algorithm.extract(pg.de).get_log()
        if logs:
            gen, fevals, best, dx, df = logs[-1]
            if (
                ctx["generation_done"] == 1
                or ctx["generation_done"] == ctx["generations_total"]
                or ctx["generation_done"] % ctx["verbosity_every"] == 0
            ):
                self.append_progress(
                    f"Gen {int(ctx['generation_done']):>3} | fevals {int(fevals):>4} | best {float(best):.6g} | dx {float(dx):.6g} | df {float(df):.6g}"
                )

        if ctx["generation_done"] >= ctx["generations_total"]:
            self._finish_optimization(cancelled=False)
            return

        self.after(1, self._optimization_step)

    def _finish_optimization(self, cancelled: bool) -> None:
        if self.optimization_context is None:
            self.optimization_running = False
            self.optimization_cancel_requested = False
            return

        ctx = self.optimization_context
        population = ctx["population"]
        champion_x = population.champion_x
        champion = ctx["evaluator"].evaluate(ctx["variables"], champion_x)

        for variable, value in zip(ctx["variables"], champion_x):
            row = self.rows[variable.surface_index]
            if variable.parameter == "Rc":
                row.rc = float(value)
            elif variable.parameter == "Thickness":
                row.thickness = float(value)

        self._sync_table()
        self.refresh_plot()
        initial = ctx["initial"]
        if cancelled:
            self.status_var.set(
                f"Optimization stopped: {initial.total:.6g} -> {champion.total:.6g}"
            )
            self.append_progress(
                f"Optimization stopped | merit {initial.total:.6g} -> {champion.total:.6g}"
            )
        else:
            self.status_var.set(
                f"Optimization finished: {initial.total:.6g} -> {champion.total:.6g}"
            )
            self.append_progress(
                f"Optimization finished | merit {initial.total:.6g} -> {champion.total:.6g}"
            )
        for operand in champion.operands:
            self.append_progress(
                f"  {operand.name}: value={operand.value:.6g} weighted={operand.weighted:.6g}"
            )

        self.optimization_context = None
        self.optimization_running = False
        self.optimization_cancel_requested = False

    def _sample_ray_heights(self, max_radius: float) -> list[float]:
        if max_radius <= 1e-9:
            return [0.0]
        count = self._current_ray_count()
        span = max_radius * self._current_ray_height_factor()
        if count == 1:
            return [0.0]
        return list(np.linspace(-span, span, count))

    def _current_ray_count(self) -> int:
        try:
            return max(1, int(self.ray_count_var.get()))
        except ValueError:
            return 5

    def _current_ray_height_factor(self) -> float:
        try:
            factor = float(self.ray_height_factor_var.get())
        except ValueError:
            factor = 0.8
        return max(min(factor, 1.5), 0.05)

    def _plot_fallback_preview(self, max_radius: float) -> None:
        positions = []
        z = 0.0
        for row in self.rows:
            positions.append(z)
            radius = max(row.diameter / 2.0, 0.5)
            color = "#4f81bd" if row.glass.upper() != "AIR" else "#7f8c8d"
            self.ax.plot([z, z], [-radius, radius], color=color, linewidth=2)
            self.ax.text(
                z,
                radius + max_radius * 0.08,
                row.name,
                rotation=45,
                ha="left",
                va="bottom",
                fontsize=8,
            )
            z += row.thickness

        total_length = max(z, 1.0)
        margin = max(total_length * 0.05, 5.0)
        self.ax.set_xlim(-margin, total_length + margin)
        self.ax.set_ylim(-(max_radius * 1.4), max_radius * 1.4)
        self.ax.axhline(0.0, color="#2c3e50", linewidth=0.8)

    def open_layout(self) -> None:
        path = filedialog.askopenfilename(
            title="Open Kraken layout",
            initialdir=str(LAYOUTS_DIR),
            filetypes=[("Python layout", "*.py")],
        )
        if not path:
            return
        self.current_layout_file = Path(path)
        try:
            info = _load_python_data(Path(path))
            self.rows = [
                SurfaceRow(
                    surface=str(item.get("surface", self._infer_surface_type(item))),
                    name=str(item.get("name", "Surface")),
                    optimize_rc=_coerce_opt_flag(item.get("optimize_rc", item.get("opt_rc", ""))),
                    optimize_rc_bounds=_coerce_bounds(item.get("optimize_rc_bounds")),
                    rc=float(item.get("rc", 0.0)),
                    optimize_thickness=_coerce_opt_flag(item.get("optimize_thickness", item.get("opt_thickness", ""))),
                    optimize_thickness_bounds=_coerce_bounds(item.get("optimize_thickness_bounds")),
                    thickness=float(item.get("thickness", 0.0)),
                    diameter=float(item.get("diameter", 25.0)),
                    glass=str(item.get("glass", "AIR")),
                )
                for item in info["surfaces"]
            ]
        except Exception:
            surfaces = self._extract_surfaces_from_example(Path(path))
            self.rows = [self._row_from_surface(surface, index, len(surfaces)) for index, surface in enumerate(surfaces)]
        self._normalize_special_rows()
        self._sync_table()
        self.refresh_plot()
        self.status_var.set(f"Opened {Path(path).name}")

    def save_layout(self) -> None:
        if self.current_layout_file is None:
            self.save_layout_as()
            return
        self._write_layout_file(self.current_layout_file)

    def save_layout_as(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Save Kraken layout",
            initialdir=str(LAYOUTS_DIR),
            defaultextension=".py",
            filetypes=[("Python layout", "*.py")],
        )
        if not path:
            return
        self.current_layout_file = Path(path)
        self._write_layout_file(self.current_layout_file)
        self.load_layouts()

    def _write_layout_file(self, path: Path) -> None:
        self._read_rows_from_table()
        title = path.stem.replace("_", " ").title()
        lines = [
            "#!/usr/bin/env python3",
            f'TITLE = "{title}"',
            "",
            "import KrakenOS as Kos",
            "",
            "",
            "def build_system():",
            "    surfaces = []",
        ]
        for index, row in enumerate(self.rows):
            var_name = f"s{index}"
            rc_bounds_repr = repr(tuple(row.optimize_rc_bounds)) if row.optimize_rc_bounds is not None else "None"
            thickness_bounds_repr = repr(tuple(row.optimize_thickness_bounds)) if row.optimize_thickness_bounds is not None else "None"
            lines.extend(
                [
                    f"    {var_name} = Kos.surf()",
                    f"    {var_name}.Name = {row.name!r}",
                    f"    {var_name}.Rc = {float(row.rc)!r}",
                    f"    {var_name}.Thickness = {float(row.thickness)!r}",
                    f"    {var_name}.Diameter = {float(row.diameter)!r}",
                    f"    {var_name}.Glass = {row.glass!r}",
                ]
            )
            if row.optimize_rc:
                lines.append(f"    {var_name}.optimize_rc = True")
            if row.optimize_rc_bounds is not None:
                lines.append(f"    {var_name}.optimize_rc_bounds = {tuple(row.optimize_rc_bounds)!r}")
            if row.optimize_thickness:
                lines.append(f"    {var_name}.optimize_thickness = True")
            if row.optimize_thickness_bounds is not None:
                lines.append(f"    {var_name}.optimize_thickness_bounds = {tuple(row.optimize_thickness_bounds)!r}")
            if row.surface == "Thin Lens":
                focal = float(row.rc) if float(row.rc) != 0.0 else 100.0
                lines.append(f"    {var_name}.Thin_Lens = {focal!r}")
                lines.append(f"    {var_name}.Rc = 0.0")
            elif row.surface == "Grating":
                lines.append(f"    {var_name}.Diff_Ord = 1.0")
                lines.append(f"    {var_name}.Grating_D = 1.0")
            lines.append(
                "    surfaces.append({"
                f"'surface': {row.surface!r}, "
                f"'name': {row.name!r}, "
                f"'rc': {float(row.rc)!r}, "
                f"'thickness': {float(row.thickness)!r}, "
                f"'diameter': {float(row.diameter)!r}, "
                f"'glass': {row.glass!r}, "
                f"'optimize_rc': {row.optimize_rc!r}, "
                f"'optimize_rc_bounds': {rc_bounds_repr}, "
                f"'optimize_thickness': {row.optimize_thickness!r}, "
                f"'optimize_thickness_bounds': {thickness_bounds_repr}"
                "})"
            )
            lines.append("")
        lines.extend(
            [
                "    return surfaces",
                "",
                "",
                "SURFACES = build_system()",
                "",
                "",
                "def build_runtime_system():",
                "    surface_dicts = SURFACES",
                "    runtime_surfaces = []",
                "    for spec in surface_dicts:",
                "        s = Kos.surf()",
                "        s.Name = spec['name']",
                "        s.Rc = spec['rc']",
                "        s.Thickness = spec['thickness']",
                "        s.Diameter = spec['diameter']",
                "        s.Glass = spec['glass']",
                "        if spec['surface'] == 'Thin Lens':",
                "            s.Thin_Lens = spec['rc'] if spec['rc'] != 0 else 100.0",
                "            s.Rc = 0.0",
                "        elif spec['surface'] == 'Grating':",
                "            s.Diff_Ord = 1.0",
                "            s.Grating_D = 1.0",
                "        runtime_surfaces.append(s)",
                "    setup = Kos.Setup()",
                "    return Kos.system(runtime_surfaces, setup)",
                "",
                "",
                "def build_rays(system):",
                "    rays = Kos.raykeeper(system)",
                "    max_radius = max((float(s.Diameter) / 2.0 for s in system.SDT), default=1.0)",
                "    ray_heights = [(-0.8 * max_radius), (-max_radius / 3.0), 0.0, (max_radius / 3.0), (0.8 * max_radius)]",
                "    for y0 in ray_heights:",
                "        system.Trace([0.0, y0, 0.0], [0.0, 0.0, 1.0], 0.55)",
                "        rays.push()",
                "    return rays",
                "",
                "",
                "if __name__ == '__main__':",
                "    system = build_runtime_system()",
                "    rays = build_rays(system)",
                "    Kos.display2d(system, rays, 0)",
                "",
            ]
        )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.status_var.set(f"Saved {path.name}")

    def _extract_surfaces_from_example(self, path: Path):
        original_system = Kos.system
        original_display2d = getattr(Kos, "display2d", None)
        original_display3d = getattr(Kos, "display3d", None)
        original_display2d_colab = getattr(Kos, "display2d_colab", None)

        def capture_system(surf_data, setup, build=1):
            raise _CapturedExample(list(surf_data))

        try:
            Kos.system = capture_system
            if original_display2d is not None:
                Kos.display2d = lambda *args, **kwargs: None
            if original_display3d is not None:
                Kos.display3d = lambda *args, **kwargs: None
            if original_display2d_colab is not None:
                Kos.display2d_colab = lambda *args, **kwargs: None

            namespace = {
                "__name__": "__main__",
                "__file__": str(path),
            }
            code = path.read_text(encoding="utf-8", errors="ignore")
            try:
                previous_cwd = os.getcwd()
                os.chdir(path.parent)
                exec(compile(code, str(path), "exec"), namespace, namespace)
            except _CapturedExample as captured:
                return captured.surfaces
            finally:
                os.chdir(previous_cwd)
        finally:
            Kos.system = original_system
            if original_display2d is not None:
                Kos.display2d = original_display2d
            if original_display3d is not None:
                Kos.display3d = original_display3d
            if original_display2d_colab is not None:
                Kos.display2d_colab = original_display2d_colab

        raise ValueError("No KrakenOS system definition was captured from the example.")

    @staticmethod
    def _row_from_surface(surface, index: int, total: int) -> SurfaceRow:
        surface_type = "Standard"
        if index == 0:
            surface_type = "Object"
        elif index == total - 1:
            surface_type = "Image"
        elif getattr(surface, "Thin_Lens", 0.0) != 0:
            surface_type = "Thin Lens"
        elif getattr(surface, "Diff_Ord", 0.0) != 0:
            surface_type = "Grating"

        rc_value = float(getattr(surface, "Rc", 0.0))
        if surface_type == "Thin Lens":
            rc_value = float(getattr(surface, "Thin_Lens", 0.0))

        return SurfaceRow(
            surface=surface_type,
            name=str(getattr(surface, "Name", "") or f"Surface {index}"),
            rc=rc_value,
            thickness=float(getattr(surface, "Thickness", 0.0)),
            diameter=float(getattr(surface, "Diameter", 25.0)),
            glass=str(getattr(surface, "Glass", "AIR")),
        )

    @staticmethod
    def _infer_surface_type(item: dict) -> str:
        if "surface" in item:
            return str(item["surface"])
        name = str(item.get("name", "")).strip().lower()
        if name == "object":
            return "Object"
        if name == "image":
            return "Image"
        return "Standard"

    def _normalize_special_rows(self) -> None:
        if not self.rows:
            return
        self.rows[0].surface = "Object"
        if not self.rows[0].name or self.rows[0].name == "Surface":
            self.rows[0].name = "Object"
        self.rows[-1].surface = "Image"
        if not self.rows[-1].name or self.rows[-1].name == "Surface":
            self.rows[-1].name = "Image"

    @staticmethod
    def _flipped_name(name: str) -> str:
        placeholder_front = "__KR_FRONT__"
        placeholder_back = "__KR_BACK__"
        placeholder_left = "__KR_LEFT__"
        placeholder_right = "__KR_RIGHT__"
        placeholder_entry = "__KR_ENTRY__"
        placeholder_exit = "__KR_EXIT__"
        value = name
        value = re.sub(r"\bFront\b", placeholder_front, value, flags=re.IGNORECASE)
        value = re.sub(r"\bBack\b", placeholder_back, value, flags=re.IGNORECASE)
        value = re.sub(r"\bLeft\b", placeholder_left, value, flags=re.IGNORECASE)
        value = re.sub(r"\bRight\b", placeholder_right, value, flags=re.IGNORECASE)
        value = re.sub(r"\bEntry\b", placeholder_entry, value, flags=re.IGNORECASE)
        value = re.sub(r"\bExit\b", placeholder_exit, value, flags=re.IGNORECASE)
        value = value.replace(placeholder_front, "Back")
        value = value.replace(placeholder_back, "Front")
        value = value.replace(placeholder_left, "Right")
        value = value.replace(placeholder_right, "Left")
        value = value.replace(placeholder_entry, "Exit")
        value = value.replace(placeholder_exit, "Entry")
        return value


def main() -> None:
    app = KrakenLayoutEditor()
    app.mainloop()


if __name__ == "__main__":
    main()
