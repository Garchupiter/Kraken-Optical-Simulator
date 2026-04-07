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
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk
import warnings

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

import KrakenOS as Kos
import pyvista as pv
from KrakenOS.Display import Plot2DRays, Plot2DSurf, plot3d, rayplot3d
from KrakenOS.Optimization import (
    OPERAND_REGISTRY,
    VARIABLE_REGISTRY,
    MeritEvaluator,
    MeritFunction,
    OpticalVariable,
)
from KrakenOS.Optimization.adapters.pygmo2_adapter import Pygmo2MeritProblem


LAYOUTS_DIR = Path(__file__).resolve().parent.parent / "common_optical_layouts"
EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "Examples"
DEFAULT_LAYOUT_TITLE = "Double Mirror Fold"
FIELDS = (
    "label",
    "surface",
    "name",
    "rc",
    "thickness",
    "diameter",
    "tilt_x",
    "tilt_y",
    "tilt_z",
    "desp_x",
    "desp_y",
    "desp_z",
    "axis_move",
    "glass",
)
COLUMN_LABELS = {
    "label": "#",
    "surface": "Surface",
    "name": "Name",
    "rc": "Rc [mm]",
    "thickness": "Thickness [mm]",
    "diameter": "Diameter [mm]",
    "tilt_x": "TiltX [deg]",
    "tilt_y": "TiltY [deg]",
    "tilt_z": "TiltZ [deg]",
    "desp_x": "DespX [mm]",
    "desp_y": "DespY [mm]",
    "desp_z": "DespZ [mm]",
    "axis_move": "AxisMove",
    "glass": "Glass",
}
NUMERIC_FIELDS = {
    "rc",
    "thickness",
    "diameter",
    "tilt_x",
    "tilt_y",
    "tilt_z",
    "desp_x",
    "desp_y",
    "desp_z",
    "axis_move",
}
SURFACE_TYPES = ("Object", "Standard", "Mirror", "Thin Lens", "Grating", "Image")


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
    tilt_x: float = 0.0
    tilt_y: float = 0.0
    tilt_z: float = 0.0
    desp_x: float = 0.0
    desp_y: float = 0.0
    desp_z: float = 0.0
    axis_move: float = 0.0
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
        self.operand_weight_vars: dict[str, tk.StringVar] = {}
        self.operand_target_vars: dict[str, tk.StringVar] = {}
        self.operand_wavelength_vars: dict[str, tk.StringVar] = {}
        self.operand_field_vars: dict[str, tk.StringVar] = {}
        self.operand_field_x_vars: dict[str, tk.StringVar] = {}
        self.operand_field_y_vars: dict[str, tk.StringVar] = {}
        self.operand_surface_vars: dict[str, tk.StringVar] = {}
        self.operand_aperture_type_vars: dict[str, tk.StringVar] = {}
        self.operand_aperture_value_vars: dict[str, tk.StringVar] = {}
        self.operand_frequency_vars: dict[str, tk.StringVar] = {}
        self.operand_mtf_mode_vars: dict[str, tk.StringVar] = {}
        self.operand_control_widgets: dict[str, dict[str, tuple[tk.Widget, ...]]] = {}
        self.operand_setup_frames: dict[str, tk.Widget] = {}
        self._spinner_phase = 0
        self._spinner_after_id: str | None = None
        self._refresh_after_id: str | None = None
        self._preview_field_ray_count = 1
        self._active_cell: tuple[str, str] | None = None
        self._cell_border_parts: list[tk.Frame] = []
        self._grid_overlays: list[tk.Frame] = []
        self._grid_after_id: str | None = None
        self._initial_layout_passes = 0
        self._last_field_type = "Angle"
        self._field_defaults_initialized = False
        self._field_type_defaults = {
            "Angle": "5.0",
            "Object Height": "5.0",
            "Paraxial Image Height": "5.0",
            "Real Image Height": "5.0",
        }

        self._build_menu()
        self._build_ui()
        self.load_layouts()
        self.load_examples()
        if self.layout_names:
            initial_layout = DEFAULT_LAYOUT_TITLE if DEFAULT_LAYOUT_TITLE in self.layout_files else self.layout_names[0]
            self.load_layout_by_name(initial_layout)

    def _build_menu(self) -> None:
        menubar = tk.Menu(self)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open", command=self.open_layout)
        file_menu.add_command(label="Save", command=self.save_layout)
        file_menu.add_command(label="Save As", command=self.save_layout_as)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        action_menu = tk.Menu(menubar, tearoff=0)
        action_menu.add_command(label="Refresh Plot", command=self.refresh_plot)
        action_menu.add_command(label="Clear Marks", command=self.clear_optimization_marks)
        menubar.add_cascade(label="Actions", menu=action_menu)

        self.config(menu=menubar)

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        style = ttk.Style(self)
        style.configure(
            "Excel.Treeview",
            background="white",
            fieldbackground="white",
            bordercolor="#d9dee7",
            borderwidth=1,
            relief="solid",
            rowheight=24,
        )
        style.configure(
            "Excel.Treeview.Heading",
            bordercolor="#d9dee7",
            borderwidth=1,
            relief="solid",
        )
        style.map(
            "Excel.Treeview",
            background=[("selected", "#eaf2ff")],
            foreground=[("selected", "black")],
        )

        main = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        main.grid(row=0, column=0, sticky="nsew")
        self.main_pane = main

        left_panel = ttk.Panedwindow(main, orient=tk.VERTICAL)
        main.add(left_panel, weight=2)

        top = ttk.Panedwindow(left_panel, orient=tk.HORIZONTAL)
        left_panel.add(top, weight=1)

        control_stack = ttk.Frame(top, padding=8)
        control_stack.columnconfigure(0, weight=1)
        control_stack.rowconfigure(0, weight=0)
        control_stack.rowconfigure(1, weight=0)
        top.add(control_stack, weight=1)

        controls = ttk.LabelFrame(control_stack, text="Display", padding=8)
        controls.grid(row=0, column=0, sticky="ew")
        controls.columnconfigure(0, weight=1)
        controls.columnconfigure(1, weight=1)

        field_panel = ttk.LabelFrame(control_stack, text="Field", padding=8)
        field_panel.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        field_panel.columnconfigure(0, weight=1)
        field_panel.columnconfigure(1, weight=1)

        table_frame = ttk.Frame(top, padding=8)
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(1, weight=1)
        top.add(table_frame, weight=4)

        plot_frame = ttk.Frame(left_panel, padding=8)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(1, weight=1)
        left_panel.add(plot_frame, weight=4)

        right_panel = ttk.Panedwindow(main, orient=tk.VERTICAL)
        main.add(right_panel, weight=1)

        info_stack = ttk.Panedwindow(right_panel, orient=tk.VERTICAL)
        right_panel.add(info_stack, weight=3)

        results = ttk.LabelFrame(info_stack, text="Information", padding=8)
        results.columnconfigure(0, weight=1)
        results.rowconfigure(0, weight=1)
        info_stack.add(results, weight=3)

        optimization = ttk.LabelFrame(info_stack, text="Optimization", padding=8)
        optimization.columnconfigure(0, weight=1)
        optimization.rowconfigure(2, weight=0)
        info_stack.add(optimization, weight=2)

        bottom = ttk.Panedwindow(right_panel, orient=tk.VERTICAL)
        right_panel.add(bottom, weight=1)

        progress_frame = ttk.LabelFrame(bottom, text="Progress", padding=8)
        progress_frame.columnconfigure(0, weight=1)
        progress_frame.rowconfigure(1, weight=1)
        bottom.add(progress_frame, weight=1)

        debug_frame = ttk.LabelFrame(bottom, text="Debug", padding=8)
        debug_frame.columnconfigure(0, weight=1)
        debug_frame.rowconfigure(0, weight=1)
        bottom.add(debug_frame, weight=1)

        table_toolbar = ttk.Frame(table_frame)
        table_toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        ttk.Button(table_toolbar, text="Add surface", command=self.add_surface).pack(side="left")
        ttk.Button(table_toolbar, text="Delete", command=self.delete_selected).pack(side="left", padx=(6, 0))
        ttk.Button(table_toolbar, text="Duplicate", command=self.duplicate_selected).pack(side="left", padx=(6, 0))
        ttk.Button(table_toolbar, text="Flip", command=self.flip_selected).pack(side="left", padx=(6, 0))
        ttk.Button(table_toolbar, text="▲", width=3, command=self.move_up).pack(side="left", padx=(10, 0))
        ttk.Button(table_toolbar, text="▼", width=3, command=self.move_down).pack(side="left", padx=(4, 0))

        self.layout_var = tk.StringVar(value="Common Optical Layout")
        self.layout_menu = ttk.Combobox(
            table_toolbar,
            textvariable=self.layout_var,
            state="readonly",
            width=28,
        )
        self.layout_menu.pack(side="left", padx=(12, 0))
        self.layout_menu.bind("<<ComboboxSelected>>", self._on_layout_selected)

        self.example_var = tk.StringVar(value="Examples")
        self.example_menu = ttk.Combobox(
            table_toolbar,
            textvariable=self.example_var,
            state="readonly",
            width=28,
        )
        self.example_menu.pack(side="left", padx=(8, 0))
        self.example_menu.bind("<<ComboboxSelected>>", self._on_example_selected)

        self.table = ttk.Treeview(
            table_frame,
            columns=FIELDS,
            show="headings",
            selectmode="extended",
            style="Excel.Treeview",
        )
        for field in FIELDS:
            self.table.heading(field, text=COLUMN_LABELS[field])
            width = (
                55 if field == "label"
                else 140 if field == "surface"
                else 160 if field == "name"
                else 95 if field in {"tilt_x", "tilt_y", "tilt_z"}
                else 95 if field in {"desp_x", "desp_y", "desp_z"}
                else 85 if field == "axis_move"
                else 110
            )
            self.table.column(field, width=width, stretch=True, anchor="center")
        self.table.grid(row=1, column=0, sticky="nsew")
        self.table.bind("<Button-1>", self._on_table_click, add="+")
        self.table.bind("<Double-1>", self.begin_edit)
        self.table.bind("<Button-3>", self.show_context_menu)
        self.table.bind("<<TreeviewSelect>>", self._update_active_cell_border, add="+")
        self.table.bind("<Configure>", self._update_active_cell_border, add="+")
        self.table.bind("<Configure>", self._schedule_table_grid_update, add="+")
        self.table.bind("<MouseWheel>", self._update_active_cell_border, add="+")
        self.table.bind("<MouseWheel>", self._schedule_table_grid_update, add="+")
        self.table.bind("<Button-4>", self._update_active_cell_border, add="+")
        self.table.bind("<Button-4>", self._schedule_table_grid_update, add="+")
        self.table.bind("<Button-5>", self._update_active_cell_border, add="+")
        self.table.bind("<Button-5>", self._schedule_table_grid_update, add="+")
        self.table.bind("<Left>", self._move_active_cell)
        self.table.bind("<Right>", self._move_active_cell)
        self.table.bind("<Up>", self._move_active_cell)
        self.table.bind("<Down>", self._move_active_cell)
        self.table.bind("<KP_Left>", self._move_active_cell)
        self.table.bind("<KP_Right>", self._move_active_cell)
        self.table.bind("<KP_Up>", self._move_active_cell)
        self.table.bind("<KP_Down>", self._move_active_cell)
        self.table.tag_configure("optimize", background="#fff4bf")

        border_color = "#4a89ff"
        self._cell_border_parts = [
            tk.Frame(self.table, bg=border_color, height=2, width=2),
            tk.Frame(self.table, bg=border_color, height=2, width=2),
            tk.Frame(self.table, bg=border_color, height=2, width=2),
            tk.Frame(self.table, bg=border_color, height=2, width=2),
        ]
        self._hide_active_cell_border()

        yscroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.table.yview)
        yscroll.grid(row=1, column=1, sticky="ns")
        self.table.configure(yscrollcommand=lambda first, last: self._on_table_scroll(yscroll, first, last))
        xscroll = ttk.Scrollbar(table_frame, orient="horizontal", command=self.table.xview)
        xscroll.grid(row=2, column=0, sticky="ew")
        self.table.configure(xscrollcommand=xscroll.set)

        self._build_controls_panel(controls)
        self._build_field_panel(field_panel)
        self._build_results_panel(results)
        self._build_optimization_panel(optimization)

        plot_toolbar = ttk.Frame(plot_frame)
        plot_toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        ttk.Button(plot_toolbar, text="Open 3D", command=self.open_3d_view).pack(side="left")
        ttk.Button(plot_toolbar, text="2D", command=lambda: self.set_analysis_mode("none")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="Native", command=lambda: self.set_analysis_mode("native_off_axis")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="Spot", command=lambda: self.set_analysis_mode("spot")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="PSF", command=lambda: self.set_analysis_mode("psf")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="RMS", command=lambda: self.set_analysis_mode("rms")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="FC/Dist", command=lambda: self.set_analysis_mode("field_curvature")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="Pupil", command=lambda: self.set_analysis_mode("pupil")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="Seidel", command=lambda: self.set_analysis_mode("seidel")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="Wavefront", command=lambda: self.set_analysis_mode("wavefront")).pack(side="left", padx=(6, 0))
        ttk.Button(plot_toolbar, text="MTF", command=lambda: self.set_analysis_mode("mtf")).pack(side="left", padx=(6, 0))
        ttk.Checkbutton(
            plot_toolbar,
            text="Show PP / EP / XP",
            variable=self.show_cardinals_var,
            command=self.refresh_plot,
        ).pack(side="left", padx=(12, 0))

        self.figure = Figure(figsize=(7, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")

        self.debug_text = tk.Text(debug_frame, wrap="word", height=8, width=24)
        self.debug_text.grid(row=0, column=0, sticky="nsew")
        debug_scroll = ttk.Scrollbar(debug_frame, orient="vertical", command=self.debug_text.yview)
        debug_scroll.grid(row=0, column=1, sticky="ns")
        self.debug_text.configure(yscrollcommand=debug_scroll.set)

        status_bar = ttk.Frame(progress_frame)
        status_bar.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        status_bar.columnconfigure(1, weight=1)
        self.progress_spinner_var = tk.StringVar(value="idle")
        self.progress_percent_var = tk.StringVar(value="0%")
        self.progress_bar_var = tk.DoubleVar(value=0.0)
        ttk.Label(status_bar, textvariable=self.progress_spinner_var, width=4).grid(row=0, column=0, sticky="w")
        self.progress_bar = ttk.Progressbar(
            status_bar,
            orient="horizontal",
            mode="determinate",
            maximum=100.0,
            variable=self.progress_bar_var,
        )
        self.progress_bar.grid(row=0, column=1, sticky="ew", padx=(0, 8))
        ttk.Label(status_bar, textvariable=self.progress_percent_var, width=14).grid(
            row=0, column=2, sticky="e"
        )

        self.progress_text = tk.Text(progress_frame, wrap="word", height=8, width=24)
        self.progress_text.grid(row=1, column=0, sticky="nsew")
        progress_scroll = ttk.Scrollbar(progress_frame, orient="vertical", command=self.progress_text.yview)
        progress_scroll.grid(row=1, column=1, sticky="ns")
        self.progress_text.configure(yscrollcommand=progress_scroll.set)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(self, textvariable=self.status_var, anchor="w").grid(
            row=1, column=0, sticky="ew", padx=8, pady=(0, 2)
        )
        self.after_idle(self._set_initial_pane_layout)
        self.bind("<Configure>", self._maybe_refresh_initial_pane_layout, add="+")

    def _build_controls_panel(self, parent) -> None:
        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=1)

        ttk.Label(parent, text="Object mode").grid(row=0, column=0, sticky="w", pady=(0, 2))
        self.object_mode_var = tk.StringVar(value="Infinity")
        self.object_mode_menu = ttk.Combobox(
            parent,
            textvariable=self.object_mode_var,
            state="readonly",
            values=["Finite", "Infinity"],
        )
        self.object_mode_menu.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self.object_mode_menu.bind("<<ComboboxSelected>>", self._on_object_mode_changed)

        ttk.Label(parent, text="Wavelength [um]").grid(row=0, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.wavelength_var = tk.StringVar(value="0.55")
        ttk.Entry(parent, textvariable=self.wavelength_var, width=12).grid(
            row=1, column=1, sticky="ew", pady=(0, 8), padx=(8, 0)
        )

        ttk.Label(parent, text="Orientation").grid(row=2, column=0, sticky="w", pady=(0, 2))
        self.display_orientation_var = tk.StringVar(value="Vertical")
        self.display_orientation_menu = ttk.Combobox(
            parent,
            textvariable=self.display_orientation_var,
            state="readonly",
            values=["Vertical", "Horizontal"],
        )
        self.display_orientation_menu.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        self.display_orientation_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())

        ttk.Label(parent, text="Ray fan count").grid(row=2, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.ray_count_var = tk.StringVar(value="5")
        ttk.Entry(parent, textvariable=self.ray_count_var, width=12).grid(
            row=3, column=1, sticky="ew", pady=(0, 8), padx=(8, 0)
        )

        ttk.Label(parent, text="Pupil factor").grid(row=4, column=0, sticky="w", pady=(0, 2))
        self.ray_height_factor_var = tk.StringVar(value="0.8")
        ttk.Entry(parent, textvariable=self.ray_height_factor_var, width=12).grid(
            row=5, column=0, sticky="ew", pady=(0, 8)
        )

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

        ttk.Label(parent, text="Aperture type").grid(row=6, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.aperture_type_var = tk.StringVar(value="STOP")
        self.aperture_type_menu = ttk.Combobox(
            parent,
            textvariable=self.aperture_type_var,
            state="readonly",
            values=["STOP", "EPD"],
        )
        self.aperture_type_menu.grid(row=7, column=1, sticky="ew", pady=(0, 8), padx=(8, 0))
        self.aperture_type_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())

        ttk.Label(parent, text="Aperture value").grid(row=8, column=0, sticky="w", pady=(0, 2))
        self.aperture_value_var = tk.StringVar(value="1.0")
        ttk.Entry(parent, textvariable=self.aperture_value_var, width=12).grid(
            row=9, column=0, sticky="ew"
        )

        self.show_cardinals_var = tk.BooleanVar(value=True)

        for variable in (
            self.wavelength_var,
            self.ray_count_var,
            self.ray_height_factor_var,
            self.aperture_value_var,
        ):
            variable.trace_add("write", self._schedule_refresh_plot)

    def _build_field_panel(self, parent) -> None:
        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=1)

        ttk.Label(parent, text="Field type").grid(row=0, column=0, sticky="w", pady=(0, 2))
        self.field_type_var = tk.StringVar(value="Angle")
        self.field_type_menu = ttk.Combobox(
            parent,
            textvariable=self.field_type_var,
            state="readonly",
            values=["Angle", "Object Height", "Paraxial Image Height", "Real Image Height"],
        )
        self.field_type_menu.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        self.field_type_menu.bind("<<ComboboxSelected>>", self._on_field_type_changed)

        self.field_mode_note_var = tk.StringVar(value="")
        ttk.Label(parent, textvariable=self.field_mode_note_var, foreground="#475569", justify="left").grid(
            row=2, column=0, columnspan=2, sticky="ew", pady=(0, 6)
        )

        self.field_value_label_var = tk.StringVar(value="Angle [deg]")
        ttk.Label(parent, textvariable=self.field_value_label_var).grid(row=3, column=0, sticky="w", pady=(0, 2))
        self.field_value_var = tk.StringVar(value="5.0")
        ttk.Entry(parent, textvariable=self.field_value_var, width=12).grid(
            row=4, column=0, sticky="ew", pady=(0, 8)
        )

        ttk.Label(parent, text="Field count").grid(row=3, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.field_count_var = tk.StringVar(value="1")
        ttk.Entry(parent, textvariable=self.field_count_var, width=12).grid(
            row=4, column=1, sticky="ew", pady=(0, 8), padx=(8, 0)
        )

        self.field_warning_var = tk.StringVar(value="")
        ttk.Label(parent, textvariable=self.field_warning_var, foreground="#b45309", justify="left").grid(
            row=5, column=0, columnspan=2, sticky="ew", pady=(0, 4)
        )
        self.field_summary_var = tk.StringVar(value="")
        ttk.Label(parent, textvariable=self.field_summary_var, justify="left").grid(
            row=6, column=0, columnspan=2, sticky="ew"
        )

        for variable in (self.field_count_var, self.field_value_var):
            variable.trace_add("write", self._schedule_refresh_plot)
        self._sync_field_mode_ui()
    def _build_optimization_panel(self, parent) -> None:
        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=1)

        button_row = ttk.Frame(parent)
        button_row.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        ttk.Button(button_row, text="Start Optimization", command=self.start_optimization).pack(side="left")
        ttk.Button(button_row, text="Stop", command=self.stop_optimization).pack(side="left", padx=(8, 0))

        ttk.Label(parent, text="Merit operands").grid(row=1, column=0, sticky="w", pady=(0, 2))
        self.merit_mode_list = tk.Listbox(
            parent,
            exportselection=False,
            selectmode="extended",
            height=min(4, max(2, len(OPERAND_REGISTRY))),
            width=14,
        )
        for spec in OPERAND_REGISTRY.values():
            self.merit_mode_list.insert("end", spec.label)
        if OPERAND_REGISTRY:
            self.merit_mode_list.selection_set(0)
        self.merit_mode_list.grid(row=2, column=0, sticky="nsew", pady=(0, 8), padx=(0, 8))
        self.merit_mode_list.bind("<<ListboxSelect>>", lambda _e: self._update_operand_setup_visibility())

        setup_holder = ttk.Frame(parent, height=320)
        setup_holder.grid(row=2, column=1, sticky="nsew", pady=(0, 8))
        setup_holder.grid_propagate(False)
        setup_holder.columnconfigure(0, weight=1)
        setup_holder.rowconfigure(0, weight=1)

        setup_frame = ttk.Frame(setup_holder)
        setup_frame.grid(row=0, column=0, sticky="nsew")
        setup_frame.columnconfigure(0, weight=1)
        ttk.Label(setup_frame, text="Operand setup").grid(row=0, column=0, sticky="w", pady=(0, 2))

        for idx, spec in enumerate(OPERAND_REGISTRY.values(), start=1):
            card = ttk.LabelFrame(setup_frame, text=spec.label, padding=6)
            card.grid(row=idx, column=0, sticky="ew", pady=(0, 8))
            card.columnconfigure(1, weight=1)
            control_widgets: dict[str, tuple[tk.Widget, ...]] = {}

            weight_var = tk.StringVar(value=f"{spec.default_weight:g}")
            self.operand_weight_vars[spec.label] = weight_var
            weight_label = ttk.Label(card, text="Weight")
            weight_label.grid(row=0, column=0, sticky="w")
            weight_entry = ttk.Entry(card, textvariable=weight_var, width=6)
            weight_entry.grid(row=0, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(weight_entry)
            control_widgets["weight"] = (weight_label, weight_entry)

            target_var = tk.StringVar(value=f"{spec.default_target:g}")
            self.operand_target_vars[spec.label] = target_var
            target_label = ttk.Label(card, text="Target")
            target_label.grid(row=1, column=0, sticky="w")
            target_entry = ttk.Entry(card, textvariable=target_var, width=6)
            target_entry.grid(row=1, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(target_entry)
            control_widgets["target"] = (target_label, target_entry)

            wavelength_var = tk.StringVar(value=self.wavelength_var.get())
            self.operand_wavelength_vars[spec.label] = wavelength_var
            wavelength_label = ttk.Label(card, text="Wvl")
            wavelength_label.grid(row=2, column=0, sticky="w")
            wavelength_entry = ttk.Entry(card, textvariable=wavelength_var, width=6)
            wavelength_entry.grid(row=2, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(wavelength_entry)
            control_widgets["wavelength"] = (wavelength_label, wavelength_entry)

            field_var = tk.StringVar(value="0")
            self.operand_field_vars[spec.label] = field_var
            field_label = ttk.Label(card, text="Field")
            field_label.grid(row=3, column=0, sticky="w")
            field_entry = ttk.Entry(card, textvariable=field_var, width=6)
            field_entry.grid(row=3, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(field_entry)
            control_widgets["field"] = (field_label, field_entry)

            surface_row = 4
            aperture_type_row = 5
            aperture_value_row = 6
            frequency_row = 7
            mode_row = 8

            if spec.label == "MTF @ freq":
                field_x_var = tk.StringVar(value="0")
                field_y_var = tk.StringVar(value="0")
                self.operand_field_x_vars[spec.label] = field_x_var
                self.operand_field_y_vars[spec.label] = field_y_var
                field_x_label = ttk.Label(card, text="Field X")
                field_x_label.grid(row=3, column=0, sticky="w")
                field_x_entry = ttk.Entry(card, textvariable=field_x_var, width=6)
                field_x_entry.grid(row=3, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
                self._bind_deferred_refresh(field_x_entry)
                field_y_label = ttk.Label(card, text="Field Y")
                field_y_label.grid(row=4, column=0, sticky="w")
                field_y_entry = ttk.Entry(card, textvariable=field_y_var, width=6)
                field_y_entry.grid(row=4, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
                self._bind_deferred_refresh(field_y_entry)
                control_widgets["field_xy"] = (field_x_label, field_x_entry, field_y_label, field_y_entry)
                field_label.grid_remove()
                field_entry.grid_remove()
                surface_row = 5
                aperture_type_row = 6
                aperture_value_row = 7
                frequency_row = 8
                mode_row = 9

            surface_var = tk.StringVar(value="Auto")
            self.operand_surface_vars[spec.label] = surface_var
            surface_label = ttk.Label(card, text="Surf")
            surface_label.grid(row=surface_row, column=0, sticky="w")
            surface_menu = ttk.Combobox(card, textvariable=surface_var, state="readonly", width=6, values=["Auto"])
            surface_menu.grid(row=surface_row, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            surface_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())
            control_widgets["surface"] = (surface_label, surface_menu)

            aperture_type_var = tk.StringVar(value=self.aperture_type_var.get())
            self.operand_aperture_type_vars[spec.label] = aperture_type_var
            aperture_type_label = ttk.Label(card, text="Aper")
            aperture_type_label.grid(row=aperture_type_row, column=0, sticky="w")
            aperture_type_menu = ttk.Combobox(
                card,
                textvariable=aperture_type_var,
                state="readonly",
                width=6,
                values=["STOP", "EPD"],
            )
            aperture_type_menu.grid(row=aperture_type_row, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            aperture_type_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())
            control_widgets["aperture"] = (aperture_type_label, aperture_type_menu)

            aperture_value_var = tk.StringVar(value=self.aperture_value_var.get())
            self.operand_aperture_value_vars[spec.label] = aperture_value_var
            aperture_value_label = ttk.Label(card, text="AVal")
            aperture_value_label.grid(row=aperture_value_row, column=0, sticky="w")
            aperture_value_entry = ttk.Entry(card, textvariable=aperture_value_var, width=6)
            aperture_value_entry.grid(row=aperture_value_row, column=1, sticky="ew", padx=(6, 0), pady=(0, 0))
            self._bind_deferred_refresh(aperture_value_entry)
            control_widgets["aperture_value"] = (aperture_value_label, aperture_value_entry)

            if spec.label == "MTF @ freq":
                frequency_var = tk.StringVar(value="50")
                self.operand_frequency_vars[spec.label] = frequency_var
                frequency_label = ttk.Label(card, text="Freq")
                frequency_label.grid(row=frequency_row, column=0, sticky="w")
                frequency_entry = ttk.Entry(card, textvariable=frequency_var, width=6)
                frequency_entry.grid(row=frequency_row, column=1, sticky="ew", padx=(6, 0), pady=(0, 0))
                self._bind_deferred_refresh(frequency_entry)
                control_widgets["frequency"] = (frequency_label, frequency_entry)

                mtf_mode_var = tk.StringVar(value="Average")
                self.operand_mtf_mode_vars[spec.label] = mtf_mode_var
                mode_label = ttk.Label(card, text="Mode")
                mode_label.grid(row=mode_row, column=0, sticky="w")
                mtf_mode_menu = ttk.Combobox(
                    card,
                    textvariable=mtf_mode_var,
                    state="readonly",
                    width=8,
                    values=["Average", "Tangential", "Sagittal"],
                )
                mtf_mode_menu.grid(row=mode_row, column=1, sticky="ew", padx=(6, 0), pady=(4, 0))
                mtf_mode_menu.bind("<<ComboboxSelected>>", lambda _e: self.refresh_plot())
                control_widgets["mtf_mode"] = (mode_label, mtf_mode_menu)

            self.operand_control_widgets[spec.label] = control_widgets
            self.operand_setup_frames[spec.label] = card
            self._apply_operand_control_visibility(spec.label)

        self._update_operand_setup_visibility()

    def _build_results_panel(self, parent) -> None:
        self.results_table = ttk.Treeview(parent, columns=("property", "value"), show="headings", selectmode="none")
        self.results_table.heading("property", text="Property")
        self.results_table.heading("value", text="Value")
        self.results_table.column("property", width=96, anchor="w", stretch=False)
        self.results_table.column("value", width=40, anchor="w", stretch=True)
        self.results_table.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(parent, orient="vertical", command=self.results_table.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.results_table.configure(yscrollcommand=scroll.set)

    def _bind_deferred_refresh(self, widget: tk.Widget) -> None:
        widget.bind("<Return>", lambda _e: self.refresh_plot())
        widget.bind("<Tab>", lambda _e: self.refresh_plot())
        widget.bind("<FocusOut>", lambda _e: self.refresh_plot())

    def _apply_operand_control_visibility(self, label: str) -> None:
        spec = self._merit_spec_for_label(label)
        if spec is None:
            return
        visible_controls = set(spec.controls)
        widget_groups = self.operand_control_widgets.get(label, {})
        for control_name, widgets in widget_groups.items():
            for widget in widgets:
                if control_name in visible_controls:
                    widget.grid()
                else:
                    widget.grid_remove()

    def _update_operand_setup_visibility(self) -> None:
        if not hasattr(self, "merit_mode_list"):
            return
        selected = {self.merit_mode_list.get(i) for i in self.merit_mode_list.curselection()}
        for label, frame in self.operand_setup_frames.items():
            visible = label in selected
            if visible:
                frame.grid()
            else:
                frame.grid_remove()

    def _set_initial_pane_layout(self) -> None:
        self.update_idletasks()
        total_width = self.main_pane.winfo_width()
        if total_width < 300:
            self.after(100, self._set_initial_pane_layout)
            return
        try:
            self.main_pane.sashpos(0, int(total_width * 0.80))
            self._initial_layout_passes += 1
        except Exception:
            self.after(100, self._set_initial_pane_layout)

    def _maybe_refresh_initial_pane_layout(self, _event=None) -> None:
        if self._initial_layout_passes >= 40:
            return
        self.after(100, self._set_initial_pane_layout)

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
        self.layout_names = sorted(
            self.layout_files,
            key=lambda name: (0 if name == DEFAULT_LAYOUT_TITLE else 1, name.lower()),
        )
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

    def open_3d_view(self) -> None:
        try:
            wavelength = self._current_wavelength()
            system = self.build_system()
            rays = Kos.raykeeper(system)
            max_radius = max((max(row.diameter / 2.0, 0.5) for row in self.rows), default=1.0)
            self._trace_preview_rays(system, rays, wavelength, max_radius)
            self.last_system = system
            self.last_rays = rays
            plotter = pv.Plotter(shape=(1, 1), title="KrakenOS 3D", notebook=False)
            plot3d(system, 0, plotter, 0.99)
            rayplot3d(rays, 0, plotter, 0.99, 0)
            plotter.add_axes(line_width=4)
            cx, cy, cz = plotter.center
            plotter.set_focus([cx, cy, cz])
            plotter.camera_position = [-1.0, 0.5, 1.0]
            plotter.set_viewup([0, 1.0, 0])
            plotter.enable_anti_aliasing()
            plotter.set_background("white", top="white")
            plotter.add_text("KrakenOS", position="upper_left", font_size=20, color="royalblue")
            plotter.show_grid(font_size=6, color="black")
            plotter.show(auto_close=True, interactive=True, interactive_update=False)
            self.status_var.set("Opened Kraken 3D view")
            self.append_debug("Opened Kraken 3D view (close with window X, Alt+F4, or q)")
        except Exception as exc:
            self.append_debug(f"3D view error: {exc}")
            self.status_var.set(f"3D view failed: {exc}")

    def load_layout_by_name(self, name: str) -> None:
        path = self.layout_files.get(name)
        if path is None:
            return
        self.current_layout_file = path
        try:
            info = _load_python_data(path)
            loaded_rows = [
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
                    tilt_x=float(item.get("tilt_x", 0.0)),
                    tilt_y=float(item.get("tilt_y", 0.0)),
                    tilt_z=float(item.get("tilt_z", 0.0)),
                    desp_x=float(item.get("desp_x", 0.0)),
                    desp_y=float(item.get("desp_y", 0.0)),
                    desp_z=float(item.get("desp_z", 0.0)),
                    axis_move=float(item.get("axis_move", 0.0)),
                    glass=str(item.get("glass", "AIR")),
                )
                for item in info["surfaces"]
            ]
        except Exception:
            surfaces = self._extract_surfaces_from_example(path)
            loaded_rows = [self._row_from_surface(surface, index, len(surfaces)) for index, surface in enumerate(surfaces)]

        loaded_rows = self._normalized_rows_copy(loaded_rows)
        if self.rows:
            self.rows = self._append_layout_rows(self.rows, loaded_rows)
        else:
            self.rows = loaded_rows
            self._apply_initial_field_defaults()
            self._apply_initial_layout_view_defaults(name)

        self._normalize_special_rows()
        self._sync_table()
        self.refresh_plot()
        self.layout_var.set(name)
        self.example_var.set("Examples")
        self.status_var.set(f"Appended {name}")

    def load_example_by_name(self, name: str) -> None:
        path = self.example_files.get(name)
        if path is None:
            return
        try:
            surfaces = self._extract_surfaces_from_example(path)
        except Exception as exc:
            self.status_var.set(f"Failed to load example {name}: {exc}")
            return
        self.current_layout_file = None
        self.rows = [self._row_from_surface(surface, index, len(surfaces)) for index, surface in enumerate(surfaces)]
        self._normalize_special_rows()
        self._apply_example_display_defaults(path)
        self._sync_table()
        self.refresh_plot()
        self.layout_var.set("Common Optical Layout")
        self.example_var.set(name)
        self.status_var.set(f"Loaded example {name}")

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
                f"{row.tilt_x:g}",
                f"{row.tilt_y:g}",
                f"{row.tilt_z:g}",
                f"{row.desp_x:g}",
                f"{row.desp_y:g}",
                f"{row.desp_z:g}",
                f"{row.axis_move:g}",
                row.glass,
            ]
            tags = ("optimize",) if self._row_has_optimization(row) else ()
            self.table.insert("", "end", values=values, tags=tags)
        self._refresh_analysis_surface_choices()
        self._refresh_operand_surface_choices()
        self._sync_object_controls()

    def _refresh_analysis_surface_choices(self) -> None:
        options = ["Auto"]
        for index, row in enumerate(self.rows):
            options.append(f"{index}: {row.name}")
        current = self.analysis_surface_var.get()
        self.analysis_surface_menu["values"] = options
        if current not in options:
            self.analysis_surface_var.set("Auto")
        self._schedule_table_grid_update()
        self.after_idle(self._update_active_cell_border)

    @staticmethod
    def _parse_numeric_display(value: str) -> float:
        return float(value.replace("*", "").strip())

    @staticmethod
    def _normalized_rows_copy(rows: list[SurfaceRow]) -> list[SurfaceRow]:
        copied = [SurfaceRow(**asdict(row)) for row in rows]
        if copied:
            copied[0].surface = "Object"
            if not copied[0].name or copied[0].name == "Surface":
                copied[0].name = "Object"
            copied[-1].surface = "Image"
            if not copied[-1].name or copied[-1].name == "Surface":
                copied[-1].name = "Image"
        return copied

    @staticmethod
    def _append_layout_rows(existing_rows: list[SurfaceRow], layout_rows: list[SurfaceRow]) -> list[SurfaceRow]:
        base = [SurfaceRow(**asdict(row)) for row in existing_rows]
        additions = [SurfaceRow(**asdict(row)) for row in layout_rows[1:-1]]
        if not additions:
            return base
        insert_at = len(base)
        if base and base[-1].surface == "Image":
            insert_at -= 1
        for offset, row in enumerate(additions):
            base.insert(insert_at + offset, row)
        return base

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
                    tilt_x=float(values[6]),
                    tilt_y=float(values[7]),
                    tilt_z=float(values[8]),
                    desp_x=float(values[9]),
                    desp_y=float(values[10]),
                    desp_z=float(values[11]),
                    axis_move=float(values[12]),
                    glass=str(values[13]),
                )
            )
        self.rows = rows
        self._sync_object_controls()

    def _on_table_click(self, event: tk.Event) -> None:
        row_id = self.table.identify_row(event.y)
        column_id = self.table.identify_column(event.x)
        self.table.focus_set()
        if not row_id or not column_id:
            self._active_cell = None
            self._hide_active_cell_border()
            return
        self._active_cell = (row_id, column_id)
        self.table.focus(row_id)
        self.table.selection_set(row_id)
        self.after_idle(self._update_active_cell_border)

    def _move_active_cell(self, event: tk.Event) -> str:
        self.table.focus_set()
        children = list(self.table.get_children())
        if not children:
            return "break"
        if self._active_cell is None:
            row_id = children[0]
            column_id = "#2"
        else:
            row_id, column_id = self._active_cell
            if row_id not in children:
                row_id = children[0]
            column_index = int(column_id.replace("#", ""))
            row_index = children.index(row_id)
            if event.keysym == "Left":
                column_index = max(2, column_index - 1)
            elif event.keysym == "Right":
                column_index = min(len(FIELDS), column_index + 1)
            elif event.keysym == "Up":
                row_index = max(0, row_index - 1)
            elif event.keysym == "Down":
                row_index = min(len(children) - 1, row_index + 1)
            row_id = children[row_index]
            column_id = f"#{column_index}"
        self._active_cell = (row_id, column_id)
        self.table.focus(row_id)
        self.table.selection_set(row_id)
        self._ensure_active_cell_visible(row_id, column_id)
        self.after_idle(self._update_active_cell_border)
        return "break"

    def _ensure_active_cell_visible(self, row_id: str, column_id: str) -> None:
        self.table.see(row_id)
        self.update_idletasks()
        columns = list(self.table["columns"])
        target_bbox = self.table.bbox(row_id, column_id)
        if target_bbox:
            x, _y, width, _height = target_bbox
            visible_width = max(self.table.winfo_width(), 1)
            if x >= 0 and (x + width) <= visible_width:
                return

        total_width = 0
        target_left = 0
        target_width = 0
        target_field = FIELDS[int(column_id.replace("#", "")) - 1]
        for field in columns:
            width = int(self.table.column(field, "width"))
            if field == target_field:
                target_left = total_width
                target_width = width
            total_width += width
        if total_width <= 0:
            return
        visible_width = max(self.table.winfo_width(), 1)
        view_left, _view_right = self.table.xview()
        visible_left = view_left * total_width
        visible_right = visible_left + visible_width
        if target_left < visible_left:
            desired_left = max(0.0, target_left - 16.0)
            self.table.xview_moveto(desired_left / total_width)
        elif target_left + target_width > visible_right:
            desired_left = max(0.0, target_left + target_width - visible_width + 16.0)
            self.table.xview_moveto(min(1.0, desired_left / total_width))
        self.update_idletasks()

    def _hide_active_cell_border(self) -> None:
        for part in self._cell_border_parts:
            part.place_forget()

    def _update_active_cell_border(self, _event: tk.Event | None = None) -> None:
        if self._active_cell is None:
            self._hide_active_cell_border()
            return
        row_id, column_id = self._active_cell
        try:
            x, y, width, height = self.table.bbox(row_id, column_id)
        except tk.TclError:
            self._hide_active_cell_border()
            return
        if width <= 0 or height <= 0:
            self._hide_active_cell_border()
            return
        top, bottom, left, right = self._cell_border_parts
        top.place(x=x, y=y, width=width, height=2)
        bottom.place(x=x, y=y + height - 2, width=width, height=2)
        left.place(x=x, y=y, width=2, height=height)
        right.place(x=x + width - 2, y=y, width=2, height=height)

    def _on_table_scroll(self, scrollbar: ttk.Scrollbar, first: str, last: str) -> None:
        scrollbar.set(first, last)
        self._schedule_table_grid_update()
        self.after_idle(self._update_active_cell_border)

    def _clear_table_grid(self) -> None:
        for part in self._grid_overlays:
            part.destroy()
        self._grid_overlays.clear()

    def _schedule_table_grid_update(self, _event: tk.Event | None = None) -> None:
        if self._grid_after_id is not None:
            self.after_cancel(self._grid_after_id)
        self._grid_after_id = self.after(30, self._update_table_grid)

    def _update_table_grid(self, _event: tk.Event | None = None) -> None:
        self._grid_after_id = None
        self._clear_table_grid()
        columns = list(self.table["columns"])
        items = self.table.get_children()
        if not columns or not items:
            return
        grid_color = "#e2e7ef"
        visible_bboxes = []
        for item in items:
            bbox = self.table.bbox(item, "#1")
            if bbox:
                visible_bboxes.append((item, bbox))
        if not visible_bboxes:
            return
        data_top = min(bbox[1] for _, bbox in visible_bboxes)
        data_bottom = max(bbox[1] + bbox[3] for _, bbox in visible_bboxes)
        data_height = max(0, data_bottom - data_top)
        if data_height <= 0:
            return

        first_item = visible_bboxes[0][0]
        for column_index in range(1, len(columns)):
            bbox = self.table.bbox(first_item, f"#{column_index}")
            if not bbox:
                continue
            x, _y, width, _height = bbox
            separator = tk.Frame(self.table, bg=grid_color, width=1)
            separator.place(x=x + width - 1, y=data_top, width=1, height=data_height)
            self._grid_overlays.append(separator)

        for item, bbox in visible_bboxes:
            _x, y, width, height = bbox
            row_line = tk.Frame(self.table, bg=grid_color, height=1)
            row_line.place(x=0, y=y + height - 1, relwidth=1.0, height=1)
            self._grid_overlays.append(row_line)

        self.after_idle(self._update_active_cell_border)

    def _refresh_operand_surface_choices(self) -> None:
        values = ["Auto"]
        for index, row in enumerate(self.rows):
            if row.surface in {"Object", "Image"}:
                continue
            values.append(f"{index}: {row.name}")
        for label, var in self.operand_surface_vars.items():
            current = var.get().strip() if var.get() else "Auto"
            if current not in values:
                var.set("Auto")
        for widget in self.winfo_children():
            self._apply_surface_values_to_descendants(widget, values)

    def _apply_surface_values_to_descendants(self, widget, values) -> None:
        if isinstance(widget, ttk.Combobox):
            textvar = widget.cget("textvariable")
            for var in self.operand_surface_vars.values():
                if str(var) == textvar:
                    widget["values"] = values
                    break
        for child in widget.winfo_children():
            self._apply_surface_values_to_descendants(child, values)

    def _sync_object_controls(self) -> None:
        if not hasattr(self, "field_summary_var"):
            return
        self._sync_field_mode_ui()
        metrics = self._field_metrics()
        self.field_summary_var.set(
            "Angle: {angle:.3g} deg\nObject: {obj:.3g} mm\nParaxial image: {parax:.3g} mm\nReal image: {real:.3g} mm".format(
                angle=metrics["angle_deg"],
                obj=metrics["object_height"],
                parax=metrics["paraxial_image_height"],
                real=metrics["real_image_height"],
            )
        )
        warning = ""
        if self.rows and self._current_object_mode() == "Finite":
            object_radius = max(float(self.rows[0].diameter) / 2.0, 0.0)
            if abs(metrics["object_height"]) > object_radius + 1e-9:
                warning = f"Field exceeds object radius ({object_radius:.3g} mm)."
        self.field_warning_var.set(warning)

    def _on_object_mode_changed(self, _event=None) -> None:
        self._sync_field_default_from_current_type()
        self._sync_field_mode_ui()
        self.refresh_plot()

    def _apply_initial_field_defaults(self) -> None:
        if self._field_defaults_initialized or not hasattr(self, "field_type_var"):
            return
        if self._current_object_mode() == "Infinity":
            self.field_type_var.set("Angle")
            self._last_field_type = "Angle"
            self._field_type_defaults["Angle"] = "0.0"
            self.field_value_var.set("0.0")
        else:
            self.field_type_var.set("Object Height")
            self._last_field_type = "Object Height"
            self._field_type_defaults["Object Height"] = "0.0"
            self.field_value_var.set("0.0")
        self._field_defaults_initialized = True
        self._sync_field_mode_ui()

    def _apply_initial_layout_view_defaults(self, name: str) -> None:
        if not hasattr(self, "display_orientation_var"):
            return
        if name == DEFAULT_LAYOUT_TITLE:
            self.display_orientation_var.set("Horizontal")
            self.object_mode_var.set("Finite")
            self.field_type_var.set("Object Height")
            self._last_field_type = "Object Height"
            self._field_type_defaults["Object Height"] = "0.0"
            self.field_value_var.set("0.0")
            self._sync_field_mode_ui()

    def _on_field_type_changed(self, _event=None) -> None:
        self._sync_field_default_from_current_type()
        self._sync_field_mode_ui()
        self.refresh_plot()

    def _sync_field_default_from_current_type(self) -> None:
        if not hasattr(self, "field_value_var"):
            return
        previous_type = getattr(self, "_last_field_type", self._current_field_type())
        field_type = self._current_field_type()
        current_text = self.field_value_var.get().strip()
        if current_text:
            self._field_type_defaults[previous_type] = current_text
        default_text = self._field_type_defaults.get(field_type, "0.0")
        self._last_field_type = field_type
        if current_text != default_text:
            self.field_value_var.set(default_text)

    def _sync_field_mode_ui(self) -> None:
        if not hasattr(self, "field_type_menu"):
            return
        if self._current_object_mode() == "Infinity":
            values = [
                "Angle",
                "Paraxial Image Height",
                "Real Image Height",
                "Object Height",
            ]
            note = "Preferred: Angle for infinity object. Image-height modes are derived targets."
        else:
            values = [
                "Object Height",
                "Paraxial Image Height",
                "Real Image Height",
                "Angle",
            ]
            note = "Preferred: Object Height for finite object. Angle remains available as a derived field."
        self.field_type_menu["values"] = values
        if hasattr(self, "field_mode_note_var"):
            self.field_mode_note_var.set(note)
        if hasattr(self, "field_value_label_var"):
            label_map = {
                "Angle": "Angle [deg]",
                "Object Height": "Object Height [mm]",
                "Paraxial Image Height": "Paraxial Image Height [mm]",
                "Real Image Height": "Real Image Height [mm]",
            }
            self.field_value_label_var.set(label_map.get(self._current_field_type(), "Field value"))

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
        spec = self._variable_spec_for_field(field)
        if spec is None:
            return
        row_index = self.table.index(row_id)
        row = self.rows[row_index]
        if row.surface in {"Object", "Image"} or not spec.is_supported(row):
            return
        if self.popup_menu is not None:
            self.popup_menu.destroy()
        self.current_menu_row_id = row_id
        self.current_menu_field = field
        marked = spec.is_enabled(row)
        bounds = spec.get_bounds(row)
        menu = tk.Menu(self, tearoff=0)
        menu.add_command(
            label=f"{'Unselect' if marked else 'Select'} {spec.label} for optimization",
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
        if field == "surface":
            index = self.table.index(row_id)
            row = self.rows[index]
            if value == "Mirror":
                row.glass = "MIRROR"
                if abs(row.tilt_x) < 1e-9 and abs(row.tilt_y) < 1e-9 and abs(row.tilt_z) < 1e-9:
                    row.tilt_x = 45.0
            elif row.surface not in {"Object", "Image"} and row.glass == "MIRROR":
                row.glass = "AIR"
        self._normalize_special_rows()
        self._sync_table()
        self.refresh_plot()
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None

    @staticmethod
    def _row_has_optimization(row: SurfaceRow) -> bool:
        return row.optimize_rc or row.optimize_thickness

    def toggle_current_optimization_cell(self) -> None:
        if self.current_menu_row_id is None or self.current_menu_field is None:
            return
        index = self.table.index(self.current_menu_row_id)
        row = self.rows[index]
        spec = self._variable_spec_for_field(self.current_menu_field)
        if spec is None:
            return
        spec.set_enabled(row, not spec.is_enabled(row))
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
        spec = self._variable_spec_for_field(self.current_menu_field)
        if spec is None:
            return
        current = spec.get_bounds(row)
        current_value = spec.value_from_row(row)
        default_bounds = current or spec.default_bounds(current_value)

        dialog = tk.Toplevel(self)
        dialog.title(f"Bounds for {row.name} {spec.label}")
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
            spec.set_bounds(row, (lower, upper))
            self.append_progress(
                f"Bounds set for row {index} {spec.label}: [{lower:g}, {upper:g}]"
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
        spec = self._variable_spec_for_field(self.current_menu_field)
        if spec is None:
            return
        spec.set_bounds(row, None)
        self.append_progress(f"Bounds cleared for row {index} {spec.label}.")
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
        clear_aperture = max(
            [max(row.diameter, 1.0) for row in self.rows if row.surface not in {"Object", "Image"}] or [100.0]
        ) * 4.0
        for row in self.rows:
            surface = Kos.surf()
            display_name = row.name if row.surface in {"Object", "Image"} else ""
            surface.Name = display_name.replace(" ", "\n")
            surface.Rc = row.rc
            surface.Thickness = row.thickness
            surface.Diameter = clear_aperture if row.surface in {"Object", "Image"} else row.diameter
            surface.Glass = row.glass
            surface.TiltX = row.tilt_x
            surface.TiltY = row.tilt_y
            surface.TiltZ = row.tilt_z
            surface.DespX = row.desp_x
            surface.DespY = row.desp_y
            surface.DespZ = row.desp_z
            surface.AxisMove = row.axis_move
            surface.Nm_Pos = self._name_offset(row)
            surface.Drawing = 0.0 if row.surface in {"Object", "Image", "Mirror"} else 1.0
            if row.surface == "Mirror":
                surface.Glass = "MIRROR"
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

    def _current_mtf_frequency(self) -> float:
        var = self.operand_frequency_vars.get("MTF @ freq")
        if var is None:
            return 50.0
        try:
            value = float(var.get())
        except ValueError:
            return 50.0
        return max(0.0, value)

    def _operand_mtf_mode(self, label: str) -> str:
        var = self.operand_mtf_mode_vars.get(label)
        if var is None:
            return "average"
        value = var.get().strip().lower()
        if value in {"tangential", "sagittal", "average"}:
            return value
        return "average"

    def _mtf_analysis_settings(self) -> dict[str, float | int | str]:
        label = "MTF @ freq"
        return {
            "wavelength": self._operand_wavelength(label),
            "surface_index": self._operand_surface_index(label),
            "aperture_type": self._operand_aperture_type(label),
            "aperture_value": self._operand_aperture_value(label),
            "field_type": self._operand_field_type(label),
            "field_x": self._operand_field_x(label),
            "field_y": self._operand_field_y(label),
        }

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
        if self.analysis_mode in {"none", "native_off_axis"}:
            self.ax = self.figure.add_subplot(111)
            analysis_ax = None
        else:
            gs = self.figure.add_gridspec(1, 2, width_ratios=[4.2, 1.35], wspace=0.22)
            self.ax = self.figure.add_subplot(gs[0])
            analysis_ax = self.figure.add_subplot(gs[1])

        if self.analysis_mode == "native_off_axis" and self._has_off_axis_geometry():
            self._plot_native_off_axis_preview(
                analysis_ax,
                max_radius,
                use_native_surfaces=True,
            )
            self.ax.grid(True, alpha=0.2)
            self.ax.set_xlabel("Native X [mm]" if self._current_display_orientation() == "Vertical" else "Native Fold X [mm]")
            self.ax.set_ylabel("Native Y [mm]" if self._current_display_orientation() == "Vertical" else "Native Fold Y [mm]")
            self.ax.set_title("")
            self.figure.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.92, wspace=0.28)
            self.figure.text(0.5, 0.035, "KrakenOS Layout", ha="center", va="center")
            self._sync_object_controls()
            self.canvas.draw_idle()
            return

        if self._is_folded_mirror_preview_mode():
            self._plot_folded_mirror_preview(analysis_ax)
            self.ax.grid(True, alpha=0.2)
            self.ax.set_xlabel("Fold X [mm]")
            self.ax.set_ylabel("Fold Y [mm]")
            self.ax.set_title("")
            self.figure.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.92, wspace=0.28)
            self.figure.text(0.5, 0.035, "KrakenOS Layout", ha="center", va="center")
            self._sync_object_controls()
            self.canvas.draw_idle()
            return

        try:
            wavelength = self._current_wavelength()
            capture = io.StringIO()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                with redirect_stdout(capture), redirect_stderr(capture):
                    system = self.build_system()
                    if getattr(system.Pr3D, "ExistSolid", 0) == 0:
                        original_build = system.BUILD
                        system.BUILD = 1
                        system.build()
                        system.BUILD = original_build
                    rays = Kos.raykeeper(system)
                    self._trace_preview_rays(system, rays, wavelength, max_radius)
            self.append_debug(capture.getvalue())
            self.last_system = system
            self.last_rays = rays
            if self._has_native_drawn_surfaces():
                Plot2DSurf(system, 0, self.ax)
            self._draw_custom_mirror_surfaces()
            self._apply_display_orientation_to_lines()
            surf_line_count = len(self.ax.lines)
            self._style_embedded_plot(surf_line_count)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                self._draw_colored_rays(rays)
            self._apply_display_orientation_to_lines(surf_line_count)
            if self._has_off_axis_geometry():
                self._set_plot_limits_from_drawn_data()
                self.ax.set_aspect("equal", adjustable="box")
            elif self._current_display_orientation() == "Horizontal":
                self._set_plot_limits_from_drawn_data()
                self.ax.set_aspect("equal", adjustable="box")
            else:
                self._set_plot_limits_from_layout(max_radius)
            self._draw_input_ray_overlay(max_radius)
            self._draw_reference_plane_labels()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                optics_info = self._collect_optics_info(system, rays, wavelength)
            self._draw_optics_markers(optics_info)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                self._plot_analysis(analysis_ax, system, rays, wavelength)
                self._update_results(system, rays, wavelength, optics_info)
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
        if self._current_display_orientation() == "Horizontal":
            self.ax.set_xlabel("-Y [mm]")
            self.ax.set_ylabel("-Z [mm]")
        else:
            self.ax.set_xlabel("Z [mm]")
            self.ax.set_ylabel("Y [mm]")
        self.ax.set_title("")
        self.figure.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.92, wspace=0.28)
        self.figure.text(0.5, 0.035, "KrakenOS Layout", ha="center", va="center")
        self._sync_object_controls()
        self.canvas.draw_idle()
        if self._initial_layout_passes < 40:
            self.after(50, self._set_initial_pane_layout)

    def _plot_analysis(self, analysis_ax, system, rays, wavelength: float) -> None:
        if analysis_ax is None:
            return
        analysis_ax.clear()
        try:
            analysis_rays = self._build_analysis_rays(system, wavelength)
            X, Y, Z, L, M, N = self._pick_image_plane_data(analysis_rays)
        except Exception:
            X = Y = Z = L = M = N = np.asarray([])

        if X.size == 0 and self.analysis_mode in {"spot", "rms"}:
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
            if np.ptp(X) < 1e-12 and np.ptp(Y) < 1e-12:
                analysis_ax.set_xlim(float(X[0]) - 1.0, float(X[0]) + 1.0)
                analysis_ax.set_ylim(float(Y[0]) - 1.0, float(Y[0]) + 1.0)
            analysis_ax.grid(True, alpha=0.2)
            return

        if self.analysis_mode == "psf":
            try:
                self._begin_analysis_progress("PSF analysis")
                field_type = "angle" if self._current_object_mode() == "Infinity" else "height"
                field_y = self._current_field_angle_deg() if field_type == "angle" else self._current_field_height()
                self._update_analysis_progress("Tracing rays", 1, 3)
                psf_rays = self._build_analysis_rays(
                    system,
                    wavelength,
                    sample_count=max(48, self._current_ray_count() * 10),
                    pattern="hexapolar",
                    surface_index=self._analysis_surface_index(),
                    aperture_type=self._current_aperture_type(),
                    aperture_value=self._current_aperture_value(),
                    field_type=field_type,
                    field_x=0.0,
                    field_y=field_y,
                )
                self._update_analysis_progress("Building PSF image", 2, 3)
                x_local, y_local, _z_local, _l_local, _m_local, _n_local = self._pick_image_plane_data(psf_rays)
                x_local = np.asarray(x_local, dtype=float)
                y_local = np.asarray(y_local, dtype=float)
                finite = np.isfinite(x_local) & np.isfinite(y_local)
                x_local = x_local[finite]
                y_local = y_local[finite]
                if x_local.size < 4:
                    raise RuntimeError("Not enough image-plane samples for PSF")
                span_x = max(float(np.ptp(x_local)), 1e-3)
                span_y = max(float(np.ptp(y_local)), 1e-3)
                span = max(span_x, span_y) * 1.25
                bins = 128
                hist, xedges, yedges = np.histogram2d(
                    x_local,
                    y_local,
                    bins=bins,
                    range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
                )
                psf = hist.T
                psf /= max(float(np.max(psf)), 1e-12)
                extent = [float(xedges[0]), float(xedges[-1]), float(yedges[0]), float(yedges[-1])]
                image = analysis_ax.imshow(psf, origin="lower", extent=extent, cmap="inferno", aspect="equal")
                analysis_ax.set_title(f"Geometric PSF  |  {field_type}={field_y:.3g}  |  {wavelength:.4g} um")
                analysis_ax.set_xlabel("X [mm]")
                analysis_ax.set_ylabel("Y [mm]")
                analysis_ax.grid(False)
                self.figure.colorbar(image, ax=analysis_ax, fraction=0.046, pad=0.04, label="Normalized intensity")
                self.append_debug(f"PSF analysis ok: rays={x_local.size}, bins={bins}")
                self._update_analysis_progress("Rendering", 3, 3)
                self._finish_analysis_progress("PSF analysis", success=True)
            except Exception as exc:
                self.append_debug(f"PSF analysis error: {exc}")
                analysis_ax.text(0.5, 0.5, "PSF analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("PSF analysis", success=False)
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
                self._begin_analysis_progress("Pupil analysis")
                self._update_analysis_progress("Building pupil", 1, 2)
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
                self._update_analysis_progress("Rendering", 2, 2)
                self._finish_analysis_progress("Pupil analysis", success=True)
            except Exception:
                analysis_ax.text(0.5, 0.5, "Pupil analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("Pupil analysis", success=False)
            return

        if self.analysis_mode == "seidel":
            try:
                self._begin_analysis_progress("Seidel analysis")
                self._update_analysis_progress("Building pupil", 1, 3)
                pupil = Kos.PupilCalc(
                    system,
                    self._analysis_surface_index(),
                    wavelength,
                    self._current_aperture_type(),
                    self._current_aperture_value(),
                )
                self._update_analysis_progress("Computing coefficients", 2, 3)
                seidel = Kos.Seidel(pupil)
                values = np.asarray(seidel.SCW_TOTAL, dtype=float)
                labels = seidel.SCW_NM
                analysis_ax.bar(np.arange(len(values)), values, color="#8e44ad")
                analysis_ax.set_xticks(np.arange(len(values)), labels, rotation=25, ha="right")
                analysis_ax.set_title("Seidel Coefficients in Waves")
                analysis_ax.set_ylabel("Waves")
                analysis_ax.grid(True, axis="y", alpha=0.2)
                self._update_analysis_progress("Rendering", 3, 3)
                self._finish_analysis_progress("Seidel analysis", success=True)
            except Exception:
                analysis_ax.text(0.5, 0.5, "Seidel analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("Seidel analysis", success=False)
            return

        if self.analysis_mode == "wavefront":
            try:
                self._begin_analysis_progress("Wavefront analysis")
                self._update_analysis_progress("Building pupil", 1, 3)
                pupil = Kos.PupilCalc(
                    system,
                    self._analysis_surface_index(),
                    wavelength,
                    self._current_aperture_type(),
                    self._current_aperture_value(),
                )
                pupil.Samp = max(6, min(16, int(np.sqrt(max(1, self._current_ray_count())) * 3)))
                self._update_analysis_progress("Computing phase", 2, 3)
                px, py, phase, p2v = Kos.Phase(pupil)
                scatter = analysis_ax.scatter(py, px, c=phase, cmap="RdBu_r", s=20)
                analysis_ax.set_title(f"Wavefront  |  P2V = {float(p2v):.4g}")
                analysis_ax.set_xlabel("X pupil")
                analysis_ax.set_ylabel("Y pupil")
                analysis_ax.set_aspect("equal", adjustable="box")
                analysis_ax.grid(True, alpha=0.2)
                self.figure.colorbar(scatter, ax=analysis_ax, fraction=0.046, pad=0.04, label="Waves")
                self._update_analysis_progress("Rendering", 3, 3)
                self._finish_analysis_progress("Wavefront analysis", success=True)
            except Exception:
                analysis_ax.text(0.5, 0.5, "Wavefront analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("Wavefront analysis", success=False)
            return

        if self.analysis_mode == "field_curvature":
            try:
                self._begin_analysis_progress("Field curvature / distortion")
                field_type = "angle" if self._current_object_mode() == "Infinity" else "height"
                field_limit = self._current_field_angle_deg() if field_type == "angle" else self._current_field_height()
                if field_limit <= 1e-9:
                    field_limit = 5.0 if field_type == "angle" else max(self._current_field_height(), 0.5)
                field_sample_count = max(11, self._current_field_count())
                field_samples = list(np.linspace(-field_limit, field_limit, field_sample_count))
                self.append_debug(
                    f"Field curvature/distortion sampling: type={field_type}, limit={field_limit:.4g}, count={field_sample_count}"
                )
                sample_count = max(18, self._current_ray_count() * 3)
                axis_results: dict[str, dict[str, np.ndarray]] = {}
                axis_defs = {
                    "X": {"color_fc": "#1f77b4", "color_dist": "#6baed6", "marker_fc": "o", "marker_dist": "s"},
                    "Y": {"color_fc": "#d62728", "color_dist": "#ff9896", "marker_fc": "^", "marker_dist": "D"},
                }
                total_steps = len(field_samples) * 2
                completed_steps = 0

                for axis_name in ("X", "Y"):
                    measured_fields: list[float] = []
                    image_heights: list[float] = []
                    focus_shifts: list[float] = []
                    for field_value in field_samples:
                        completed_steps += 1
                        self._update_analysis_progress(
                            f"Sampling {axis_name}-field",
                            completed_steps,
                            total_steps,
                        )
                        field_x = field_value if axis_name == "X" else 0.0
                        field_y = field_value if axis_name == "Y" else 0.0
                        sample_rays = self._build_analysis_rays(
                            system,
                            wavelength,
                            sample_count=sample_count,
                            pattern="hexapolar",
                            surface_index=self._analysis_surface_index(),
                            aperture_type=self._current_aperture_type(),
                            aperture_value=self._current_aperture_value(),
                            field_type=field_type,
                            field_x=field_x,
                            field_y=field_y,
                        )
                        x_local, y_local, _z_local, l_local, m_local, n_local = self._pick_image_plane_data(sample_rays)
                        x_local = np.asarray(x_local, dtype=float)
                        y_local = np.asarray(y_local, dtype=float)
                        l_local = np.asarray(l_local, dtype=float)
                        m_local = np.asarray(m_local, dtype=float)
                        n_local = np.asarray(n_local, dtype=float)
                        finite = (
                            np.isfinite(x_local)
                            & np.isfinite(y_local)
                            & np.isfinite(l_local)
                            & np.isfinite(m_local)
                            & np.isfinite(n_local)
                        )
                        x_local = x_local[finite]
                        y_local = y_local[finite]
                        l_local = l_local[finite]
                        m_local = m_local[finite]
                        n_local = n_local[finite]
                        if x_local.size < 4:
                            continue

                        slopes_x = l_local / np.where(np.abs(n_local) < 1e-9, np.sign(n_local) * 1e-9 + 1e-9, n_local)
                        slopes_y = m_local / np.where(np.abs(n_local) < 1e-9, np.sign(n_local) * 1e-9 + 1e-9, n_local)

                        if axis_name == "X":
                            image_height = float(np.mean(x_local))
                            denom = float(np.sum(slopes_x**2))
                            focus_shift = 0.0 if denom <= 1e-12 else -float(np.sum(x_local * slopes_x) / denom)
                        else:
                            image_height = float(np.mean(y_local))
                            denom = float(np.sum(slopes_y**2))
                            focus_shift = 0.0 if denom <= 1e-12 else -float(np.sum(y_local * slopes_y) / denom)

                        measured_fields.append(field_value)
                        image_heights.append(image_height)
                        focus_shifts.append(focus_shift)

                    if len(measured_fields) < 2:
                        continue

                    fields = np.asarray(measured_fields, dtype=float)
                    heights = np.asarray(image_heights, dtype=float)
                    focus = np.asarray(focus_shifts, dtype=float)
                    abs_fields = np.abs(fields)
                    abs_heights = np.abs(heights)
                    mask = abs_fields > 1e-9
                    if np.count_nonzero(mask) >= 2:
                        slope = float(np.dot(abs_fields[mask], abs_heights[mask]) / max(np.dot(abs_fields[mask], abs_fields[mask]), 1e-12))
                        ideal = slope * abs_fields
                        distortion = np.zeros_like(abs_heights)
                        valid = ideal > 1e-12
                        distortion[valid] = (abs_heights[valid] - ideal[valid]) / ideal[valid] * 100.0
                    else:
                        distortion = np.zeros_like(abs_heights)

                    axis_results[axis_name] = {
                        "fields": abs_fields,
                        "focus": focus,
                        "distortion": distortion,
                    }

                if not axis_results:
                    raise RuntimeError("Not enough field samples for field-curvature/distortion")

                ax2 = analysis_ax.twinx()
                legend_lines = []
                legend_labels = []
                for axis_name, data in axis_results.items():
                    style = axis_defs[axis_name]
                    line_fc, = analysis_ax.plot(
                        data["fields"],
                        data["focus"],
                        color=style["color_fc"],
                        marker=style["marker_fc"],
                        linewidth=0.0,
                        markersize=4.5,
                        label=f"{axis_name} focus",
                    )
                    line_dist, = ax2.plot(
                        data["fields"],
                        data["distortion"],
                        color=style["color_dist"],
                        marker=style["marker_dist"],
                        linewidth=0.0,
                        markersize=4.0,
                        label=f"{axis_name} distortion",
                    )
                    if len(data["fields"]) >= 4:
                        smooth_fields = np.linspace(float(np.min(data["fields"])), float(np.max(data["fields"])), 200)
                        focus_degree = min(3, len(data["fields"]) - 1)
                        dist_degree = min(3, len(data["distortion"]) - 1)
                        try:
                            focus_poly = np.poly1d(np.polyfit(data["fields"], data["focus"], deg=focus_degree))
                            dist_poly = np.poly1d(np.polyfit(data["fields"], data["distortion"], deg=dist_degree))
                            analysis_ax.plot(
                                smooth_fields,
                                focus_poly(smooth_fields),
                                color=style["color_fc"],
                                linewidth=2.0,
                                alpha=0.9,
                            )
                            ax2.plot(
                                smooth_fields,
                                dist_poly(smooth_fields),
                                color=style["color_dist"],
                                linewidth=1.8,
                                linestyle="--",
                                alpha=0.9,
                            )
                        except Exception:
                            analysis_ax.plot(
                                data["fields"],
                                data["focus"],
                                color=style["color_fc"],
                                linewidth=2.0,
                                alpha=0.9,
                            )
                            ax2.plot(
                                data["fields"],
                                data["distortion"],
                                color=style["color_dist"],
                                linewidth=1.8,
                                linestyle="--",
                                alpha=0.9,
                            )
                    else:
                        analysis_ax.plot(
                            data["fields"],
                            data["focus"],
                            color=style["color_fc"],
                            linewidth=2.0,
                            alpha=0.9,
                        )
                        ax2.plot(
                            data["fields"],
                            data["distortion"],
                            color=style["color_dist"],
                            linewidth=1.8,
                            linestyle="--",
                            alpha=0.9,
                        )
                    legend_lines.extend([line_fc, line_dist])
                    legend_labels.extend([f"{axis_name} focus", f"{axis_name} distortion"])

                analysis_ax.set_title(f"Field Curvature / Distortion  |  {field_type}")
                analysis_ax.set_xlabel("Field")
                analysis_ax.set_ylabel("Best-focus shift [mm]", color="#2c3e50")
                ax2.set_ylabel("Distortion [%]", color="#2c3e50")
                analysis_ax.grid(True, alpha=0.2)
                analysis_ax.legend(legend_lines, legend_labels, loc="best", fontsize=8)
                self.append_debug(
                    "Field curvature/distortion ok: "
                    + ", ".join(f"{axis}={len(data['fields'])}" for axis, data in axis_results.items())
                )
                self._finish_analysis_progress("Field curvature / distortion", success=True)
            except Exception as exc:
                self.append_debug(f"Field curvature/distortion error: {exc}")
                analysis_ax.text(0.5, 0.5, "Field curvature/distortion unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("Field curvature / distortion", success=False)
            return

        if self.analysis_mode == "mtf":
            try:
                self._begin_analysis_progress("MTF analysis")
                mtf_settings = self._mtf_analysis_settings()
                wavelength = float(mtf_settings["wavelength"])
                self._update_analysis_progress("Preparing pupil", 1, 4)
                pupil = Kos.PupilCalc(
                    system,
                    int(mtf_settings["surface_index"]),
                    wavelength,
                    str(mtf_settings["aperture_type"]),
                    float(mtf_settings["aperture_value"]),
                )
                pupil.FieldType = str(mtf_settings["field_type"])
                pupil.FieldX = float(mtf_settings["field_x"])
                pupil.FieldY = float(mtf_settings["field_y"])
                px, py, phase, _p2v = Kos.Phase(pupil)
                px = np.asarray(px, dtype=float)
                py = np.asarray(py, dtype=float)
                phase = np.asarray(phase, dtype=float)
                finite = np.isfinite(px) & np.isfinite(py) & np.isfinite(phase)
                px = px[finite]
                py = py[finite]
                phase = phase[finite]
                if px.size < 6:
                    raise RuntimeError("Not enough finite pupil samples for MTF fitting")

                zcoef = None
                used_terms = None
                last_error = None
                self._update_analysis_progress("Fitting Zernike terms", 2, 4)
                for term_count in (15, 10, 6, 4):
                    if px.size <= term_count:
                        continue
                    try:
                        coef_seed = np.ones(term_count)
                        zcoef, _mat, _rms_chief, _rms_centroid, _fitting_error = Kos.Zernike_Fitting(
                            px,
                            py,
                            phase,
                            coef_seed,
                        )
                        used_terms = term_count
                        break
                    except Exception as exc:
                        last_error = exc
                if zcoef is None:
                    raise RuntimeError(f"Zernike fit failed: {last_error}")

                focal = max(0.01, float(getattr(pupil, "EFFL", 0.0)))
                diameter = max(0.01, 2.0 * float(getattr(pupil, "RadPupInp", 1.0)))
                mtf = Kos.calculate_mtf(
                    np.asarray(zcoef, dtype=float),
                    focal,
                    diameter,
                    wavelength,
                    pixels=256,
                    PupilSample=4,
                )
                samples = int(mtf.shape[0])
                freq_max = diameter / max(wavelength * 1e-3, 1e-9)
                freq = np.linspace(0.0, freq_max, samples // 2) / 10.0
                tangential = np.asarray(mtf[samples // 2, samples // 2:], dtype=float)
                sagittal = np.asarray(mtf[samples // 2:, samples // 2], dtype=float)
                count = min(len(freq), len(tangential), len(sagittal))
                mtf_mode = self._operand_mtf_mode("MTF @ freq")
                tan_label = "Tangential"
                sag_label = "Sagittal"
                tan_style = {"linewidth": 2.2, "alpha": 1.0}
                sag_style = {"linewidth": 2.2, "alpha": 1.0}
                if mtf_mode == "tangential":
                    sag_style = {"linewidth": 1.2, "alpha": 0.35}
                    tan_label = "Tangential (selected)"
                elif mtf_mode == "sagittal":
                    tan_style = {"linewidth": 1.2, "alpha": 0.35}
                    sag_label = "Sagittal (selected)"
                analysis_ax.plot(
                    freq[:count],
                    tangential[:count],
                    label=tan_label,
                    color="#1f77b4",
                    **tan_style,
                )
                analysis_ax.plot(
                    freq[:count],
                    sagittal[:count],
                    label=sag_label,
                    color="#d62728",
                    **sag_style,
                )
                analysis_ax.set_title(
                    f"Diffraction MTF  |  {pupil.FieldType}=({float(mtf_settings['field_x']):.3g}, {float(mtf_settings['field_y']):.3g})  |  {wavelength:.4g} um"
                )
                analysis_ax.set_xlabel("Spatial frequency [cycles/mm]")
                analysis_ax.set_ylabel("MTF")
                analysis_ax.set_ylim(0.0, 1.05)
                analysis_ax.grid(True, alpha=0.2)
                analysis_ax.legend(loc="upper right", fontsize=8)
                self.append_debug(f"MTF analysis ok: terms={used_terms}, samples={px.size}")
                self._update_analysis_progress("Rendering diffraction MTF", 4, 4)
                self._finish_analysis_progress("MTF analysis", success=True)
            except Exception as exc:
                self.append_debug(f"MTF diffraction path failed: {exc}")
                try:
                    self._update_analysis_progress("Building geometric fallback", 3, 4)
                    dense_count = max(24, self._current_ray_count() * 6)
                    mtf_rays = self._build_analysis_rays(
                        system,
                        wavelength,
                        sample_count=dense_count,
                        pattern="hexapolar",
                        surface_index=int(mtf_settings["surface_index"]),
                        aperture_type=str(mtf_settings["aperture_type"]),
                        aperture_value=float(mtf_settings["aperture_value"]),
                        field_type=str(mtf_settings["field_type"]),
                        field_x=float(mtf_settings["field_x"]),
                        field_y=float(mtf_settings["field_y"]),
                    )
                    x_local, y_local, _z_local, _l_local, _m_local, _n_local = self._pick_image_plane_data(mtf_rays)
                    x_local = np.asarray(x_local, dtype=float)
                    y_local = np.asarray(y_local, dtype=float)
                    finite = np.isfinite(x_local) & np.isfinite(y_local)
                    x_local = x_local[finite]
                    y_local = y_local[finite]
                    if x_local.size < 4:
                        raise RuntimeError("Not enough image-plane ray samples for geometric MTF")

                    span_x = max(float(np.ptp(x_local)), 1e-3)
                    span_y = max(float(np.ptp(y_local)), 1e-3)
                    span = max(span_x, span_y) * 1.25
                    if span <= 0:
                        span = 1.0
                    bins = 128
                    hist, xedges, yedges = np.histogram2d(
                        x_local,
                        y_local,
                        bins=bins,
                        range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
                    )
                    psf = hist / max(np.sum(hist), 1.0)
                    otf = np.fft.fftshift(np.fft.fft2(psf))
                    mtf = np.abs(otf)
                    mtf /= max(float(np.max(mtf)), 1e-12)

                    dx = float(xedges[1] - xedges[0])
                    freq = np.fft.fftshift(np.fft.fftfreq(bins, d=dx))
                    center = bins // 2
                    positive = freq[center:]
                    tangential = mtf[center, center:]
                    sagittal = mtf[center:, center]
                    count = min(len(positive), len(tangential), len(sagittal))

                    plot_freq = positive[:count]
                    plot_tan = tangential[:count]
                    plot_sag = sagittal[:count]
                    target_freq = self._current_mtf_frequency()
                    mtf_mode = self._operand_mtf_mode("MTF @ freq")
                    tan_value = float(np.interp(target_freq, plot_freq, plot_tan, left=plot_tan[0], right=plot_tan[-1]))
                    sag_value = float(np.interp(target_freq, plot_freq, plot_sag, left=plot_sag[0], right=plot_sag[-1]))
                    if mtf_mode == "tangential":
                        selected_value = tan_value
                        selected_color = "#1f77b4"
                        selected_label = "Tangential"
                    elif mtf_mode == "sagittal":
                        selected_value = sag_value
                        selected_color = "#d62728"
                        selected_label = "Sagittal"
                    else:
                        selected_value = 0.5 * (tan_value + sag_value)
                        selected_color = "#2c3e50"
                        selected_label = "Average"
                    tan_label = "Tangential"
                    sag_label = "Sagittal"
                    tan_style = {"linewidth": 2.2, "alpha": 1.0}
                    sag_style = {"linewidth": 2.2, "alpha": 1.0}
                    if mtf_mode == "tangential":
                        sag_style = {"linewidth": 1.2, "alpha": 0.35}
                        tan_label = "Tangential (selected)"
                    elif mtf_mode == "sagittal":
                        tan_style = {"linewidth": 1.2, "alpha": 0.35}
                        sag_label = "Sagittal (selected)"
                    analysis_ax.plot(
                        plot_freq,
                        plot_tan,
                        label=tan_label,
                        color="#1f77b4",
                        **tan_style,
                    )
                    analysis_ax.plot(
                        plot_freq,
                        plot_sag,
                        label=sag_label,
                        color="#d62728",
                        **sag_style,
                    )
                    analysis_ax.axvline(target_freq, color="#2c3e50", linewidth=1.0, linestyle="--", alpha=0.8)
                    analysis_ax.scatter(
                        [target_freq],
                        [tan_value],
                        color="#1f77b4",
                        s=32 if mtf_mode == "tangential" else 20,
                        alpha=1.0 if mtf_mode != "sagittal" else 0.45,
                        zorder=4,
                    )
                    analysis_ax.scatter(
                        [target_freq],
                        [sag_value],
                        color="#d62728",
                        s=32 if mtf_mode == "sagittal" else 20,
                        alpha=1.0 if mtf_mode != "tangential" else 0.45,
                        zorder=4,
                    )
                    analysis_ax.scatter(
                        [target_freq],
                        [selected_value],
                        color=selected_color,
                        edgecolors="white",
                        linewidths=0.8,
                        s=48,
                        zorder=5,
                    )
                    analysis_ax.set_title(
                        f"Geometric MTF  |  {mtf_settings['field_type']}=({float(mtf_settings['field_x']):.3g}, {float(mtf_settings['field_y']):.3g})  |  {wavelength:.4g} um"
                    )
                    analysis_ax.set_xlabel("Spatial frequency [cycles/mm]")
                    analysis_ax.set_ylabel("MTF")
                    analysis_ax.set_ylim(0.0, 1.05)
                    xmax = min(float(plot_freq[-1]), max(100.0, target_freq * 2.5))
                    analysis_ax.set_xlim(0.0, xmax)
                    analysis_ax.grid(True, alpha=0.2)
                    analysis_ax.legend(loc="upper right", fontsize=8)
                    analysis_ax.text(
                        0.02,
                        0.02,
                        f"{target_freq:.1f} cy/mm\\n{selected_label}={selected_value:.3f}\\nT={tan_value:.3f}  S={sag_value:.3f}",
                        transform=analysis_ax.transAxes,
                        ha="left",
                        va="bottom",
                        fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
                    )
                    self.append_debug(f"MTF fallback ok: rays={x_local.size}, bins={bins}, pupil_samp={dense_count}")
                    self._update_analysis_progress("Rendering geometric MTF", 4, 4)
                    self._finish_analysis_progress("MTF analysis", success=True)
                except Exception as fallback_exc:
                    self.append_debug(f"MTF analysis error: {fallback_exc}")
                    analysis_ax.text(0.5, 0.5, "MTF analysis unavailable", ha="center", va="center")
                    analysis_ax.set_axis_off()
                    self._finish_analysis_progress("MTF analysis", success=False)
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

    def _current_object_mode(self) -> str:
        mode = self.object_mode_var.get().strip()
        return mode if mode in {"Finite", "Infinity"} else "Finite"

    def _current_object_distance(self) -> float:
        if self.rows:
            distance = float(self.rows[0].thickness)
        else:
            distance = 100.0
        return max(distance, 1e-6)

    def _current_field_count(self) -> int:
        try:
            return max(1, int(self.field_count_var.get()))
        except ValueError:
            return 1

    def _current_field_type(self) -> str:
        value = self.field_type_var.get().strip()
        if value in {"Angle", "Object Height", "Paraxial Image Height", "Real Image Height"}:
            return value
        return "Angle"

    def _current_field_value(self) -> float:
        try:
            return float(self.field_value_var.get())
        except ValueError:
            return 0.0

    def _current_field_angle_deg(self) -> float:
        return float(self._field_metrics().get("angle_deg", 0.0))

    def _current_field_height(self) -> float:
        return float(self._field_metrics().get("object_height", 0.0))

    def _field_metrics(self) -> dict[str, float]:
        field_type = self._current_field_type()
        raw_value = self._current_field_value()
        object_distance = self._current_object_distance()
        effl = self._current_effl_estimate()
        image_distance = self._current_image_distance()

        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            if field_type == "Angle":
                angle_deg = raw_value
                object_height = object_distance * np.tan(np.deg2rad(angle_deg))
            elif field_type == "Object Height":
                object_height = raw_value
                angle_deg = np.rad2deg(np.arctan2(object_height, object_distance))
            elif field_type == "Paraxial Image Height":
                paraxial_image_height = raw_value
                angle_deg = np.rad2deg(np.arctan2(paraxial_image_height, max(effl, 1e-6)))
                object_height = object_distance * np.tan(np.deg2rad(angle_deg))
            else:
                real_image_height = raw_value
                angle_deg = np.rad2deg(np.arctan2(real_image_height, max(image_distance, 1e-6)))
                object_height = object_distance * np.tan(np.deg2rad(angle_deg))

            paraxial_image_height = effl * np.tan(np.deg2rad(angle_deg))
            real_image_height = image_distance * np.tan(np.deg2rad(angle_deg))

        if not np.isfinite(angle_deg):
            angle_deg = 0.0
        if not np.isfinite(object_height):
            object_height = 0.0
        if not np.isfinite(paraxial_image_height):
            paraxial_image_height = 0.0
        if not np.isfinite(real_image_height):
            real_image_height = 0.0
        return {
            "angle_deg": float(angle_deg),
            "object_height": float(object_height),
            "paraxial_image_height": float(paraxial_image_height),
            "real_image_height": float(real_image_height),
        }

    def _current_effl_estimate(self) -> float:
        if self.last_system is not None:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    _a, _b, _c, _d, effl, *_rest = self.last_system.EFL(self._current_wavelength())  # type: ignore[misc]
                return max(abs(float(effl)), 1e-6)
            except Exception:
                pass
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    _, _, _, _, _, _, _, effl, *_rest = self.last_system.Parax(self._current_wavelength())
                return max(abs(float(effl)), 1e-6)
            except Exception:
                pass
        return 100.0

    def _current_image_distance(self) -> float:
        if len(self.rows) >= 2:
            return max(float(self.rows[-2].thickness), 1e-6)
        return 100.0

    def _schedule_refresh_plot(self, *_args) -> None:
        if not self.winfo_exists():
            return
        if hasattr(self, "_refresh_after_id") and self._refresh_after_id is not None:
            self.after_cancel(self._refresh_after_id)
        self._refresh_after_id = self.after(120, self._refresh_plot_from_controls)

    def _refresh_plot_from_controls(self) -> None:
        self._refresh_after_id = None
        if self.optimization_running:
            return
        self.refresh_plot()

    def _style_embedded_plot(self, surf_line_count: int) -> None:
        for index, line in enumerate(self.ax.lines):
            if index < surf_line_count:
                line.set_linewidth(max(line.get_linewidth(), 1.25))

    def _field_colors(self, count: int) -> list[str]:
        if count <= 1:
            return ["#39FF14"]
        cmap = [
            "#39FF14",
            "#00E5FF",
            "#FF9F1C",
            "#FF4D6D",
            "#9B5DE5",
            "#FFD166",
            "#2EC4B6",
            "#E71D36",
        ]
        return [cmap[i % len(cmap)] for i in range(count)]

    def _draw_colored_rays(self, rays) -> None:
        if self._has_off_axis_geometry():
            before = len(self.ax.lines)
            try:
                Plot2DRays(rays, 0, 0, self.ax, 0)
                for line in self.ax.lines[before:]:
                    line.set_linewidth(1.8)
                    line.set_alpha(0.95)
            except Exception:
                for ray in rays.CC:
                    points = np.asarray(ray, dtype=float)
                    if points.shape[0] < 2:
                        continue
                    self.ax.plot(points[:, 2], points[:, 1], color="#39FF14", linewidth=1.8, alpha=0.95)
            return
        ray_count = max(1, self._preview_field_ray_count)
        field_count = max(1, self._current_field_count())
        colors = self._field_colors(field_count)
        for index, ray in enumerate(rays.CC):
            points = np.asarray(ray)
            if points.shape[0] < 2:
                continue
            field_index = min(index // ray_count, field_count - 1)
            color = colors[field_index]
            self.ax.plot(points[:, 2], points[:, 1], color=color, linewidth=1.8, alpha=0.95)

    def _current_display_orientation(self) -> str:
        value = getattr(self, "display_orientation_var", None)
        if value is None:
            return "Vertical"
        mode = value.get().strip()
        return mode if mode in {"Vertical", "Horizontal"} else "Vertical"

    def _project_xy(self, z, y):
        z_arr = np.asarray(z, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        if self._current_display_orientation() == "Horizontal":
            return -y_arr, -z_arr
        return z_arr, y_arr

    def _apply_display_orientation_to_lines(self, start_index: int = 0) -> None:
        if self._current_display_orientation() == "Vertical":
            return
        for line in self.ax.lines[start_index:]:
            xdata = np.asarray(line.get_xdata(orig=False), dtype=float)
            ydata = np.asarray(line.get_ydata(orig=False), dtype=float)
            if xdata.size == 0 or ydata.size == 0:
                continue
            proj_x, proj_y = self._project_xy(xdata, ydata)
            line.set_xdata(proj_x)
            line.set_ydata(proj_y)

    def _has_off_axis_geometry(self) -> bool:
        for row in self.rows:
            if row.surface == "Mirror":
                return True
            if any(
                abs(value) > 1e-9
                for value in (row.tilt_x, row.tilt_y, row.tilt_z, row.desp_x, row.desp_y, row.desp_z, row.axis_move)
            ):
                return True
        return False

    def _has_native_drawn_surfaces(self) -> bool:
        for row in self.rows:
            if row.surface not in {"Object", "Image", "Mirror"}:
                return True
        return False

    def _is_folded_mirror_preview_mode(self) -> bool:
        if self._current_display_orientation() != "Horizontal":
            return False
        mirror_count = 0
        for row in self.rows:
            if row.surface == "Mirror":
                mirror_count += 1
            elif row.surface not in {"Object", "Image", "Standard"}:
                return False
        return mirror_count >= 1

    @staticmethod
    def _reflect_2d(direction: np.ndarray, line_angle_deg: float) -> np.ndarray:
        theta = np.deg2rad(float(line_angle_deg))
        tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
        tangent_norm = np.linalg.norm(tangent)
        if tangent_norm <= 1e-12:
            return direction
        tangent /= tangent_norm
        normal = np.array([-tangent[1], tangent[0]], dtype=float)
        reflected = direction - 2.0 * np.dot(direction, normal) * normal
        norm = np.linalg.norm(reflected)
        if norm <= 1e-12:
            return direction
        return reflected / norm

    @staticmethod
    def _intersect_ray_with_line(
        origin: np.ndarray,
        direction: np.ndarray,
        center: np.ndarray,
        line_angle_deg: float,
    ) -> tuple[np.ndarray | None, float | None]:
        theta = np.deg2rad(float(line_angle_deg))
        tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
        matrix = np.column_stack((direction, -tangent))
        try:
            t_ray, t_line = np.linalg.solve(matrix, center - origin)
        except np.linalg.LinAlgError:
            return None, None
        if t_ray < 0:
            return None, None
        point = origin + direction * t_ray
        return point, float(t_line)

    @staticmethod
    def _glass_index_for_preview(name: str) -> float:
        glass = str(name).strip().upper()
        if glass in {"", "AIR", "NULL"}:
            return 1.0
        if glass == "BK7":
            return 1.5168
        if glass == "F2":
            return 1.6200
        if glass in {"FS", "F_SILICA", "SILICA"}:
            return 1.4585
        return 1.5

    @staticmethod
    def _intersect_ray_with_spherical_surface(
        origin: np.ndarray,
        direction: np.ndarray,
        vertex: np.ndarray,
        axis_dir: np.ndarray,
        radius: float,
    ) -> tuple[np.ndarray | None, float | None]:
        if abs(radius) <= 1e-9:
            return None, None
        axis = np.asarray(axis_dir, dtype=float)
        axis /= max(np.linalg.norm(axis), 1e-12)
        tangent = np.array([-axis[1], axis[0]], dtype=float)
        center = vertex + axis * float(radius)
        oc = origin - center
        a = float(np.dot(direction, direction))
        b = 2.0 * float(np.dot(direction, oc))
        c = float(np.dot(oc, oc) - radius * radius)
        disc = b * b - 4.0 * a * c
        if disc < 0.0:
            return None, None
        root = np.sqrt(disc)
        candidates = [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]
        candidates = [t for t in candidates if t >= 1e-9]
        if not candidates:
            return None, None
        t_ray = min(candidates)
        point = origin + direction * t_ray
        local = point - vertex
        return point, float(np.dot(local, tangent))

    @staticmethod
    def _refract_ray_2d(direction: np.ndarray, normal: np.ndarray, n_before: float, n_after: float) -> np.ndarray:
        d = np.asarray(direction, dtype=float)
        d /= max(np.linalg.norm(d), 1e-12)
        n = np.asarray(normal, dtype=float)
        n /= max(np.linalg.norm(n), 1e-12)
        if np.dot(d, n) > 0.0:
            n = -n
        eta = float(n_before) / float(n_after)
        cos_i = -float(np.dot(n, d))
        k = 1.0 - eta * eta * (1.0 - cos_i * cos_i)
        if k < 0.0:
            reflected = d + 2.0 * cos_i * n
            return reflected / max(np.linalg.norm(reflected), 1e-12)
        refracted = eta * d + (eta * cos_i - np.sqrt(k)) * n
        return refracted / max(np.linalg.norm(refracted), 1e-12)

    def _plot_folded_mirror_preview(self, analysis_ax, draw_optical_surfaces: bool = True, draw_rays: bool = True) -> None:
        self.ax.clear()
        if analysis_ax is not None:
            analysis_ax.clear()
            analysis_ax.text(0.5, 0.5, "Analysis unavailable\nfor folded preview", ha="center", va="center")
            analysis_ax.set_axis_off()

        point, direction, max_half, extent_points, elements = self._draw_folded_layout_geometry()

        colors = self._field_colors(max(1, self._current_field_count()))
        field_values = self._sample_field_values(self._current_field_height())
        pupil_samples = self._sample_ray_heights(self._entrance_radius(max_half))
        fan_angles = self._sample_fan_angles_deg() if self._current_object_mode() == "Finite" else [0.0]
        ray_paths: list[list[np.ndarray]] = []
        for field_index, field_value in enumerate(field_values):
            if self._current_object_mode() == "Infinity":
                for pupil_y in pupil_samples:
                    d = direction.copy()
                    origin = point + np.array([float(field_value) + float(pupil_y), 0.0], dtype=float)
                    p = origin.copy()
                    path = [origin.copy()]
                    current_dir = d
                    current_medium = 1.0
                    for surface_type, center, row, branch_dir in elements:
                        if surface_type == "Mirror":
                            hit, along = self._intersect_ray_with_line(p, current_dir, center, float(row.tilt_x))
                            if hit is None:
                                break
                            half = max(row.diameter / 2.0, 0.5)
                            if along is not None and abs(along) > half:
                                break
                            if np.linalg.norm(hit - path[-1]) > 1e-9:
                                path.append(hit.copy())
                            p = hit
                            current_dir = self._reflect_2d(current_dir, float(row.tilt_x))
                        elif surface_type == "Standard":
                            hit, along = self._intersect_ray_with_spherical_surface(
                                p, current_dir, center, branch_dir, float(row.rc)
                            )
                            if hit is None:
                                break
                            half = max(row.diameter / 2.0, 0.5)
                            if along is not None and abs(along) > half:
                                break
                            if np.linalg.norm(hit - path[-1]) > 1e-9:
                                path.append(hit.copy())
                            axis = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
                            sphere_center = center + axis * float(row.rc)
                            normal = hit - sphere_center
                            next_medium = self._glass_index_for_preview(row.glass)
                            current_dir = self._refract_ray_2d(current_dir, normal, current_medium, next_medium)
                            current_medium = next_medium
                            p = hit
                        elif surface_type == "Image":
                            tangent = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
                            angle = np.rad2deg(np.arctan2(tangent[1], tangent[0]))
                            hit, along = self._intersect_ray_with_line(p, current_dir, center, angle)
                            if hit is None:
                                break
                            half = max(row.diameter / 2.0, 0.5)
                            if along is not None and abs(along) > half:
                                break
                            if np.linalg.norm(hit - path[-1]) > 1e-9:
                                path.append(hit.copy())
                            p = hit
                            break
                    ray_paths.append(path)
                    extent_points.extend(path)
            else:
                for fan_angle in fan_angles:
                    angle = np.deg2rad(float(fan_angle))
                    d = np.array([np.sin(angle), np.cos(angle)], dtype=float)
                    norm = np.linalg.norm(d)
                    if norm > 1e-12:
                        d /= norm
                    origin = point + np.array([float(field_value), 0.0], dtype=float)
                    p = origin.copy()
                    path = [origin.copy()]
                    current_dir = d
                    current_medium = 1.0
                    for surface_type, center, row, branch_dir in elements:
                        if surface_type == "Mirror":
                            hit, along = self._intersect_ray_with_line(p, current_dir, center, float(row.tilt_x))
                            if hit is None:
                                break
                            half = max(row.diameter / 2.0, 0.5)
                            if along is not None and abs(along) > half:
                                break
                            if np.linalg.norm(hit - path[-1]) > 1e-9:
                                path.append(hit.copy())
                            p = hit
                            current_dir = self._reflect_2d(current_dir, float(row.tilt_x))
                        elif surface_type == "Standard":
                            hit, along = self._intersect_ray_with_spherical_surface(
                                p, current_dir, center, branch_dir, float(row.rc)
                            )
                            if hit is None:
                                break
                            half = max(row.diameter / 2.0, 0.5)
                            if along is not None and abs(along) > half:
                                break
                            if np.linalg.norm(hit - path[-1]) > 1e-9:
                                path.append(hit.copy())
                            axis = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
                            sphere_center = center + axis * float(row.rc)
                            normal = hit - sphere_center
                            next_medium = self._glass_index_for_preview(row.glass)
                            current_dir = self._refract_ray_2d(current_dir, normal, current_medium, next_medium)
                            current_medium = next_medium
                            p = hit
                        elif surface_type == "Image":
                            tangent = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
                            angle = np.rad2deg(np.arctan2(tangent[1], tangent[0]))
                            hit, along = self._intersect_ray_with_line(p, current_dir, center, angle)
                            if hit is None:
                                break
                            half = max(row.diameter / 2.0, 0.5)
                            if along is not None and abs(along) > half:
                                break
                            if np.linalg.norm(hit - path[-1]) > 1e-9:
                                path.append(hit.copy())
                            p = hit
                            break
                    ray_paths.append(path)
                    extent_points.extend(path)
        for surface_type, center, row, branch_dir in elements:
            if surface_type in {"Mirror", "Standard"} and not draw_optical_surfaces:
                continue
            if surface_type == "Mirror":
                half = max(row.diameter / 2.0, 0.5)
                theta = np.deg2rad(float(row.tilt_x))
                tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
                tangent /= max(np.linalg.norm(tangent), 1e-12)
                p0 = center - tangent * half
                p1 = center + tangent * half
                self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#202020", linewidth=2.2, alpha=0.95)
                extent_points.extend([p0, p1])
            elif surface_type == "Standard":
                axis = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
                tangent = np.array([-axis[1], axis[0]], dtype=float)
                half = max(row.diameter / 2.0, 0.5)
                yy = np.linspace(-half, half, 128)
                if abs(float(row.rc)) <= half + 1e-9:
                    xx = np.zeros_like(yy)
                else:
                    rr = abs(float(row.rc))
                    sign = 1.0 if float(row.rc) >= 0.0 else -1.0
                    xx = float(row.rc) - sign * np.sqrt(np.maximum(rr * rr - yy * yy, 0.0))
                pts = center[None, :] + np.outer(xx, axis) + np.outer(yy, tangent)
                self.ax.plot(pts[:, 0], pts[:, 1], color="#2563eb", linewidth=1.8, alpha=0.95)
                extent_points.extend(pts.tolist())
            elif surface_type == "Image":
                self._draw_folded_plane(center, row.diameter, branch_dir, row.name or "Image", max_half, extent_points)

        if draw_rays:
            rays_per_field = max(1, len(ray_paths) // max(1, len(field_values)))
            for index, path in enumerate(ray_paths):
                pts = np.asarray(path, dtype=float)
                if pts.shape[0] < 2:
                    continue
                field_index = min(index // rays_per_field, len(colors) - 1)
                self.ax.plot(pts[:, 0], pts[:, 1], color=colors[field_index], linewidth=1.8, alpha=0.95)

        ext = np.asarray(extent_points, dtype=float)
        if ext.size:
            x_min, y_min = np.min(ext[:, 0]), np.min(ext[:, 1])
            x_max, y_max = np.max(ext[:, 0]), np.max(ext[:, 1])
            span_x = max(x_max - x_min, 1.0)
            span_y = max(y_max - y_min, 1.0)
            self.ax.set_xlim(x_min - 0.1 * span_x, x_max + 0.1 * span_x)
            self.ax.set_ylim(y_max + 0.1 * span_y, y_min - 0.1 * span_y)
            self.ax.set_aspect("equal", adjustable="box")

        self._set_results(
            [
                ("Mode", "Folded mirror preview"),
                ("Surface count", str(len(self.rows))),
                ("Object mode", self._current_object_mode()),
                ("Field type", self._current_field_type()),
                ("Field count", str(self._current_field_count())),
                ("Ray fan count", str(self._current_ray_count())),
            ]
        )
        self.status_var.set("Folded mirror preview")

    def _draw_folded_layout_geometry(self):
        point = np.array([0.0, 0.0], dtype=float)
        direction = np.array([0.0, 1.0], dtype=float)
        max_half = max((max(row.diameter / 2.0, 0.5) for row in self.rows), default=1.0)
        extent_points = [point.copy()]

        if self.rows:
            first = self.rows[0]
            self._draw_folded_plane(point, first.diameter, direction, first.name or "Object", max_half, extent_points)

        elements: list[tuple[str, np.ndarray, SurfaceRow, np.ndarray]] = []
        current_dir = direction.copy()
        object_thickness = max(float(self.rows[0].thickness), 0.0) if self.rows else 0.0
        current_point = point + current_dir * object_thickness
        extent_points.append(current_point.copy())
        for row in self.rows[1:]:
            elements.append((row.surface, current_point.copy(), row, current_dir.copy()))
            if row.surface == "Mirror":
                current_dir = self._reflect_2d(current_dir, float(row.tilt_x))
            travel = max(float(row.thickness), 0.0)
            current_point = current_point + current_dir * travel
            extent_points.append(current_point.copy())

        for surface_type, center, row, branch_dir in elements:
            if surface_type == "Mirror":
                half = max(row.diameter / 2.0, 0.5)
                theta = np.deg2rad(float(row.tilt_x))
                tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
                tangent /= max(np.linalg.norm(tangent), 1e-12)
                p0 = center - tangent * half
                p1 = center + tangent * half
                self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#202020", linewidth=2.2, alpha=0.95)
                extent_points.extend([p0, p1])
            elif surface_type == "Standard":
                axis = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
                tangent = np.array([-axis[1], axis[0]], dtype=float)
                half = max(row.diameter / 2.0, 0.5)
                yy = np.linspace(-half, half, 128)
                if abs(float(row.rc)) <= half + 1e-9:
                    xx = np.zeros_like(yy)
                else:
                    rr = abs(float(row.rc))
                    sign = 1.0 if float(row.rc) >= 0.0 else -1.0
                    xx = float(row.rc) - sign * np.sqrt(np.maximum(rr * rr - yy * yy, 0.0))
                pts = center[None, :] + np.outer(xx, axis) + np.outer(yy, tangent)
                self.ax.plot(pts[:, 0], pts[:, 1], color="#2563eb", linewidth=1.8, alpha=0.95)
                extent_points.extend(pts.tolist())
            elif surface_type == "Image":
                self._draw_folded_plane(center, row.diameter, branch_dir, row.name or "Image", max_half, extent_points)

        return point, direction, max_half, extent_points, elements

    def _draw_folded_plane(
        self,
        center: np.ndarray,
        diameter: float,
        along: np.ndarray,
        label: str,
        max_half: float,
        extent_points: list[np.ndarray],
    ) -> None:
        tangent = np.array([-along[1], along[0]], dtype=float)
        norm = np.linalg.norm(tangent)
        if norm <= 1e-12:
            tangent = np.array([1.0, 0.0], dtype=float)
        else:
            tangent /= norm
        half = max(diameter / 2.0, 0.5)
        p0 = center - tangent * half
        p1 = center + tangent * half
        self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#202020", linewidth=1.2, alpha=0.9)
        self.ax.text(center[0], center[1] - max_half * 0.15, label, ha="center", va="bottom", fontsize=9, color="#202020")
        extent_points.extend([p0, p1])

    def _plot_native_off_axis_preview(self, analysis_ax, max_radius: float, use_native_surfaces: bool = True) -> None:
        folded_visual_mode = self._is_folded_mirror_preview_mode()
        if folded_visual_mode:
            self._plot_folded_mirror_preview(analysis_ax, draw_optical_surfaces=False, draw_rays=False)
        else:
            self.ax.clear()
            if analysis_ax is not None:
                analysis_ax.clear()
                analysis_ax.text(0.5, 0.5, "Native off-axis ray preview\nAnalysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()

        wavelength = self._current_wavelength()
        capture = io.StringIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with redirect_stdout(capture), redirect_stderr(capture):
                system = self.build_system()
                if getattr(system.Pr3D, "ExistSolid", 0) == 0:
                    original_build = system.BUILD
                    system.BUILD = 1
                    system.build()
                    system.BUILD = original_build
                rays = Kos.raykeeper(system)
                self._trace_preview_rays(system, rays, wavelength, max_radius)
        self.append_debug(capture.getvalue())
        self.last_system = system
        self.last_rays = rays

        extent_points: list[np.ndarray] = []
        folded_elements = None
        folded_max_half = max_radius
        if folded_visual_mode:
            _, _, folded_max_half, extent_points, folded_elements = self._draw_folded_layout_geometry()
            if use_native_surfaces:
                self._overlay_native_folded_surfaces(system, folded_elements, extent_points)
        else:
            native_surface_ok = False
            try:
                native_surface_ok = self._plot_native_off_axis_surfaces(system)
            except Exception as exc:
                self.append_debug(f"Native off-axis surface draw fallback: {exc}")
            if not native_surface_ok:
                if self._has_native_drawn_surfaces():
                    try:
                        Plot2DSurf(system, 0, self.ax)
                    except Exception as exc:
                        self.append_debug(f"Native off-axis Plot2DSurf fallback: {exc}")
                self._draw_custom_mirror_surfaces()
            self._draw_reference_plane_labels()

        field_count = max(1, self._current_field_count())
        colors = self._field_colors(field_count)
        ray_count = max(1, self._preview_field_ray_count)
        ray_lengths: list[int] = []
        last_surface_counts: dict[int, int] = {}
        hit_sequences: list[list[int]] = []
        surface_hit_counts: dict[int, int] = {}
        final_surface_index = max(0, len(self.rows) - 1)
        rays_reaching_image = 0
        native_fold_paths = (
            self._native_folded_display_paths(system, rays, folded_elements, folded_max_half)
            if folded_elements is not None
            else None
        )
        for index, ray in enumerate(rays.CC):
            pts = np.asarray(ray, dtype=float)
            if pts.shape[0] < 2:
                continue
            ray_lengths.append(int(pts.shape[0]))
            try:
                surf_ids = np.asarray(rays.SURFACE[index], dtype=int).ravel()
                if surf_ids.size:
                    sequence = [int(v) for v in surf_ids.tolist()]
                    hit_sequences.append(sequence)
                    seen_in_ray: set[int] = set()
                    for surface_index in sequence:
                        if surface_index in seen_in_ray:
                            continue
                        seen_in_ray.add(surface_index)
                        surface_hit_counts[surface_index] = surface_hit_counts.get(surface_index, 0) + 1
                    last_surface = int(surf_ids[-1])
                    last_surface_counts[last_surface] = last_surface_counts.get(last_surface, 0) + 1
                    if last_surface == final_surface_index:
                        rays_reaching_image += 1
            except Exception:
                pass
            field_index = min(index // ray_count, field_count - 1)
            color = colors[field_index]
            if native_fold_paths is not None and index < len(native_fold_paths):
                folded_pts = native_fold_paths[index]
                if folded_pts.shape[0] >= 2:
                    self.ax.plot(folded_pts[:, 0], folded_pts[:, 1], color=color, linewidth=1.8, alpha=0.95)
                    extent_points.extend(folded_pts)
                    continue
                if folded_visual_mode:
                    continue
            x_vals, y_vals = self._project_xy(pts[:, 2], pts[:, 1])
            self.ax.plot(x_vals, y_vals, color=color, linewidth=1.8, alpha=0.95)
            extent_points.extend(np.column_stack((x_vals, y_vals)))

        if folded_visual_mode and folded_elements is not None:
            self._overlay_native_folded_surface_hits(folded_elements, surface_hit_counts, folded_max_half, extent_points)

        if self._current_display_orientation() == "Horizontal":
            if extent_points:
                ext = np.asarray(extent_points, dtype=float)
                x_min, y_min = np.min(ext[:, 0]), np.min(ext[:, 1])
                x_max, y_max = np.max(ext[:, 0]), np.max(ext[:, 1])
                span_x = max(x_max - x_min, 1.0)
                span_y = max(y_max - y_min, 1.0)
                self.ax.set_xlim(x_min - 0.1 * span_x, x_max + 0.1 * span_x)
                self.ax.set_ylim(y_max + 0.1 * span_y, y_min - 0.1 * span_y)
            else:
                self._set_plot_limits_from_drawn_data()
            self.ax.set_aspect("equal", adjustable="box")
        elif self.ax.lines:
            self._set_plot_limits_from_drawn_data()
        else:
            self._set_plot_limits_from_layout(max_radius)

        self._set_results(
            [
                ("Mode", "Native folded preview" if folded_visual_mode else "Native off-axis ray preview"),
                ("Surface count", str(len(self.rows))),
                ("Object mode", self._current_object_mode()),
                ("Field type", self._current_field_type()),
                ("Field count", str(self._current_field_count())),
                ("Ray fan count", str(self._current_ray_count())),
                ("Native rays", str(len(rays.CC))),
                ("Rays to image", str(rays_reaching_image)),
                ("Avg ray hits", f"{np.mean(ray_lengths):.3g}" if ray_lengths else "0"),
                ("Max ray hits", str(max(ray_lengths) if ray_lengths else 0)),
                ("Surface hits", self._native_surface_count_summary(surface_hit_counts)),
                ("Last-hit surfaces", self._native_surface_count_summary(last_surface_counts)),
                *self._native_folded_prescription_summary(folded_elements),
            ]
        )
        if ray_lengths:
            self.append_debug(
                "Native preview rays: count={count}, avg_hits={avg:.3g}, max_hits={mx}".format(
                    count=len(ray_lengths),
                    avg=float(np.mean(ray_lengths)),
                    mx=max(ray_lengths),
                )
            )
        if last_surface_counts:
            self.append_debug(f"Native preview last-hit surfaces: {self._native_surface_count_summary(last_surface_counts)}")
        if surface_hit_counts:
            self.append_debug(f"Native preview surface hits: {self._native_surface_count_summary(surface_hit_counts)}")
        for line in self._native_hit_sequence_lines(hit_sequences):
            self.append_debug(line)
        self.status_var.set("Native folded preview" if folded_visual_mode else "Native off-axis ray preview")

    def _overlay_native_folded_surface_hits(
        self,
        elements,
        surface_hit_counts: dict[int, int],
        max_half: float,
        extent_points: list[np.ndarray],
    ) -> None:
        for surface_index, count in sorted(surface_hit_counts.items()):
            if surface_index <= 0 or surface_index > len(elements):
                continue
            surface_type, center, row, branch_dir = elements[surface_index - 1]
            accent = "#d97706" if surface_type == "Mirror" else "#1d4ed8"
            if surface_type == "Mirror":
                half = max(row.diameter / 2.0, 0.5)
                theta = np.deg2rad(float(row.tilt_x))
                tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
                tangent /= max(np.linalg.norm(tangent), 1e-12)
                p0 = center - tangent * half
                p1 = center + tangent * half
                self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=accent, linewidth=3.0, alpha=0.35, linestyle="--")
                extent_points.extend([p0, p1])
            elif surface_type == "Standard":
                axis = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
                tangent = np.array([-axis[1], axis[0]], dtype=float)
                half = max(row.diameter / 2.0, 0.5)
                yy = np.linspace(-half, half, 128)
                if abs(float(row.rc)) <= half + 1e-9:
                    xx = np.zeros_like(yy)
                else:
                    rr = abs(float(row.rc))
                    sign = 1.0 if float(row.rc) >= 0.0 else -1.0
                    xx = float(row.rc) - sign * np.sqrt(np.maximum(rr * rr - yy * yy, 0.0))
                pts = center[None, :] + np.outer(xx, axis) + np.outer(yy, tangent)
                self.ax.plot(pts[:, 0], pts[:, 1], color=accent, linewidth=2.4, alpha=0.28, linestyle="--")
                extent_points.extend(pts.tolist())
            elif surface_type == "Image":
                tangent = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
                tangent /= max(np.linalg.norm(tangent), 1e-12)
                half = max(row.diameter / 2.0, 0.5)
                p0 = center - tangent * half
                p1 = center + tangent * half
                self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#059669", linewidth=2.6, alpha=0.25, linestyle="--")
                extent_points.extend([p0, p1])
            self.ax.text(
                center[0],
                center[1] + max_half * 0.12,
                f"x{count}",
                ha="center",
                va="bottom",
                fontsize=8,
                color=accent,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 0.4},
            )

    def _overlay_native_folded_surfaces(self, system, elements, extent_points: list[np.ndarray]) -> None:
        transforms = getattr(system, "TRANS_2A", None)
        surfaces = getattr(system, "AAA", None)
        if transforms is None or surfaces is None or not elements:
            return
        block_count = min(len(self.rows), getattr(surfaces, "n_blocks", 0), len(transforms))
        for surface_index in range(1, block_count):
            if surface_index > len(elements):
                break
            surface_type, center2, row, branch_dir = elements[surface_index - 1]
            if surface_type not in {"Mirror", "Standard"}:
                continue
            try:
                poly = surfaces[surface_index]
                pts = np.asarray(poly.points, dtype=float)
            except Exception:
                continue
            if pts.size == 0:
                continue
            try:
                transform = np.asarray(transforms[surface_index], dtype=float)
                if transform.shape == (4, 4):
                    pts_h = np.c_[pts, np.ones(len(pts))]
                    pts = (pts_h @ transform.T)[:, :3]
            except Exception:
                pass
            display_pts = self._map_native_surface_to_folded(surface_type, center2, row, branch_dir, pts)
            if display_pts is None or display_pts.shape[0] < 8:
                continue
            overlay_color = "#7c3aed" if surface_type == "Standard" else "#b45309"
            self.ax.plot(display_pts[:, 0], display_pts[:, 1], color="white", linewidth=4.6, alpha=0.96)
            self.ax.plot(display_pts[:, 0], display_pts[:, 1], color=overlay_color, linewidth=2.2, alpha=0.96)
            extent_points.extend(display_pts.tolist())

    def _map_native_surface_to_folded(
        self, surface_type: str, center2: np.ndarray, row: SurfaceRow, branch_dir: np.ndarray, pts3: np.ndarray
    ) -> np.ndarray | None:
        pts3 = np.asarray(pts3, dtype=float)
        if pts3.size == 0:
            return None
        center_x = float(np.median(pts3[:, 0]))
        tolerance = max(0.08, 0.015 * max(np.ptp(pts3[:, 1]), np.ptp(pts3[:, 2]), 1.0))
        mask = np.abs(pts3[:, 0] - center_x) <= tolerance
        sliced = pts3[mask]
        if sliced.shape[0] < 16:
            tolerance = max(0.2, tolerance * 2.0)
            mask = np.abs(pts3[:, 0] - center_x) <= tolerance
            sliced = pts3[mask]
        if sliced.shape[0] < 8:
            return None
        center3 = np.mean(sliced, axis=0)
        centered = sliced - center3
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            return None
        vectors = [np.asarray(v, dtype=float) for v in vh]
        yz_ranked = [
            (float(np.linalg.norm(vec[1:])), np.asarray(vec, dtype=float) / max(np.linalg.norm(vec), 1e-12))
            for vec in vectors
        ]
        yz_ranked.sort(key=lambda item: item[0], reverse=True)
        if len(yz_ranked) < 2:
            return None
        basis_a = yz_ranked[0][1]
        basis_b = yz_ranked[1][1]
        u_a = centered @ basis_a
        u_b = centered @ basis_b
        if np.ptp(u_a) >= np.ptp(u_b):
            aperture_component = u_a
            sag_component = u_b
        else:
            aperture_component = u_b
            sag_component = u_a
        edge_threshold = 0.7 * max(np.max(np.abs(aperture_component)), 1e-9)
        edge_mask = np.abs(aperture_component) >= edge_threshold
        center_mask = np.abs(aperture_component) <= 0.15 * max(np.max(np.abs(aperture_component)), 1e-9)
        if np.any(edge_mask) and np.any(center_mask) and abs(float(row.rc)) > 1e-9:
            edge_sag = float(np.mean(sag_component[edge_mask]))
            center_sag = float(np.mean(sag_component[center_mask]))
            expected_sign = 1.0 if float(row.rc) > 0.0 else -1.0
            observed_sign = np.sign(edge_sag - center_sag)
            if observed_sign != 0.0 and observed_sign != expected_sign:
                sag_component = -sag_component
        tangent2 = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
        tangent2 /= max(np.linalg.norm(tangent2), 1e-12)
        axis2 = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
        if surface_type == "Mirror":
            display_pts = center2[None, :] + np.outer(aperture_component, tangent2)
        else:
            display_pts = center2[None, :] + np.outer(aperture_component, tangent2) + np.outer(sag_component, axis2)
        order = np.argsort(display_pts[:, 1])
        return display_pts[order]

    def _native_folded_prescription_summary(self, elements) -> list[tuple[str, str]]:
        if not elements:
            return []
        summary: list[tuple[str, str]] = []
        summary.append(("Obj -> first [mm]", f"{max(float(self.rows[0].thickness), 0.0):.3g}"))
        for surface_index, (surface_type, _, row, _) in enumerate(elements, start=1):
            if surface_type == "Image":
                continue
            label = f"{surface_index}:{row.name or row.surface} -> next"
            summary.append((f"{label} [mm]", f"{max(float(row.thickness), 0.0):.3g}"))
            if surface_type == "Mirror":
                summary.append((f"{surface_index}:{row.name or row.surface} TiltX [deg]", f"{float(row.tilt_x):.3g}"))
        return summary

    def _preview_ray_start_specs(self, max_half: float) -> list[tuple[np.ndarray, np.ndarray]]:
        point = np.array([0.0, 0.0], dtype=float)
        direction = np.array([0.0, 1.0], dtype=float)
        starts: list[tuple[np.ndarray, np.ndarray]] = []
        field_values = self._sample_field_values(
            self._current_field_angle_deg() if self._current_object_mode() == "Infinity" else self._current_field_height()
        )
        pupil_samples = self._sample_ray_heights(self._entrance_radius(max_half))
        if self._current_object_mode() == "Infinity":
            for field_value in field_values:
                angle = np.deg2rad(float(field_value))
                d = np.array([np.sin(angle), np.cos(angle)], dtype=float)
                d /= max(np.linalg.norm(d), 1e-12)
                for pupil_y in pupil_samples:
                    origin = point + np.array([float(pupil_y), 0.0], dtype=float)
                    starts.append((origin, d.copy()))
        else:
            object_distance = max(float(self.rows[0].thickness), 1e-9) if self.rows else 1.0
            for field_value in field_values:
                origin = point + np.array([float(field_value), 0.0], dtype=float)
                for pupil_y in pupil_samples:
                    target = np.array([float(pupil_y), object_distance], dtype=float)
                    d = target - origin
                    d /= max(np.linalg.norm(d), 1e-12)
                    starts.append((origin.copy(), d))
        return starts

    def _native_folded_display_paths(self, system, rays, elements, max_half: float) -> list[np.ndarray]:
        if elements is None:
            return []
        starts = self._preview_ray_start_specs(max_half)
        element_map = {index + 1: item for index, item in enumerate(elements)}
        paths: list[np.ndarray] = []
        for ray_index, surface_ids_raw in enumerate(rays.SURFACE):
            if ray_index >= len(starts):
                break
            origin, current_dir = starts[ray_index]
            current_point = origin.copy()
            current_medium = 1.0
            path = [origin.copy()]
            surface_ids = [int(v) for v in np.asarray(surface_ids_raw, dtype=int).ravel().tolist()]
            last_id: int | None = None
            for surface_index in surface_ids:
                if surface_index == last_id:
                    continue
                last_id = surface_index
                element = element_map.get(surface_index)
                if element is None:
                    continue
                surface_type, center, row, branch_dir = element
                if surface_type == "Mirror":
                    hit, along = self._intersect_ray_with_line(current_point, current_dir, center, float(row.tilt_x))
                    if hit is None:
                        break
                    half = max(row.diameter / 2.0, 0.5)
                    if along is not None and abs(along) > half:
                        break
                    if np.linalg.norm(hit - path[-1]) > 1e-9:
                        path.append(hit.copy())
                    current_point = hit
                    current_dir = self._reflect_2d(current_dir, float(row.tilt_x))
                elif surface_type == "Standard":
                    hit, along = self._intersect_ray_with_spherical_surface(
                        current_point, current_dir, center, branch_dir, float(row.rc)
                    )
                    if hit is None:
                        break
                    half = max(row.diameter / 2.0, 0.5)
                    if along is not None and abs(along) > half:
                        break
                    if np.linalg.norm(hit - path[-1]) > 1e-9:
                        path.append(hit.copy())
                    axis = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
                    sphere_center = center + axis * float(row.rc)
                    normal = hit - sphere_center
                    next_medium = self._glass_index_for_preview(row.glass)
                    current_dir = self._refract_ray_2d(current_dir, normal, current_medium, next_medium)
                    current_medium = next_medium
                    current_point = hit
                elif surface_type == "Image":
                    tangent = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
                    angle = np.rad2deg(np.arctan2(tangent[1], tangent[0]))
                    hit, along = self._intersect_ray_with_line(current_point, current_dir, center, angle)
                    if hit is None:
                        break
                    half = max(row.diameter / 2.0, 0.5)
                    if along is not None and abs(along) > half:
                        break
                    if np.linalg.norm(hit - path[-1]) > 1e-9:
                        path.append(hit.copy())
                    current_point = hit
            paths.append(np.asarray(path, dtype=float))
        return paths

    def _native_surface_count_summary(self, counts: dict[int, int]) -> str:
        if not counts:
            return "none"
        parts = []
        for surface_index in sorted(counts):
            label = str(surface_index)
            if 0 <= surface_index < len(self.rows):
                row = self.rows[surface_index]
                label = f"{surface_index}:{row.name or row.surface}"
            parts.append(f"{label} x{counts[surface_index]}")
        return ", ".join(parts)

    def _native_hit_sequence_lines(self, sequences: list[list[int]]) -> list[str]:
        if not sequences:
            return []
        lines: list[str] = []
        preview_count = min(8, len(sequences))
        for ray_index in range(preview_count):
            seq = sequences[ray_index]
            labels: list[str] = []
            last_value: int | None = None
            for surface_index in seq:
                if surface_index == last_value:
                    continue
                last_value = surface_index
                if 0 <= surface_index < len(self.rows):
                    row = self.rows[surface_index]
                    labels.append(f"{surface_index}:{row.name or row.surface}")
                else:
                    labels.append(str(surface_index))
            lines.append(f"Native ray {ray_index} hits: {' -> '.join(labels) if labels else 'none'}")
        if len(sequences) > preview_count:
            lines.append(f"Native hit sequences truncated: {len(sequences) - preview_count} more rays")
        return lines

    def _plot_native_off_axis_surfaces(self, system) -> bool:
        transforms = getattr(system, "TRANS_2A", None)
        surfaces = getattr(system, "AAA", None)
        if transforms is None or surfaces is None:
            return False
        drew_any = False
        block_count = min(len(self.rows), getattr(surfaces, "n_blocks", 0), len(transforms))
        for index in range(block_count):
            row = self.rows[index]
            if row.surface in {"Object", "Image"}:
                continue
            try:
                poly = surfaces[index]
                pts = np.asarray(poly.points, dtype=float)
            except Exception:
                continue
            if pts.size == 0:
                continue
            try:
                transform = np.asarray(transforms[index], dtype=float)
                if transform.shape == (4, 4):
                    pts_h = np.c_[pts, np.ones(len(pts))]
                    pts = (pts_h @ transform.T)[:, :3]
            except Exception:
                pass
            proj_x, proj_y = self._native_surface_projection(pts, row)
            color = "#202020" if row.surface == "Mirror" else "#2563eb"
            if proj_x is None or proj_y is None or len(proj_x) == 0:
                continue
            if len(proj_x) >= 8:
                self.ax.plot(proj_x, proj_y, color=color, linewidth=1.4, alpha=0.85)
            else:
                self.ax.scatter(proj_x, proj_y, s=4, c=color, alpha=0.35, linewidths=0)
            drew_any = True
        return drew_any

    def _native_surface_projection(self, pts: np.ndarray, row: SurfaceRow):
        if pts.size == 0:
            return None, None
        pts = np.asarray(pts, dtype=float)
        if row.surface == "Mirror":
            return self._native_mirror_projection(pts)
        return self._native_refractive_projection(pts, row)

    def _native_mirror_projection(self, pts: np.ndarray):
        proj_x, proj_y = self._project_xy(pts[:, 2], pts[:, 1])
        proj_x = np.asarray(proj_x, dtype=float)
        proj_y = np.asarray(proj_y, dtype=float)
        if proj_x.size == 0:
            return None, None
        coords = np.column_stack((proj_x, proj_y))
        center = np.mean(coords, axis=0)
        centered = coords - center
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
            axis = vh[0]
        except np.linalg.LinAlgError:
            order = np.argsort(proj_x)
            return proj_x[order], proj_y[order]
        t = centered @ axis
        order = np.argsort(t)
        return coords[order, 0], coords[order, 1]

    def _native_refractive_projection(self, pts: np.ndarray, row: SurfaceRow):
        center_x = float(np.median(pts[:, 0]))
        tolerance = max(float(row.diameter) * 0.015, 0.08)
        mask = np.abs(pts[:, 0] - center_x) <= tolerance
        sliced = pts[mask]
        if sliced.shape[0] < 16:
            tolerance = max(float(row.diameter) * 0.03, 0.2)
            mask = np.abs(pts[:, 0] - center_x) <= tolerance
            sliced = pts[mask]
        if sliced.shape[0] < 8:
            stride = max(1, len(pts) // 300)
            sliced = pts[::stride]
        proj_x, proj_y = self._project_xy(sliced[:, 2], sliced[:, 1])
        proj_x = np.asarray(proj_x, dtype=float)
        proj_y = np.asarray(proj_y, dtype=float)
        if proj_x.size == 0:
            return None, None
        order = np.argsort(proj_x)
        return proj_x[order], proj_y[order]

    def _draw_reference_plane_labels(self) -> None:
        if not self.rows:
            return
        y0, y1 = self.ax.get_ylim()
        y_text = y1 - 0.08 * (y1 - y0)
        z_pos = 0.0
        for row in self.rows:
            if row.surface in {"Object", "Image"} and row.name:
                half_height = max(row.diameter / 2.0, 0.5)
                center_z = z_pos + float(row.desp_z)
                center_y = float(row.desp_y)
                angle = np.deg2rad(float(row.tilt_x))
                dz = np.cos(angle) * 0.0
                dy = np.sin(angle) * 0.0
                x_vals, y_vals = self._project_xy(
                    [center_z - dz, center_z + dz],
                    [center_y - half_height - dy, center_y + half_height + dy],
                )
                self.ax.plot(
                    x_vals,
                    y_vals,
                    color="#202020",
                    linewidth=1.2,
                    alpha=0.9,
                )
                text_x = float(x_vals[0])
                if self._current_display_orientation() == "Horizontal":
                    self.ax.text(
                        float(np.mean(x_vals)),
                        float(y_vals[0]),
                        row.name,
                        ha="center",
                        va="bottom",
                        fontsize=9,
                        color="#202020",
                    )
                else:
                    self.ax.text(
                        text_x,
                        y_text,
                        row.name,
                        ha="center",
                        va="top",
                        fontsize=9,
                        color="#202020",
                    )
            z_pos += row.thickness

    def _draw_custom_mirror_surfaces(self) -> None:
        z_pos = 0.0
        for row in self.rows:
            if row.surface == "Mirror":
                half_length = max(row.diameter / 2.0, 0.5)
                angle = np.deg2rad(float(row.tilt_x))
                dz = np.cos(angle) * half_length
                dy = np.sin(angle) * half_length
                center_z = z_pos + float(row.desp_z)
                center_y = float(row.desp_y)
                self.ax.plot(
                    [center_z - dz, center_z + dz],
                    [center_y - dy, center_y + dy],
                    color="#202020",
                    linewidth=2.2,
                    alpha=0.95,
                    solid_capstyle="round",
                )
            z_pos += row.thickness

    def _draw_optics_markers(self, optics_info: dict) -> None:
        if not self.show_cardinals_var.get():
            return
        y0, y1 = self.ax.get_ylim()
        span = y1 - y0
        y_top = y1 - 0.18 * span
        marker_specs = [
            ("Front PP", optics_info.get("ppa"), "#ff9f1c"),
            ("Back PP", optics_info.get("ppp"), "#ff9f1c"),
            ("EP", optics_info.get("ep_z"), "#00bcd4"),
            ("XP", optics_info.get("xp_z"), "#e91e63"),
        ]
        for label, z_pos, color in marker_specs:
            if z_pos is None:
                continue
            z_val = float(z_pos)
            if self._current_display_orientation() == "Horizontal":
                _, y_vals = self._project_xy([z_val, z_val], [0.0, 0.0])
                y_mark = float(y_vals[0])
                self.ax.axhline(y_mark, color=color, linewidth=1.0, linestyle=":", alpha=0.9)
                x0, x1 = self.ax.get_xlim()
                self.ax.text(
                    x0 + 0.04 * (x1 - x0),
                    y_mark,
                    label,
                    color=color,
                    fontsize=8,
                    ha="left",
                    va="bottom",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.6},
                )
            else:
                self.ax.axvline(z_val, color=color, linewidth=1.0, linestyle=":", alpha=0.9)
                self.ax.text(
                    z_val,
                    y_top,
                    label,
                    color=color,
                    fontsize=8,
                    ha="center",
                    va="top",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.6},
                )

    def _set_plot_limits_from_layout(self, max_radius: float) -> None:
        total_length = sum(float(row.thickness) for row in self.rows)
        total_length = max(total_length, 1.0)
        margin_x = max(total_length * 0.05, 5.0)
        margin_y = max(max_radius * 0.2, 2.0)
        self.ax.set_xlim(-margin_x, total_length + margin_x)
        self.ax.set_ylim(-(max_radius + margin_y), max_radius + margin_y)

    def _set_plot_limits_from_drawn_data(self) -> None:
        x_values: list[float] = []
        y_values: list[float] = []
        for line in self.ax.lines:
            xdata = np.asarray(line.get_xdata(orig=False), dtype=float)
            ydata = np.asarray(line.get_ydata(orig=False), dtype=float)
            finite = np.isfinite(xdata) & np.isfinite(ydata)
            if np.any(finite):
                x_values.extend(xdata[finite].tolist())
                y_values.extend(ydata[finite].tolist())
        if not x_values or not y_values:
            return
        x_min = min(x_values)
        x_max = max(x_values)
        y_min = min(y_values)
        y_max = max(y_values)
        span_x = max(x_max - x_min, 1.0)
        span_y = max(y_max - y_min, 1.0)
        self.ax.set_xlim(x_min - 0.08 * span_x, x_max + 0.08 * span_x)
        self.ax.set_ylim(y_min - 0.12 * span_y, y_max + 0.12 * span_y)

    def _draw_input_ray_overlay(self, max_radius: float) -> None:
        if not self.rows:
            return
        if self._current_object_mode() == "Infinity":
            return
        object_distance = self._current_object_distance()
        if object_distance <= 1e-9:
            return
        field_samples = self._sample_field_values(self._current_field_height())
        angle_samples = self._sample_fan_angles_deg()
        colors = self._field_colors(len(field_samples))
        for field_index, field_height in enumerate(field_samples):
            color = colors[field_index]
            for angle_deg in angle_samples:
                angle_rad = np.deg2rad(angle_deg)
                pupil_y = float(field_height) + float(np.tan(angle_rad) * object_distance)
                x_vals, y_vals = self._project_xy([0.0, object_distance], [float(field_height), float(pupil_y)])
                self.ax.plot(
                    x_vals,
                    y_vals,
                    color=color,
                    linewidth=1.8,
                    alpha=0.95,
                )

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
        measure = tkfont.nametofont("TkDefaultFont").measure
        property_width = measure("Property") + 18
        for key, value in items:
            self.results_table.insert("", "end", values=(key, value))
            property_width = max(property_width, measure(str(key)) + 18)
        self.results_table.column("property", width=min(property_width, 150), anchor="w", stretch=False)

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

    def _begin_analysis_progress(self, label: str) -> None:
        if self.optimization_running:
            return
        self.progress_spinner_var.set("...")
        self.progress_percent_var.set("working")
        self.progress_bar_var.set(0.0)
        self.append_progress(f"{label} started.")

    def _update_analysis_progress(self, label: str, done: int | None = None, total: int | None = None) -> None:
        if self.optimization_running:
            return
        frames = ("|", "/", "-", "\\")
        self.progress_spinner_var.set(frames[self._spinner_phase % len(frames)])
        self._spinner_phase += 1
        if done is not None and total is not None and total > 0:
            percent = max(0.0, min(100.0, (done / total) * 100.0))
            self.progress_percent_var.set(f"{int(percent)}% ({done}/{total})")
            self.progress_bar_var.set(percent)
        else:
            self.progress_percent_var.set(label)
        self.update_idletasks()

    def _finish_analysis_progress(self, label: str, success: bool = True) -> None:
        if self.optimization_running:
            return
        self.progress_spinner_var.set("ok" if success else "err")
        self.progress_percent_var.set("100%" if success else "failed")
        self.progress_bar_var.set(100.0 if success else 0.0)
        self.append_progress(f"{label} {'completed' if success else 'failed'}.")

    def _update_progress_indicators(self) -> None:
        if not self.optimization_running or self.optimization_context is None:
            self.progress_spinner_var.set("idle")
            self.progress_percent_var.set("0%")
            self.progress_bar_var.set(0.0)
            return
        done = int(self.optimization_context.get("generation_done", 0))
        total = max(1, int(self.optimization_context.get("generations_total", 1)))
        percent = max(0.0, min(100.0, (done / total) * 100.0))
        self.progress_percent_var.set(f"{int(percent)}% ({done}/{total})")
        self.progress_bar_var.set(percent)

    @staticmethod
    def _pick_image_plane_data(rays):
        try:
            X, Y, Z, L, M, N = rays.pick(-1, coordinates="local")
            if np.asarray(X).size:
                return X, Y, Z, L, M, N
        except Exception:
            pass
        return rays.pick(-1)

    def _build_analysis_rays(
        self,
        system,
        wavelength: float,
        sample_count: int | None = None,
        pattern: str = "hexapolar",
        *,
        surface_index: int | None = None,
        aperture_type: str | None = None,
        aperture_value: float | None = None,
        field_type: str | None = None,
        field_x: float | None = None,
        field_y: float | None = None,
    ):
        rays = Kos.raykeeper(system)
        pupil = Kos.PupilCalc(
            system,
            self._analysis_surface_index() if surface_index is None else int(surface_index),
            wavelength,
            self._current_aperture_type() if aperture_type is None else str(aperture_type),
            self._current_aperture_value() if aperture_value is None else float(aperture_value),
        )
        pupil.Samp = max(2, int(sample_count if sample_count is not None else self._current_ray_count()))
        pupil.Ptype = pattern

        clean = 1
        resolved_field_type = field_type or ("angle" if self._current_object_mode() == "Infinity" else "height")
        resolved_field_x = 0.0 if field_x is None else float(field_x)
        resolved_field_y = field_y
        if resolved_field_type == "angle":
            pupil.FieldType = "angle"
            field_values = (
                [float(resolved_field_y)]
                if resolved_field_y is not None
                else self._sample_field_values(self._current_field_angle_deg())
            )
            for value in field_values:
                pupil.FieldX = resolved_field_x
                pupil.FieldY = float(value)
                x, y, z, L, M, N = pupil.Pattern2Field()
                Kos.TraceLoop(x, y, z, L, M, N, wavelength, rays, clean=clean)
                clean = 0
        else:
            pupil.FieldType = "height"
            field_values = (
                [float(resolved_field_y)]
                if resolved_field_y is not None
                else self._sample_field_values(self._current_field_height())
            )
            for value in field_values:
                pupil.FieldX = resolved_field_x
                pupil.FieldY = float(value)
                x, y, z, L, M, N = pupil.Pattern2Field()
                Kos.TraceLoop(x, y, z, L, M, N, wavelength, rays, clean=clean)
                clean = 0
        return rays

    def _collect_optics_info(self, system, rays, wavelength: float) -> dict:
        info: dict[str, float | None | str] = {
            "effl": None,
            "magnification": None,
            "ppa": None,
            "ppp": None,
            "spot_rms": None,
            "spot_cen_x": None,
            "spot_cen_y": None,
            "ep_radius": None,
            "ep_z": None,
            "xp_radius": None,
            "xp_z": None,
            "airy_radius": None,
        }
        try:
            _, _, _, a, b, c, d, effl, ppa, ppp, _, _, _ = system.Parax(wavelength)
            info.update(
                {
                    "magnification": float(a),
                    "effl": float(effl),
                    "ppa": float(ppa),
                    "ppp": float(ppp),
                    "parax_a": float(a),
                    "parax_b": float(b),
                    "parax_c": float(c),
                    "parax_d": float(d),
                }
            )
        except Exception:
            pass
        try:
            analysis_rays = self._build_analysis_rays(system, wavelength)
            X, Y, Z, L, M, N = self._pick_image_plane_data(analysis_rays)
            if X.size:
                rms, cenX, cenY = Kos.RMS(X, Y, Z, L, M, N)
                info.update(
                    {
                        "spot_rms": float(rms),
                        "spot_cen_x": float(cenX),
                        "spot_cen_y": float(cenY),
                    }
                )
        except Exception:
            pass
        try:
            pupil = Kos.PupilCalc(
                system,
                self._analysis_surface_index(),
                wavelength,
                self._current_aperture_type(),
                self._current_aperture_value(),
            )
            info.update(
                {
                    "ep_radius": float(pupil.RadPupInp),
                    "ep_z": float(pupil.PosPupInp[2]),
                    "xp_radius": float(pupil.RadPupOut),
                    "xp_z": float(pupil.PosPupOut[2]),
                    "airy_radius": float(pupil.FocusAiryRadius),
                }
            )
        except Exception:
            pass
        return info

    def _start_progress_spinner(self) -> None:
        if self._spinner_after_id is not None:
            self.after_cancel(self._spinner_after_id)
            self._spinner_after_id = None
        self._spinner_phase = 0
        self._animate_progress_spinner()

    def _stop_progress_spinner(self) -> None:
        if self._spinner_after_id is not None:
            self.after_cancel(self._spinner_after_id)
            self._spinner_after_id = None
        if not self.optimization_running:
            self.progress_spinner_var.set("idle")

    def _animate_progress_spinner(self) -> None:
        if not self.optimization_running:
            self._spinner_after_id = None
            self.progress_spinner_var.set("idle")
            return
        frames = ("|", "/", "-", "\\")
        self.progress_spinner_var.set(frames[self._spinner_phase % len(frames)])
        self._spinner_phase += 1
        self._spinner_after_id = self.after(120, self._animate_progress_spinner)

    def _update_results(self, system, rays, wavelength: float, optics_info: dict | None = None) -> None:
        optics_info = optics_info or self._collect_optics_info(system, rays, wavelength)
        items = []
        items.append(("Surface count", str(len(self.rows))))
        items.append(("Optimized vars", str(len(self._build_optimization_variables()))))
        items.append(("Object mode", self._current_object_mode()))
        items.append(("Wavelength [um]", f"{wavelength:.4g}"))
        items.append(("Analysis surface", str(self._analysis_surface_index())))
        items.append(("Aperture type", self._current_aperture_type()))
        items.append(("Aperture value", f"{self._current_aperture_value():.4g}"))
        field_metrics = self._field_metrics()
        items.append(("Field type", self._current_field_type()))
        items.append(("Field angle [deg]", f"{field_metrics['angle_deg']:.4g}"))
        items.append(("Object height [mm]", f"{field_metrics['object_height']:.4g}"))
        items.append(("Paraxial image height [mm]", f"{field_metrics['paraxial_image_height']:.4g}"))
        items.append(("Real image height [mm]", f"{field_metrics['real_image_height']:.4g}"))

        total_length = sum(max(float(row.thickness), 0.0) for row in self.rows)
        items.append(("Total length [mm]", f"{total_length:.4g}"))

        if optics_info.get("effl") is not None:
            items.append(("Imaging", ""))
            items.append(("EFFL [mm]", f"{float(optics_info['effl']):.4g}"))
            items.append(("Magnification", f"{float(optics_info['magnification']):.4g}"))
            items.append(("Principal Planes", ""))
            items.append(("Front principal plane [mm]", f"{float(optics_info['ppa']):.4g}"))
            items.append(("Back principal plane [mm]", f"{float(optics_info['ppp']):.4g}"))
        else:
            items.append(("Paraxial data", "Unavailable"))

        if optics_info.get("ep_radius") is not None:
            items.append(("Pupils", ""))
            items.append(("Entrance pupil radius [mm]", f"{float(optics_info['ep_radius']):.4g}"))
            items.append(("Entrance pupil diameter [mm]", f"{2.0 * float(optics_info['ep_radius']):.4g}"))
            items.append(("Entrance pupil z [mm]", f"{float(optics_info['ep_z']):.4g}"))
            items.append(("Exit pupil radius [mm]", f"{float(optics_info['xp_radius']):.4g}"))
            items.append(("Exit pupil diameter [mm]", f"{2.0 * float(optics_info['xp_radius']):.4g}"))
            items.append(("Exit pupil z [mm]", f"{float(optics_info['xp_z']):.4g}"))
            items.append(("Airy radius [mm]", f"{float(optics_info['airy_radius']):.4g}"))
        else:
            items.append(("Pupil data", "Unavailable"))

        items.append(("Spot", ""))
        if optics_info.get("spot_rms") is not None:
            items.append(("Spot RMS [mm]", f"{float(optics_info['spot_rms']):.4g}"))
            items.append(("Spot centroid X [mm]", f"{float(optics_info['spot_cen_x']):.4g}"))
            items.append(("Spot centroid Y [mm]", f"{float(optics_info['spot_cen_y']):.4g}"))
        else:
            items.append(("Spot RMS [mm]", "Unavailable"))

        self._set_results(items)

    @staticmethod
    def _optimization_bounds(parameter: str, value: float) -> tuple[float, float]:
        for spec in VARIABLE_REGISTRY.values():
            if spec.parameter == parameter:
                return spec.default_bounds(value)
        raise ValueError(f"Unsupported optimization parameter: {parameter}")

    @staticmethod
    def _variable_spec_for_field(field: str):
        return VARIABLE_REGISTRY.get(field)

    @staticmethod
    def _merit_spec_for_label(label: str):
        for spec in OPERAND_REGISTRY.values():
            if spec.label == label:
                return spec
        return None

    def _selected_operand_specs(self) -> list:
        if not hasattr(self, "merit_mode_list"):
            return []
        labels = [self.merit_mode_list.get(i) for i in self.merit_mode_list.curselection()]
        specs = []
        for label in labels:
            spec = self._merit_spec_for_label(label)
            if spec is not None:
                specs.append(spec)
        return specs

    def _operand_weight(self, label: str) -> float:
        var = self.operand_weight_vars.get(label)
        if var is None:
            spec = self._merit_spec_for_label(label)
            return 1.0 if spec is None else spec.default_weight
        try:
            return float(var.get())
        except ValueError:
            spec = self._merit_spec_for_label(label)
            return 1.0 if spec is None else spec.default_weight

    def _operand_target(self, label: str) -> float:
        var = self.operand_target_vars.get(label)
        if var is None:
            spec = self._merit_spec_for_label(label)
            return 0.0 if spec is None else spec.default_target
        try:
            return float(var.get())
        except ValueError:
            spec = self._merit_spec_for_label(label)
            return 0.0 if spec is None else spec.default_target

    def _operand_wavelength(self, label: str) -> float:
        var = self.operand_wavelength_vars.get(label)
        if var is None:
            return self._current_wavelength()
        try:
            return float(var.get())
        except ValueError:
            return self._current_wavelength()

    def _operand_field(self, label: str) -> float:
        var = self.operand_field_vars.get(label)
        if var is None:
            return 0.0
        try:
            return float(var.get())
        except ValueError:
            return 0.0

    def _operand_field_x(self, label: str) -> float:
        var = self.operand_field_x_vars.get(label)
        if var is None:
            return 0.0
        try:
            return float(var.get())
        except ValueError:
            return 0.0

    def _operand_field_y(self, label: str) -> float:
        var = self.operand_field_y_vars.get(label)
        if var is not None:
            try:
                return float(var.get())
            except ValueError:
                return 0.0
        return self._operand_field(label)

    def _operand_field_type(self, label: str) -> str:
        if self._current_field_type() == "Angle":
            return "angle"
        return "height"

    def _operand_surface_index(self, label: str) -> int:
        var = self.operand_surface_vars.get(label)
        if var is None:
            return self._analysis_surface_index()
        value = var.get().strip()
        if not value or value == "Auto":
            return self._analysis_surface_index()
        try:
            return int(value.split(":", 1)[0].strip())
        except ValueError:
            return self._analysis_surface_index()

    def _operand_aperture_type(self, label: str) -> str:
        var = self.operand_aperture_type_vars.get(label)
        if var is None:
            return self._current_aperture_type()
        value = var.get().strip().upper()
        return value if value in {"STOP", "EPD"} else self._current_aperture_type()

    def _operand_aperture_value(self, label: str) -> float:
        var = self.operand_aperture_value_vars.get(label)
        if var is None:
            return self._current_aperture_value()
        try:
            value = float(var.get())
        except ValueError:
            return self._current_aperture_value()
        return value if value != 0.0 else self._current_aperture_value()

    def _build_optimization_variables(self) -> list[OpticalVariable]:
        variables: list[OpticalVariable] = []
        for index, row in enumerate(self.rows):
            if row.surface in {"Object", "Image"}:
                continue
            for spec in VARIABLE_REGISTRY.values():
                if not spec.is_supported(row) or not spec.is_enabled(row):
                    continue
                value = spec.value_from_row(row)
                lower, upper = spec.get_bounds(row) or spec.default_bounds(value)
                variables.append(
                    OpticalVariable(
                        index,
                        spec.parameter,
                        lower,
                        upper,
                        name=f"{row.name} {spec.label}",
                    )
                )
        return variables

    def _build_merit_function(self) -> MeritFunction:
        selected_specs = self._selected_operand_specs()
        if not selected_specs:
            return MeritFunction(operands=[])
        operands = []
        for spec in selected_specs:
            operands.extend(spec.build_merit_function(self).operands)
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

        merit_specs = self._selected_operand_specs()
        merit = self._build_merit_function()
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
        self.append_progress(
            "Optimization start | operands: "
            + ", ".join(spec.label for spec in merit_specs)
        )
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
        self._update_progress_indicators()
        self._start_progress_spinner()
        self.after(0, self._optimization_step)

    def stop_optimization(self) -> None:
        if not self.optimization_running:
            self.append_progress("Stop ignored: no optimization is running.")
            return
        self.optimization_cancel_requested = True
        self.append_progress("Stop requested. Optimization will stop after the current generation.")
        self._update_progress_indicators()

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
        self._update_progress_indicators()
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
            self._stop_progress_spinner()
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
        self._stop_progress_spinner()
        self._update_progress_indicators()

    def _sample_ray_heights(self, max_radius: float) -> list[float]:
        if max_radius <= 1e-9:
            return [0.0]
        count = self._current_ray_count()
        span = max_radius * self._current_ray_height_factor()
        if count == 1:
            return [0.0]
        return list(np.linspace(-span, span, count))

    def _sample_field_values(self, maximum: float) -> list[float]:
        count = self._current_field_count()
        if count == 1 or maximum <= 1e-9:
            return [0.0]
        return list(np.linspace(-maximum, maximum, count))

    def _sample_fan_angles_deg(self) -> list[float]:
        maximum = self._current_field_angle_deg()
        count = self._current_ray_count()
        if count == 1 or maximum <= 1e-9:
            return [0.0]
        return list(np.linspace(-maximum, maximum, count))

    def _entrance_radius(self, fallback_radius: float) -> float:
        object_radius = None
        if self.rows:
            object_radius = max(float(self.rows[0].diameter) / 2.0, 0.5)
        for row in self.rows[1:]:
            if row.surface not in {"Object", "Image"}:
                radius = max(row.diameter / 2.0, 0.5)
                if object_radius is not None:
                    return min(radius, object_radius)
                return radius
        if object_radius is not None:
            return min(fallback_radius, object_radius)
        return fallback_radius

    def _trace_preview_rays(self, system, rays, wavelength: float, max_radius: float) -> None:
        system.IgnoreVignetting(0)
        if self._has_off_axis_geometry():
            rays.clean()
            if self._current_object_mode() == "Infinity":
                field_values = self._sample_field_values(self._current_field_angle_deg())
                pupil_samples = self._sample_ray_heights(self._entrance_radius(max_radius))
                for field_angle in field_values:
                    angle_rad = np.deg2rad(float(field_angle))
                    direction = np.array([0.0, np.sin(angle_rad), np.cos(angle_rad)], dtype=float)
                    norm = np.linalg.norm(direction)
                    if norm <= 1e-12:
                        continue
                    direction /= norm
                    for pupil_y in pupil_samples:
                        origin = [0.0, float(pupil_y), 0.0]
                        system.Trace(origin, direction.tolist(), wavelength)
                        rays.push()
                self._preview_field_ray_count = len(pupil_samples)
            else:
                field_values = self._sample_field_values(self._current_field_height())
                pupil_samples = self._sample_ray_heights(self._entrance_radius(max_radius))
                object_distance = self._current_object_distance()
                for field_value in field_values:
                    origin = np.array([0.0, float(field_value), 0.0], dtype=float)
                    for pupil_y in pupil_samples:
                        target = np.array([0.0, float(pupil_y), object_distance], dtype=float)
                        direction = target - origin
                        norm = np.linalg.norm(direction)
                        if norm <= 1e-12:
                            continue
                        direction /= norm
                        system.Trace(origin.tolist(), direction.tolist(), wavelength)
                        rays.push()
                self._preview_field_ray_count = len(pupil_samples)
        elif self._current_object_mode() == "Infinity":
            pupil = Kos.PupilCalc(
                system,
                self._analysis_surface_index(),
                wavelength,
                self._current_aperture_type(),
                self._current_aperture_value(),
            )
            pupil.Samp = max(2, self._current_ray_count())
            pupil.Ptype = "fany"
            clean = 1
            last_bundle = 1
            pupil.FieldType = "angle"
            field_values = self._sample_field_values(self._current_field_angle_deg())
            for field_value in field_values:
                pupil.FieldX = 0.0
                pupil.FieldY = float(field_value)
                x, y, z, L, M, N = pupil.Pattern2Field()
                last_bundle = max(1, len(np.asarray(x)))
                Kos.TraceLoop(x, y, z, L, M, N, wavelength, rays, clean=clean)
                clean = 0
            self._preview_field_ray_count = last_bundle
        else:
            rays.clean()
            field_values = self._sample_field_values(self._current_field_height())
            pupil_samples = self._sample_ray_heights(self._entrance_radius(max_radius))
            object_distance = self._current_object_distance()
            for field_value in field_values:
                origin = np.array([0.0, float(field_value), 0.0], dtype=float)
                for pupil_y in pupil_samples:
                    target = np.array([0.0, float(pupil_y), object_distance], dtype=float)
                    direction = target - origin
                    norm = np.linalg.norm(direction)
                    if norm <= 1e-12:
                        continue
                    direction /= norm
                    system.Trace(origin.tolist(), direction.tolist(), wavelength)
                    rays.push()
            self._preview_field_ray_count = len(pupil_samples)
        system.Vignetting(0)

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
            x_vals, y_vals = self._project_xy([z, z], [-radius, radius])
            self.ax.plot(x_vals, y_vals, color=color, linewidth=2)
            self.ax.text(
                float(x_vals[0]),
                float(np.max(y_vals) + max_radius * 0.08),
                row.name,
                rotation=45,
                ha="left",
                va="bottom",
                fontsize=8,
            )
            z += row.thickness

        total_length = max(z, 1.0)
        margin = max(total_length * 0.05, 5.0)
        if self._current_display_orientation() == "Horizontal":
            self._set_plot_limits_from_drawn_data()
        else:
            self.ax.set_xlim(-margin, total_length + margin)
            self.ax.set_ylim(-(max_radius * 1.4), max_radius * 1.4)
        axis_x, axis_y = self._project_xy([0.0, total_length], [0.0, 0.0])
        self.ax.plot(axis_x, axis_y, color="#2c3e50", linewidth=0.8)

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
                    f"    {var_name}.TiltX = {float(row.tilt_x)!r}",
                    f"    {var_name}.TiltY = {float(row.tilt_y)!r}",
                    f"    {var_name}.TiltZ = {float(row.tilt_z)!r}",
                    f"    {var_name}.DespX = {float(row.desp_x)!r}",
                    f"    {var_name}.DespY = {float(row.desp_y)!r}",
                    f"    {var_name}.DespZ = {float(row.desp_z)!r}",
                    f"    {var_name}.AxisMove = {float(row.axis_move)!r}",
                    f"    {var_name}.Glass = {row.glass!r}",
                ]
            )
            if row.surface == "Mirror":
                lines.append(f"    {var_name}.Glass = 'MIRROR'")
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
                f"'tilt_x': {float(row.tilt_x)!r}, "
                f"'tilt_y': {float(row.tilt_y)!r}, "
                f"'tilt_z': {float(row.tilt_z)!r}, "
                f"'desp_x': {float(row.desp_x)!r}, "
                f"'desp_y': {float(row.desp_y)!r}, "
                f"'desp_z': {float(row.desp_z)!r}, "
                f"'axis_move': {float(row.axis_move)!r}, "
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
                "    clear_aperture = max((max(float(spec['diameter']), 1.0) for spec in surface_dicts if spec['surface'] not in {'Object', 'Image'}), default=100.0) * 4.0",
                "    for spec in surface_dicts:",
                "        s = Kos.surf()",
                "        s.Name = spec['name']",
                "        s.Rc = spec['rc']",
                "        s.Thickness = spec['thickness']",
                "        s.Diameter = clear_aperture if spec['surface'] in {'Object', 'Image'} else spec['diameter']",
                "        s.TiltX = spec.get('tilt_x', 0.0)",
                "        s.TiltY = spec.get('tilt_y', 0.0)",
                "        s.TiltZ = spec.get('tilt_z', 0.0)",
                "        s.DespX = spec.get('desp_x', 0.0)",
                "        s.DespY = spec.get('desp_y', 0.0)",
                "        s.DespZ = spec.get('desp_z', 0.0)",
                "        s.AxisMove = spec.get('axis_move', 0.0)",
                "        s.Glass = spec['glass']",
                "        if spec['surface'] == 'Mirror':",
                "            s.Glass = 'MIRROR'",
                "        s.Drawing = 0.0 if spec['surface'] in {'Object', 'Image'} else 1.0",
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
        elif str(getattr(surface, "Glass", "AIR")).upper() == "MIRROR":
            surface_type = "Mirror"

        rc_value = float(getattr(surface, "Rc", 0.0))
        if surface_type == "Thin Lens":
            rc_value = float(getattr(surface, "Thin_Lens", 0.0))

        return SurfaceRow(
            surface=surface_type,
            name=str(getattr(surface, "Name", "") or f"Surface {index}"),
            rc=rc_value,
            thickness=float(getattr(surface, "Thickness", 0.0)),
            diameter=float(getattr(surface, "Diameter", 25.0)),
            tilt_x=float(getattr(surface, "TiltX", 0.0)),
            tilt_y=float(getattr(surface, "TiltY", 0.0)),
            tilt_z=float(getattr(surface, "TiltZ", 0.0)),
            desp_x=float(getattr(surface, "DespX", 0.0)),
            desp_y=float(getattr(surface, "DespY", 0.0)),
            desp_z=float(getattr(surface, "DespZ", 0.0)),
            axis_move=float(getattr(surface, "AxisMove", 0.0)),
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
        glass = str(item.get("glass", "AIR")).strip().upper()
        if glass == "MIRROR":
            return "Mirror"
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
