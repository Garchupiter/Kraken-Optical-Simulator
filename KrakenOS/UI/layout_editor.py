"""Simple KrakenOS layout editor.

This is an initial editor scaffold that mirrors the RayTracing workflow:
- file-backed starter layouts
- editable surface table
- embedded axial sketch with a small traced ray fan
"""

from __future__ import annotations

import importlib.util
import io
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stderr, redirect_stdout
import ctypes
from dataclasses import dataclass, asdict
import html
import multiprocessing as mp
import os
from pathlib import Path
from pprint import pformat
from queue import Empty
import re
import signal
import shutil
import subprocess
import sys
import time
import traceback
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox, ttk
import warnings
import webbrowser

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox
import numpy as np

import KrakenOS as Kos
import pyvista as pv
from KrakenOS.Display import Plot2DRays, Plot2DSurf, edge_3d, filter_face_2dplot, wavelength_to_rgb
from KrakenOS.Optimization import (
    OPERAND_REGISTRY,
    VARIABLE_REGISTRY,
    MeritEvaluator,
    MeritFunction,
    MTFAtFrequencyOperand,
    OpticalVariable,
)
from KrakenOS.Optimization.adapters.pygmo2_adapter import Pygmo2MeritProblem

try:
    from vtkmodules.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
except Exception:
    vtkTkRenderWindowInteractor = None

try:
    from vtkmodules.vtkFiltersCore import vtkTubeFilter
    from vtkmodules.vtkInteractionWidgets import vtkOrientationMarkerWidget
    from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
    from vtkmodules.vtkRenderingCore import vtkActor, vtkCellPicker, vtkDataSetMapper, vtkRenderer
except Exception:
    vtkTubeFilter = None
    vtkOrientationMarkerWidget = None
    vtkAxesActor = None
    vtkActor = None
    vtkCellPicker = None
    vtkDataSetMapper = None
    vtkRenderer = None


LAYOUTS_DIR = Path(__file__).resolve().parent.parent / "common_optical_layouts"
EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "Examples"
DEFAULT_LAYOUT_TITLE = "Doublet Lens"
FOLDED_STARTER_LAYOUT_TITLE = "Double Mirror Fold"
AUTO_PLOT_PATH = Path.home() / "Pictures" / "kraken_layout_latest.jpg"
DEBUG_LOG_PATH = Path.home() / "Pictures" / "kraken_debug_latest.log"
FIELDS = (
    "label",
    "surface",
    "name",
    "glass",
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
)
COLUMN_LABELS = {
    "label": "#",
    "surface": "Surface",
    "name": "Name",
    "glass": "Material",
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
    settings = getattr(module, "SETTINGS", {})
    if not isinstance(settings, dict):
        settings = {}
    return {"title": title, "surfaces": surfaces, "settings": settings}


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


_CUPY_IMPORT_ATTEMPTED = False
_CUPY_MODULE = None
_TORCH_IMPORT_ATTEMPTED = False
_TORCH_MODULE = None
_CUDA_LIBS_PRELOADED = False
_WORKER_SYSTEM_CACHE_SIGNATURE = None
_WORKER_SYSTEM_CACHE_SYSTEM = None


def _preload_cuda_libraries():
    global _CUDA_LIBS_PRELOADED
    if _CUDA_LIBS_PRELOADED:
        return
    _CUDA_LIBS_PRELOADED = True

    driver_candidates = (
        "/run/opengl-driver/lib/libcuda.so.1",
        "/run/opengl-driver/lib/libcuda.so",
        "/run/opengl-driver-32/lib/libcuda.so.1",
        "/run/opengl-driver-32/lib/libcuda.so",
    )
    for candidate in driver_candidates:
        if not os.path.exists(candidate):
            continue
        try:
            ctypes.CDLL(candidate, mode=ctypes.RTLD_GLOBAL)
            break
        except Exception:
            continue

    # Best-effort preload of CUDA runtime/NVRTC from pip wheels.
    package_libs = (
        ("nvidia.cuda_nvrtc", ("libnvrtc.so",)),
        ("nvidia.cuda_runtime", ("libcudart.so",)),
        ("nvidia.cu13", ("libnvrtc-builtins.so", "libnvrtc.so", "libcudart.so", "libcufft.so")),
    )
    for module_name, lib_prefixes in package_libs:
        try:
            spec = importlib.util.find_spec(module_name)
            if spec is None or not spec.submodule_search_locations:
                continue
            for search_path in spec.submodule_search_locations:
                lib_dir = Path(search_path) / "lib"
                if not lib_dir.exists():
                    continue
                for prefix in lib_prefixes:
                    for lib_path in sorted(lib_dir.glob(f"{prefix}*")):
                        try:
                            ctypes.CDLL(str(lib_path), mode=ctypes.RTLD_GLOBAL)
                            break
                        except Exception:
                            continue
        except Exception:
            continue


def _short_error_message(exc: Exception, limit: int = 220) -> str:
    text = str(exc).strip()
    if not text:
        return exc.__class__.__name__
    first = text.splitlines()[0].strip()
    if len(first) > limit:
        return first[:limit] + "..."
    return first


class Kraken3DInspector(tk.Toplevel):
    def __init__(self, editor: "KrakenLayoutEditor") -> None:
        super().__init__(editor)
        self.editor = editor
        self.available = False
        self.unavailable_reason = ""
        self.title("KrakenOS 3D Inspector")
        self.geometry("1100x780")
        self.minsize(720, 520)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        self._renderer = None
        self._vtk_widget = None
        self._vtk_interactor = None
        self._orientation_widget = None
        self._picker = None
        self._picked_row_index: int | None = None
        self._actor_row_map: dict[str, int] = {}
        self._row_actor_map: dict[int, list[str]] = {}
        self._camera_preset = "iso"
        self.show_rays_var = tk.BooleanVar(value=True)
        self.status_var = tk.StringVar(value="3D inspector ready")

        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        host = ttk.Frame(self, padding=8)
        host.grid(row=1, column=0, sticky="nsew")
        host.columnconfigure(0, weight=1)
        host.rowconfigure(0, weight=1)

        if vtkTkRenderWindowInteractor is None or vtkRenderer is None:
            self.unavailable_reason = "Embedded VTK/Tk viewer unavailable."
            self.status_var.set(self.unavailable_reason)
            return

        try:
            toolbar = ttk.Frame(self, padding=(8, 8, 8, 0))
            toolbar.grid(row=0, column=0, sticky="ew")
            ttk.Button(toolbar, text="Refresh", command=self.refresh_from_editor).pack(side="left")
            ttk.Button(toolbar, text="Iso", command=lambda: self.set_camera_preset("iso")).pack(side="left", padx=(8, 0))
            ttk.Button(toolbar, text="ZY", command=lambda: self.set_camera_preset("zy")).pack(side="left", padx=(4, 0))
            ttk.Button(toolbar, text="XY", command=lambda: self.set_camera_preset("xy")).pack(side="left", padx=(4, 0))
            ttk.Button(toolbar, text="XZ", command=lambda: self.set_camera_preset("xz")).pack(side="left", padx=(4, 0))
            ttk.Checkbutton(
                toolbar,
                text="Show rays",
                variable=self.show_rays_var,
                command=self.refresh_from_editor,
            ).pack(side="left", padx=(12, 0))
            ttk.Label(
                toolbar,
                text="Click a surface in 3D to select its table row",
                foreground="#4b5563",
            ).pack(side="right")

            self._vtk_widget = vtkTkRenderWindowInteractor(host, width=1100, height=720)
            self._vtk_widget.grid(row=0, column=0, sticky="nsew")
            render_window = self._vtk_widget.GetRenderWindow()
            self._renderer = vtkRenderer()
            render_window.AddRenderer(self._renderer)
            self._renderer.SetBackground(1.0, 1.0, 1.0)

            self._vtk_interactor = render_window.GetInteractor()
            if self._vtk_interactor is not None:
                self._vtk_interactor.AddObserver("LeftButtonPressEvent", self._on_left_button_press)

            if vtkCellPicker is not None:
                self._picker = vtkCellPicker()
                self._picker.SetTolerance(0.0005)

            if vtkOrientationMarkerWidget is not None and vtkAxesActor is not None and self._vtk_interactor is not None:
                axes = vtkAxesActor()
                self._orientation_widget = vtkOrientationMarkerWidget()
                self._orientation_widget.SetOrientationMarker(axes)
                self._orientation_widget.SetInteractor(self._vtk_interactor)
                self._orientation_widget.SetViewport(0.0, 0.0, 0.18, 0.18)
                self._orientation_widget.SetEnabled(1)
                self._orientation_widget.InteractiveOff()

            self._vtk_widget.Initialize()
            ttk.Label(self, textvariable=self.status_var, padding=(8, 0, 8, 8)).grid(row=2, column=0, sticky="ew")
            self.available = True
        except Exception as exc:
            self.unavailable_reason = _short_error_message(exc)
            self.status_var.set(f"Embedded 3D unavailable: {self.unavailable_reason}")
            try:
                host.destroy()
            except Exception:
                pass

    @staticmethod
    def _actor_key(actor) -> str | None:
        if actor is None:
            return None
        try:
            return str(actor.GetAddressAsString(""))
        except Exception:
            return str(id(actor))

    @staticmethod
    def _surface_color(surface) -> tuple[float, float, float]:
        absorb_color = (10 / 256.0, 23 / 256.0, 24 / 256.0)
        mirror_color = (189 / 256.0, 189 / 256.0, 189 / 256.0)
        glass_color = (12 / 256.0, 238 / 256.0, 246 / 256.0)
        try:
            color = tuple(float(v) for v in surface.Color)
            if len(color) == 3 and any(abs(v) > 1e-9 for v in color):
                return color
        except Exception:
            pass
        glass = str(getattr(surface, "Glass", "") or "").upper()
        if glass == "MIRROR":
            return mirror_color
        if glass == "ABSORB":
            return absorb_color
        return glass_color

    @staticmethod
    def _mesh_with_transform(poly, transform) -> pv.DataSet | None:
        try:
            mesh = pv.wrap(poly)
        except Exception:
            return None
        try:
            mesh = mesh.extract_surface()
        except Exception:
            pass
        try:
            mesh = mesh.copy(deep=True)
        except Exception:
            return None
        try:
            pts = np.asarray(mesh.points, dtype=float)
        except Exception:
            return None
        if pts.size == 0:
            return None
        try:
            transform_arr = np.asarray(transform, dtype=float)
            if transform_arr.shape == (4, 4):
                pts_h = np.c_[pts, np.ones(len(pts))]
                mesh.points = (pts_h @ transform_arr.T)[:, :3]
        except Exception:
            pass
        return mesh

    def _set_row_highlight(self, row_index: int | None) -> None:
        if row_index == self._picked_row_index:
            return
        if self._renderer is None:
            self._picked_row_index = row_index
            return
        collection = self._renderer.GetActors()
        collection.InitTraversal()
        for _ in range(collection.GetNumberOfItems()):
            actor = collection.GetNextActor()
            key = self._actor_key(actor)
            prop = actor.GetProperty()
            if prop is None:
                continue
            actor_row_index = self._actor_row_map.get(key) if key is not None else None
            if row_index is not None and actor_row_index == row_index:
                prop.SetEdgeVisibility(1)
                prop.SetEdgeColor(1.0, 0.55, 0.05)
                prop.SetLineWidth(2.0)
            else:
                prop.SetEdgeVisibility(0)
                prop.SetLineWidth(1.0)
        self._picked_row_index = row_index

    def highlight_row(self, row_index: int | None) -> None:
        self._set_row_highlight(row_index)
        self.render()

    def _add_mesh_actor(
        self,
        mesh,
        *,
        color: tuple[float, float, float],
        opacity: float = 1.0,
        pick_row_index: int | None = None,
        line_width: float = 1.0,
        wireframe: bool = False,
    ) -> None:
        if self._renderer is None or vtkActor is None or vtkDataSetMapper is None:
            return
        mapper = vtkDataSetMapper()
        mapper.SetInputData(mesh)
        actor = vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        prop.SetColor(*color)
        prop.SetOpacity(opacity)
        prop.SetLineWidth(line_width)
        if wireframe:
            prop.SetRepresentationToWireframe()
        else:
            prop.SetInterpolationToPhong()
            prop.SetSpecular(0.18)
            prop.SetSpecularPower(12.0)
        if pick_row_index is None:
            actor.PickableOff()
        else:
            actor_key = self._actor_key(actor)
            if actor_key is not None:
                self._actor_row_map[actor_key] = pick_row_index
                self._row_actor_map.setdefault(pick_row_index, []).append(actor_key)
        self._renderer.AddActor(actor)

    def _add_ray_actor(self, mesh, *, radius: float, color: tuple[float, float, float]) -> None:
        if self._renderer is None or vtkActor is None or vtkDataSetMapper is None:
            return
        actor = vtkActor()
        if vtkTubeFilter is not None:
            tube = vtkTubeFilter()
            tube.SetInputData(mesh)
            tube.SetRadius(radius)
            tube.SetNumberOfSides(10)
            tube.CappingOn()
            mapper = vtkDataSetMapper()
            mapper.SetInputConnection(tube.GetOutputPort())
            actor.SetMapper(mapper)
        else:
            mapper = vtkDataSetMapper()
            mapper.SetInputData(mesh)
            actor.SetMapper(mapper)
            actor.GetProperty().SetLineWidth(2.0)
        actor.GetProperty().SetColor(*color)
        actor.GetProperty().SetOpacity(0.95)
        actor.PickableOff()
        self._renderer.AddActor(actor)

    def _scene_bounds(self) -> tuple[np.ndarray, float]:
        if self._renderer is None:
            return np.zeros(3, dtype=float), 1.0
        bounds = np.asarray(self._renderer.ComputeVisiblePropBounds(), dtype=float)
        if bounds.size != 6 or not np.all(np.isfinite(bounds)) or bounds[0] > bounds[1]:
            return np.zeros(3, dtype=float), 1.0
        center = np.array(
            [
                0.5 * (bounds[0] + bounds[1]),
                0.5 * (bounds[2] + bounds[3]),
                0.5 * (bounds[4] + bounds[5]),
            ],
            dtype=float,
        )
        radius = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4], 1.0)
        return center, radius

    def set_camera_preset(self, preset: str) -> None:
        self._camera_preset = preset
        if self._renderer is None:
            return
        camera = self._renderer.GetActiveCamera()
        if camera is None:
            return
        center, radius = self._scene_bounds()
        distance = max(radius * 2.2, 50.0)
        if preset == "zy":
            position = center + np.array([distance, 0.0, 0.0], dtype=float)
            view_up = (0.0, 1.0, 0.0)
        elif preset == "xy":
            position = center + np.array([0.0, 0.0, distance], dtype=float)
            view_up = (0.0, 1.0, 0.0)
        elif preset == "xz":
            position = center + np.array([0.0, distance, 0.0], dtype=float)
            view_up = (0.0, 0.0, 1.0)
        else:
            position = center + np.array([-distance * 0.95, distance * 0.55, distance * 0.8], dtype=float)
            view_up = (0.0, 1.0, 0.0)
        camera.SetPosition(*position.tolist())
        camera.SetFocalPoint(*center.tolist())
        camera.SetViewUp(*view_up)
        self._renderer.ResetCameraClippingRange()
        self.render()

    def render(self) -> None:
        if self._vtk_widget is None:
            return
        try:
            self._vtk_widget.GetRenderWindow().Render()
        except Exception:
            pass

    def refresh_scene(self, system, rays, row_names: list[str], *, reset_camera: bool = False) -> None:
        if self._renderer is None:
            raise RuntimeError("Embedded VTK/Tk viewer unavailable")

        self._renderer.RemoveAllViewProps()
        self._actor_row_map.clear()
        self._row_actor_map.clear()
        self._picked_row_index = None

        transforms = getattr(system, "TRANS_2A", None)
        surfaces = getattr(system, "AAA", None)
        block_count = min(len(row_names), getattr(surfaces, "n_blocks", 0), len(transforms) if transforms is not None else 0)
        drew_surfaces = 0
        for index in range(block_count):
            surface = system.SDT_0[index]
            mesh = self._mesh_with_transform(surfaces[index], transforms[index])
            if mesh is None or int(getattr(mesh, "n_points", 0)) == 0:
                continue
            color = self._surface_color(surface)
            self._add_mesh_actor(mesh, color=color, opacity=0.68, pick_row_index=index)
            try:
                edges = mesh.extract_feature_edges(
                    feature_angle=10,
                    boundary_edges=True,
                    feature_edges=False,
                    manifold_edges=False,
                )
                if int(getattr(edges, "n_points", 0)) > 0:
                    self._add_mesh_actor(edges, color=(0.15, 0.15, 0.15), opacity=1.0, line_width=1.0)
            except Exception:
                pass
            drew_surfaces += 1

        side_index = 0
        for row_index in getattr(system, "side_number", []):
            if row_index >= len(row_names):
                side_index += 1
                continue
            try:
                mesh = pv.wrap(system.BBB[side_index]).extract_surface().copy(deep=True)
            except Exception:
                side_index += 1
                continue
            side_index += 1
            if int(getattr(mesh, "n_points", 0)) == 0:
                continue
            color = self._surface_color(system.SDT_0[row_index])
            self._add_mesh_actor(mesh, color=color, opacity=0.18, pick_row_index=row_index)

        if self.show_rays_var.get():
            center, radius = self._scene_bounds()
            ray_radius = max(radius * 0.0015, 0.08)
            for wave, ray_pts in zip(getattr(rays, "RayWave", []), getattr(rays, "CC", [])):
                try:
                    ray_mesh = pv.lines_from_points(ray_pts)
                except Exception:
                    continue
                if int(getattr(ray_mesh, "n_points", 0)) < 2:
                    continue
                color = tuple(wavelength_to_rgb(float(wave) * 1000.0))
                self._add_ray_actor(ray_mesh, radius=ray_radius, color=color)

        self._renderer.ResetCamera()
        self.set_camera_preset(self._camera_preset)
        self.highlight_row(self.editor._current_selected_row_index())
        self.status_var.set(f"3D scene ready | surfaces={drew_surfaces} | rays={len(getattr(rays, 'CC', []))}")
        self.render()

    def refresh_from_editor(self) -> None:
        try:
            system, rays = self.editor._build_preview_system_and_rays()
            row_names = [row.name for row in self.editor.rows]
            self.refresh_scene(system, rays, row_names, reset_camera=False)
            self.editor.status_var.set("3D inspector updated")
        except Exception as exc:
            self.status_var.set(f"3D refresh failed: {_short_error_message(exc)}")
            self.editor.append_debug(f"3D inspector refresh error: {exc}")

    def _on_left_button_press(self, obj, _event) -> None:
        if self._picker is None or self._renderer is None or self._vtk_interactor is None:
            return
        x, y = self._vtk_interactor.GetEventPosition()
        self._picker.Pick(x, y, 0.0, self._renderer)
        actor = self._picker.GetActor()
        actor_key = self._actor_key(actor)
        row_index = self._actor_row_map.get(actor_key) if actor_key is not None else None
        if row_index is None:
            self._set_row_highlight(None)
            self.status_var.set("3D scene ready")
            return
        self._set_row_highlight(row_index)
        self.editor._select_table_row(row_index)
        row_name = self.editor.rows[row_index].name if 0 <= row_index < len(self.editor.rows) else "Surface"
        self.status_var.set(f"Selected row {row_index}: {row_name}")
        self.render()

    def _on_close(self) -> None:
        self.editor._three_d_inspector = None
        try:
            self.destroy()
        except Exception:
            pass


def _optional_cupy():
    global _CUPY_IMPORT_ATTEMPTED, _CUPY_MODULE
    if not _CUPY_IMPORT_ATTEMPTED:
        _CUPY_IMPORT_ATTEMPTED = True
        _preload_cuda_libraries()
        try:
            import cupy as cp  # type: ignore
            _CUPY_MODULE = cp
        except Exception:
            _CUPY_MODULE = None
    return _CUPY_MODULE


def _optional_torch():
    global _TORCH_IMPORT_ATTEMPTED, _TORCH_MODULE
    if not _TORCH_IMPORT_ATTEMPTED:
        _TORCH_IMPORT_ATTEMPTED = True
        _preload_cuda_libraries()
        try:
            import torch  # type: ignore
            _TORCH_MODULE = torch
        except Exception:
            _TORCH_MODULE = None
    return _TORCH_MODULE


def _build_system_from_specs(row_specs: list[dict]) -> object:
    surfaces = []
    clear_aperture = max(
        [max(float(spec["diameter"]), 1.0) for spec in row_specs if spec["surface"] not in {"Object", "Image"}] or [100.0]
    ) * 4.0
    for spec in row_specs:
        surface = Kos.surf()
        surface.Name = ""
        surface.Rc = float(spec["rc"])
        surface.Thickness = float(spec["thickness"])
        surface.Diameter = clear_aperture if spec["surface"] in {"Object", "Image"} else float(spec["diameter"])
        surface.Glass = str(spec["glass"])
        surface.TiltX = float(spec.get("tilt_x", 0.0))
        surface.TiltY = float(spec.get("tilt_y", 0.0))
        surface.TiltZ = float(spec.get("tilt_z", 0.0))
        surface.DespX = float(spec.get("desp_x", 0.0))
        surface.DespY = float(spec.get("desp_y", 0.0))
        surface.DespZ = float(spec.get("desp_z", 0.0))
        surface.AxisMove = float(spec.get("axis_move", 0.0))
        surface.Drawing = 0.0 if spec["surface"] in {"Object", "Image", "Mirror"} else 1.0
        if spec["surface"] == "Mirror":
            surface.Glass = "MIRROR"
        if spec["surface"] == "Thin Lens":
            focal = float(spec["rc"])
            surface.Thin_Lens = focal if focal != 0.0 else 100.0
            surface.Rc = 0.0
        elif spec["surface"] == "Grating":
            surface.Diff_Ord = 1.0
            surface.Grating_D = 1.0
        surfaces.append(surface)
    return Kos.system(surfaces, Kos.Setup(), build=1)


def _row_specs_signature(row_specs: list[dict]):
    signature = []
    for spec in row_specs:
        signature.append(
            (
                str(spec.get("surface", "")),
                str(spec.get("name", "")),
                float(spec.get("rc", 0.0)),
                float(spec.get("thickness", 0.0)),
                float(spec.get("diameter", 0.0)),
                str(spec.get("glass", "AIR")),
                float(spec.get("tilt_x", 0.0)),
                float(spec.get("tilt_y", 0.0)),
                float(spec.get("tilt_z", 0.0)),
                float(spec.get("desp_x", 0.0)),
                float(spec.get("desp_y", 0.0)),
                float(spec.get("desp_z", 0.0)),
                float(spec.get("axis_move", 0.0)),
            )
        )
    return tuple(signature)


def _build_cached_system_from_specs(row_specs: list[dict]) -> object:
    global _WORKER_SYSTEM_CACHE_SIGNATURE, _WORKER_SYSTEM_CACHE_SYSTEM
    signature = _row_specs_signature(row_specs)
    if _WORKER_SYSTEM_CACHE_SYSTEM is None or _WORKER_SYSTEM_CACHE_SIGNATURE != signature:
        _WORKER_SYSTEM_CACHE_SYSTEM = _build_system_from_specs(row_specs)
        _WORKER_SYSTEM_CACHE_SIGNATURE = signature
    return _WORKER_SYSTEM_CACHE_SYSTEM


def _pick_image_plane_data_static(rays):
    try:
        X, Y, Z, L, M, N = rays.pick(-1, coordinates="local")
        if np.asarray(X).size:
            return X, Y, Z, L, M, N
    except Exception:
        pass
    return rays.pick(-1)


def _trace_analysis_chunk(
    row_specs: list[dict],
    wavelength: float,
    x_bundle,
    y_bundle,
    z_bundle,
    l_bundle,
    m_bundle,
    n_bundle,
):
    system = _build_cached_system_from_specs(row_specs)
    rays = Kos.raykeeper(system)
    Kos.TraceLoop(
        np.asarray(x_bundle, dtype=float),
        np.asarray(y_bundle, dtype=float),
        np.asarray(z_bundle, dtype=float),
        np.asarray(l_bundle, dtype=float),
        np.asarray(m_bundle, dtype=float),
        np.asarray(n_bundle, dtype=float),
        float(wavelength),
        rays,
        clean=1,
    )
    x_local, y_local, _z_local, _l_local, _m_local, _n_local = _pick_image_plane_data_static(rays)
    return np.asarray(x_local, dtype=float), np.asarray(y_local, dtype=float)


def _trace_analysis_chunk_full(
    row_specs: list[dict],
    wavelength: float,
    x_bundle,
    y_bundle,
    z_bundle,
    l_bundle,
    m_bundle,
    n_bundle,
):
    system = _build_cached_system_from_specs(row_specs)
    rays = Kos.raykeeper(system)
    Kos.TraceLoop(
        np.asarray(x_bundle, dtype=float),
        np.asarray(y_bundle, dtype=float),
        np.asarray(z_bundle, dtype=float),
        np.asarray(l_bundle, dtype=float),
        np.asarray(m_bundle, dtype=float),
        np.asarray(n_bundle, dtype=float),
        float(wavelength),
        rays,
        clean=1,
    )
    x_local, y_local, z_local, l_local, m_local, n_local = _pick_image_plane_data_static(rays)
    return (
        np.asarray(x_local, dtype=float),
        np.asarray(y_local, dtype=float),
        np.asarray(z_local, dtype=float),
        np.asarray(l_local, dtype=float),
        np.asarray(m_local, dtype=float),
        np.asarray(n_local, dtype=float),
    )


def _serialize_operand_results(operands) -> list[dict]:
    serialized = []
    for operand in operands:
        serialized.append(
            {
                "name": str(getattr(operand, "name", "")),
                "value": float(getattr(operand, "value", 0.0)),
                "weighted": float(getattr(operand, "weighted", 0.0)),
                "target": float(getattr(operand, "target", 0.0)),
            }
        )
    return serialized


def _run_optimization_job(
    progress_queue,
    stop_event,
    row_specs: list[dict],
    merit_function: MeritFunction,
    variables: list[OpticalVariable],
    x0: list[float],
    generations_total: int,
    verbosity_every: int,
    population_size: int,
    optimization_workers: int,
    parallel_enabled: bool,
):
    try:
        if os.name == "posix":
            try:
                os.setsid()
            except Exception:
                pass
        pagmo_lib = Path(os.path.expanduser("~/Projects/pagmo2/_install/lib64/libpagmo.so"))
        if pagmo_lib.exists():
            try:
                ctypes.CDLL(str(pagmo_lib), mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass
        import pygmo as pg  # type: ignore

        system = _build_system_from_specs(row_specs)
        has_mtf_operand = any(isinstance(operand, MTFAtFrequencyOperand) for operand in merit_function.operands)
        evaluator = MeritEvaluator(
            system.SDT,
            setup=system.SETUP,
            merit_function=merit_function,
            mtf_worker_count=max(1, int(optimization_workers)) if has_mtf_operand else 1,
        )
        try:
            initial = evaluator.evaluate(variables, x0)

            udp = Pygmo2MeritProblem(evaluator=evaluator, variables=variables)
            problem = pg.problem(udp)
            workers = 1
            backend = "sequential"
            population_kwargs: dict[str, object] = {"size": int(population_size), "seed": 42}
            debug_messages: list[str] = []
            if has_mtf_operand and int(optimization_workers) > 1:
                workers = max(1, int(optimization_workers))
                backend = f"mtf_chunks ({workers} workers)"
                debug_messages.append("Optimization uses internal MTF chunk tracing instead of pygmo mp_bfe.")
            elif parallel_enabled and int(optimization_workers) > 1:
                workers = max(1, int(optimization_workers))
                try:
                    pg.mp_bfe.resize_pool(workers)
                    population_kwargs["b"] = pg.bfe(pg.mp_bfe())
                    backend = f"mp_bfe ({workers} workers)"
                except Exception as exc:
                    workers = 1
                    backend = "sequential"
                    debug_messages.append(f"Optimization parallel backend disabled: {exc}")

            progress_queue.put(
                {
                    "type": "bootstrap",
                    "initial_total": float(initial.total),
                    "compute_backend": backend,
                    "workers": workers,
                    "debug_messages": debug_messages,
                }
            )

            try:
                population = pg.population(problem, **population_kwargs)
            except Exception as exc:
                if "b" not in population_kwargs:
                    raise RuntimeError(f"failed to initialize population: {exc}") from exc
                debug_messages.append(f"Optimization population batch evaluator failed: {exc}")
                workers = 1
                backend = "sequential"
                progress_queue.put(
                    {
                        "type": "bootstrap",
                        "initial_total": float(initial.total),
                        "compute_backend": backend,
                        "workers": workers,
                        "debug_messages": [debug_messages[-1]],
                    }
                )
                population = pg.population(problem, size=int(population_size), seed=42)
            population.push_back(x0)

            for generation_done in range(int(generations_total)):
                if stop_event.is_set():
                    break
                algorithm = pg.algorithm(pg.de(gen=1, seed=42 + int(generation_done)))
                algorithm.set_verbosity(1)
                capture = io.StringIO()
                with redirect_stdout(capture), redirect_stderr(capture):
                    population = algorithm.evolve(population)
                logs = algorithm.extract(pg.de).get_log()
                payload = {
                    "type": "generation",
                    "generation_done": int(generation_done) + 1,
                    "debug": capture.getvalue(),
                    "champion_x": [float(value) for value in population.champion_x],
                }
                if logs:
                    gen, fevals, best, dx, df = logs[-1]
                    payload.update(
                        {
                            "log_gen": int(gen),
                            "log_fevals": int(fevals),
                            "log_best": float(best),
                            "log_dx": float(dx),
                            "log_df": float(df),
                            "verbosity_every": int(verbosity_every),
                            "generations_total": int(generations_total),
                        }
                    )
                progress_queue.put(payload)

            champion_x = [float(value) for value in population.champion_x]
            champion = evaluator.evaluate(variables, champion_x)
            progress_queue.put(
                {
                    "type": "complete",
                    "cancelled": bool(stop_event.is_set()),
                    "champion_x": champion_x,
                    "initial_total": float(initial.total),
                    "final_total": float(champion.total),
                    "compute_backend": backend,
                    "workers": workers,
                    "operands": _serialize_operand_results(champion.operands),
                }
            )
        finally:
            try:
                evaluator._shutdown_mtf_executor()
            except Exception:
                pass
    except Exception as exc:
        progress_queue.put(
            {
                "type": "error",
                "message": str(exc),
                "traceback": traceback.format_exc(),
            }
        )


class KrakenLayoutEditor(tk.Tk):
    def __init__(self, *, headless: bool = False) -> None:
        super().__init__()
        self.headless = headless
        self.title("KrakenOS Layout Editor")
        self.geometry("1400x850")
        self.minsize(1100, 720)
        if not self.headless:
            self.after(50, self._maximize_window)

        self.current_layout_file: Path | None = None
        self.layout_files: dict[str, Path] = {}
        self.example_files: dict[str, Path] = {}
        self.rows: list[SurfaceRow] = []
        self.editor: tk.Widget | None = None
        self.popup_menu: tk.Menu | None = None
        self.current_menu_row_id: str | None = None
        self.current_menu_field: str | None = None
        self._text_popup_menu: tk.Menu | None = None
        self._formula_help_path: Path | None = None
        self._undo_button: ttk.Button | None = None
        self._redo_button: ttk.Button | None = None
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
        self.operand_mtf_algorithm_vars: dict[str, tk.StringVar] = {}
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
        self._autosave_after_id: str | None = None
        self._initial_layout_passes = 0
        self._last_field_type = "Angle"
        self._field_defaults_initialized = False
        self._field_type_defaults = {
            "Angle": "0.0",
            "Object Height": "0.0",
            "Paraxial Image Height": "0.0",
            "Real Image Height": "0.0",
        }
        self.auto_save_plot_var = tk.BooleanVar(value=not self.headless)
        self.show_native_overlays_var = tk.BooleanVar(value=True)
        self.show_native_active_spans_var = tk.BooleanVar(value=False)
        self.show_native_hit_labels_var = tk.BooleanVar(value=False)
        self.show_clipped_rays_var = tk.BooleanVar(value=True)
        self._last_analysis_label = "2D"
        self._last_analysis_workers = 1
        self._last_analysis_parallel_capable = False
        self._last_analysis_accelerator = "CPU"
        self._gpu_backend_reported = False
        self._analysis_executor: ProcessPoolExecutor | None = None
        self._analysis_executor_workers = 0
        self._optimization_process = None
        self._optimization_queue = None
        self._optimization_stop_event = None
        self._last_optics_info: dict | None = None
        self._cardinal_marker_artists: list = []
        self._analysis_ax = None
        self._hover_hint_artists: dict = {}
        self._hover_axis = None
        self._last_viewer_open_time = 0.0
        self._three_d_inspector: Kraken3DInspector | None = None
        self._legacy_3d_plotter = None
        self._legacy_3d_after_id = None
        self._layout_pick_regions: dict[int, np.ndarray] = {}
        self._layout_selection_artists: list = []
        self._undo_stack: list[dict[str, object]] = []
        self._redo_stack: list[dict[str, object]] = []
        self._history_pending_state: dict[str, object] | None = None
        self._history_restoring = False
        self._history_limit = 80

        self._build_menu()
        self._build_ui()
        self._bind_global_copy_shortcuts()
        self.bind_all("<Control-z>", self._undo_event, add="+")
        self.bind_all("<Control-y>", self._redo_event, add="+")
        self.bind_all("<Control-Shift-Z>", self._redo_event, add="+")
        self._reset_debug_log()
        self.load_layouts()
        self.load_examples()
        if self.layout_names:
            initial_layout = DEFAULT_LAYOUT_TITLE if DEFAULT_LAYOUT_TITLE in self.layout_files else self.layout_names[0]
            self.load_layout_by_name(initial_layout)
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._history_pending_state = None
        self._update_undo_redo_buttons()
        self.after(0, self._report_compute_backends)

    def _maximize_window(self) -> None:
        # Prefer maximize/zoom over fullscreen so copy/paste and WM behavior remain normal.
        try:
            self.state("zoomed")
            return
        except Exception:
            pass
        try:
            self.attributes("-zoomed", True)
            return
        except Exception:
            pass
        try:
            width = max(1200, int(self.winfo_screenwidth() * 0.96))
            height = max(800, int(self.winfo_screenheight() * 0.95))
            self.geometry(f"{width}x{height}+0+0")
        except Exception:
            pass

    def _build_menu(self) -> None:
        menubar = tk.Menu(self)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open", command=self.open_layout)
        file_menu.add_command(label="Save", command=self.save_layout)
        file_menu.add_command(label="Save As", command=self.save_layout_as)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        edit_menu = tk.Menu(menubar, tearoff=0)
        edit_menu.add_command(label="Undo", command=self.undo, accelerator="Ctrl+Z")
        edit_menu.add_command(label="Redo", command=self.redo, accelerator="Ctrl+Y")
        menubar.add_cascade(label="Edit", menu=edit_menu)

        action_menu = tk.Menu(menubar, tearoff=0)
        action_menu.add_command(label="Refresh Plot", command=self.refresh_plot)
        action_menu.add_command(label="Benchmark PSF/MTF", command=self.benchmark_psf_mtf)
        action_menu.add_command(label="Copy Debug", command=self.copy_debug_to_clipboard)
        action_menu.add_command(label="Clear Marks", command=self.clear_optimization_marks)
        action_menu.add_checkbutton(label="Auto-save JPG", variable=self.auto_save_plot_var)
        action_menu.add_separator()
        action_menu.add_checkbutton(label="Show Native Overlays", variable=self.show_native_overlays_var, command=self.refresh_plot)
        action_menu.add_checkbutton(label="Show Native Active Spans", variable=self.show_native_active_spans_var, command=self.refresh_plot)
        action_menu.add_checkbutton(label="Show Native Hit Labels", variable=self.show_native_hit_labels_var, command=self.refresh_plot)
        menubar.add_cascade(label="Actions", menu=action_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="Paraxial Calculator", command=self.open_paraxial_calculator)
        help_menu.add_command(label="Optics Formula Sheet", command=self.show_formula_help)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menubar)

    def destroy(self) -> None:
        if self._three_d_inspector is not None:
            try:
                self._three_d_inspector.destroy()
            except Exception:
                pass
            self._three_d_inspector = None
        self._close_legacy_3d_plotter()
        self._shutdown_analysis_executor()
        self._shutdown_optimization_worker(force=True)
        super().destroy()

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

        control_host = ttk.Frame(top, padding=(8, 8, 4, 8))
        control_host.columnconfigure(0, weight=1)
        control_host.rowconfigure(0, weight=1)
        top.add(control_host, weight=1)

        self.control_canvas = tk.Canvas(control_host, highlightthickness=0, borderwidth=0)
        self.control_canvas.grid(row=0, column=0, sticky="nsew")
        control_scroll = ttk.Scrollbar(control_host, orient="vertical", command=self.control_canvas.yview)
        control_scroll.grid(row=0, column=1, sticky="ns", padx=(6, 0))
        self.control_canvas.configure(yscrollcommand=control_scroll.set)

        control_stack = ttk.Frame(self.control_canvas)
        control_stack.columnconfigure(0, weight=1)
        self.control_stack_window = self.control_canvas.create_window((0, 0), window=control_stack, anchor="nw")
        control_stack.bind("<Configure>", self._on_control_stack_configure)
        self.control_canvas.bind("<Configure>", self._on_control_canvas_configure)

        controls = ttk.LabelFrame(control_stack, text="Display", padding=8)
        controls.grid(row=0, column=0, sticky="ew")
        for column in range(2):
            controls.columnconfigure(column, weight=1, uniform="display_cols")

        field_panel = ttk.LabelFrame(control_stack, text="Field", padding=8)
        field_panel.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        for column in range(2):
            field_panel.columnconfigure(column, weight=1, uniform="field_cols")

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
        self._undo_button = ttk.Button(table_toolbar, text="Undo", command=self.undo)
        self._undo_button.pack(side="left")
        self._redo_button = ttk.Button(table_toolbar, text="Redo", command=self.redo)
        self._redo_button.pack(side="left", padx=(6, 6))
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
                else 120 if field == "glass"
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
        self.table.bind("<<TreeviewSelect>>", self._on_table_selection_changed, add="+")
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
        self._selection_anchor_row: str | None = None
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
        self.analysis_mode_button_var = tk.StringVar(value=self.analysis_mode)
        mode_buttons = (
            ("2D", "none"),
            ("Native", "native_off_axis"),
            ("Spot", "spot"),
            ("PSF", "psf"),
            ("RMS", "rms"),
            ("FC/Dist", "field_curvature"),
            ("Pupil", "pupil"),
            ("Seidel", "seidel"),
            ("Wavefront", "wavefront"),
            ("MTF", "mtf"),
        )
        for text, mode in mode_buttons:
            ttk.Radiobutton(
                plot_toolbar,
                text=text,
                style="Toolbutton",
                variable=self.analysis_mode_button_var,
                value=mode,
                command=lambda m=mode: self.set_analysis_mode(m),
            ).pack(side="left", padx=(6, 0))
        ttk.Checkbutton(
            plot_toolbar,
            text="Show PP / EP / XP",
            variable=self.show_cardinals_var,
            command=self._on_toggle_cardinal_markers,
        ).pack(side="left", padx=(12, 0))
        ttk.Button(plot_toolbar, text="Update", command=self._manual_update_plot).pack(side="right")

        self.figure = Figure(figsize=(7, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.canvas.mpl_connect("motion_notify_event", self._on_plot_canvas_motion)
        self.canvas.mpl_connect("figure_leave_event", self._on_plot_canvas_leave)
        self.canvas.get_tk_widget().bind("<Button-1>", self._on_plot_widget_click, add="+")

        self.debug_text = tk.Text(debug_frame, wrap="word", height=8, width=24)
        self.debug_text.grid(row=0, column=0, sticky="nsew")
        debug_scroll = ttk.Scrollbar(debug_frame, orient="vertical", command=self.debug_text.yview)
        debug_scroll.grid(row=0, column=1, sticky="ns")
        self.debug_text.configure(yscrollcommand=debug_scroll.set)
        self._bind_text_copy_shortcuts(self.debug_text)
        debug_actions = ttk.Frame(debug_frame)
        debug_actions.grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Button(debug_actions, text="Copy Selected", command=lambda: self._copy_selection_from_text_widget(self.debug_text)).pack(side="left")
        ttk.Button(debug_actions, text="Copy All", command=lambda: self._copy_all_from_text_widget(self.debug_text)).pack(side="left", padx=(6, 0))

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
        self._bind_text_copy_shortcuts(self.progress_text)
        self._bind_text_context_menu(self.debug_text)
        self._bind_text_context_menu(self.progress_text)

        status_bar = ttk.Frame(self)
        status_bar.grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 2))
        status_bar.columnconfigure(0, weight=0)
        status_bar.columnconfigure(1, weight=1)
        self.status_var = tk.StringVar(value="Ready")
        self.status_hint_var = tk.StringVar(value="")
        ttk.Label(status_bar, textvariable=self.status_var, anchor="w").grid(row=0, column=0, sticky="w")
        ttk.Label(
            status_bar,
            textvariable=self.status_hint_var,
            anchor="e",
            justify="right",
            foreground="#475569",
        ).grid(row=0, column=1, sticky="ew", padx=(12, 0))
        self.after_idle(self._set_initial_pane_layout)
        self.bind("<Configure>", self._maybe_refresh_initial_pane_layout, add="+")

    def _build_controls_panel(self, parent) -> None:
        for column in range(2):
            parent.columnconfigure(column, weight=1)

        ttk.Label(parent, text="Object mode").grid(row=0, column=0, sticky="w", pady=(0, 2))
        self.object_mode_var = tk.StringVar(value="Infinity")
        self.object_mode_menu = ttk.Combobox(
            parent,
            textvariable=self.object_mode_var,
            state="readonly",
            width=12,
            values=["Finite", "Infinity"],
        )
        self.object_mode_menu.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self.object_mode_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.object_mode_menu.bind("<<ComboboxSelected>>", self._on_object_mode_changed)

        ttk.Label(parent, text="Wavelength [um]").grid(row=0, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.wavelength_var = tk.StringVar(value="0.55")
        wavelength_entry = ttk.Entry(parent, textvariable=self.wavelength_var, width=12)
        wavelength_entry.grid(
            row=1, column=1, sticky="ew", pady=(0, 8), padx=(8, 0)
        )

        ttk.Label(parent, text="Orientation").grid(row=2, column=0, sticky="w", pady=(0, 2))
        self.display_orientation_var = tk.StringVar(value="Vertical")
        self.display_orientation_menu = ttk.Combobox(
            parent,
            textvariable=self.display_orientation_var,
            state="readonly",
            width=12,
            values=["Vertical", "Horizontal"],
        )
        self.display_orientation_menu.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        self.display_orientation_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.display_orientation_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)

        ttk.Label(parent, text="Ray fan count").grid(row=2, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.ray_count_var = tk.StringVar(value="31")
        ray_count_entry = ttk.Entry(parent, textvariable=self.ray_count_var, width=12)
        ray_count_entry.grid(
            row=3, column=1, sticky="ew", pady=(0, 8), padx=(8, 0)
        )

        ttk.Label(parent, text="Pupil factor").grid(row=4, column=0, sticky="w", pady=(0, 2))
        self.ray_height_factor_var = tk.StringVar(value="0.8")
        ray_height_entry = ttk.Entry(parent, textvariable=self.ray_height_factor_var, width=12)
        ray_height_entry.grid(
            row=5, column=0, sticky="ew", pady=(0, 8)
        )

        ttk.Label(parent, text="Analysis stop surface").grid(row=4, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.analysis_surface_var = tk.StringVar(value="Auto")
        self.analysis_surface_menu = ttk.Combobox(
            parent,
            textvariable=self.analysis_surface_var,
            state="readonly",
            width=12,
            values=["Auto"],
        )
        self.analysis_surface_menu.grid(row=5, column=1, sticky="ew", pady=(0, 8), padx=(8, 0))
        self.analysis_surface_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.analysis_surface_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)

        ttk.Label(parent, text="Aperture type").grid(row=6, column=0, sticky="w", pady=(0, 2))
        self.aperture_type_var = tk.StringVar(value="EPD")
        self.aperture_type_menu = ttk.Combobox(
            parent,
            textvariable=self.aperture_type_var,
            state="readonly",
            width=12,
            values=["STOP", "EPD"],
        )
        self.aperture_type_menu.grid(row=7, column=0, sticky="ew")
        self.aperture_type_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.aperture_type_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)

        ttk.Label(parent, text="Aperture value").grid(row=6, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.aperture_value_var = tk.StringVar(value="4.0")
        aperture_value_entry = ttk.Entry(parent, textvariable=self.aperture_value_var, width=12)
        aperture_value_entry.grid(
            row=7, column=1, sticky="ew", padx=(8, 0)
        )

        ttk.Label(parent, text="Spot view").grid(row=8, column=0, sticky="w", pady=(8, 2))
        self.spot_view_mode_var = tk.StringVar(value="Grid")
        self.spot_view_mode_menu = ttk.Combobox(
            parent,
            textvariable=self.spot_view_mode_var,
            state="readonly",
            width=12,
            values=["Grid", "Absolute", "Centroid"],
        )
        self.spot_view_mode_menu.grid(row=9, column=0, sticky="ew")
        self.spot_view_mode_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.spot_view_mode_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)

        clipped_check = ttk.Checkbutton(
            parent,
            text="Show clipped rays",
            variable=self.show_clipped_rays_var,
            command=self._mark_plot_update_pending,
        )
        clipped_check.grid(row=9, column=1, sticky="w", padx=(8, 0))
        clipped_check.bind("<ButtonPress-1>", self._begin_history_capture, add="+")

        self.show_cardinals_var = tk.BooleanVar(value=True)

        self._bind_deferred_manual_update(wavelength_entry)
        self._bind_deferred_manual_update(ray_count_entry)
        self._bind_deferred_manual_update(ray_height_entry)
        self._bind_deferred_manual_update(aperture_value_entry)

    def _build_field_panel(self, parent) -> None:
        for column in range(2):
            parent.columnconfigure(column, weight=1)

        ttk.Label(parent, text="Field type").grid(row=0, column=0, sticky="w", pady=(0, 2))
        self.field_type_var = tk.StringVar(value="Angle")
        self.field_type_menu = ttk.Combobox(
            parent,
            textvariable=self.field_type_var,
            state="readonly",
            width=12,
            values=["Angle", "Object Height", "Paraxial Image Height", "Real Image Height"],
        )
        self.field_type_menu.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self.field_type_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.field_type_menu.bind("<<ComboboxSelected>>", self._on_field_type_changed)

        self.field_mode_note_var = tk.StringVar(value="")

        self.field_value_label_var = tk.StringVar(value="Angle [deg]")
        ttk.Label(parent, textvariable=self.field_value_label_var).grid(row=0, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.field_value_var = tk.StringVar(value="5.0")
        field_value_entry = ttk.Entry(parent, textvariable=self.field_value_var, width=12)
        field_value_entry.grid(
            row=1, column=1, sticky="ew", pady=(0, 8), padx=(8, 0)
        )

        ttk.Label(parent, text="Field count").grid(row=2, column=0, sticky="w", pady=(0, 2))
        self.field_count_var = tk.StringVar(value="1")
        field_count_entry = ttk.Entry(parent, textvariable=self.field_count_var, width=12)
        field_count_entry.grid(
            row=3, column=0, sticky="ew", pady=(0, 8)
        )

        ttk.Label(parent, text="Image size").grid(row=2, column=1, sticky="w", pady=(0, 2), padx=(8, 0))
        self.image_diameter_mode_var = tk.StringVar(value="Auto")
        self.image_diameter_mode_menu = ttk.Combobox(
            parent,
            textvariable=self.image_diameter_mode_var,
            state="readonly",
            width=12,
            values=["Auto", "Manual"],
        )
        self.image_diameter_mode_menu.grid(row=3, column=1, sticky="ew", pady=(0, 8), padx=(8, 0))
        self.image_diameter_mode_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
        self.image_diameter_mode_menu.bind("<<ComboboxSelected>>", self._on_image_diameter_mode_changed)

        self.field_warning_var = tk.StringVar(value="")
        self.field_summary_var = tk.StringVar(value="")

        self._bind_deferred_manual_update(field_value_entry, sync_fields=True)
        self._bind_deferred_manual_update(field_count_entry, sync_fields=True)
        self._sync_field_mode_ui()

    def _on_control_stack_configure(self, _event=None) -> None:
        if not hasattr(self, "control_canvas"):
            return
        self.control_canvas.configure(scrollregion=self.control_canvas.bbox("all"))

    def _on_control_canvas_configure(self, event=None) -> None:
        if not hasattr(self, "control_canvas") or not hasattr(self, "control_stack_window"):
            return
        width = self.control_canvas.winfo_width() if event is None else int(event.width)
        self.control_canvas.itemconfigure(self.control_stack_window, width=max(width, 1))

    def _update_field_status_hint(self) -> None:
        if not hasattr(self, "status_hint_var"):
            return
        note = self.field_mode_note_var.get().strip() if hasattr(self, "field_mode_note_var") else ""
        warning = self.field_warning_var.get().strip() if hasattr(self, "field_warning_var") else ""
        summary = self.field_summary_var.get().strip() if hasattr(self, "field_summary_var") else ""
        summary = summary.replace("\n", " | ")
        parts = [part for part in (note, warning, summary) if part]
        self.status_hint_var.set("  ||  ".join(parts))
    def _build_optimization_panel(self, parent) -> None:
        operand_list_width = max((len(spec.label) for spec in OPERAND_REGISTRY.values()), default=14) + 2
        operand_list_minsize = max(150, operand_list_width * 8)
        parent.columnconfigure(0, weight=0, minsize=operand_list_minsize)
        parent.columnconfigure(1, weight=1, minsize=220)

        button_row = ttk.Frame(parent)
        button_row.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        ttk.Button(button_row, text="Start Optimization", command=self.start_optimization).pack(side="left")
        ttk.Button(button_row, text="Stop", command=self.stop_optimization).pack(side="left", padx=(8, 0))
        cpu_total = max(1, int(os.cpu_count() or 1))
        worker_choices = ["Auto", "1"]
        for candidate in (2, 4, 6, 8, 12, 16, cpu_total):
            candidate = max(1, min(cpu_total, int(candidate)))
            text = str(candidate)
            if text not in worker_choices:
                worker_choices.append(text)
        ttk.Label(button_row, text="Workers").pack(side="left", padx=(14, 0))
        self.optimization_workers_var = tk.StringVar(value="Auto")
        ttk.Combobox(
            button_row,
            textvariable=self.optimization_workers_var,
            state="readonly",
            width=6,
            values=worker_choices,
        ).pack(side="left", padx=(6, 0))

        ttk.Label(parent, text="Merit operands").grid(row=1, column=0, sticky="w", pady=(0, 2))
        self.merit_mode_list = tk.Listbox(
            parent,
            exportselection=False,
            selectmode="extended",
            height=min(4, max(2, len(OPERAND_REGISTRY))),
            width=operand_list_width,
        )
        for spec in OPERAND_REGISTRY.values():
            self.merit_mode_list.insert("end", spec.label)
        if OPERAND_REGISTRY:
            self.merit_mode_list.selection_set(0)
        self.merit_mode_list.grid(row=2, column=0, sticky="nsw", pady=(0, 8), padx=(0, 8))
        self.merit_mode_list.bind("<ButtonPress-1>", self._begin_history_capture, add="+")
        self.merit_mode_list.bind("<<ListboxSelect>>", lambda _e: self._update_operand_setup_visibility())

        setup_holder = ttk.Frame(parent, height=320)
        setup_holder.grid(row=2, column=1, sticky="nsew", pady=(0, 8), padx=(4, 0))
        setup_holder.grid_propagate(False)
        setup_holder.columnconfigure(0, weight=1)
        setup_holder.rowconfigure(0, weight=1)

        setup_frame = ttk.Frame(setup_holder)
        setup_frame.grid(row=0, column=0, sticky="nsew")
        setup_frame.columnconfigure(0, weight=1, minsize=220)
        ttk.Label(setup_frame, text="Operand setup").grid(row=0, column=0, sticky="w", pady=(0, 2))

        for idx, spec in enumerate(OPERAND_REGISTRY.values(), start=1):
            card = ttk.LabelFrame(setup_frame, text=spec.label, padding=6)
            card.grid(row=idx, column=0, sticky="ew", pady=(0, 8))
            card.columnconfigure(1, weight=1, minsize=120)
            control_widgets: dict[str, tuple[tk.Widget, ...]] = {}

            weight_var = tk.StringVar(value=f"{spec.default_weight:g}")
            self.operand_weight_vars[spec.label] = weight_var
            weight_label = ttk.Label(card, text="Weight")
            weight_label.grid(row=0, column=0, sticky="w")
            weight_entry = ttk.Entry(card, textvariable=weight_var, width=12)
            weight_entry.grid(row=0, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(weight_entry)
            control_widgets["weight"] = (weight_label, weight_entry)

            target_var = tk.StringVar(value=f"{spec.default_target:g}")
            self.operand_target_vars[spec.label] = target_var
            target_label = ttk.Label(card, text="Target")
            target_label.grid(row=1, column=0, sticky="w")
            target_entry = ttk.Entry(card, textvariable=target_var, width=12)
            target_entry.grid(row=1, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(target_entry)
            control_widgets["target"] = (target_label, target_entry)

            wavelength_var = tk.StringVar(value=self.wavelength_var.get())
            self.operand_wavelength_vars[spec.label] = wavelength_var
            wavelength_label = ttk.Label(card, text="Wvl")
            wavelength_label.grid(row=2, column=0, sticky="w")
            wavelength_entry = ttk.Entry(card, textvariable=wavelength_var, width=12)
            wavelength_entry.grid(row=2, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(wavelength_entry)
            control_widgets["wavelength"] = (wavelength_label, wavelength_entry)

            field_var = tk.StringVar(value="0")
            self.operand_field_vars[spec.label] = field_var
            field_label = ttk.Label(card, text="Field")
            field_label.grid(row=3, column=0, sticky="w")
            field_entry = ttk.Entry(card, textvariable=field_var, width=12)
            field_entry.grid(row=3, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            self._bind_deferred_refresh(field_entry)
            control_widgets["field"] = (field_label, field_entry)

            surface_row = 4
            frequency_row = 6
            mode_row = 7
            algorithm_row = 8

            if spec.label == "MTF @ freq":
                field_x_var = tk.StringVar(value="0")
                field_y_var = tk.StringVar(value="0")
                self.operand_field_x_vars[spec.label] = field_x_var
                self.operand_field_y_vars[spec.label] = field_y_var
                field_x_label = ttk.Label(card, text="Field X")
                field_x_label.grid(row=3, column=0, sticky="w")
                field_x_entry = ttk.Entry(card, textvariable=field_x_var, width=12)
                field_x_entry.grid(row=3, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
                self._bind_deferred_refresh(field_x_entry)
                field_y_label = ttk.Label(card, text="Field Y(s)")
                field_y_label.grid(row=4, column=0, sticky="w")
                field_y_entry = ttk.Entry(card, textvariable=field_y_var, width=12)
                field_y_entry.grid(row=4, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
                self._bind_deferred_refresh(field_y_entry)
                control_widgets["field_xy"] = (field_x_label, field_x_entry, field_y_label, field_y_entry)
                field_label.grid_remove()
                field_entry.grid_remove()
                surface_row = 5
                frequency_row = 6
                mode_row = 7
                algorithm_row = 8

            surface_var = tk.StringVar(value="Auto")
            self.operand_surface_vars[spec.label] = surface_var
            surface_label = ttk.Label(card, text="Surf")
            surface_label.grid(row=surface_row, column=0, sticky="w")
            surface_menu = ttk.Combobox(card, textvariable=surface_var, state="readonly", width=12, values=["Auto"])
            surface_menu.grid(row=surface_row, column=1, sticky="ew", padx=(6, 0), pady=(0, 4))
            surface_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
            surface_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)
            control_widgets["surface"] = (surface_label, surface_menu)

            if spec.label == "MTF @ freq":
                frequency_var = tk.StringVar(value="5")
                self.operand_frequency_vars[spec.label] = frequency_var
                frequency_label = ttk.Label(card, text="Freq")
                frequency_label.grid(row=frequency_row, column=0, sticky="w")
                frequency_entry = ttk.Entry(card, textvariable=frequency_var, width=12)
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
                    width=12,
                    values=["Average", "Tangential", "Sagittal"],
                )
                mtf_mode_menu.grid(row=mode_row, column=1, sticky="ew", padx=(6, 0), pady=(4, 0))
                mtf_mode_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
                mtf_mode_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)
                control_widgets["mtf_mode"] = (mode_label, mtf_mode_menu)

                mtf_algorithm_var = tk.StringVar(value="PSF FFT")
                self.operand_mtf_algorithm_vars[spec.label] = mtf_algorithm_var
                algorithm_label = ttk.Label(card, text="Alg")
                algorithm_label.grid(row=algorithm_row, column=0, sticky="w")
                mtf_algorithm_menu = ttk.Combobox(
                    card,
                    textvariable=mtf_algorithm_var,
                    state="readonly",
                    width=12,
                    values=["PSF FFT", "LSF FFT"],
                )
                mtf_algorithm_menu.grid(row=algorithm_row, column=1, sticky="ew", padx=(6, 0), pady=(4, 0))
                mtf_algorithm_menu.bind("<FocusIn>", self._begin_history_capture, add="+")
                mtf_algorithm_menu.bind("<<ComboboxSelected>>", self._mark_plot_update_pending)
                control_widgets["mtf_algorithm"] = (algorithm_label, mtf_algorithm_menu)

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
        widget.bind("<FocusIn>", self._begin_history_capture, add="+")
        widget.bind("<Return>", self._mark_plot_update_pending)
        widget.bind("<Tab>", self._mark_plot_update_pending)
        widget.bind("<FocusOut>", self._mark_plot_update_pending)

    def _bind_deferred_manual_update(self, widget: tk.Widget, *, sync_fields: bool = False) -> None:
        def _on_commit(_event=None):
            if sync_fields:
                self._sync_object_controls()
            self._mark_plot_update_pending()

        widget.bind("<FocusIn>", self._begin_history_capture, add="+")
        widget.bind("<Return>", _on_commit)
        widget.bind("<Tab>", _on_commit)
        widget.bind("<FocusOut>", _on_commit)

    def _mark_plot_update_pending(self, _event=None) -> None:
        self._commit_history_capture()
        if hasattr(self, "status_var"):
            self.status_var.set("Display settings changed. Click Update.")

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
        self._commit_history_capture()
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
        if hasattr(self, "analysis_mode_button_var"):
            self.analysis_mode_button_var.set(mode)
        mode_label_map = {
            "none": "2D",
            "native_off_axis": "Native",
            "spot": "Spot",
            "psf": "PSF",
            "rms": "RMS",
            "field_curvature": "FC/Dist",
            "pupil": "Pupil",
            "seidel": "Seidel",
            "wavefront": "Wavefront",
            "mtf": "MTF",
        }
        mode_label = mode_label_map.get(mode, mode or "2D")
        if hasattr(self, "status_var"):
            self.status_var.set(f"Analysis mode set to {mode_label}. Click Update.")
        self.append_progress(f"Mode selected: {mode_label} (pending update).")

    def _manual_update_plot(self) -> None:
        mode = (self.analysis_mode or "none").strip()
        mode_label_map = {
            "none": "2D",
            "native_off_axis": "Native",
            "spot": "Spot",
            "psf": "PSF",
            "rms": "RMS",
            "field_curvature": "FC/Dist",
            "pupil": "Pupil",
            "seidel": "Seidel",
            "wavefront": "Wavefront",
            "mtf": "MTF",
        }
        mode_label = mode_label_map.get(mode, mode or "2D")
        modes_with_internal_progress = {"psf", "pupil", "seidel", "wavefront", "field_curvature", "mtf"}
        if mode in modes_with_internal_progress:
            self.append_progress(f"Display update requested ({mode_label}).")
            self.refresh_plot()
            self.append_progress(f"Display update completed ({mode_label}).")
            return
        self._begin_analysis_progress("Display update")
        self._update_analysis_progress(f"Refreshing {mode_label}", 1, 2)
        self.refresh_plot()
        self._update_analysis_progress("Rendering", 2, 2)
        self._finish_analysis_progress("Display update", success=True)

    def _open_plot_axis_once(self, target_ax) -> None:
        if target_ax not in {self.ax, self._analysis_ax}:
            return
        now = time.monotonic()
        if now - self._last_viewer_open_time < 0.4:
            return
        self._last_viewer_open_time = now
        self._open_high_res_plot_in_system_viewer(target_ax)

    def _on_plot_canvas_motion(self, event) -> None:
        target_ax = getattr(event, "inaxes", None)
        if target_ax not in self._hover_hint_artists:
            target_ax = None
        if target_ax is self._hover_axis:
            return
        self._set_hover_axis(target_ax)

    def _on_plot_canvas_leave(self, _event=None) -> None:
        self._set_hover_axis(None)

    def _on_plot_widget_click(self, event) -> str | None:
        try:
            self.canvas.draw()
            renderer = self.figure.canvas.get_renderer()
            widget = self.canvas.get_tk_widget()
            x_display = float(event.x)
            y_display = float(widget.winfo_height() - event.y)
            if self.ax is not None and self.ax in self.figure.axes:
                if self.ax.get_window_extent(renderer).contains(x_display, y_display):
                    row_index = self._find_layout_pick_row(x_display, y_display)
                    if row_index is not None:
                        self._select_table_row(row_index)
                        return "break"
                    self._open_plot_axis_once(self.ax)
                    return "break"
            if self._analysis_ax is not None and self._analysis_ax in self.figure.axes:
                if self._analysis_ax.get_window_extent(renderer).contains(x_display, y_display):
                    self._open_plot_axis_once(self._analysis_ax)
                    return "break"
        except Exception as exc:
            self.append_debug(f"Plot viewer dispatch failed: {exc}")
        return None

    @staticmethod
    def _convex_hull_2d(points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        if pts.ndim != 2 or pts.shape[1] != 2:
            return np.empty((0, 2), dtype=float)
        pts = pts[np.all(np.isfinite(pts), axis=1)]
        if pts.shape[0] <= 2:
            return pts
        pts = np.unique(np.round(pts, decimals=6), axis=0)
        if pts.shape[0] <= 2:
            return pts

        def cross(o, a, b) -> float:
            return float((a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]))

        ordered = sorted((float(x), float(y)) for x, y in pts)
        lower: list[tuple[float, float]] = []
        for p in ordered:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0.0:
                lower.pop()
            lower.append(p)
        upper: list[tuple[float, float]] = []
        for p in reversed(ordered):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0.0:
                upper.pop()
            upper.append(p)
        hull = lower[:-1] + upper[:-1]
        return np.asarray(hull, dtype=float)

    def _clear_layout_selection_overlay(self) -> None:
        for artist in self._layout_selection_artists:
            try:
                artist.remove()
            except Exception:
                pass
        self._layout_selection_artists = []

    def _project_layout_points(self, points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        if pts.ndim != 2 or pts.shape[1] < 3:
            return np.empty((0, 2), dtype=float)
        proj_x, proj_y = self._project_xy(pts[:, 2], pts[:, 1])
        return np.column_stack((proj_x, proj_y))

    def _project_layout_polyline(self, z_values, y_values) -> np.ndarray:
        z_arr = np.asarray(z_values, dtype=float).reshape(-1)
        y_arr = np.asarray(y_values, dtype=float).reshape(-1)
        if z_arr.size == 0 or y_arr.size == 0 or z_arr.size != y_arr.size:
            return np.empty((0, 2), dtype=float)
        mask = np.isfinite(z_arr) & np.isfinite(y_arr)
        if not np.any(mask):
            return np.empty((0, 2), dtype=float)
        proj_x, proj_y = self._project_xy(z_arr[mask], y_arr[mask])
        return np.column_stack((proj_x, proj_y))

    def _row_layout_polylines(self, system, row_index: int, z_pos: float) -> list[np.ndarray]:
        if not (0 <= row_index < len(self.rows)):
            return []
        row = self.rows[row_index]
        polylines: list[np.ndarray] = []
        if row.surface == "Mirror":
            half_length = max(float(row.diameter) / 2.0, 0.5)
            angle = np.deg2rad(float(row.tilt_x))
            dz = np.cos(angle) * half_length
            dy = np.sin(angle) * half_length
            center_z = z_pos + float(row.desp_z)
            center_y = float(row.desp_y)
            poly = self._project_layout_polyline(
                [center_z - dz, center_z + dz],
                [center_y - dy, center_y + dy],
            )
            if poly.size > 0:
                polylines.append(poly)
            return polylines
        if row.surface in {"Object", "Image"}:
            half_height = max(float(row.diameter) / 2.0, 0.5)
            center_z = z_pos + float(row.desp_z)
            center_y = float(row.desp_y)
            poly = self._project_layout_polyline(
                [center_z, center_z],
                [center_y - half_height, center_y + half_height],
            )
            if poly.size > 0:
                polylines.append(poly)
            return polylines
        surfaces = getattr(system, "AAA", None)
        surface_data = getattr(system, "SDT_0", None)
        if surfaces is None or surface_data is None:
            return polylines
        if row_index >= getattr(surfaces, "n_blocks", 0) or row_index >= len(surface_data):
            return polylines
        surface = surface_data[row_index]
        if int(getattr(surface, "Drawing", 1)) != 1:
            return polylines
        if str(getattr(surface, "Glass", "") or "").upper() == "NULL":
            return polylines
        solid = 1 if getattr(surface, "Solid_3d_stl", "None") != "None" else 0
        mesh = surfaces[row_index]
        for direction in (1, -1):
            try:
                _ax, ay, az = edge_3d(mesh, direction, 0, 0, solid)
                az, ay = filter_face_2dplot(np.asarray(az, dtype=float), np.asarray(ay, dtype=float), solid)
            except Exception:
                continue
            poly = self._project_layout_polyline(az, ay)
            if int(poly.shape[0]) >= 2:
                polylines.append(poly)
        return polylines

    def _rebuild_layout_pick_regions(self, system) -> None:
        self._layout_pick_regions = {}
        z_pos = 0.0
        for row_index, row in enumerate(self.rows):
            polylines = self._row_layout_polylines(system, row_index, z_pos)
            if polylines:
                self._layout_pick_regions[row_index] = polylines
            z_pos += float(row.thickness)

    @staticmethod
    def _distance_to_polyline(point_xy: np.ndarray, polyline_xy: np.ndarray) -> float:
        pts = np.asarray(polyline_xy, dtype=float)
        point = np.asarray(point_xy, dtype=float)
        if pts.ndim != 2 or pts.shape[0] == 0:
            return float("inf")
        if pts.shape[0] == 1:
            return float(np.linalg.norm(point - pts[0]))
        best = float("inf")
        for start, end in zip(pts[:-1], pts[1:]):
            seg = end - start
            denom = float(np.dot(seg, seg))
            if denom <= 1e-12:
                dist = float(np.linalg.norm(point - start))
            else:
                t = float(np.clip(np.dot(point - start, seg) / denom, 0.0, 1.0))
                proj = start + t * seg
                dist = float(np.linalg.norm(point - proj))
            if dist < best:
                best = dist
        return best

    def _find_layout_pick_row(self, x_display: float, y_display: float) -> int | None:
        if self.ax is None or not self._layout_pick_regions:
            return None
        best_row = None
        best_distance = float("inf")
        threshold = 14.0
        click_xy = np.array([x_display, y_display], dtype=float)
        for row_index, polylines in self._layout_pick_regions.items():
            row_distance = float("inf")
            for polyline in polylines:
                try:
                    display_pts = self.ax.transData.transform(np.asarray(polyline, dtype=float))
                except Exception:
                    continue
                if display_pts.size == 0:
                    continue
                row_distance = min(row_distance, self._distance_to_polyline(click_xy, display_pts))
            if row_distance < best_distance:
                best_distance = row_distance
                best_row = int(row_index)
        if best_distance <= threshold:
            return best_row
        return None

    def _update_layout_selection_overlay(self, row_index: int | None = None) -> None:
        self._clear_layout_selection_overlay()
        if self.ax is None:
            return
        if row_index is None:
            row_index = self._current_selected_row_index()
        if row_index is None:
            return
        polylines = self._layout_pick_regions.get(int(row_index))
        if not polylines:
            return
        artists: list = []
        for polyline in polylines:
            pts = np.asarray(polyline, dtype=float)
            if pts.ndim != 2 or pts.shape[0] == 0:
                continue
            if pts.shape[0] == 1:
                artists.append(
                    self.ax.scatter(
                        pts[:, 0],
                        pts[:, 1],
                        s=55,
                        c="#f97316",
                        edgecolors="white",
                        linewidths=1.4,
                        zorder=950,
                    )
                )
                continue
            underlay, = self.ax.plot(
                pts[:, 0],
                pts[:, 1],
                color="white",
                linewidth=5.0,
                alpha=0.92,
                zorder=940,
            )
            overlay, = self.ax.plot(
                pts[:, 0],
                pts[:, 1],
                color="#f97316",
                linewidth=2.2,
                alpha=0.98,
                zorder=941,
            )
            artists.extend([underlay, overlay])
        self._layout_selection_artists = artists
        self.canvas.draw_idle()

    def _configure_plot_hover_hints(self) -> None:
        self._hover_hint_artists = {}
        self._hover_axis = None
        if hasattr(self, "canvas"):
            self.canvas.get_tk_widget().configure(cursor="")
        candidate_axes = [self.ax]
        if self._analysis_ax is not None:
            candidate_axes.append(self._analysis_ax)
        for axis in candidate_axes:
            if axis is None:
                continue
            highlight = Rectangle(
                (0.0, 0.0),
                1.0,
                1.0,
                transform=axis.transAxes,
                facecolor="#60a5fa",
                edgecolor="#2563eb",
                linewidth=1.0,
                alpha=0.06,
                visible=False,
                zorder=1000,
            )
            axis.add_patch(highlight)
            hint = axis.text(
                0.5,
                0.985,
                "Click to open in viewer",
                transform=axis.transAxes,
                ha="center",
                va="top",
                fontsize=8,
                color="#334155",
                visible=False,
                zorder=1001,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.8},
            )
            self._hover_hint_artists[axis] = (highlight, hint)

    def _set_hover_axis(self, axis) -> None:
        self._hover_axis = axis
        for current_ax, artists in self._hover_hint_artists.items():
            active = current_ax is axis
            for artist in artists:
                artist.set_visible(active)
        if hasattr(self, "canvas"):
            cursor = "hand2" if axis is not None else ""
            self.canvas.get_tk_widget().configure(cursor=cursor)
            self.canvas.draw_idle()

    @staticmethod
    def _viewer_open_command(image_path: Path) -> list[str] | None:
        preferred = os.getenv("KRAKEN_IMAGE_VIEWER", "").strip()
        if preferred:
            parts = preferred.split()
            binary = parts[0]
            if shutil.which(binary):
                return [*parts, str(image_path)]
        for binary in ("nomacs-x11", "nomacs"):
            if shutil.which(binary):
                return [binary, str(image_path)]
        if sys.platform == "darwin":
            return ["open", str(image_path)]
        if os.name == "nt":
            return None
        if shutil.which("xdg-open"):
            return ["xdg-open", str(image_path)]
        if shutil.which("gio"):
            return ["gio", "open", str(image_path)]
        for binary in ("imv", "feh", "eog", "gwenview", "ristretto", "pqiv", "sxiv", "nsxiv"):
            if shutil.which(binary):
                return [binary, str(image_path)]
        return None

    def _open_image_with_system_viewer(self, image_path: Path) -> None:
        if os.name == "nt":
            os.startfile(str(image_path))  # type: ignore[attr-defined]
            return
        command = self._viewer_open_command(image_path)
        if command is None:
            raise RuntimeError("No system image viewer command found.")
        subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)

    def _open_high_res_plot_in_system_viewer(self, target_ax=None) -> None:
        previous_hover_axis = self._hover_axis if self._hover_axis in self._hover_hint_artists else None
        try:
            # Hide hover hint overlays so exported images only contain plot content.
            self._set_hover_axis(None)
            out_dir = Path.home() / "Pictures"
            out_dir.mkdir(parents=True, exist_ok=True)
            stamp = time.strftime("%Y%m%d_%H%M%S")
            axis_label = "analysis" if target_ax is self._analysis_ax else "layout"
            image_path = out_dir / f"kraken_plot_{axis_label}_{stamp}.png"

            self.canvas.draw()
            if target_ax is not None and target_ax in self.figure.axes:
                renderer = self.figure.canvas.get_renderer()
                tight_bbox = target_ax.get_tightbbox(renderer)
                if tight_bbox is not None:
                    # savefig expects bbox_inches in inches, convert from display pixels
                    bbox = tight_bbox.transformed(self.figure.dpi_scale_trans.inverted()).padded(0.08)
                else:
                    fig_w, fig_h = self.figure.get_size_inches()
                    pos = target_ax.get_position()
                    bbox = Bbox.from_extents(
                        float(pos.x0) * fig_w,
                        float(pos.y0) * fig_h,
                        float(pos.x1) * fig_w,
                        float(pos.y1) * fig_h,
                    ).expanded(1.08, 1.12)
                self.figure.savefig(image_path, dpi=320, bbox_inches=bbox)
            else:
                self.figure.savefig(image_path, dpi=320)

            self._open_image_with_system_viewer(image_path)
            self.status_var.set(f"Opened image in system viewer: {image_path.name}")
            self.append_progress(f"Opened high-res image: {image_path}")
        except Exception as exc:
            self.append_debug(f"High-resolution viewer launch failed: {exc}")
        finally:
            if previous_hover_axis is not None:
                self._set_hover_axis(previous_hover_axis)

    def open_3d_view(self) -> None:
        try:
            if vtkTkRenderWindowInteractor is not None:
                try:
                    if self._three_d_inspector is None or not self._three_d_inspector.winfo_exists():
                        self._three_d_inspector = Kraken3DInspector(self)
                    if self._three_d_inspector.available:
                        self._three_d_inspector.deiconify()
                        self._three_d_inspector.lift()
                        self._three_d_inspector.focus_force()
                        self._three_d_inspector.refresh_from_editor()
                        self.status_var.set("Opened Kraken 3D inspector")
                        self.append_debug("Opened Kraken 3D inspector")
                        return
                    reason = self._three_d_inspector.unavailable_reason or "VTK/Tk unavailable"
                    try:
                        self._three_d_inspector.destroy()
                    except Exception:
                        pass
                    self._three_d_inspector = None
                    self.append_debug(f"Embedded 3D inspector unavailable: {reason}. Falling back to legacy PyVista viewer.")
                except Exception as exc:
                    self.append_debug(f"Embedded 3D inspector failed: {exc}. Falling back to legacy PyVista viewer.")
                    if self._three_d_inspector is not None:
                        try:
                            self._three_d_inspector.destroy()
                        except Exception:
                            pass
                        self._three_d_inspector = None

            self._close_legacy_3d_plotter()
            system, rays = self._build_preview_system_and_rays()
            plotter = self._build_legacy_3d_plotter(system, rays)
            plotter.show(auto_close=False, interactive=True, interactive_update=True)
            self._legacy_3d_plotter = plotter
            self._schedule_legacy_3d_poll()
            self.status_var.set("Opened legacy Kraken 3D view")
            self.append_debug("Opened legacy Kraken 3D view")
        except Exception as exc:
            self.append_debug(f"3D view error: {exc}")
            self.status_var.set(f"3D view failed: {exc}")

    def _schedule_legacy_3d_poll(self) -> None:
        if self._legacy_3d_after_id is not None:
            return
        try:
            self._legacy_3d_after_id = self.after(20, self._poll_legacy_3d_plotter)
        except Exception:
            self._legacy_3d_after_id = None

    def _poll_legacy_3d_plotter(self) -> None:
        self._legacy_3d_after_id = None
        plotter = self._legacy_3d_plotter
        if plotter is None:
            return
        if bool(getattr(plotter, "_closed", False)):
            self._legacy_3d_plotter = None
            self.status_var.set("Closed legacy Kraken 3D view")
            return
        try:
            plotter.update(stime=1, force_redraw=True)
        except Exception as exc:
            self.append_debug(f"Legacy 3D update failed: {exc}")
            self._close_legacy_3d_plotter()
            self.status_var.set(f"3D view closed: {_short_error_message(exc)}")
            return
        self._schedule_legacy_3d_poll()

    def _close_legacy_3d_plotter(self) -> None:
        if self._legacy_3d_after_id is not None:
            try:
                self.after_cancel(self._legacy_3d_after_id)
            except Exception:
                pass
            self._legacy_3d_after_id = None
        plotter = self._legacy_3d_plotter
        self._legacy_3d_plotter = None
        if plotter is None:
            return
        try:
            plotter.close()
        except Exception:
            pass

    def _build_preview_system_and_rays(self):
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
                max_radius = max((max(row.diameter / 2.0, 0.5) for row in self.rows), default=1.0)
                self._trace_preview_rays(system, rays, wavelength, max_radius)
        self.append_debug(capture.getvalue())
        self.last_system = system
        self.last_rays = rays
        return system, rays

    def _build_legacy_3d_plotter(self, system, rays):
        plotter = pv.Plotter(shape=(1, 1), title="KrakenOS 3D", notebook=False)
        plotter.set_background("white", top="white")
        plotter.enable_anti_aliasing()
        try:
            plotter.render_window.SetStereoTypeToFake()
            plotter.render_window.StereoRenderOff()
            plotter.render_window.StereoCapableWindowOff()
        except Exception:
            pass
        plotter.add_axes(line_width=3)
        plotter.show_grid(font_size=6, color="black", n_xlabels=2, n_ylabels=2, n_zlabels=2, fmt="%.0f", bold=False)
        scene_info = self._populate_legacy_3d_plotter_scene(
            plotter,
            system,
            rays,
            add_clip_plane=True,
            add_labels=True,
        )
        self._configure_legacy_3d_plotter(
            plotter,
            ray_actors=scene_info["ray_actors"],
            mirror_actors=scene_info["mirror_actors"],
            lens_actors=scene_info["lens_actors"],
            helper_actors=scene_info["helper_actors"],
        )
        self._enable_legacy_3d_close_handling(plotter)
        setattr(plotter, "_kraken_scene", scene_info)
        setattr(plotter, "_kraken_system", system)
        setattr(plotter, "_kraken_rays", rays)
        return plotter

    def _populate_legacy_3d_plotter_scene(
        self,
        plotter,
        system,
        rays,
        *,
        add_clip_plane: bool,
        add_labels: bool,
    ) -> dict[str, list]:
        label_points: list[np.ndarray] = []
        label_text: list[str] = []
        merged_shell = None
        merged_bodies = None
        ray_actors = []
        mirror_actors = []
        lens_actors = []
        helper_actors = []
        actor_row_map: dict[str, int] = {}
        row_actor_map: dict[int, list] = {}

        def register_actor(actor, row_index: int | None = None, *, pickable: bool = False):
            if actor is None:
                return None
            try:
                actor.SetPickable(bool(pickable))
            except Exception:
                pass
            if row_index is None:
                return actor
            actor_key = Kraken3DInspector._actor_key(actor)
            if actor_key is None:
                return actor
            actor_row_map[actor_key] = row_index
            row_actor_map.setdefault(row_index, []).append(actor)
            try:
                actor._kraken_row_index = row_index
                prop = actor.GetProperty()
                if prop is not None:
                    actor._kraken_base_style = {
                        "edge_visibility": int(prop.GetEdgeVisibility()),
                        "line_width": float(prop.GetLineWidth()),
                        "edge_color": tuple(float(v) for v in prop.GetEdgeColor()),
                        "opacity": float(prop.GetOpacity()),
                        "ambient": float(prop.GetAmbient()),
                        "diffuse": float(prop.GetDiffuse()),
                    }
            except Exception:
                pass
            return actor

        transforms = getattr(system, "TRANS_2A", None)
        surfaces = getattr(system, "AAA", None)
        if transforms is not None and surfaces is not None:
            block_count = min(len(self.rows), getattr(surfaces, "n_blocks", 0), len(transforms))
            for index in range(block_count):
                row = self.rows[index]
                if row.surface in {"Object", "Image"}:
                    continue
                mesh = Kraken3DInspector._mesh_with_transform(surfaces[index], transforms[index])
                if mesh is None or int(getattr(mesh, "n_points", 0)) == 0:
                    continue
                if self._legacy_3d_is_stop_plane(row):
                    ring_mesh = self._legacy_3d_stop_ring_mesh(mesh, row)
                    if ring_mesh is not None and int(getattr(ring_mesh, "n_points", 0)) > 0:
                        actor = register_actor(
                            plotter.add_mesh(
                                ring_mesh,
                                color="#f59e0b",
                                opacity=0.95,
                                smooth_shading=False,
                                show_edges=False,
                                pickable=True,
                            ),
                            index,
                            pickable=True,
                        )
                        if actor is not None:
                            helper_actors.append(actor)
                    if add_labels:
                        try:
                            label_points.append(np.mean(np.asarray(mesh.points, dtype=float), axis=0))
                            label_text.append(f"{index}: {row.name}")
                        except Exception:
                            pass
                    continue
                color = Kraken3DInspector._surface_color(system.SDT_0[index])
                opacity = 0.88 if row.surface == "Mirror" else 0.68
                actor = register_actor(
                    plotter.add_mesh(
                        mesh,
                        color=color,
                        opacity=opacity,
                        smooth_shading=True,
                        show_edges=False,
                        pickable=True,
                    ),
                    index,
                    pickable=True,
                )
                if row.surface == "Mirror":
                    mirror_actors.append(actor)
                else:
                    lens_actors.append(actor)
                try:
                    edges = mesh.extract_feature_edges(
                        feature_angle=10,
                        boundary_edges=True,
                        feature_edges=False,
                        manifold_edges=False,
                    )
                    if int(getattr(edges, "n_points", 0)) > 0:
                        plotter.add_mesh(edges, color="#1f2937", line_width=1.0, pickable=False)
                except Exception:
                    pass
                try:
                    merged_shell = mesh.copy(deep=True) if merged_shell is None else merged_shell.merge(mesh)
                except Exception:
                    pass
                if add_labels:
                    try:
                        label_points.append(np.mean(np.asarray(mesh.points, dtype=float), axis=0))
                        label_text.append(f"{index}: {row.name}")
                    except Exception:
                        pass

        side_index = 0
        for row_index in getattr(system, "side_number", []):
            try:
                body = pv.wrap(system.BBB[side_index]).extract_surface().copy(deep=True)
            except Exception:
                side_index += 1
                continue
            side_index += 1
            if int(getattr(body, "n_points", 0)) == 0:
                continue
            color = Kraken3DInspector._surface_color(system.SDT_0[row_index])
            actor = register_actor(
                plotter.add_mesh(
                    body,
                    color=color,
                    opacity=0.18,
                    smooth_shading=False,
                    show_edges=False,
                    pickable=True,
                ),
                row_index,
                pickable=True,
            )
            if self.rows[row_index].surface == "Mirror":
                mirror_actors.append(actor)
            else:
                lens_actors.append(actor)
            try:
                merged_bodies = body.copy(deep=True) if merged_bodies is None else merged_bodies.merge(body)
            except Exception:
                pass

        merged_scene = merged_shell
        if merged_bodies is not None:
            try:
                merged_scene = merged_bodies.copy(deep=True) if merged_scene is None else merged_scene.merge(merged_bodies)
            except Exception:
                pass

        ray_radius = self._legacy_3d_ray_radius(system, rays)
        for wave, ray_pts in zip(getattr(rays, "RayWave", []), getattr(rays, "CC", [])):
            try:
                line = pv.lines_from_points(ray_pts)
            except Exception:
                continue
            if int(getattr(line, "n_points", 0)) < 2:
                continue
            try:
                ray_mesh = line.tube(radius=ray_radius, n_sides=10, capping=True)
            except Exception:
                ray_mesh = line
            actor = plotter.add_mesh(
                ray_mesh,
                color=tuple(wavelength_to_rgb(float(wave) * 1000.0)),
                opacity=0.96,
                pickable=False,
            )
            try:
                actor.SetPickable(False)
            except Exception:
                pass
            ray_actors.append(actor)

        if merged_scene is not None and int(getattr(merged_scene, "n_points", 0)) > 0:
            if add_clip_plane:
                try:
                    actor = plotter.add_mesh_clip_plane(
                        merged_scene,
                        normal=(1.0, 0.0, 0.0),
                        normal_rotation=True,
                        color="#f59e0b",
                        opacity=0.45,
                        name="folded_clip",
                    )
                    try:
                        actor.SetPickable(False)
                    except Exception:
                        pass
                    helper_actors.append(actor)
                except Exception as exc:
                    self.append_debug(f"Legacy 3D clip-plane setup failed: {exc}")
            try:
                slice_mesh = merged_scene.slice(normal=(1.0, 0.0, 0.0), origin=np.asarray(merged_scene.center))
                if int(getattr(slice_mesh, "n_points", 0)) > 0:
                    actor = plotter.add_mesh(
                        slice_mesh,
                        color="#dc2626",
                        line_width=3,
                        name="folded_section",
                        pickable=False,
                    )
                    try:
                        actor.SetPickable(False)
                    except Exception:
                        pass
                    helper_actors.append(actor)
            except Exception as exc:
                self.append_debug(f"Legacy 3D section setup failed: {exc}")

        if add_labels and label_points and hasattr(plotter, "add_point_labels"):
            try:
                plotter.add_point_labels(
                    np.asarray(label_points, dtype=float),
                    label_text,
                    font_size=10,
                    point_size=0,
                    text_color="black",
                    fill_shape=False,
                    shape_opacity=0.0,
                    margin=2,
                    always_visible=True,
                )
            except Exception:
                pass

        return {
            "ray_actors": ray_actors,
            "mirror_actors": mirror_actors,
            "lens_actors": lens_actors,
            "helper_actors": helper_actors,
            "actor_row_map": actor_row_map,
            "row_actor_map": row_actor_map,
        }

    def _build_clean_legacy_3d_plotter(self, system, rays):
        plotter = pv.Plotter(off_screen=True, window_size=(2200, 1400), notebook=False)
        plotter.set_background("white", top="white")
        plotter.enable_anti_aliasing()
        try:
            plotter.render_window.SetStereoTypeToFake()
            plotter.render_window.StereoRenderOff()
            plotter.render_window.StereoCapableWindowOff()
        except Exception:
            pass
        plotter.add_axes(line_width=3)
        plotter.show_grid(font_size=6, color="black", n_xlabels=2, n_ylabels=2, n_zlabels=2, fmt="%.0f", bold=False)
        scene_info = self._populate_legacy_3d_plotter_scene(
            plotter,
            system,
            rays,
            add_clip_plane=False,
            add_labels=False,
        )
        setattr(plotter, "_kraken_scene", scene_info)
        return plotter

    @staticmethod
    def _legacy_3d_is_stop_plane(row: SurfaceRow) -> bool:
        name = (row.name or "").strip().lower()
        return "stop" in name or "aperture" in name

    @staticmethod
    def _legacy_3d_stop_ring_mesh(mesh, row: SurfaceRow):
        try:
            pts = np.asarray(mesh.points, dtype=float)
        except Exception:
            return None
        if pts.size == 0:
            return None
        center = np.mean(pts, axis=0)
        centered = pts - center
        try:
            _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
            normal = np.asarray(vh[-1], dtype=float)
        except Exception:
            normal = np.array([0.0, 0.0, 1.0], dtype=float)
        norm = max(float(np.linalg.norm(normal)), 1e-12)
        normal = normal / norm
        outer = max(float(row.diameter) * 0.5, 0.5)
        inner = max(outer * 0.82, outer - 0.8)
        if inner >= outer:
            inner = outer * 0.9
        try:
            return pv.Disc(center=center, inner=inner, outer=outer, normal=normal, r_res=1, c_res=96)
        except Exception:
            return None

    def _legacy_3d_ray_radius(self, system, rays) -> float:
        try:
            bounds = np.asarray(system.AAA.bounds, dtype=float)
        except Exception:
            bounds = np.array([0.0, 100.0, -25.0, 25.0, -25.0, 25.0], dtype=float)
        span = max(
            float(bounds[1] - bounds[0]) if bounds.size >= 2 else 0.0,
            float(bounds[3] - bounds[2]) if bounds.size >= 4 else 0.0,
            float(bounds[5] - bounds[4]) if bounds.size >= 6 else 0.0,
            1.0,
        )
        ray_count = max(len(getattr(rays, "CC", [])), 1)
        return max(span * 0.0009 / max(np.sqrt(ray_count), 1.0), 0.05)

    def _configure_legacy_3d_plotter(
        self,
        plotter,
        *,
        ray_actors=None,
        mirror_actors=None,
        lens_actors=None,
        helper_actors=None,
    ) -> None:
        help_lines = [
            "KrakenOS 3D",
            "Click a surface to select its row in the editor",
            "Drag the orange clipping plane through the folded system",
            "Keys: I Iso  Y YZ  T Top  X XZ  H Home  K Save PNG  Q Close",
        ]
        plotter.add_text("\n".join(help_lines), position="upper_left", font_size=12, color="royalblue")
        self._set_legacy_3d_camera(plotter, "iso")
        plotter.add_key_event("i", lambda: self._set_legacy_3d_camera(plotter, "iso"))
        plotter.add_key_event("I", lambda: self._set_legacy_3d_camera(plotter, "iso"))
        plotter.add_key_event("y", lambda: self._set_legacy_3d_camera(plotter, "yz"))
        plotter.add_key_event("Y", lambda: self._set_legacy_3d_camera(plotter, "yz"))
        plotter.add_key_event("t", lambda: self._set_legacy_3d_camera(plotter, "top"))
        plotter.add_key_event("T", lambda: self._set_legacy_3d_camera(plotter, "top"))
        plotter.add_key_event("x", lambda: self._set_legacy_3d_camera(plotter, "xz"))
        plotter.add_key_event("X", lambda: self._set_legacy_3d_camera(plotter, "xz"))
        plotter.add_key_event("h", lambda: self._set_legacy_3d_camera(plotter, "reset"))
        plotter.add_key_event("H", lambda: self._set_legacy_3d_camera(plotter, "reset"))
        plotter.add_key_event("k", lambda: self._save_legacy_3d_screenshot(plotter))
        plotter.add_key_event("K", lambda: self._save_legacy_3d_screenshot(plotter))
        plotter.add_key_event("q", self._close_legacy_3d_plotter)
        plotter.add_key_event("Q", self._close_legacy_3d_plotter)

        ray_actors = list(ray_actors or [])
        mirror_actors = list(mirror_actors or [])
        lens_actors = list(lens_actors or [])
        helper_actors = list(helper_actors or [])
        setattr(
            plotter,
            "_kraken_visibility",
            {
                "rays": True,
                "mirrors": True,
                "lenses": True,
                "helpers": True,
            },
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Save PNG",
            position=(10, 10),
            callback=lambda _state: self._save_legacy_3d_screenshot(plotter),
            color="#0f766e",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Close",
            position=(10, 50),
            callback=lambda _state: self._close_legacy_3d_plotter(),
            color="#991b1b",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Iso",
            position=(10, 90),
            callback=lambda _state: self._set_legacy_3d_camera(plotter, "iso"),
            color="#475569",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="YZ",
            position=(10, 130),
            callback=lambda _state: self._set_legacy_3d_camera(plotter, "yz"),
            color="#475569",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Top",
            position=(10, 170),
            callback=lambda _state: self._set_legacy_3d_camera(plotter, "top"),
            color="#475569",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="XZ",
            position=(10, 210),
            callback=lambda _state: self._set_legacy_3d_camera(plotter, "xz"),
            color="#475569",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Home",
            position=(10, 250),
            callback=lambda _state: self._set_legacy_3d_camera(plotter, "reset"),
            color="#475569",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Rays",
            position=(10, 310),
            callback=lambda state: self._legacy_3d_set_actor_visibility(ray_actors, state, plotter, "rays"),
            value=True,
            color="#2563eb",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Mirrors",
            position=(10, 350),
            callback=lambda state: self._legacy_3d_set_actor_visibility(mirror_actors, state, plotter, "mirrors"),
            value=True,
            color="#6b7280",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Lenses",
            position=(10, 390),
            callback=lambda state: self._legacy_3d_set_actor_visibility(lens_actors, state, plotter, "lenses"),
            value=True,
            color="#0284c7",
        )
        self._add_legacy_3d_action_button(
            plotter,
            label="Helpers",
            position=(10, 430),
            callback=lambda state: self._legacy_3d_set_actor_visibility(helper_actors, state, plotter, "helpers"),
            value=True,
            color="#f59e0b",
        )
        self._enable_legacy_3d_picking(plotter)

    def _enable_legacy_3d_close_handling(self, plotter) -> None:
        def request_close(*_args):
            try:
                self.after(0, self._close_legacy_3d_plotter)
            except Exception:
                self._close_legacy_3d_plotter()

        try:
            if getattr(plotter, "iren", None) is not None:
                plotter.iren.add_observer("ExitEvent", request_close)
        except Exception as exc:
            self.append_debug(f"Legacy 3D close observer unavailable: {exc}")
        try:
            if getattr(plotter, "render_window", None) is not None:
                plotter.render_window.AddObserver("DeleteEvent", request_close)
        except Exception:
            pass

    def _enable_legacy_3d_picking(self, plotter) -> None:
        if vtkCellPicker is None or getattr(plotter, "iren", None) is None:
            return
        try:
            picker = vtkCellPicker()
            picker.SetTolerance(0.0005)
            setattr(plotter, "_kraken_picker", picker)
            setattr(plotter, "_kraken_selected_row", None)
            plotter.iren.add_observer(
                "LeftButtonPressEvent",
                lambda *_args: self._legacy_3d_pick_click(plotter),
            )
        except Exception as exc:
            self.append_debug(f"Legacy 3D picking unavailable: {exc}")

    def _legacy_3d_pick_click(self, plotter) -> None:
        picker = getattr(plotter, "_kraken_picker", None)
        if picker is None or getattr(plotter, "iren", None) is None:
            return
        try:
            x, y = plotter.iren.get_event_position()
            renderer = plotter.iren.get_poked_renderer(x, y)
            if renderer is None:
                renderer = getattr(plotter, "renderer", None)
            if renderer is None:
                return
            picker.Pick(x, y, 0.0, renderer)
            actor = picker.GetActor()
        except Exception as exc:
            self.append_debug(f"Legacy 3D pick failed: {exc}")
            return
        actor_key = Kraken3DInspector._actor_key(actor)
        scene_info = dict(getattr(plotter, "_kraken_scene", {}) or {})
        row_index = scene_info.get("actor_row_map", {}).get(actor_key) if actor_key is not None else None
        if row_index is None:
            self._legacy_3d_set_selected_row(plotter, None)
            self.status_var.set("3D view ready")
            return
        self._legacy_3d_set_selected_row(plotter, int(row_index))
        self._select_table_row(int(row_index))
        row_name = self.rows[int(row_index)].name if 0 <= int(row_index) < len(self.rows) else "Surface"
        self.status_var.set(f"3D selected row {int(row_index)}: {row_name}")

    def _legacy_3d_set_selected_row(self, plotter, row_index: int | None) -> None:
        current = getattr(plotter, "_kraken_selected_row", None)
        if current == row_index:
            return
        scene_info = dict(getattr(plotter, "_kraken_scene", {}) or {})
        row_actor_map = dict(scene_info.get("row_actor_map", {}) or {})
        if current is not None:
            for actor in row_actor_map.get(int(current), []):
                self._legacy_3d_set_actor_highlight(actor, False)
        if row_index is not None:
            for actor in row_actor_map.get(int(row_index), []):
                self._legacy_3d_set_actor_highlight(actor, True)
        setattr(plotter, "_kraken_selected_row", row_index)
        try:
            plotter.render()
        except Exception:
            pass

    @staticmethod
    def _legacy_3d_set_actor_highlight(actor, selected: bool) -> None:
        if actor is None:
            return
        try:
            prop = actor.GetProperty()
        except Exception:
            return
        if prop is None:
            return
        base = getattr(actor, "_kraken_base_style", None)
        if not isinstance(base, dict):
            try:
                base = {
                    "edge_visibility": int(prop.GetEdgeVisibility()),
                    "line_width": float(prop.GetLineWidth()),
                    "edge_color": tuple(float(v) for v in prop.GetEdgeColor()),
                    "opacity": float(prop.GetOpacity()),
                    "ambient": float(prop.GetAmbient()),
                    "diffuse": float(prop.GetDiffuse()),
                }
                actor._kraken_base_style = base
            except Exception:
                base = {}
        if selected:
            try:
                prop.SetEdgeVisibility(1)
                prop.SetEdgeColor(1.0, 0.45, 0.05)
                prop.SetLineWidth(max(float(base.get("line_width", 1.0)), 2.5))
                prop.SetOpacity(min(max(float(base.get("opacity", 1.0)), 0.3) + 0.08, 1.0))
                prop.SetAmbient(max(float(base.get("ambient", 0.0)), 0.18))
            except Exception:
                return
            return
        try:
            prop.SetEdgeVisibility(int(base.get("edge_visibility", 0)))
            edge_color = tuple(base.get("edge_color", (0.0, 0.0, 0.0)))
            if len(edge_color) == 3:
                prop.SetEdgeColor(*edge_color)
            prop.SetLineWidth(float(base.get("line_width", 1.0)))
            prop.SetOpacity(float(base.get("opacity", 1.0)))
            prop.SetAmbient(float(base.get("ambient", 0.0)))
            prop.SetDiffuse(float(base.get("diffuse", 1.0)))
        except Exception:
            pass

    def _add_legacy_3d_action_button(
        self,
        plotter,
        *,
        label: str,
        position: tuple[int, int],
        callback,
        value: bool = False,
        color: str = "#2563eb",
    ) -> None:
        plotter.add_checkbox_button_widget(
            callback,
            value=value,
            position=position,
            size=26,
            border_size=4,
            color_on=color,
            color_off=color,
            background_color="white",
        )
        plotter.add_text(
            label,
            position=(position[0] + 34, position[1] + 2),
            font_size=10,
            color="black",
        )

    def _legacy_3d_set_actor_visibility(self, actors, visible: bool, plotter, group: str | None = None) -> None:
        if group:
            state = dict(getattr(plotter, "_kraken_visibility", {}) or {})
            state[group] = bool(visible)
            setattr(plotter, "_kraken_visibility", state)
        for actor in actors:
            try:
                actor.SetVisibility(bool(visible))
            except Exception:
                continue
        plotter.render()

    def _save_legacy_3d_screenshot(self, plotter) -> None:
        try:
            default_name = f"kraken_3d_{time.strftime('%Y%m%d_%H%M%S')}.png"
            selected_path = filedialog.asksaveasfilename(
                parent=self,
                title="Save 3D view as PNG",
                initialdir=str(Path.home() / "Pictures"),
                initialfile=default_name,
                defaultextension=".png",
                filetypes=[("PNG image", "*.png")],
            )
            if not selected_path:
                self.status_var.set("3D PNG save cancelled")
                return
            image_path = Path(selected_path)
            system = getattr(plotter, "_kraken_system", None)
            rays = getattr(plotter, "_kraken_rays", None)
            if system is None or rays is None:
                raise RuntimeError("3D scene data unavailable for clean screenshot")
            clean_plotter = self._build_clean_legacy_3d_plotter(system, rays)
            clean_scene = dict(getattr(clean_plotter, "_kraken_scene", {}) or {})
            visibility = dict(getattr(plotter, "_kraken_visibility", {}) or {})
            if not visibility.get("rays", True):
                self._legacy_3d_set_actor_visibility(clean_scene.get("ray_actors", []), False, clean_plotter)
            if not visibility.get("mirrors", True):
                self._legacy_3d_set_actor_visibility(clean_scene.get("mirror_actors", []), False, clean_plotter)
            if not visibility.get("lenses", True):
                self._legacy_3d_set_actor_visibility(clean_scene.get("lens_actors", []), False, clean_plotter)
            if not visibility.get("helpers", True):
                self._legacy_3d_set_actor_visibility(clean_scene.get("helper_actors", []), False, clean_plotter)
            try:
                clean_plotter.camera_position = plotter.camera_position
            except Exception:
                self._set_legacy_3d_camera(clean_plotter, "iso")
            clean_plotter.screenshot(str(image_path))
            try:
                clean_plotter.close()
            except Exception:
                pass
            self.status_var.set(f"Saved 3D PNG: {image_path.name}")
            self.append_progress(f"Saved 3D PNG: {image_path}")
        except Exception as exc:
            self.append_debug(f"3D screenshot failed: {exc}")
            self.status_var.set(f"3D screenshot failed: {_short_error_message(exc)}")

    @staticmethod
    def _set_legacy_3d_camera(plotter, preset: str) -> None:
        cx, cy, cz = plotter.center
        bounds = np.asarray(plotter.bounds, dtype=float)
        span_x = float(bounds[1] - bounds[0]) if bounds.size >= 2 else 0.0
        span_y = float(bounds[3] - bounds[2]) if bounds.size >= 4 else 0.0
        span_z = float(bounds[5] - bounds[4]) if bounds.size >= 6 else 0.0
        span = max(span_x, span_y, span_z, 1.0)
        distance = max(span * 2.6, 180.0)
        if preset == "yz":
            position = (cx + distance, cy, cz)
            view_up = (0.0, 0.0, 1.0)
            parallel_scale = max(span_z, span_y, 1.0) * 0.55
        elif preset == "top":
            # Exact top view. This used to conflict only because VTK binds the numeric `3` key to stereo.
            position = (cx, cy, cz + distance)
            view_up = (0.0, 1.0, 0.0)
            parallel_scale = max(span_x, span_y, 1.0) * 0.55
        elif preset == "xz":
            position = (cx, cy + distance, cz)
            view_up = (0.0, 0.0, 1.0)
            parallel_scale = max(span_x, span_z, 1.0) * 0.55
        else:
            position = (cx - distance, cy + distance * 0.55, cz + distance * 0.75)
            view_up = (0.0, 1.0, 0.0)
            parallel_scale = None
        plotter.camera_position = [position, (cx, cy, cz), view_up]
        try:
            plotter.camera.parallel_projection = parallel_scale is not None
            if parallel_scale is not None:
                plotter.camera.parallel_scale = float(parallel_scale)
        except Exception:
            pass
        plotter.reset_camera_clipping_range()
        plotter.set_background("white", top="white")
        plotter.render()

    def _refresh_3d_inspector_if_open(self) -> None:
        if self._three_d_inspector is None:
            return
        try:
            if not self._three_d_inspector.winfo_exists():
                self._three_d_inspector = None
                return
        except Exception:
            self._three_d_inspector = None
            return
        try:
            if self.last_system is None or self.last_rays is None:
                return
            row_names = [row.name for row in self.rows]
            self._three_d_inspector.refresh_scene(self.last_system, self.last_rays, row_names, reset_camera=False)
        except Exception as exc:
            self.append_debug(f"3D inspector sync failed: {exc}")

    def _current_selected_row_index(self) -> int | None:
        items = self.table.selection()
        if not items:
            return None
        try:
            return int(self.table.index(items[0]))
        except Exception:
            return None

    def _on_table_selection_changed(self, _event: tk.Event | None = None) -> None:
        self._sync_surface_selection(self._current_selected_row_index(), from_table=True)

    def _sync_surface_selection(self, row_index: int | None, *, from_table: bool = False) -> None:
        if self._three_d_inspector is not None:
            try:
                if self._three_d_inspector.winfo_exists() and self._three_d_inspector.available:
                    self._three_d_inspector.highlight_row(row_index)
            except Exception:
                pass
        if self._legacy_3d_plotter is not None:
            try:
                self._legacy_3d_set_selected_row(self._legacy_3d_plotter, row_index)
            except Exception:
                pass
        self._update_layout_selection_overlay(row_index)
        if from_table and row_index is not None and 0 <= row_index < len(self.rows):
            self.status_var.set(f"Selected row {row_index}: {self.rows[row_index].name}")

    def _select_table_row(self, index: int) -> None:
        items = self.table.get_children()
        if not (0 <= index < len(items)):
            return
        row_id = items[index]
        self.table.selection_set(row_id)
        self.table.focus(row_id)
        self.table.see(row_id)
        self._active_cell = (row_id, "#1")
        self._update_active_cell_border()
        self._sync_surface_selection(index)

    def load_layout_by_name(self, name: str) -> None:
        path = self.layout_files.get(name)
        if path is None:
            return
        if self.rows:
            self._begin_history_capture()
        self.current_layout_file = path
        had_existing_rows = bool(self.rows)
        info: dict[str, object] = {"surfaces": [], "settings": {}}
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
        insert_after = self._selected_insert_index()
        if had_existing_rows:
            self.rows = self._append_layout_rows(self.rows, loaded_rows, insert_after=insert_after)
        else:
            self.rows = loaded_rows
            self._apply_initial_field_defaults()
            self._apply_initial_layout_view_defaults(name)

        self._normalize_special_rows()
        self._sync_table()
        if not had_existing_rows:
            self._apply_layout_settings(info.get("settings", {}))
        self._select_inserted_layout_rows(loaded_rows, insert_after=insert_after)
        if had_existing_rows:
            self._commit_history_capture()
        self.refresh_plot()
        self.layout_var.set(name)
        self.example_var.set("Examples")
        self.status_var.set(f"Appended {name}")

    def _selected_operand_labels(self) -> list[str]:
        if not hasattr(self, "merit_mode_list"):
            return []
        return [self.merit_mode_list.get(i) for i in self.merit_mode_list.curselection()]

    def _set_selected_operand_labels(self, labels: list[str]) -> None:
        if not hasattr(self, "merit_mode_list"):
            return
        self.merit_mode_list.selection_clear(0, "end")
        wanted = {str(label) for label in labels}
        for index in range(self.merit_mode_list.size()):
            label = self.merit_mode_list.get(index)
            if label in wanted:
                self.merit_mode_list.selection_set(index)
        self._update_operand_setup_visibility()

    def _capture_editor_state(self) -> dict[str, object]:
        selected_indices = []
        if hasattr(self, "table"):
            try:
                selected_indices = [int(self.table.index(item)) for item in self.table.selection()]
            except Exception:
                selected_indices = []
        active_cell = None
        if self._active_cell is not None:
            row_id, field = self._active_cell
            try:
                active_cell = {"row": int(self.table.index(row_id)), "field": str(field)}
            except Exception:
                active_cell = None
        layout_path = str(self.current_layout_file) if self.current_layout_file is not None else None
        return {
            "rows": [asdict(row) for row in self.rows],
            "settings": self._collect_layout_settings(),
            "selected_indices": selected_indices,
            "active_cell": active_cell,
            "current_layout_file": layout_path,
        }

    def _begin_history_capture(self, _event: tk.Event | None = None) -> None:
        if self._history_restoring or self._history_pending_state is not None:
            return
        self._history_pending_state = self._capture_editor_state()

    def _commit_history_capture(self) -> None:
        if self._history_restoring:
            self._history_pending_state = None
            return
        snapshot = self._history_pending_state
        self._history_pending_state = None
        if snapshot is None:
            return
        current = self._capture_editor_state()
        if snapshot == current:
            return
        self._undo_stack.append(snapshot)
        if len(self._undo_stack) > self._history_limit:
            self._undo_stack = self._undo_stack[-self._history_limit :]
        self._redo_stack.clear()
        self._update_undo_redo_buttons()

    def _push_history_snapshot(self) -> None:
        if self._history_restoring:
            return
        self._history_pending_state = self._capture_editor_state()
        self._commit_history_capture()

    def _restore_history_state(self, state: dict[str, object]) -> None:
        self._history_restoring = True
        try:
            rows = state.get("rows", [])
            restored_rows = [SurfaceRow(**dict(item)) for item in rows if isinstance(item, dict)]
            self.rows = self._normalized_rows_copy(restored_rows)
            layout_path = state.get("current_layout_file")
            self.current_layout_file = Path(layout_path) if isinstance(layout_path, str) and layout_path else None
            self._sync_table()
            self._apply_layout_settings(state.get("settings", {}))
            self._normalize_special_rows()
            self._sync_table()
            selected_indices = [int(index) for index in state.get("selected_indices", []) if isinstance(index, int)]
            items = list(self.table.get_children())
            selected_items = [items[index] for index in selected_indices if 0 <= index < len(items)]
            if selected_items:
                self.table.selection_set(selected_items)
                self.table.focus(selected_items[0])
                self.table.see(selected_items[0])
            else:
                self.table.selection_remove(*items)
            active_cell = state.get("active_cell")
            self._active_cell = None
            if isinstance(active_cell, dict):
                row_index = int(active_cell.get("row", -1))
                field = str(active_cell.get("field", ""))
                if 0 <= row_index < len(items) and field in FIELDS:
                    self._active_cell = (items[row_index], field)
            self._update_active_cell_border()
            self._refresh_analysis_surface_choices()
            self._refresh_operand_surface_choices()
        finally:
            self._history_restoring = False
            self._history_pending_state = None
        self.refresh_plot()
        self._update_undo_redo_buttons()

    def _update_undo_redo_buttons(self) -> None:
        if self._undo_button is not None:
            self._undo_button.configure(state=("normal" if self._undo_stack else "disabled"))
        if self._redo_button is not None:
            self._redo_button.configure(state=("normal" if self._redo_stack else "disabled"))

    def undo(self) -> None:
        if not self._undo_stack:
            return
        current = self._capture_editor_state()
        state = self._undo_stack.pop()
        self._redo_stack.append(current)
        self._restore_history_state(state)
        self.status_var.set("Undo applied.")

    def redo(self) -> None:
        if not self._redo_stack:
            return
        current = self._capture_editor_state()
        state = self._redo_stack.pop()
        self._undo_stack.append(current)
        self._restore_history_state(state)
        self.status_var.set("Redo applied.")

    def _undo_event(self, _event=None) -> str:
        self.undo()
        return "break"

    def _redo_event(self, _event=None) -> str:
        self.redo()
        return "break"

    def _collect_layout_settings(self) -> dict[str, object]:
        operand_settings: dict[str, dict[str, object]] = {}
        all_labels = {spec.label for spec in OPERAND_REGISTRY.values()}
        for label in all_labels:
            payload: dict[str, object] = {}
            var = self.operand_weight_vars.get(label)
            if var is not None:
                payload["weight"] = var.get()
            var = self.operand_target_vars.get(label)
            if var is not None:
                payload["target"] = var.get()
            var = self.operand_wavelength_vars.get(label)
            if var is not None:
                payload["wavelength"] = var.get()
            var = self.operand_field_vars.get(label)
            if var is not None:
                payload["field"] = var.get()
            var = self.operand_field_x_vars.get(label)
            if var is not None:
                payload["field_x"] = var.get()
            var = self.operand_field_y_vars.get(label)
            if var is not None:
                payload["field_y"] = var.get()
            var = self.operand_surface_vars.get(label)
            if var is not None:
                payload["surface"] = var.get()
            var = self.operand_aperture_type_vars.get(label)
            if var is not None:
                payload["aperture_type"] = var.get()
            var = self.operand_aperture_value_vars.get(label)
            if var is not None:
                payload["aperture_value"] = var.get()
            var = self.operand_frequency_vars.get(label)
            if var is not None:
                payload["frequency"] = var.get()
            var = self.operand_mtf_mode_vars.get(label)
            if var is not None:
                payload["mtf_mode"] = var.get()
            var = self.operand_mtf_algorithm_vars.get(label)
            if var is not None:
                payload["mtf_algorithm"] = var.get()
            if payload:
                operand_settings[label] = payload

        return {
            "object_mode": self.object_mode_var.get().strip(),
            "display_orientation": self.display_orientation_var.get().strip(),
            "wavelength": self.wavelength_var.get().strip(),
            "ray_count": self.ray_count_var.get().strip(),
            "ray_height_factor": self.ray_height_factor_var.get().strip(),
            "analysis_surface": self.analysis_surface_var.get().strip(),
            "aperture_type": self.aperture_type_var.get().strip(),
            "aperture_value": self.aperture_value_var.get().strip(),
            "spot_view_mode": self.spot_view_mode_var.get().strip(),
            "show_clipped_rays": bool(self.show_clipped_rays_var.get()),
            "show_cardinals": bool(self.show_cardinals_var.get()),
            "field_type": self.field_type_var.get().strip(),
            "field_value": self.field_value_var.get().strip(),
            "field_count": self.field_count_var.get().strip(),
            "image_diameter_mode": self.image_diameter_mode_var.get().strip() if hasattr(self, "image_diameter_mode_var") else "Auto",
            "analysis_mode": str(self.analysis_mode or "none").strip(),
            "auto_save_plot": bool(self.auto_save_plot_var.get()),
            "show_native_overlays": bool(self.show_native_overlays_var.get()),
            "show_native_active_spans": bool(self.show_native_active_spans_var.get()),
            "show_native_hit_labels": bool(self.show_native_hit_labels_var.get()),
            "optimization_workers": self.optimization_workers_var.get().strip() if hasattr(self, "optimization_workers_var") else "Auto",
            "selected_operands": self._selected_operand_labels(),
            "operands": operand_settings,
        }

    def _apply_layout_settings(self, settings: object) -> None:
        if not isinstance(settings, dict):
            return

        def _parse_bool(value) -> bool:
            if isinstance(value, str):
                return value.strip().lower() in {"1", "true", "yes", "on"}
            return bool(value)

        def _set_text(var, key: str) -> None:
            value = settings.get(key)
            if value is not None:
                var.set(str(value).strip())

        display_orientation = str(settings.get("display_orientation", "")).strip()
        if display_orientation in {"Vertical", "Horizontal"}:
            self.display_orientation_var.set(display_orientation)

        object_mode = str(settings.get("object_mode", "")).strip()
        if object_mode in {"Finite", "Infinity"}:
            self.object_mode_var.set(object_mode)

        _set_text(self.wavelength_var, "wavelength")
        _set_text(self.ray_count_var, "ray_count")
        _set_text(self.ray_height_factor_var, "ray_height_factor")

        aperture_type = str(settings.get("aperture_type", "")).strip().upper()
        if aperture_type in {"STOP", "EPD"}:
            self.aperture_type_var.set(aperture_type)
        _set_text(self.aperture_value_var, "aperture_value")

        if "spot_view_mode" in settings:
            spot_view_mode = str(settings.get("spot_view_mode", "")).strip()
            if spot_view_mode in {"Grid", "Absolute", "Centroid"}:
                self.spot_view_mode_var.set(spot_view_mode)

        field_count = settings.get("field_count")
        if field_count is not None:
            self.field_count_var.set(str(field_count).strip())

        image_diameter_mode = str(settings.get("image_diameter_mode", "")).strip()
        if image_diameter_mode in {"Auto", "Manual"} and hasattr(self, "image_diameter_mode_var"):
            self.image_diameter_mode_var.set(image_diameter_mode)

        self._sync_field_mode_ui()

        field_type = str(settings.get("field_type", "")).strip()
        if field_type in {"Angle", "Object Height", "Paraxial Image Height", "Real Image Height"}:
            self.field_type_var.set(field_type)
            self._last_field_type = field_type

        field_value = settings.get("field_value")
        if field_value is not None:
            field_value_text = str(field_value).strip()
            self.field_value_var.set(field_value_text)
            self._field_type_defaults[self._current_field_type()] = field_value_text

        analysis_surface = settings.get("analysis_surface")
        if analysis_surface is not None:
            analysis_surface_text = str(analysis_surface).strip()
            if analysis_surface_text == "Auto":
                self.analysis_surface_var.set("Auto")
            elif analysis_surface_text in set(self.analysis_surface_menu["values"]):
                self.analysis_surface_var.set(analysis_surface_text)
            else:
                try:
                    index = int(float(analysis_surface_text))
                except (TypeError, ValueError):
                    index = -1
                if 0 <= index < len(self.rows):
                    self.analysis_surface_var.set(f"{index}: {self.rows[index].name}")

        bool_vars = {
            "show_clipped_rays": self.show_clipped_rays_var,
            "show_cardinals": self.show_cardinals_var,
            "auto_save_plot": self.auto_save_plot_var,
            "show_native_overlays": self.show_native_overlays_var,
            "show_native_active_spans": self.show_native_active_spans_var,
            "show_native_hit_labels": self.show_native_hit_labels_var,
        }
        for key, var in bool_vars.items():
            if key in settings:
                var.set(_parse_bool(settings.get(key)))

        if hasattr(self, "optimization_workers_var") and "optimization_workers" in settings:
            worker_text = str(settings.get("optimization_workers", "Auto")).strip() or "Auto"
            self.optimization_workers_var.set(worker_text)

        selected_operands = settings.get("selected_operands")
        if isinstance(selected_operands, (list, tuple)):
            self._set_selected_operand_labels([str(label) for label in selected_operands])

        operand_settings = settings.get("operands")
        if isinstance(operand_settings, dict):
            for label, payload in operand_settings.items():
                if not isinstance(payload, dict):
                    continue
                if label in self.operand_weight_vars and "weight" in payload:
                    self.operand_weight_vars[label].set(str(payload["weight"]).strip())
                if label in self.operand_target_vars and "target" in payload:
                    self.operand_target_vars[label].set(str(payload["target"]).strip())
                if label in self.operand_wavelength_vars and "wavelength" in payload:
                    self.operand_wavelength_vars[label].set(str(payload["wavelength"]).strip())
                if label in self.operand_field_vars and "field" in payload:
                    self.operand_field_vars[label].set(str(payload["field"]).strip())
                if label in self.operand_field_x_vars and "field_x" in payload:
                    self.operand_field_x_vars[label].set(str(payload["field_x"]).strip())
                if label in self.operand_field_y_vars and "field_y" in payload:
                    self.operand_field_y_vars[label].set(str(payload["field_y"]).strip())
                if label in self.operand_surface_vars and "surface" in payload:
                    surface_text = str(payload["surface"]).strip()
                    if surface_text:
                        self.operand_surface_vars[label].set(surface_text)
                if label in self.operand_aperture_type_vars and "aperture_type" in payload:
                    aperture_text = str(payload["aperture_type"]).strip().upper()
                    if aperture_text in {"STOP", "EPD"}:
                        self.operand_aperture_type_vars[label].set(aperture_text)
                if label in self.operand_aperture_value_vars and "aperture_value" in payload:
                    self.operand_aperture_value_vars[label].set(str(payload["aperture_value"]).strip())
                if label in self.operand_frequency_vars and "frequency" in payload:
                    self.operand_frequency_vars[label].set(str(payload["frequency"]).strip())
                if label in self.operand_mtf_mode_vars and "mtf_mode" in payload:
                    mode_text = str(payload["mtf_mode"]).strip()
                    if mode_text in {"Average", "Tangential", "Sagittal"}:
                        self.operand_mtf_mode_vars[label].set(mode_text)
                if label in self.operand_mtf_algorithm_vars and "mtf_algorithm" in payload:
                    algorithm_text = str(payload["mtf_algorithm"]).strip()
                    if algorithm_text in {"PSF FFT", "LSF FFT"}:
                        self.operand_mtf_algorithm_vars[label].set(algorithm_text)

        analysis_mode = str(settings.get("analysis_mode", "")).strip()
        if analysis_mode in {"none", "native_off_axis", "spot", "psf", "rms", "field_curvature", "pupil", "seidel", "wavefront", "mtf"}:
            self.analysis_mode = analysis_mode
            if hasattr(self, "analysis_mode_button_var"):
                self.analysis_mode_button_var.set(analysis_mode)

        self._sync_object_controls()
        self._update_field_status_hint()

    def load_example_by_name(self, name: str) -> None:
        path = self.example_files.get(name)
        if path is None:
            return
        if self.rows:
            self._begin_history_capture()
        try:
            surfaces = self._extract_surfaces_from_example(path)
        except Exception as exc:
            self._history_pending_state = None
            self.status_var.set(f"Failed to load example {name}: {exc}")
            return
        self.current_layout_file = None
        self.rows = [self._row_from_surface(surface, index, len(surfaces)) for index, surface in enumerate(surfaces)]
        self._normalize_special_rows()
        self._apply_example_display_defaults(path)
        self._sync_table()
        self._commit_history_capture()
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
                row.glass,
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

    def _selected_insert_index(self) -> int | None:
        selected = self.table.selection()
        if not selected:
            return None
        indices = sorted(self.table.index(item) for item in selected)
        if not indices:
            return None
        return indices[-1]

    def _select_inserted_layout_rows(self, layout_rows: list[SurfaceRow], insert_after: int | None) -> None:
        additions = max(len(layout_rows) - 2, 0)
        if additions <= 0:
            return
        if insert_after is None:
            insert_at = len(self.rows) - additions
            if self.rows and self.rows[-1].surface == "Image":
                insert_at -= 1
        else:
            insert_at = min(insert_after + 1, len(self.rows) - additions)
        items = self.table.get_children()
        new_items = items[insert_at : insert_at + additions]
        if new_items:
            self.table.selection_set(new_items)
            self.table.focus(new_items[0])
            self.table.see(new_items[0])

    @staticmethod
    def _append_layout_rows(
        existing_rows: list[SurfaceRow], layout_rows: list[SurfaceRow], insert_after: int | None = None
    ) -> list[SurfaceRow]:
        base = [SurfaceRow(**asdict(row)) for row in existing_rows]
        additions = [SurfaceRow(**asdict(row)) for row in layout_rows[1:-1]]
        if not additions:
            return base
        if insert_after is None:
            insert_at = len(base)
            if base and base[-1].surface == "Image":
                insert_at -= 1
        else:
            insert_at = max(0, min(insert_after + 1, len(base)))
            if base and base[-1].surface == "Image":
                insert_at = min(insert_at, len(base) - 1)
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
                    glass=str(values[3]),
                    optimize_rc=self.rows[len(rows)].optimize_rc if len(rows) < len(self.rows) else False,
                    optimize_rc_bounds=self.rows[len(rows)].optimize_rc_bounds if len(rows) < len(self.rows) else None,
                    rc=self._parse_numeric_display(str(values[4])),
                    optimize_thickness=self.rows[len(rows)].optimize_thickness if len(rows) < len(self.rows) else False,
                    optimize_thickness_bounds=self.rows[len(rows)].optimize_thickness_bounds if len(rows) < len(self.rows) else None,
                    thickness=self._parse_numeric_display(str(values[5])),
                    diameter=float(values[6]),
                    tilt_x=float(values[7]),
                    tilt_y=float(values[8]),
                    tilt_z=float(values[9]),
                    desp_x=float(values[10]),
                    desp_y=float(values[11]),
                    desp_z=float(values[12]),
                    axis_move=float(values[13]),
                )
            )
        self.rows = rows
        self._sync_object_controls()

    def _on_table_click(self, event: tk.Event) -> str | None:
        row_id = self.table.identify_row(event.y)
        column_id = self.table.identify_column(event.x)
        self.table.focus_set()
        if not row_id or not column_id:
            self._active_cell = None
            self._hide_active_cell_border()
            return None
        self._active_cell = (row_id, column_id)
        children = list(self.table.get_children())
        shift_pressed = bool(event.state & 0x0001)
        control_pressed = bool(event.state & 0x0004)
        if shift_pressed and children:
            anchor = self._selection_anchor_row
            if anchor not in children:
                anchor = self.table.focus() or row_id
            if anchor not in children:
                anchor = row_id
            start = children.index(anchor)
            end = children.index(row_id)
            if start <= end:
                selected = children[start : end + 1]
            else:
                selected = children[end : start + 1]
            self.table.selection_set(selected)
        elif control_pressed:
            selected = set(self.table.selection())
            if row_id in selected:
                selected.remove(row_id)
            else:
                selected.add(row_id)
            self.table.selection_set(list(selected))
            self._selection_anchor_row = row_id
        else:
            self.table.selection_set(row_id)
            self._selection_anchor_row = row_id
        self.table.focus(row_id)
        self.after_idle(self._update_active_cell_border)
        return "break"

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
        self.after_idle(self._schedule_table_grid_update)
        return "break"

    def _ensure_active_cell_visible(self, row_id: str, column_id: str) -> None:
        self.table.see(row_id)
        self.update_idletasks()
        columns = list(self.table["columns"])
        if column_id == "#2":
            self.table.xview_moveto(0.0)
            self.update_idletasks()
        target_bbox = self.table.bbox(row_id, column_id)
        if target_bbox:
            x, _y, width, _height = target_bbox
            visible_width = max(self.table.winfo_width(), 1)
            if x >= 0 and (x + width) <= visible_width:
                self.after_idle(self._update_active_cell_border)
                self.after_idle(self._schedule_table_grid_update)
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
        self.after_idle(self._update_active_cell_border)
        self.after_idle(self._schedule_table_grid_update)

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
        self._apply_image_diameter_mode()
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
        self._update_field_status_hint()

    def _on_object_mode_changed(self, _event=None) -> None:
        self._sync_field_default_from_current_type()
        self._sync_field_mode_ui()
        self._sync_object_controls()
        self._mark_plot_update_pending()

    def _on_image_diameter_mode_changed(self, _event=None) -> None:
        self._apply_image_diameter_mode()
        self._sync_table()
        self._sync_object_controls()
        self._mark_plot_update_pending()

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
        if name == FOLDED_STARTER_LAYOUT_TITLE:
            self.display_orientation_var.set("Horizontal")
            self.object_mode_var.set("Finite")
            self.field_type_var.set("Object Height")
            self._last_field_type = "Object Height"
            self._field_type_defaults["Object Height"] = "0.0"
            self.field_value_var.set("0.0")
            self._sync_field_mode_ui()
        elif name == "Doublet Lens":
            self.display_orientation_var.set("Vertical")
            self.object_mode_var.set("Infinity")
            self.field_type_var.set("Angle")
            self._last_field_type = "Angle"
            self._field_type_defaults["Angle"] = "0.0"
            self.field_value_var.set("0.0")
            self._sync_field_mode_ui()
        else:
            self.display_orientation_var.set("Vertical")
            self.object_mode_var.set("Finite")
            self.field_type_var.set("Object Height")
            self._last_field_type = "Object Height"
            self._field_type_defaults["Object Height"] = "0.0"
            self.field_value_var.set("0.0")
            self._sync_field_mode_ui()

    def _on_field_type_changed(self, _event=None) -> None:
        self._sync_field_default_from_current_type()
        self._sync_field_mode_ui()
        self._sync_object_controls()
        self._mark_plot_update_pending()

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
        self._update_field_status_hint()

    def add_surface(self) -> None:
        self._begin_history_capture()
        insert_at = len(self.rows)
        if self.rows and self.rows[-1].surface == "Image":
            insert_at -= 1
        self.rows.insert(insert_at, SurfaceRow())
        self._sync_table()
        self._commit_history_capture()
        self.refresh_plot()

    def delete_selected(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        self._begin_history_capture()
        indices = sorted(self.table.index(item) for item in selected)
        for index in reversed(indices):
            del self.rows[index]
        self._sync_table()
        self._commit_history_capture()
        self.refresh_plot()

    def duplicate_selected(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        self._begin_history_capture()
        indices = sorted(self.table.index(item) for item in selected)
        insert_at = indices[-1] + 1
        duplicates = [SurfaceRow(**asdict(self.rows[index])) for index in indices]
        for offset, row in enumerate(duplicates):
            self.rows.insert(insert_at + offset, row)
        self._normalize_special_rows()
        self._sync_table()
        new_items = self.table.get_children()[insert_at:insert_at + len(duplicates)]
        self.table.selection_set(new_items)
        self._commit_history_capture()
        self.refresh_plot()

    def flip_selected(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        self._begin_history_capture()
        indices = sorted(self.table.index(item) for item in selected)
        if len(indices) < 2:
            self._history_pending_state = None
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
        self._commit_history_capture()
        self.refresh_plot()

    def move_up(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        self._begin_history_capture()
        index = min(self.table.index(item) for item in selected)
        if index == 0:
            self._history_pending_state = None
            return
        self.rows[index - 1], self.rows[index] = self.rows[index], self.rows[index - 1]
        self._sync_table()
        self.table.selection_set(self.table.get_children()[index - 1])
        self._commit_history_capture()
        self.refresh_plot()

    def move_down(self) -> None:
        selected = self.table.selection()
        if not selected:
            return
        self._begin_history_capture()
        index = max(self.table.index(item) for item in selected)
        if index >= len(self.rows) - 1:
            self._history_pending_state = None
            return
        self.rows[index + 1], self.rows[index] = self.rows[index], self.rows[index + 1]
        self._sync_table()
        self.table.selection_set(self.table.get_children()[index + 1])
        self._commit_history_capture()
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
        row_index = self.table.index(row_id)
        paraxial_target = self._paraxial_solve_target_for_cell(row_index, field)
        spec = self._variable_spec_for_field(field)
        row = self.rows[row_index]
        supports_optimization = False
        bounds = None
        if spec is not None and row.surface != "Image" and spec.is_supported(row):
            supports_optimization = True
            bounds = spec.get_bounds(row)
        if not supports_optimization and paraxial_target is None:
            return
        if self.popup_menu is not None:
            self.popup_menu.destroy()
        self.current_menu_row_id = row_id
        self.current_menu_field = field
        menu = tk.Menu(self, tearoff=0)
        if supports_optimization and spec is not None:
            marked = spec.is_enabled(row)
            menu.add_command(
                label=f"{'Unselect' if marked else 'Select'} {spec.label} for optimization",
                command=self.toggle_current_optimization_cell,
            )
            menu.add_separator()
            menu.add_command(label="Set bounds...", command=self.edit_current_bounds)
            menu.add_command(
                label="Clear bounds",
                command=self.clear_current_bounds,
                state=("normal" if bounds else "disabled"),
            )
        if paraxial_target is not None:
            if supports_optimization:
                menu.add_separator()
            menu.add_command(
                label=(
                    "Paraxial Solve Object Distance"
                    if paraxial_target == "object"
                    else "Paraxial Solve Image Distance"
                ),
                command=self.solve_current_paraxial_distance,
            )
            if paraxial_target == "object":
                menu.add_command(label="Set Object to 2F", command=self.set_current_object_to_two_f)
                menu.add_command(label="Set 2F <-> 2F", command=self.set_current_two_f_pair)
            else:
                menu.add_command(label="Set Image to 2F", command=self.set_current_image_to_two_f)
                menu.add_command(label="Set 2F <-> 2F", command=self.set_current_two_f_pair)
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
        row_index = self.table.index(row_id)
        self._begin_history_capture()
        if (
            field == "diameter"
            and hasattr(self, "image_diameter_mode_var")
            and 0 <= row_index < len(self.rows)
            and self.rows[row_index].surface == "Image"
        ):
            self.image_diameter_mode_var.set("Manual")
        self.table.set(row_id, field, value)
        self._read_rows_from_table()
        self._normalize_special_rows()
        self._sync_table()
        self._commit_history_capture()
        self._mark_plot_update_pending()

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
        self._begin_history_capture()
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
        self._commit_history_capture()
        self._mark_plot_update_pending()
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
        self._begin_history_capture()
        spec.set_enabled(row, not spec.is_enabled(row))
        self._sync_table()
        self._commit_history_capture()
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
        dialog.withdraw()
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
            self._begin_history_capture()
            spec.set_bounds(row, (lower, upper))
            self._commit_history_capture()
            self.append_progress(
                f"Bounds set for row {index} {spec.label}: [{lower:g}, {upper:g}]"
            )
            dialog.destroy()

        buttons = ttk.Frame(dialog)
        buttons.grid(row=4, column=0, padx=12, pady=(0, 12), sticky="w")
        ttk.Button(buttons, text="Save", command=accept).pack(side="left")
        ttk.Button(buttons, text="Cancel", command=dialog.destroy).pack(side="left", padx=(8, 0))

        self._show_centered_dialog(dialog)
        self.wait_window(dialog)
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None
        self.current_menu_row_id = None
        self.current_menu_field = None

    def _show_centered_dialog(self, dialog: tk.Toplevel) -> None:
        def place_dialog() -> None:
            if not dialog.winfo_exists():
                return
            dialog.update_idletasks()
            dialog_width = max(dialog.winfo_reqwidth(), dialog.winfo_width(), 1)
            dialog_height = max(dialog.winfo_reqheight(), dialog.winfo_height(), 1)
            screen_width = max(dialog.winfo_screenwidth(), 1)
            screen_height = max(dialog.winfo_screenheight(), 1)
            pos_x = max((screen_width - dialog_width) // 2, 0)
            pos_y = max((screen_height - dialog_height) // 2, 0)
            dialog.geometry(f"{dialog_width}x{dialog_height}+{pos_x}+{pos_y}")

        place_dialog()
        dialog.deiconify()
        dialog.lift()
        dialog.focus_force()
        dialog.after_idle(place_dialog)
        dialog.after(80, place_dialog)

    def clear_current_bounds(self) -> None:
        if self.current_menu_row_id is None or self.current_menu_field is None:
            return
        index = self.table.index(self.current_menu_row_id)
        row = self.rows[index]
        spec = self._variable_spec_for_field(self.current_menu_field)
        if spec is None:
            return
        self._begin_history_capture()
        spec.set_bounds(row, None)
        self._commit_history_capture()
        self.append_progress(f"Bounds cleared for row {index} {spec.label}.")
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None
        self.current_menu_row_id = None
        self.current_menu_field = None

    def _paraxial_solve_target_for_cell(self, row_index: int, field: str) -> str | None:
        if field != "thickness" or not self.rows:
            return None
        if row_index == 0:
            return "object"
        if row_index in {len(self.rows) - 2, len(self.rows) - 1}:
            return "image"
        return None

    def _cleanup_current_popup_menu(self) -> None:
        if self.popup_menu is not None:
            self.popup_menu.destroy()
            self.popup_menu = None
        self.current_menu_row_id = None
        self.current_menu_field = None

    def _center_dialog_over_main_window(self, dialog: tk.Toplevel) -> None:
        dialog.update_idletasks()
        parent_x = self.winfo_rootx()
        parent_y = self.winfo_rooty()
        parent_w = max(self.winfo_width(), 1)
        parent_h = max(self.winfo_height(), 1)
        dialog_w = max(dialog.winfo_width(), 1)
        dialog_h = max(dialog.winfo_height(), 1)
        pos_x = parent_x + max((parent_w - dialog_w) // 2, 0)
        pos_y = parent_y + max((parent_h - dialog_h) // 2, 0)
        dialog.geometry(f"+{pos_x}+{pos_y}")

    @staticmethod
    def _center_dialog_on_screen(dialog: tk.Toplevel) -> None:
        def place_dialog() -> None:
            if not dialog.winfo_exists():
                return
            dialog.update_idletasks()
            dialog_w = max(dialog.winfo_width(), dialog.winfo_reqwidth(), 1)
            dialog_h = max(dialog.winfo_height(), dialog.winfo_reqheight(), 1)
            screen_w = max(dialog.winfo_screenwidth(), 1)
            screen_h = max(dialog.winfo_screenheight(), 1)
            pos_x = max((screen_w - dialog_w) // 2, 0)
            pos_y = max((screen_h - dialog_h) // 2, 0)
            dialog.geometry(f"+{pos_x}+{pos_y}")

        place_dialog()
        dialog.after_idle(place_dialog)
        dialog.after(80, place_dialog)

    def show_formula_help(self) -> None:
        try:
            html_doc = self._build_formula_help_html()
            cache_dir = Path.home() / ".cache" / "krakenos"
            cache_dir.mkdir(parents=True, exist_ok=True)
            help_path = cache_dir / "optics_formula_sheet.html"
            help_path.write_text(html_doc, encoding="utf-8")
            self._formula_help_path = help_path
            opened = webbrowser.open_new_tab(help_path.as_uri())
            if not opened:
                fallback = self._open_document_with_system_viewer(help_path)
                if not fallback:
                    raise RuntimeError("No browser opener available.")
            self.status_var.set(f"Opened optics help: {help_path}")
        except Exception as exc:
            self.append_debug(f"Help page open failed: {exc}")
            self.status_var.set(f"Help page failed: {_short_error_message(exc)}")

    def open_paraxial_calculator(self) -> None:
        dialog = tk.Toplevel(self)
        dialog.withdraw()
        dialog.title("Paraxial Calculator")
        dialog.transient(self)
        dialog.grab_set()
        dialog.resizable(False, False)
        dialog.columnconfigure(1, weight=1)

        object_default = float(self.rows[0].thickness) if self.rows else 0.0
        image_row = max(0, len(self.rows) - 2)
        image_default = float(self.rows[image_row].thickness) if self.rows else 0.0
        object_mode_default = self._current_object_mode()

        effl_var = tk.StringVar(value=f"{self._current_effl_estimate():.6g}")
        ppa_var = tk.StringVar(value="0")
        ppp_var = tk.StringVar(value="0")
        ep_z_var = tk.StringVar(value="n/a")
        xp_z_var = tk.StringVar(value="n/a")
        magnification_var = tk.StringVar(value="0")
        solve_for_var = tk.StringVar(value="Image distance")
        object_mode_var = tk.StringVar(value=object_mode_default)
        object_distance_var = tk.StringVar(value=f"{object_default:.6g}")
        image_distance_var = tk.StringVar(value=f"{image_default:.6g}")
        load_note_var = tk.StringVar(value="Set known values, then click Solve.")
        result_var = tk.StringVar(value="Set known values, then click Solve.")
        detail_var = tk.StringVar(value="")
        solved_payload: dict[str, object] = {}

        def _format_calc(value: float) -> str:
            if not np.isfinite(value):
                return "Infinity"
            return f"{float(value):.6g}"

        def _try_load_from_layout() -> None:
            note_parts: list[str] = []
            try:
                effl, ppa, ppp = self._exact_paraxial_cardinals()
                effl_var.set(f"{float(effl):.6g}")
                ppa_var.set(f"{float(ppa):.6g}")
                ppp_var.set(f"{float(ppp):.6g}")
                note_parts.append("Loaded EFL/H1/H2 from layout.")
            except Exception as exc:
                note_parts.append(f"Cardinal extraction unavailable ({_short_error_message(exc)}).")
            try:
                system = self.build_system()
                pupil = Kos.PupilCalc(
                    system,
                    self._analysis_surface_index(),
                    self._current_wavelength(),
                    self._current_aperture_type(),
                    self._current_aperture_value(),
                )
                ep_z_var.set(_format_calc(float(pupil.PosPupInp[2])))
                xp_z_var.set(_format_calc(float(pupil.PosPupOut[2])))
                note_parts.append("Loaded EP/XP from current aperture settings.")
            except Exception:
                ep_z_var.set("n/a")
                xp_z_var.set("n/a")
            load_note_var.set(" ".join(note_parts) if note_parts else "Using manual values.")

        ttk.Label(dialog, text="Solve for").grid(row=0, column=0, padx=(12, 8), pady=(12, 4), sticky="w")
        solve_for_menu = ttk.Combobox(
            dialog,
            textvariable=solve_for_var,
            state="readonly",
            width=26,
            values=[
                "Image distance",
                "Object distance",
                "Magnification",
                "Distances from magnification",
            ],
        )
        solve_for_menu.grid(row=0, column=1, padx=(0, 12), pady=(12, 4), sticky="ew")

        ttk.Label(dialog, text="EFL / EFFL [mm]").grid(row=1, column=0, padx=(12, 8), pady=2, sticky="w")
        effl_entry = ttk.Entry(dialog, textvariable=effl_var, width=22)
        effl_entry.grid(row=1, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="H1 offset PPA [mm]").grid(row=2, column=0, padx=(12, 8), pady=2, sticky="w")
        ppa_entry = ttk.Entry(dialog, textvariable=ppa_var, width=22)
        ppa_entry.grid(row=2, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="H2 offset PPP [mm]").grid(row=3, column=0, padx=(12, 8), pady=2, sticky="w")
        ppp_entry = ttk.Entry(dialog, textvariable=ppp_var, width=22)
        ppp_entry.grid(row=3, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="Object mode").grid(row=4, column=0, padx=(12, 8), pady=(8, 2), sticky="w")
        object_mode_menu = ttk.Combobox(
            dialog,
            textvariable=object_mode_var,
            state="readonly",
            width=22,
            values=["Finite", "Infinity"],
        )
        object_mode_menu.grid(row=4, column=1, padx=(0, 12), pady=(8, 2), sticky="ew")

        ttk.Label(dialog, text="Object distance [mm]").grid(row=5, column=0, padx=(12, 8), pady=2, sticky="w")
        object_distance_entry = ttk.Entry(dialog, textvariable=object_distance_var, width=22)
        object_distance_entry.grid(row=5, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="Image distance [mm]").grid(row=6, column=0, padx=(12, 8), pady=2, sticky="w")
        image_distance_entry = ttk.Entry(dialog, textvariable=image_distance_var, width=22)
        image_distance_entry.grid(row=6, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="Magnification m").grid(row=7, column=0, padx=(12, 8), pady=2, sticky="w")
        magnification_entry = ttk.Entry(dialog, textvariable=magnification_var, width=22)
        magnification_entry.grid(row=7, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="EP z [mm]").grid(row=8, column=0, padx=(12, 8), pady=2, sticky="w")
        ep_z_entry = ttk.Entry(dialog, textvariable=ep_z_var, width=22, state="readonly")
        ep_z_entry.grid(row=8, column=1, padx=(0, 12), pady=2, sticky="ew")

        ttk.Label(dialog, text="XP z [mm]").grid(row=9, column=0, padx=(12, 8), pady=2, sticky="w")
        xp_z_entry = ttk.Entry(dialog, textvariable=xp_z_var, width=22, state="readonly")
        xp_z_entry.grid(row=9, column=1, padx=(0, 12), pady=2, sticky="ew")

        note_label = ttk.Label(dialog, textvariable=load_note_var, foreground="#475569", wraplength=500, justify="left")
        note_label.grid(row=10, column=0, columnspan=2, padx=12, pady=(8, 2), sticky="w")

        ttk.Label(dialog, textvariable=result_var, font=("TkDefaultFont", 10, "bold")).grid(
            row=11, column=0, columnspan=2, padx=12, pady=(4, 0), sticky="w"
        )
        ttk.Label(dialog, textvariable=detail_var, foreground="#475569", wraplength=500, justify="left").grid(
            row=12, column=0, columnspan=2, padx=12, pady=(2, 0), sticky="w"
        )

        def _read_float(var: tk.StringVar, label: str) -> float:
            text = var.get().strip()
            if not text:
                raise RuntimeError(f"{label} is required")
            try:
                value = float(text)
            except ValueError as exc:
                raise RuntimeError(f"{label} must be numeric") from exc
            if not np.isfinite(value):
                raise RuntimeError(f"{label} must be finite")
            return float(value)

        def _refresh_mode_state(_event=None) -> None:
            target = solve_for_var.get().strip()
            mode = object_mode_var.get().strip()
            if target == "Image distance":
                if mode == "Infinity":
                    object_distance_entry.configure(state="disabled")
                else:
                    object_distance_entry.configure(state="normal")
                image_distance_entry.configure(state="disabled")
                magnification_entry.configure(state="readonly")
            elif target == "Object distance":
                object_distance_entry.configure(state="disabled")
                image_distance_entry.configure(state="normal")
                magnification_entry.configure(state="readonly")
            elif target == "Magnification":
                if mode == "Infinity":
                    object_distance_entry.configure(state="disabled")
                else:
                    object_distance_entry.configure(state="normal")
                image_distance_entry.configure(state="normal")
                magnification_entry.configure(state="disabled")
            else:
                object_distance_entry.configure(state="disabled")
                image_distance_entry.configure(state="disabled")
                magnification_entry.configure(state="normal")
            solved_payload.clear()

        def _solve(_event=None) -> None:
            try:
                f = _read_float(effl_var, "EFL")
                if abs(f) <= 1e-12:
                    raise RuntimeError("EFL must be non-zero")
                h1 = _read_float(ppa_var, "H1 offset")
                h2 = _read_float(ppp_var, "H2 offset")
                target = solve_for_var.get().strip()
                mode = object_mode_var.get().strip()
                solved_payload.clear()

                if target == "Image distance":
                    if mode == "Infinity":
                        image_distance = f + h2
                        object_principal = float("inf")
                        image_principal = float(f)
                        magnification = 0.0
                    else:
                        object_distance = _read_float(object_distance_var, "Object distance")
                        object_principal = object_distance + h1
                        if abs(object_principal) <= 1e-12:
                            raise RuntimeError("Object is on H1; cannot solve image distance")
                        balance = (1.0 / f) - (1.0 / object_principal)
                        if abs(balance) <= 1e-12:
                            image_distance = float("inf")
                            image_principal = float("inf")
                            magnification = float("inf")
                        else:
                            image_principal = 1.0 / balance
                            image_distance = image_principal + h2
                            magnification = image_principal / object_principal
                    solved_payload.update(
                        {
                            "target": "image",
                            "value": image_distance,
                            "object_mode_after": mode,
                        }
                    )
                    magnification_var.set(_format_calc(magnification))
                    result_var.set(f"Image distance = {self._format_paraxial_value(image_distance)} mm")
                    detail_var.set(
                        "s={obj}, s'={img}, m={mag}".format(
                            obj=self._format_paraxial_value(object_principal),
                            img=self._format_paraxial_value(image_principal),
                            mag=self._format_paraxial_value(magnification),
                        )
                    )
                elif target == "Object distance":
                    image_distance = _read_float(image_distance_var, "Image distance")
                    image_principal = image_distance - h2
                    if abs(image_principal) <= 1e-12:
                        raise RuntimeError("Image is on H2; cannot solve object distance")
                    balance = (1.0 / f) - (1.0 / image_principal)
                    if abs(balance) <= 1e-12:
                        object_principal = float("inf")
                        object_distance = float("inf")
                        mode_after = "Infinity"
                    else:
                        object_principal = 1.0 / balance
                        object_distance = object_principal - h1
                        mode_after = "Infinity" if (not np.isfinite(object_distance) or abs(object_distance) > 1e9) else "Finite"
                    magnification = image_principal / object_principal if np.isfinite(object_principal) and abs(object_principal) > 1e-12 else float("inf")
                    solved_payload.update(
                        {
                            "target": "object",
                            "value": object_distance,
                            "object_mode_after": mode_after,
                        }
                    )
                    magnification_var.set(_format_calc(magnification))
                    result_var.set(f"Object distance = {self._format_paraxial_value(object_distance)} mm")
                    detail_var.set(
                        "s={obj}, s'={img}, m={mag}".format(
                            obj=self._format_paraxial_value(object_principal),
                            img=self._format_paraxial_value(image_principal),
                            mag=self._format_paraxial_value(magnification),
                        )
                    )
                elif target == "Magnification":
                    if mode == "Infinity":
                        object_principal = float("inf")
                        image_principal = _read_float(image_distance_var, "Image distance") - h2
                        magnification = 0.0
                    else:
                        object_distance = _read_float(object_distance_var, "Object distance")
                        image_distance = _read_float(image_distance_var, "Image distance")
                        object_principal = object_distance + h1
                        image_principal = image_distance - h2
                        if abs(object_principal) <= 1e-12:
                            raise RuntimeError("Object is on H1; cannot solve magnification")
                        magnification = image_principal / object_principal
                    solved_payload.update({"target": "magnification", "value": magnification, "object_mode_after": mode})
                    magnification_var.set(_format_calc(magnification))
                    result_var.set(f"Magnification = {self._format_paraxial_value(magnification)}")
                    detail_var.set(
                        "s={obj}, s'={img} from H1/H2".format(
                            obj=self._format_paraxial_value(object_principal),
                            img=self._format_paraxial_value(image_principal),
                        )
                    )
                else:
                    magnification = _read_float(magnification_var, "Magnification")
                    if abs(magnification) <= 1e-12:
                        raise RuntimeError("Magnification too close to zero; object distance goes to infinity")
                    if abs(1.0 + magnification) <= 1e-12:
                        raise RuntimeError("Magnification of -1 makes object/image distance singular")
                    object_principal = f * (1.0 + (1.0 / magnification))
                    image_principal = f * (1.0 + magnification)
                    object_distance = object_principal - h1
                    image_distance = image_principal + h2
                    mode_after = "Infinity" if (not np.isfinite(object_distance) or abs(object_distance) > 1e9) else "Finite"
                    object_distance_var.set(_format_calc(object_distance))
                    image_distance_var.set(_format_calc(image_distance))
                    solved_payload.update(
                        {
                            "target": "pair",
                            "object_value": object_distance,
                            "image_value": image_distance,
                            "object_mode_after": mode_after,
                        }
                    )
                    result_var.set(
                        f"Object={self._format_paraxial_value(object_distance)} mm, Image={self._format_paraxial_value(image_distance)} mm"
                    )
                    detail_var.set(
                        "From m={mag}: s={obj}, s'={img}".format(
                            mag=self._format_paraxial_value(magnification),
                            obj=self._format_paraxial_value(object_principal),
                            img=self._format_paraxial_value(image_principal),
                        )
                    )
            except Exception as exc:
                solved_payload.clear()
                result_var.set(f"Solve failed: {_short_error_message(exc)}")
                detail_var.set("")

        def _apply_to_layout() -> bool:
            try:
                if not solved_payload:
                    _solve()
                    if not solved_payload:
                        return False
                target = str(solved_payload.get("target", ""))
                solved_value = float(solved_payload.get("value", 0.0))
                mode_after = str(solved_payload.get("object_mode_after", self._current_object_mode()))

                if target == "image":
                    if not np.isfinite(solved_value):
                        raise RuntimeError("Solved image distance is infinity and cannot be applied")
                    row_index = max(0, len(self.rows) - 2)
                    self.rows[row_index].thickness = solved_value
                    self._select_table_row(row_index)
                elif target == "object":
                    self.object_mode_var.set(mode_after)
                    if mode_after == "Finite":
                        if not np.isfinite(solved_value):
                            raise RuntimeError("Solved object distance is infinity and cannot be applied in Finite mode")
                        self.rows[0].thickness = solved_value
                    self._select_table_row(0)
                elif target == "pair":
                    object_value = float(solved_payload.get("object_value", float("nan")))
                    image_value = float(solved_payload.get("image_value", float("nan")))
                    self.object_mode_var.set(mode_after)
                    if mode_after == "Finite":
                        if not np.isfinite(object_value):
                            raise RuntimeError("Solved object distance is infinity and cannot be applied in Finite mode")
                        self.rows[0].thickness = object_value
                    if not np.isfinite(image_value):
                        raise RuntimeError("Solved image distance is infinity and cannot be applied")
                    row_index = max(0, len(self.rows) - 2)
                    self.rows[row_index].thickness = image_value
                    self._select_table_row(row_index)
                elif target == "magnification":
                    self.status_var.set("Magnification computed. No layout cell to apply.")
                    return False
                else:
                    raise RuntimeError("No solved target to apply")

                self._normalize_special_rows()
                self._sync_table()
                self._sync_object_controls()
                self._mark_plot_update_pending()
                self.append_progress(f"Paraxial calculator applied: {result_var.get()}")
                self.status_var.set(f"{result_var.get()}  |  Click Update.")
                return True
            except Exception as exc:
                message = _short_error_message(exc)
                self.append_debug(f"Paraxial calculator apply failed: {exc}")
                messagebox.showerror("Paraxial Calculator", message)
                self.status_var.set(f"Paraxial calculator apply failed: {message}")
                return False

        def _apply_and_close() -> None:
            if _apply_to_layout():
                dialog.destroy()

        buttons = ttk.Frame(dialog)
        buttons.grid(row=13, column=0, columnspan=2, padx=12, pady=(10, 12), sticky="e")
        ttk.Button(buttons, text="Use Current Layout", command=lambda: (
            object_mode_var.set(self._current_object_mode()),
            object_distance_var.set(f"{(float(self.rows[0].thickness) if self.rows else 0.0):.6g}"),
            image_distance_var.set(f"{(float(self.rows[max(0, len(self.rows) - 2)].thickness) if self.rows else 0.0):.6g}"),
            _try_load_from_layout(),
            _refresh_mode_state(),
            _solve(),
        )).pack(side="left")
        ttk.Button(buttons, text="Solve", command=_solve).pack(side="left", padx=(8, 0))
        ttk.Button(buttons, text="Apply to Layout", command=_apply_and_close).pack(side="left", padx=(8, 0))
        ttk.Button(buttons, text="Close", command=dialog.destroy).pack(side="left", padx=(8, 0))

        solve_for_menu.bind("<<ComboboxSelected>>", lambda _e: (_refresh_mode_state(), _solve()))
        object_mode_menu.bind("<<ComboboxSelected>>", lambda _e: (_refresh_mode_state(), _solve()))
        for entry in (effl_entry, ppa_entry, ppp_entry, object_distance_entry, image_distance_entry, magnification_entry):
            entry.bind("<Return>", _solve)
            entry.bind("<Tab>", _solve, add="+")
            entry.bind("<FocusOut>", _solve)

        _try_load_from_layout()
        _refresh_mode_state()
        _solve()
        self._center_dialog_on_screen(dialog)
        dialog.deiconify()
        dialog.lift()
        dialog.focus_force()

    @staticmethod
    def _open_document_with_system_viewer(document_path: Path) -> bool:
        try:
            if os.name == "nt":
                os.startfile(str(document_path))  # type: ignore[attr-defined]
                return True
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(document_path)])
                return True
            if shutil.which("xdg-open"):
                subprocess.Popen(["xdg-open", str(document_path)])
                return True
            if shutil.which("gio"):
                subprocess.Popen(["gio", "open", str(document_path)])
                return True
        except Exception:
            return False
        return False

    def _build_formula_help_html(self) -> str:
        object_mode = html.escape(self._current_object_mode())
        wavelength = self._current_wavelength()
        object_gap = float(self.rows[0].thickness) if self.rows else float("nan")
        image_gap = self._current_image_distance()
        object_size = float(self.rows[0].diameter) if self.rows else float("nan")
        sensor_size = float(self.rows[-1].diameter) if self.rows else float("nan")
        field_type = html.escape(self._current_field_type())
        field_value = self._current_field_value()
        effl_text = "Unavailable"
        ppa_text = "Unavailable"
        ppp_text = "Unavailable"
        image_size_text = "Unavailable"
        fill_text = "Unavailable"
        try:
            effl, ppa, ppp = self._exact_paraxial_cardinals(wavelength)
            effl_text = f"{effl:.6g} mm"
            ppa_text = f"{ppa:.6g} mm"
            ppp_text = f"{ppp:.6g} mm"
            if self._current_object_mode() == "Finite" and object_gap > 0:
                s = object_gap + ppa
                if abs(s) > 1e-12:
                    denom = (1.0 / effl) - (1.0 / s)
                    if abs(denom) > 1e-12:
                        sp = 1.0 / denom
                        magnification = sp / s
                        image_size = abs(magnification) * max(object_size, 0.0)
                        image_size_text = f"{image_size:.6g} mm"
                        if sensor_size > 1e-12:
                            fill_text = f"{100.0 * image_size / sensor_size:.4g}%"
        except Exception:
            pass

        return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>KrakenOS Optics Formula Sheet</title>
  <style>
    :root {{
      --bg: #f6f8fc;
      --panel: #ffffff;
      --ink: #1f2937;
      --muted: #4b5563;
      --accent: #0f766e;
      --line: #d1d5db;
    }}
    body {{
      margin: 0;
      background: linear-gradient(180deg, #eef2ff 0%, var(--bg) 40%);
      color: var(--ink);
      font-family: \"Iosevka Aile\", \"Source Sans 3\", \"Noto Sans\", sans-serif;
      line-height: 1.55;
    }}
    .wrap {{
      max-width: 980px;
      margin: 0 auto;
      padding: 28px 20px 40px;
    }}
    .hero {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 14px;
      padding: 18px 20px;
      box-shadow: 0 8px 22px rgba(15, 23, 42, 0.08);
      margin-bottom: 16px;
    }}
    h1 {{
      margin: 0 0 8px 0;
      font-size: 1.42rem;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 8px 16px;
      font-size: 0.96rem;
    }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-left: 4px solid var(--accent);
      border-radius: 12px;
      padding: 14px 16px;
      margin: 12px 0;
      box-shadow: 0 6px 16px rgba(15, 23, 42, 0.05);
    }}
    h2 {{
      margin: 0 0 8px 0;
      font-size: 1.08rem;
    }}
    .note {{
      color: var(--muted);
      font-size: 0.94rem;
    }}
    ul {{
      margin-top: 8px;
      padding-left: 18px;
    }}
    code {{
      background: #eef2ff;
      border: 1px solid #c7d2fe;
      border-radius: 5px;
      padding: 0 5px;
    }}
  </style>
  <script>
    window.MathJax = {{
      tex: {{ inlineMath: [['\\\\(','\\\\)'], ['$', '$']], displayMath: [['\\\\[','\\\\]']] }},
      svg: {{ fontCache: 'global' }}
    }};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
</head>
<body>
  <div class="wrap">
    <section class="hero">
      <h1>KrakenOS Formula Sheet</h1>
      <p class="note">This page is generated from your current UI state. It uses the same centered paraxial model as the <code>Paraxial Solve</code> tool.</p>
      <div class="grid">
        <div><strong>Object mode:</strong> {object_mode}</div>
        <div><strong>Wavelength:</strong> {wavelength:.6g} um</div>
        <div><strong>Object gap:</strong> {object_gap:.6g} mm</div>
        <div><strong>Image gap:</strong> {image_gap:.6g} mm</div>
        <div><strong>Object size:</strong> {object_size:.6g} mm</div>
        <div><strong>Sensor size:</strong> {sensor_size:.6g} mm</div>
        <div><strong>Field:</strong> {field_type} = {field_value:.6g}</div>
        <div><strong>EFFL / PPA / PPP:</strong> {html.escape(effl_text)} / {html.escape(ppa_text)} / {html.escape(ppp_text)}</div>
        <div><strong>Paraxial image size:</strong> {html.escape(image_size_text)}</div>
        <div><strong>Sensor fill:</strong> {html.escape(fill_text)}</div>
      </div>
    </section>

    <section class="card">
      <h2>Paraxial imaging</h2>
      <p>\\[\\frac{{1}}{{f}} = \\frac{{1}}{{s}} + \\frac{{1}}{{s'}}\\]</p>
      <p>\\[m = \\frac{{s'}}{{s}} = \\frac{{y'}}{{y}}\\]</p>
      <p class="note">Solve in principal-plane space, then map to table thickness values.</p>
    </section>

    <section class="card">
      <h2>UI thickness conversion</h2>
      <p>\\[s = z_{{\\mathrm{{obj}}}} + \\mathrm{{PPA}}\\]</p>
      <p>\\[z_{{\\mathrm{{img}}}} = s' + \\mathrm{{PPP}}\\]</p>
      <p class="note"><code>Object Thickness</code> is \\(z_\\mathrm{{obj}}\\). Image solve writes to the last optical row thickness before <code>Image</code>.</p>
    </section>

    <section class="card">
      <h2>2F rule for thick lenses</h2>
      <p>\\[s = 2f \\Rightarrow s' = 2f\\]</p>
      <p>\\[z_{{\\mathrm{{obj}},2F}} = 2f - \\mathrm{{PPA}},\\qquad z_{{\\mathrm{{img}},2F}} = 2f + \\mathrm{{PPP}}\\]</p>
      <p class="note">The symmetry is around principal planes H1/H2, not lens vertices.</p>
    </section>

    <section class="card">
      <h2>Image size and sensor fill</h2>
      <p>\\[y' = \\left|\\frac{{s'}}{{s}}\\right|\\,y\\]</p>
      <p>\\[\\mathrm{{fill}} = \\frac{{y'}}{{y_{{\\mathrm{{sensor}}}}}}\\]</p>
      <p class="note">Changing <code>Image Diameter</code> does not change focus distance. It changes framing/fill.</p>
    </section>

    <section class="card">
      <h2>Aperture quick rule</h2>
      <p>\\[N = \\frac{{f}}{{D_{{\\mathrm{{EP}}}}}},\\qquad D_{{\\mathrm{{EP}}}} \\approx \\frac{{f}}{{N}}\\]</p>
      <p class="note">Keep <code>STOP</code> and <code>EPD</code> choices consistent between analysis and optimization.</p>
    </section>

    <section class="card">
      <h2>Practical UI reminders</h2>
      <ul>
        <li><code>Object Diameter</code> and <code>Image Diameter</code> are full sizes.</li>
        <li><code>Field type = Object Height</code> uses semi-field values.</li>
        <li>Paraxial solve is intended for centered refractive layouts. Mirror/tilt/decenter cases still need full trace validation.</li>
      </ul>
    </section>
  </div>
</body>
</html>
"""

    def _exact_paraxial_cardinals(self, wavelength: float | None = None) -> tuple[float, float, float]:
        if len(self.rows) < 3:
            raise RuntimeError("Not enough surfaces for paraxial solve")
        solve_rows = [SurfaceRow(**asdict(row)) for row in self.rows]
        optical_rows = solve_rows[1:-1]
        if not optical_rows:
            raise RuntimeError("No optical block available for paraxial solve")
        unsupported: list[str] = []
        for row in optical_rows:
            if row.surface not in {"Standard", "Thin Lens"}:
                unsupported.append(row.name or row.surface)
                continue
            if any(
                abs(value) > 1e-9
                for value in (
                    row.tilt_x,
                    row.tilt_y,
                    row.tilt_z,
                    row.desp_x,
                    row.desp_y,
                    row.desp_z,
                    row.axis_move,
                )
            ):
                unsupported.append(row.name or row.surface)
        if unsupported:
            raise RuntimeError("Paraxial solve supports centered refractive systems only")
        solve_rows[0].thickness = 0.0
        solve_rows[-2].thickness = 0.0
        solve_rows[-1].thickness = 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            solve_system = _build_system_from_specs([asdict(row) for row in solve_rows])
            _, _, _, _a, _b, _c, _d, effl, ppa, ppp, *_rest = solve_system.Parax(
                self._current_wavelength() if wavelength is None else float(wavelength)
            )
        return float(effl), float(ppa), float(ppp)

    def _paraxial_two_f_gaps(self) -> tuple[float, float, float, float, float]:
        effl, ppa, ppp = self._exact_paraxial_cardinals()
        if not np.isfinite(effl) or abs(effl) <= 1e-12:
            raise RuntimeError("Paraxial solve failed: invalid effective focal length")
        two_f = 2.0 * effl
        object_gap = two_f - ppa
        image_gap = two_f + ppp
        return float(effl), float(ppa), float(ppp), float(object_gap), float(image_gap)

    @staticmethod
    def _format_paraxial_value(value: float | str) -> str:
        if isinstance(value, str):
            return value
        if not np.isfinite(value):
            return "Infinity"
        return f"{float(value):.6g}"

    def _compute_paraxial_solve_result(self, target: str) -> dict[str, float | str]:
        effl, ppa, ppp = self._exact_paraxial_cardinals()
        if not np.isfinite(effl) or abs(effl) <= 1e-12:
            raise RuntimeError("Paraxial solve failed: invalid effective focal length")
        result: dict[str, float | str] = {
            "target": target,
            "effl": float(effl),
            "ppa": float(ppa),
            "ppp": float(ppp),
            "object_mode_before": self._current_object_mode(),
            "object_distance_before": float(self.rows[0].thickness) if self.rows else 0.0,
            "image_distance_before": float(self._current_image_distance()),
        }
        if target == "image":
            if self._current_object_mode() == "Infinity":
                solved_distance = effl + ppp
                result["object_principal"] = "Infinity"
                result["image_principal"] = float(effl)
            else:
                object_distance = float(self.rows[0].thickness)
                object_principal = object_distance + ppa
                if abs(object_principal) <= 1e-12:
                    raise RuntimeError("Object is located on the front principal plane")
                power_balance = (1.0 / effl) - (1.0 / object_principal)
                if abs(power_balance) <= 1e-12:
                    raise RuntimeError("Paraxial image is at infinity for the current object distance")
                image_principal = 1.0 / power_balance
                solved_distance = image_principal + ppp
                result["object_principal"] = float(object_principal)
                result["image_principal"] = float(image_principal)
            result["solved_distance"] = float(solved_distance)
            result["message"] = f"Paraxial solve: image distance -> {float(solved_distance):.6g} mm | EFFL={effl:.6g} mm"
            result["selected_row"] = max(0, len(self.rows) - 2)
            result["object_mode_after"] = result["object_mode_before"]
            return result

        image_distance = self._current_image_distance()
        image_principal = image_distance - ppp
        if abs(image_principal) <= 1e-12:
            raise RuntimeError("Image is located on the back principal plane")
        power_balance = (1.0 / effl) - (1.0 / image_principal)
        if abs(power_balance) <= 1e-12:
            solved_distance = float("inf")
            object_principal = float("inf")
            object_mode_after = "Infinity"
            message = f"Paraxial solve: object at infinity | EFFL={effl:.6g} mm"
        else:
            object_principal = 1.0 / power_balance
            solved_distance = object_principal - ppa
            object_mode_after = "Infinity" if (not np.isfinite(solved_distance) or abs(solved_distance) > 1e9) else "Finite"
            if object_mode_after == "Infinity":
                message = f"Paraxial solve: object at infinity | EFFL={effl:.6g} mm"
            else:
                message = f"Paraxial solve: object distance -> {float(solved_distance):.6g} mm | EFFL={effl:.6g} mm"
        result["image_principal"] = float(image_principal)
        result["object_principal"] = object_principal
        result["solved_distance"] = solved_distance
        result["message"] = message
        result["selected_row"] = 0
        result["object_mode_after"] = object_mode_after
        return result

    def _show_paraxial_solve_dialog(self, result: dict[str, float | str]) -> bool:
        dialog = tk.Toplevel(self)
        dialog.title("Paraxial Solve")
        dialog.transient(self)
        dialog.grab_set()
        dialog.resizable(False, False)

        target = str(result["target"])
        intro = (
            "Review the paraxial solve before applying it."
            if target == "image"
            else "Review the paraxial object-distance solve before applying it."
        )
        ttk.Label(dialog, text=intro, padding=(12, 12, 12, 4)).grid(row=0, column=0, columnspan=2, sticky="w")

        rows = [
            ("EFFL [mm]", self._format_paraxial_value(result["effl"])),
            ("Front principal plane PPA [mm]", self._format_paraxial_value(result["ppa"])),
            ("Back principal plane PPP [mm]", self._format_paraxial_value(result["ppp"])),
            ("Object mode", str(result["object_mode_before"])),
            ("Object distance before [mm]", self._format_paraxial_value(result["object_distance_before"])),
            ("Image distance before [mm]", self._format_paraxial_value(result["image_distance_before"])),
            ("Object distance from H1 [mm]", self._format_paraxial_value(result["object_principal"])),
            ("Image distance from H2 [mm]", self._format_paraxial_value(result["image_principal"])),
        ]
        if target == "image":
            rows.append(("Solved image gap [mm]", self._format_paraxial_value(result["solved_distance"])))
            rows.append(("Apply to row", str(int(result["selected_row"]))))
        else:
            rows.append(("Solved object gap [mm]", self._format_paraxial_value(result["solved_distance"])))
            rows.append(("Object mode after", str(result["object_mode_after"])))

        for row_index, (label, value) in enumerate(rows, start=1):
            ttk.Label(dialog, text=label).grid(row=row_index, column=0, padx=(12, 12), pady=2, sticky="w")
            ttk.Label(dialog, text=value, font=("TkDefaultFont", 10, "bold")).grid(
                row=row_index,
                column=1,
                padx=(0, 12),
                pady=2,
                sticky="e",
            )

        formula = "Thin-lens with principal planes: 1/f = 1/s + 1/s'"
        ttk.Label(dialog, text=formula, foreground="#4b5563", padding=(12, 8, 12, 4)).grid(
            row=len(rows) + 1,
            column=0,
            columnspan=2,
            sticky="w",
        )

        decision = {"apply": False}

        def accept() -> None:
            decision["apply"] = True
            dialog.destroy()

        def cancel() -> None:
            dialog.destroy()

        dialog.protocol("WM_DELETE_WINDOW", cancel)
        buttons = ttk.Frame(dialog, padding=(12, 4, 12, 12))
        buttons.grid(row=len(rows) + 2, column=0, columnspan=2, sticky="e")
        ttk.Button(buttons, text="Apply", command=accept).pack(side="left")
        ttk.Button(buttons, text="Cancel", command=cancel).pack(side="left", padx=(8, 0))
        self._center_dialog_over_main_window(dialog)
        self.wait_window(dialog)
        return bool(decision["apply"])

    def _apply_paraxial_two_f(self, mode: str) -> None:
        effl, _ppa, _ppp, object_gap, image_gap = self._paraxial_two_f_gaps()
        self._begin_history_capture()
        if mode in {"object", "pair"}:
            self.object_mode_var.set("Finite")
            self.rows[0].thickness = object_gap
        if mode in {"image", "pair"}:
            self.rows[max(0, len(self.rows) - 2)].thickness = image_gap
        self._normalize_special_rows()
        self._sync_table()
        selected_row = 0 if mode == "object" else max(0, len(self.rows) - 2)
        self._select_table_row(selected_row)
        self._commit_history_capture()
        self.refresh_plot()
        if mode == "pair":
            message = (
                f"Paraxial 2F setup applied | EFFL={effl:.6g} mm | "
                f"object={object_gap:.6g} mm | image={image_gap:.6g} mm"
            )
        elif mode == "object":
            message = f"Paraxial 2F object gap -> {object_gap:.6g} mm | EFFL={effl:.6g} mm"
        else:
            message = f"Paraxial 2F image gap -> {image_gap:.6g} mm | EFFL={effl:.6g} mm"
        self.status_var.set(message)
        self.append_progress(message)

    def set_current_object_to_two_f(self) -> None:
        try:
            self._apply_paraxial_two_f("object")
        except Exception as exc:
            error = _short_error_message(exc)
            messagebox.showerror("Paraxial 2F", error)
            self.append_debug(f"Paraxial 2F object failed: {exc}")
            self.status_var.set(f"Paraxial 2F object failed: {error}")
        finally:
            self._cleanup_current_popup_menu()

    def set_current_image_to_two_f(self) -> None:
        try:
            self._apply_paraxial_two_f("image")
        except Exception as exc:
            error = _short_error_message(exc)
            messagebox.showerror("Paraxial 2F", error)
            self.append_debug(f"Paraxial 2F image failed: {exc}")
            self.status_var.set(f"Paraxial 2F image failed: {error}")
        finally:
            self._cleanup_current_popup_menu()

    def set_current_two_f_pair(self) -> None:
        try:
            self._apply_paraxial_two_f("pair")
        except Exception as exc:
            error = _short_error_message(exc)
            messagebox.showerror("Paraxial 2F", error)
            self.append_debug(f"Paraxial 2F pair failed: {exc}")
            self.status_var.set(f"Paraxial 2F pair failed: {error}")
        finally:
            self._cleanup_current_popup_menu()

    def solve_current_paraxial_distance(self) -> None:
        if self.current_menu_row_id is None or self.current_menu_field is None:
            return
        row_index = self.table.index(self.current_menu_row_id)
        target = self._paraxial_solve_target_for_cell(row_index, self.current_menu_field)
        try:
            if target is None:
                raise RuntimeError("Paraxial solve is only available for object/image distance cells")
            result = self._compute_paraxial_solve_result(target)
            if not self._show_paraxial_solve_dialog(result):
                self.status_var.set("Paraxial solve cancelled")
                return
            selected_row = int(result["selected_row"])
            solved_distance = result["solved_distance"]
            self._begin_history_capture()
            if target == "image":
                self.rows[selected_row].thickness = float(solved_distance)
            else:
                object_mode_after = str(result["object_mode_after"])
                self.object_mode_var.set(object_mode_after)
                if object_mode_after == "Finite":
                    self.rows[0].thickness = float(solved_distance)
            message = str(result["message"])
            self._normalize_special_rows()
            self._sync_table()
            self._select_table_row(selected_row)
            self._commit_history_capture()
            self.refresh_plot()
            self.status_var.set(message)
            self.append_progress(message)
        except Exception as exc:
            error = _short_error_message(exc)
            messagebox.showerror("Paraxial Solve", error)
            self.append_debug(f"Paraxial solve failed: {exc}")
            self.status_var.set(f"Paraxial solve failed: {error}")
        finally:
            self._cleanup_current_popup_menu()

    def clear_optimization_marks(self) -> None:
        for row in self.rows:
            row.optimize_rc = False
            row.optimize_thickness = False
        self._sync_table()

    def benchmark_psf_mtf(self) -> None:
        self.append_progress("Benchmark PSF/MTF started.")
        try:
            self._read_rows_from_table()
            system = self.build_system()
            wavelength = self._current_wavelength()
            field_type = "angle" if self._current_object_mode() == "Infinity" else "height"
            field_y = self._current_field_angle_deg() if field_type == "angle" else self._current_field_height()
            sample_count = max(64, self._current_ray_count() * 12)
            self.append_progress(f"Tracing benchmark rays: sample_count={sample_count}")
            x_local, y_local, workers = self._build_geometric_image_samples(
                system,
                wavelength,
                sample_count=sample_count,
                pattern="hexapolar",
                surface_index=self._analysis_surface_index(),
                aperture_type=self._current_aperture_type(),
                aperture_value=self._current_aperture_value(),
                field_type=field_type,
                field_x=0.0,
                field_y=field_y,
            )
            if x_local.size < 4:
                raise RuntimeError("Not enough traced image-plane samples for benchmark")

            span_x = max(float(np.ptp(x_local)), 1e-3)
            span_y = max(float(np.ptp(y_local)), 1e-3)
            span = max(span_x, span_y) * 1.25
            bins = 256

            t0 = time.perf_counter()
            hist_cpu, xedges_cpu, _yedges_cpu = np.histogram2d(
                x_local,
                y_local,
                bins=bins,
                range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
            )
            psf_cpu = hist_cpu / max(np.sum(hist_cpu), 1.0)
            otf_cpu = np.fft.fftshift(np.fft.fft2(psf_cpu))
            mtf_cpu = np.abs(otf_cpu)
            mtf_cpu /= max(float(np.max(mtf_cpu)), 1e-12)
            _freq_cpu = np.fft.fftshift(np.fft.fftfreq(bins, d=float(xedges_cpu[1] - xedges_cpu[0])))
            cpu_sec = time.perf_counter() - t0

            gpu_results: list[tuple[str, float]] = []

            cp = _optional_cupy()
            if cp is not None:
                try:
                    if int(cp.cuda.runtime.getDeviceCount()) > 0:
                        _ = cp.zeros((1,), dtype=cp.float32)
                        cp.cuda.Stream.null.synchronize()
                        t1 = time.perf_counter()
                        x_gpu = cp.asarray(x_local, dtype=cp.float64)
                        y_gpu = cp.asarray(y_local, dtype=cp.float64)
                        hist_gpu, xedges_gpu, _yedges_gpu = cp.histogram2d(
                            x_gpu,
                            y_gpu,
                            bins=bins,
                            range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
                        )
                        psf_gpu = hist_gpu / cp.maximum(cp.sum(hist_gpu), 1.0)
                        otf_gpu = cp.fft.fftshift(cp.fft.fft2(psf_gpu))
                        mtf_gpu = cp.abs(otf_gpu)
                        mtf_gpu /= cp.maximum(cp.max(mtf_gpu), 1e-12)
                        _freq_gpu = cp.fft.fftshift(
                            cp.fft.fftfreq(bins, d=float(cp.asnumpy(xedges_gpu[1] - xedges_gpu[0])))
                        )
                        cp.cuda.Stream.null.synchronize()
                        gpu_results.append(("CuPy", time.perf_counter() - t1))
                except Exception as exc:
                    self.append_debug(f"Benchmark CuPy path failed: {_short_error_message(exc)}")

            torch = _optional_torch()
            if torch is not None:
                try:
                    if bool(torch.cuda.is_available()):
                        device = torch.device("cuda")
                        _ = torch.zeros((1,), dtype=torch.float32, device=device)
                        if hasattr(torch.cuda, "synchronize"):
                            torch.cuda.synchronize()
                        t2 = time.perf_counter()
                        lower = -span / 2.0
                        upper = span / 2.0
                        step = (upper - lower) / float(bins)
                        x_t = torch.as_tensor(x_local, dtype=torch.float64, device=device)
                        y_t = torch.as_tensor(y_local, dtype=torch.float64, device=device)
                        ix = torch.floor((x_t - lower) / step).to(torch.int64)
                        iy = torch.floor((y_t - lower) / step).to(torch.int64)
                        valid = (ix >= 0) & (ix < bins) & (iy >= 0) & (iy < bins)
                        ix = ix[valid]
                        iy = iy[valid]
                        lin = ix * bins + iy
                        hist_t = torch.zeros(bins * bins, dtype=torch.float64, device=device)
                        hist_t.scatter_add_(0, lin, torch.ones_like(lin, dtype=torch.float64))
                        hist_t = hist_t.view(bins, bins)
                        psf_t = hist_t / torch.clamp(torch.sum(hist_t), min=1.0)
                        otf_t = torch.fft.fftshift(torch.fft.fft2(psf_t))
                        mtf_t = torch.abs(otf_t)
                        mtf_t = mtf_t / torch.clamp(torch.max(mtf_t), min=1e-12)
                        _freq_t = torch.fft.fftshift(torch.fft.fftfreq(bins, d=step, device=device))
                        if hasattr(torch.cuda, "synchronize"):
                            torch.cuda.synchronize()
                        gpu_results.append(("Torch", time.perf_counter() - t2))
                except Exception as exc:
                    self.append_debug(f"Benchmark Torch path failed: {_short_error_message(exc)}")

            self.append_progress(
                f"Benchmark traced rays={x_local.size} | trace workers={workers} | bins={bins} | CPU post={cpu_sec:.6f}s"
            )
            if gpu_results:
                gpu_results.sort(key=lambda item: item[1])
                best_name, best_sec = gpu_results[0]
                speedup = cpu_sec / max(best_sec, 1e-12)
                for name, timing in gpu_results:
                    self.append_progress(f"Benchmark {name} post={timing:.6f}s")
                self.append_progress(
                    f"Benchmark best GPU={best_name} {best_sec:.6f}s | speedup={speedup:.2f}x"
                )
                gpu_summary = ", ".join(f"{name}={timing:.6f}s" for name, timing in gpu_results)
                self.append_debug(
                    f"PSF/MTF benchmark: rays={x_local.size}, workers={workers}, cpu={cpu_sec:.6f}s, {gpu_summary}, best={best_name}, speedup={speedup:.2f}x"
                )
            else:
                self.append_progress("Benchmark GPU post=unavailable")
                self.append_debug(
                    f"PSF/MTF benchmark: rays={x_local.size}, workers={workers}, cpu={cpu_sec:.6f}s, gpu=unavailable"
                )
            self.status_var.set("Benchmark PSF/MTF completed")
        except Exception as exc:
            self.append_progress(f"Benchmark PSF/MTF failed: {exc}")
            self.append_debug(f"Benchmark PSF/MTF failed: {exc}")
            self.status_var.set("Benchmark PSF/MTF failed")

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
            return 5.0
        try:
            value = float(var.get())
        except ValueError:
            return 5.0
        return max(0.0, value)

    def _operand_mtf_mode(self, label: str) -> str:
        var = self.operand_mtf_mode_vars.get(label)
        if var is None:
            return "average"
        value = var.get().strip().lower()
        if value in {"tangential", "sagittal", "average"}:
            return value
        return "average"

    def _operand_mtf_algorithm(self, label: str) -> str:
        var = self.operand_mtf_algorithm_vars.get(label)
        if var is None:
            return "psf_fft"
        value = var.get().strip().lower()
        if value == "lsf fft":
            return "lsf_fft"
        return "psf_fft"

    def _mtf_analysis_settings(self) -> dict[str, float | int | str]:
        return {
            "wavelength": self._current_wavelength(),
            "surface_index": self._analysis_surface_index(),
            "aperture_type": self._current_aperture_type(),
            "aperture_value": self._current_aperture_value(),
            "field_type": ("angle" if self._current_object_mode() == "Infinity" else "height"),
            "field_x": 0.0,
            "field_y": (self._current_field_angle_deg() if self._current_object_mode() == "Infinity" else self._current_field_height()),
            "algorithm": self._operand_mtf_algorithm("MTF @ freq"),
        }

    @staticmethod
    def _field_type_unit(field_type: str) -> str:
        units = {
            "Angle": "deg",
            "Object Height": "mm",
            "Paraxial Image Height": "mm",
            "Real Image Height": "mm",
        }
        return units.get(str(field_type), "")

    @staticmethod
    def _format_field_sample_value(value: float) -> str:
        return f"{float(value):.3f}".rstrip("0").rstrip(".")

    def _parse_numeric_series(self, value: str) -> list[float]:
        text = str(value or "").strip()
        if not text:
            return []
        samples: list[float] = []
        invalid: list[str] = []
        for token in re.split(r"[\s,;]+", text):
            if not token:
                continue
            try:
                samples.append(float(token))
            except ValueError:
                invalid.append(token)
        if invalid:
            self.append_debug(f"Ignoring invalid numeric samples: {', '.join(invalid)}")
        return samples

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
        self._set_analysis_parallel_status(self.analysis_mode or "2D", 1, False)
        self._clear_cardinal_marker_artists()
        self._clear_layout_selection_overlay()
        self._layout_pick_regions = {}
        self._last_optics_info = None
        self._analysis_ax = None
        if not self.rows:
            self.ax.clear()
            self.ax.set_title("Axial Layout")
            self._configure_plot_hover_hints()
            self.canvas.draw_idle()
            self._autosave_plot()
            return

        max_radius = 1.0
        for row in self.rows:
            radius = max(row.diameter / 2.0, 0.5)
            max_radius = max(max_radius, radius)

        self.figure.clear()
        if self.analysis_mode in {"none", "native_off_axis"}:
            self.ax = self.figure.add_subplot(111)
            analysis_ax = None
            self._analysis_ax = None
        else:
            gs = self.figure.add_gridspec(1, 2, width_ratios=[4.2, 1.35], wspace=0.22)
            self.ax = self.figure.add_subplot(gs[0])
            analysis_ax = self.figure.add_subplot(gs[1])
            self._analysis_ax = analysis_ax

        if self.analysis_mode == "native_off_axis" and self._has_off_axis_geometry():
            self._set_analysis_parallel_status("Native", 1, False)
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
            self._configure_plot_hover_hints()
            self._update_layout_selection_overlay()
            self.canvas.draw_idle()
            self._autosave_plot()
            return

        if self._is_folded_mirror_preview_mode():
            self._set_analysis_parallel_status("Folded preview", 1, False)
            self._plot_folded_mirror_preview(analysis_ax)
            self.ax.grid(True, alpha=0.2)
            self.ax.set_xlabel("Fold X [mm]")
            self.ax.set_ylabel("Fold Y [mm]")
            self.ax.set_title("")
            self.figure.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.92, wspace=0.28)
            self.figure.text(0.5, 0.035, "KrakenOS Layout", ha="center", va="center")
            self._sync_object_controls()
            self._configure_plot_hover_hints()
            self._update_layout_selection_overlay()
            self.canvas.draw_idle()
            self._autosave_plot()
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
            self._refresh_3d_inspector_if_open()
            self._rebuild_layout_pick_regions(system)
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
            self._draw_reference_plane_labels()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                optics_info = self._collect_optics_info(system, rays, wavelength)
            self._last_optics_info = dict(optics_info)
            self._draw_optics_markers(optics_info)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                self._plot_analysis(analysis_ax, system, rays, wavelength)
                self._update_results(system, rays, wavelength, optics_info)
            self.status_var.set(f"Plot refreshed | {self._last_analysis_label} | {self._analysis_compute_summary()}")
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
        self._configure_plot_hover_hints()
        self._update_layout_selection_overlay()
        self.canvas.draw_idle()
        self._autosave_plot()
        if self._initial_layout_passes < 40:
            self.after(50, self._set_initial_pane_layout)

    def _autosave_plot(self) -> None:
        if not self.auto_save_plot_var.get():
            return
        if self._autosave_after_id is not None:
            try:
                self.after_cancel(self._autosave_after_id)
            except Exception:
                pass
        self._autosave_after_id = self.after(400, self._do_autosave_plot)

    def _do_autosave_plot(self) -> None:
        self._autosave_after_id = None
        if not self.auto_save_plot_var.get():
            return
        if self.winfo_width() < 1200 or self.winfo_height() < 700:
            self._autosave_after_id = self.after(400, self._do_autosave_plot)
            return
        try:
            AUTO_PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
            self.update_idletasks()
            self.canvas.draw()
            self.figure.savefig(AUTO_PLOT_PATH, dpi=150)
        except Exception as exc:
            self.append_debug(f"Auto-save plot failed: {exc}")

    def _plot_analysis(self, analysis_ax, system, rays, wavelength: float) -> None:
        if analysis_ax is None:
            return
        analysis_ax.clear()
        spot_field_series: list[tuple[np.ndarray, np.ndarray, float]] = []
        try:
            if self.analysis_mode in {"spot", "rms"}:
                field_type = "angle" if self._current_object_mode() == "Infinity" else "height"
                center_field = self._current_field_angle_deg() if field_type == "angle" else self._current_field_height()
                sampled_fields = [center_field]
                if self.analysis_mode == "spot":
                    sampled_fields = self._sample_field_values(center_field)
                    if not sampled_fields:
                        sampled_fields = [center_field]

                field_results = []
                analysis_workers = 1
                for field_value in sampled_fields:
                    Xi, Yi, Zi, Li, Mi, Ni, worker_count = self._build_geometric_image_samples_full(
                        system,
                        wavelength,
                        sample_count=max(24, self._current_ray_count() * 6),
                        pattern="hexapolar",
                        surface_index=self._analysis_surface_index(),
                        aperture_type=self._current_aperture_type(),
                        aperture_value=self._current_aperture_value(),
                        field_type=field_type,
                        field_x=0.0,
                        field_y=float(field_value),
                    )
                    analysis_workers = max(analysis_workers, int(worker_count))
                    if Xi.size == 0:
                        continue
                    field_results.append((Xi, Yi, Zi, Li, Mi, Ni, float(field_value)))

                if not field_results:
                    X = Y = Z = L = M = N = np.asarray([])
                else:
                    X = np.concatenate([result[0] for result in field_results])
                    Y = np.concatenate([result[1] for result in field_results])
                    Z = np.concatenate([result[2] for result in field_results])
                    L = np.concatenate([result[3] for result in field_results])
                    M = np.concatenate([result[4] for result in field_results])
                    N = np.concatenate([result[5] for result in field_results])
                    if self.analysis_mode == "spot":
                        spot_field_series = [
                            (result[0], result[1], result[6]) for result in field_results
                        ]
            else:
                analysis_rays = self._build_analysis_rays(system, wavelength)
                X, Y, Z, L, M, N = self._pick_image_plane_data(analysis_rays)
                analysis_workers = 1
        except Exception as exc:
            self.append_debug(f"{self.analysis_mode.upper()} analysis trace error: {exc}")
            X = Y = Z = L = M = N = np.asarray([])
            analysis_workers = 1

        if X.size == 0 and self.analysis_mode in {"spot", "rms"}:
            analysis_ax.text(0.5, 0.5, "No ray data", ha="center", va="center")
            analysis_ax.set_axis_off()
            return

        if self.analysis_mode == "spot":
            self._set_analysis_parallel_status("Spot", analysis_workers, True)
            spot_mode = self._current_spot_view_mode()
            if spot_field_series:
                colors = self._field_colors(len(spot_field_series))
                field_unit = "deg" if self._current_object_mode() == "Infinity" else "mm"
                prepared_series: list[tuple[np.ndarray, np.ndarray, float]] = []
                if spot_mode == "Absolute":
                    for x_field, y_field, field_value in spot_field_series:
                        prepared_series.append((x_field, y_field, field_value))
                else:
                    centered_series = []
                    max_span = 1e-3
                    for x_field, y_field, field_value in spot_field_series:
                        cx = float(np.mean(x_field))
                        cy = float(np.mean(y_field))
                        sx = x_field - cx
                        sy = y_field - cy
                        max_span = max(max_span, float(max(np.ptp(sx), np.ptp(sy), 1e-6)))
                        centered_series.append((sx, sy, field_value))
                    if spot_mode == "Grid" and centered_series:
                        cols = max(1, int(np.ceil(np.sqrt(len(centered_series)))))
                        rows = int(np.ceil(len(centered_series) / cols))
                        spacing = max(max_span * 1.8, 2e-3)
                        for index, (sx, sy, field_value) in enumerate(centered_series):
                            row = index // cols
                            col = index % cols
                            offset_x = (col - (cols - 1) / 2.0) * spacing
                            offset_y = ((rows - 1) / 2.0 - row) * spacing
                            prepared_series.append((sx + offset_x, sy + offset_y, field_value))
                    else:
                        prepared_series = centered_series
                X_plot = np.concatenate([item[0] for item in prepared_series]) if prepared_series else X
                Y_plot = np.concatenate([item[1] for item in prepared_series]) if prepared_series else Y
                draw_order = sorted(
                    range(len(prepared_series)),
                    key=lambda idx: abs(float(prepared_series[idx][2])),
                    reverse=True,
                )
                for rank, index in enumerate(draw_order):
                    _x_field, _y_field, field_value = prepared_series[index]
                    analysis_ax.scatter(
                        _x_field,
                        _y_field,
                        s=8,
                        c=colors[index],
                        alpha=0.45,
                        label=f"{field_value:.3g} {field_unit}",
                        zorder=3 + rank,
                    )
                if len(spot_field_series) > 1:
                    analysis_ax.legend(loc="upper right", fontsize=7, title="Field")
            else:
                if spot_mode == "Absolute":
                    X_plot = X
                    Y_plot = Y
                else:
                    X_plot = X - float(np.mean(X))
                    Y_plot = Y - float(np.mean(Y))
                analysis_ax.scatter(X_plot, Y_plot, s=18, c="#c0392b", alpha=0.8)
            analysis_ax.axhline(0.0, color="#2c3e50", linewidth=0.6, alpha=0.5)
            analysis_ax.axvline(0.0, color="#2c3e50", linewidth=0.6, alpha=0.5)
            title_suffix = {
                "Absolute": "Absolute",
                "Centroid": "Centroid Referenced",
                "Grid": "Grid View",
            }.get(spot_mode, "Grid View")
            analysis_ax.set_title(f"Spot Diagram ({title_suffix})")
            analysis_ax.set_xlabel("X [mm]")
            analysis_ax.set_ylabel("Y [mm]")
            if np.ptp(X_plot) < 1e-12 and np.ptp(Y_plot) < 1e-12:
                analysis_ax.set_xlim(float(X_plot[0]) - 1.0, float(X_plot[0]) + 1.0)
                analysis_ax.set_ylim(float(Y_plot[0]) - 1.0, float(Y_plot[0]) + 1.0)
            elif np.ptp(X_plot) < 1e-12:
                center_x = float(np.mean(X_plot))
                span_y = max(float(np.ptp(Y_plot)), 1e-6)
                half_width = max(span_y * 0.35, 1e-3)
                analysis_ax.set_xlim(center_x - half_width, center_x + half_width)
            if spot_mode == "Absolute":
                x_min = float(np.min(X_plot))
                x_max = float(np.max(X_plot))
                y_min = float(np.min(Y_plot))
                y_max = float(np.max(Y_plot))
                x_span = max(x_max - x_min, 1e-9)
                y_span = max(y_max - y_min, 1e-9)
                center_x = 0.5 * (x_min + x_max)
                center_y = 0.5 * (y_min + y_max)
                half_x = max(x_span * 0.7, y_span * 0.55, 0.25)
                half_y = max(y_span * 0.55, x_span * 2.0, 0.25)
                analysis_ax.set_xlim(center_x - half_x, center_x + half_x)
                analysis_ax.set_ylim(center_y - half_y, center_y + half_y)
                analysis_ax.set_aspect("auto")
                analysis_ax.set_box_aspect(1.0)
                analysis_ax.xaxis.set_major_locator(MaxNLocator(5))
                analysis_ax.yaxis.set_major_locator(MaxNLocator(7))
                analysis_ax.tick_params(axis="x", labelrotation=0)
            else:
                analysis_ax.set_aspect("equal", adjustable="box")
            analysis_ax.grid(True, alpha=0.2)
            self.append_debug(f"Spot analysis ok: rays={len(X)}, workers={analysis_workers}")
            return

        if self.analysis_mode == "psf":
            try:
                self._begin_analysis_progress("PSF analysis")
                field_type = "angle" if self._current_object_mode() == "Infinity" else "height"
                field_y = self._current_field_angle_deg() if field_type == "angle" else self._current_field_height()
                self._update_analysis_progress("Tracing rays", 1, 3)
                x_local, y_local, worker_count = self._build_geometric_image_samples(
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
                if x_local.size < 4:
                    raise RuntimeError("Not enough image-plane samples for PSF")
                span_x = max(float(np.ptp(x_local)), 1e-3)
                span_y = max(float(np.ptp(y_local)), 1e-3)
                span = max(span_x, span_y) * 1.25
                bins = 128
                hist, xedges, yedges, accelerator = self._compute_psf_histogram(x_local, y_local, bins, span)
                psf = hist.T
                psf /= max(float(np.max(psf)), 1e-12)
                extent = [float(xedges[0]), float(xedges[-1]), float(yedges[0]), float(yedges[-1])]
                image = analysis_ax.imshow(psf, origin="lower", extent=extent, cmap="inferno", aspect="equal")
                self._set_analysis_parallel_status("PSF", worker_count, True)
                self._set_analysis_accelerator(accelerator)
                analysis_ax.set_title(f"Geometric PSF  |  {field_type}={field_y:.3g}  |  {wavelength:.4g} um")
                analysis_ax.set_xlabel("X [mm]")
                analysis_ax.set_ylabel("Y [mm]")
                analysis_ax.grid(False)
                self.figure.colorbar(image, ax=analysis_ax, fraction=0.046, pad=0.04, label="Normalized intensity")
                self.append_debug(
                    f"PSF analysis ok: rays={x_local.size}, bins={bins}, workers={worker_count}, accel={accelerator}"
                )
                self._update_analysis_progress("Rendering", 3, 3)
                self._finish_analysis_progress("PSF analysis", success=True)
            except Exception as exc:
                self._set_analysis_parallel_status("PSF", 1, True)
                self.append_debug(f"PSF analysis error: {exc}")
                analysis_ax.text(0.5, 0.5, "PSF analysis unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("PSF analysis", success=False)
            return

        if self.analysis_mode == "rms":
            self._set_analysis_parallel_status("RMS", analysis_workers, True)
            rms, cenX, cenY = Kos.RMS(X, Y, Z, L, M, N)
            radii = np.sqrt((X - cenX) ** 2 + (Y - cenY) ** 2)
            bins = min(max(5, int(np.sqrt(max(len(radii), 1)))), 20)
            analysis_ax.hist(radii, bins=bins, color="#4f81bd", edgecolor="white")
            analysis_ax.set_title(f"Spot Radius Histogram  |  RMS = {float(rms):.4g} mm")
            analysis_ax.set_xlabel("Radius [mm]")
            analysis_ax.set_ylabel("Count")
            analysis_ax.grid(True, axis="y", alpha=0.2)
            self.append_debug(f"RMS analysis ok: rays={len(X)}, workers={analysis_workers}")
            return

        if self.analysis_mode == "pupil":
            try:
                self._set_analysis_parallel_status("Pupil", 1, False)
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
                self._set_analysis_parallel_status("Seidel", 1, False)
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
                self._set_analysis_parallel_status("Wavefront", 1, False)
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
                self._set_analysis_parallel_status("Field curvature / distortion", 1, True)
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
                    worker_counts: list[int] = []
                    for field_value in field_samples:
                        completed_steps += 1
                        self._update_analysis_progress(
                            f"Sampling {axis_name}-field",
                            completed_steps,
                            total_steps,
                        )
                        field_x = field_value if axis_name == "X" else 0.0
                        field_y = field_value if axis_name == "Y" else 0.0
                        x_local, y_local, _z_local, l_local, m_local, n_local, worker_count = self._build_geometric_image_samples_full(
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
                        if x_local.size < 4:
                            continue
                        worker_counts.append(worker_count)

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
                        "workers": np.asarray([max(worker_counts) if worker_counts else 1], dtype=float),
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
                    + ", ".join(
                        f"{axis}={len(data['fields'])},workers={int(data['workers'][0])}"
                        for axis, data in axis_results.items()
                    )
                )
                self._set_analysis_parallel_status(
                    "Field curvature / distortion",
                    max(int(data["workers"][0]) for data in axis_results.values()),
                    True,
                )
                self._finish_analysis_progress("Field curvature / distortion", success=True)
            except Exception as exc:
                self._set_analysis_parallel_status("Field curvature / distortion", 1, True)
                self.append_debug(f"Field curvature/distortion error: {exc}")
                analysis_ax.text(0.5, 0.5, "Field curvature/distortion unavailable", ha="center", va="center")
                analysis_ax.set_axis_off()
                self._finish_analysis_progress("Field curvature / distortion", success=False)
            return

        if self.analysis_mode == "mtf":
            try:
                self._set_analysis_parallel_status("MTF", 1, True)
                self._begin_analysis_progress("MTF analysis")
                mtf_settings = self._mtf_analysis_settings()
                wavelength = float(mtf_settings["wavelength"])
                mtf_mode = self._operand_mtf_mode("MTF @ freq")
                mtf_algorithm = str(mtf_settings.get("algorithm", "psf_fft"))
                target_freq = self._current_mtf_frequency()
                field_samples = self._resolved_mtf_field_samples("MTF @ freq")
                if not field_samples:
                    raise RuntimeError("No valid MTF field samples")

                sample_results: list[dict[str, object]] = []
                max_workers = 1
                accelerators: set[str] = set()
                total_steps = max(2, len(field_samples) + 1)
                for index, sample in enumerate(field_samples, start=1):
                    legend = str(sample["legend"])
                    self._update_analysis_progress(
                        f"MTF field {index}/{len(field_samples)}: {legend}",
                        index,
                        total_steps,
                    )
                    try:
                        result = self._compute_diffraction_mtf_sample(
                            system,
                            wavelength=wavelength,
                            surface_index=int(mtf_settings["surface_index"]),
                            aperture_type=str(mtf_settings["aperture_type"]),
                            aperture_value=float(mtf_settings["aperture_value"]),
                            field_type=str(sample["field_type"]),
                            field_x=float(sample["field_x"]),
                            field_y=float(sample["field_y"]),
                        )
                        self.append_debug(
                            f"MTF sample {legend}: diffraction ok: method={result.get('phase_method', 'Phase')}, terms={result['used_terms']}, samples={result['sample_count']}"
                        )
                    except Exception as diff_exc:
                        self.append_debug(f"MTF sample {legend}: diffraction failed: {diff_exc}")
                        try:
                            result = self._compute_geometric_mtf_sample(
                                system,
                                wavelength=wavelength,
                                surface_index=int(mtf_settings["surface_index"]),
                                aperture_type=str(mtf_settings["aperture_type"]),
                                aperture_value=float(mtf_settings["aperture_value"]),
                                field_type=str(sample["field_type"]),
                                field_x=float(sample["field_x"]),
                                field_y=float(sample["field_y"]),
                                algorithm=mtf_algorithm,
                            )
                            self.append_debug(
                                "MTF sample {legend}: geometric ok: rays={rays}, pupil_samp={pupil_samp}, workers={workers}, accel={accel}".format(
                                    legend=legend,
                                    rays=int(result["sample_count"]),
                                    pupil_samp=int(result.get("pupil_samp", 0)),
                                    workers=int(result["worker_count"]),
                                    accel=str(result["accelerator"]),
                                )
                            )
                        except Exception as geom_exc:
                            self.append_debug(f"MTF sample {legend}: geometric failed: {geom_exc}")
                            continue

                    plot_freq = np.asarray(result["plot_freq"], dtype=float)
                    plot_tan = np.asarray(result["plot_tan"], dtype=float)
                    plot_sag = np.asarray(result["plot_sag"], dtype=float)
                    if plot_freq.size == 0 or plot_tan.size == 0 or plot_sag.size == 0:
                        continue

                    tan_value = float(np.interp(target_freq, plot_freq, plot_tan, left=plot_tan[0], right=plot_tan[-1]))
                    sag_value = float(np.interp(target_freq, plot_freq, plot_sag, left=plot_sag[0], right=plot_sag[-1]))
                    if mtf_mode == "tangential":
                        selected_value = tan_value
                        selected_label = "Tangential"
                    elif mtf_mode == "sagittal":
                        selected_value = sag_value
                        selected_label = "Sagittal"
                    else:
                        selected_value = 0.5 * (tan_value + sag_value)
                        selected_label = "Average"

                    result.update(
                        {
                            "legend": legend,
                            "basis": str(sample["basis"]),
                            "unit": str(sample["unit"]),
                            "display_x": float(sample["display_x"]),
                            "display_y": float(sample["display_y"]),
                            "tan_value": tan_value,
                            "sag_value": sag_value,
                            "selected_value": float(selected_value),
                            "selected_label": selected_label,
                        }
                    )
                    sample_results.append(result)
                    max_workers = max(max_workers, int(result.get("worker_count", 1)))
                    accelerators.add(str(result.get("accelerator", "CPU")))

                if not sample_results:
                    raise RuntimeError("MTF analysis unavailable for all selected field samples")

                colors = self._field_colors(len(sample_results))
                max_plot_freq = 0.0
                label_specs: list[dict[str, object]] = []
                for color, result in zip(colors, sample_results):
                    plot_freq = np.asarray(result["plot_freq"], dtype=float)
                    plot_tan = np.asarray(result["plot_tan"], dtype=float)
                    plot_sag = np.asarray(result["plot_sag"], dtype=float)
                    max_plot_freq = max(max_plot_freq, float(plot_freq[-1]))
                    legend = str(result["legend"])
                    overlap = bool(
                        np.allclose(
                            plot_tan,
                            plot_sag,
                            rtol=1e-3,
                            atol=max(1e-4, 5e-3 * float(max(np.max(np.abs(plot_tan)), np.max(np.abs(plot_sag)), 1e-9))),
                        )
                    )
                    result["ts_overlap"] = overlap
                    if overlap:
                        analysis_ax.plot(
                            plot_freq,
                            plot_tan,
                            label=f"T=S {legend}",
                            color=color,
                            linewidth=1.1,
                            alpha=1.0,
                            linestyle="-",
                            zorder=3,
                        )
                        label_specs.append(
                            {
                                "label": f"T=S {legend}",
                                "curve_x": plot_freq,
                                "curve_y": plot_tan,
                                "color": color,
                                "linestyle": (0, (2, 2)),
                            }
                        )
                    else:
                        analysis_ax.plot(
                            plot_freq,
                            plot_tan,
                            label=f"T {legend}",
                            color=color,
                            linewidth=1.1,
                            alpha=1.0,
                            linestyle="-",
                            zorder=3,
                        )
                        analysis_ax.plot(
                            plot_freq,
                            plot_sag,
                            label=f"S {legend}",
                            color=color,
                            linewidth=1.0,
                            alpha=1.0,
                            linestyle=(0, (6, 3)),
                            zorder=4,
                        )
                        label_specs.extend(
                            [
                                {
                                    "label": f"T {legend}",
                                    "curve_x": plot_freq,
                                    "curve_y": plot_tan,
                                    "color": color,
                                    "linestyle": (0, (2, 2)),
                                },
                                {
                                    "label": f"S {legend}",
                                    "curve_x": plot_freq,
                                    "curve_y": plot_sag,
                                    "color": color,
                                    "linestyle": (0, (1, 2)),
                                },
                            ]
                        )
                analysis_ax.plot(
                    [target_freq, target_freq],
                    [0.0, 0.08],
                    color="#2c3e50",
                    linewidth=0.9,
                    linestyle="-",
                    alpha=0.9,
                    zorder=1.8,
                )
                analysis_ax.text(
                    target_freq,
                    0.085,
                    f"ref {target_freq:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                    color="#2c3e50",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.45, "pad": 0.2},
                )
                basis = str(sample_results[0]["basis"])
                unit = str(sample_results[0]["unit"])
                x_text = self._format_field_sample_value(float(sample_results[0]["display_x"]))
                method_label = str(sample_results[0].get("method", "MTF"))
                dl_fc = None
                try:
                    effl, _ppa, _ppp = self._exact_paraxial_cardinals(wavelength)
                    pupil_ref = Kos.PupilCalc(
                        system,
                        int(mtf_settings["surface_index"]),
                        wavelength,
                        str(mtf_settings["aperture_type"]),
                        float(mtf_settings["aperture_value"]),
                    )
                    ep_diameter = max(2.0 * abs(float(getattr(pupil_ref, "RadPupInp", 0.0))), 1e-9)
                    f_number = abs(float(effl)) / ep_diameter
                    if np.isfinite(f_number) and f_number > 1e-12:
                        dl_fc = 1.0 / (max(wavelength, 1e-12) * 1e-3 * f_number)
                        if np.isfinite(dl_fc) and dl_fc > 0.0:
                            dl_freq = np.linspace(0.0, min(max_plot_freq, max(100.0, target_freq * 2.5)), 512)
                            nu = np.clip(dl_freq / dl_fc, 0.0, 1.0)
                            dl_curve = (2.0 / np.pi) * (
                                np.arccos(nu) - nu * np.sqrt(np.clip(1.0 - nu * nu, 0.0, 1.0))
                            )
                            dl_curve = np.where(dl_freq <= dl_fc, dl_curve, 0.0)
                            analysis_ax.plot(
                                dl_freq,
                                dl_curve,
                                color="#475569",
                                linewidth=1.3,
                                linestyle=(0, (4, 2)),
                                alpha=0.9,
                                label="DL ref",
                                zorder=2,
                            )
                except Exception:
                    dl_fc = None
                analysis_ax.set_title(
                    f"MTF ({method_label})  |  {basis} samples  |  ref {target_freq:.1f} cy/mm  |  {wavelength:.4g} um"
                )
                analysis_ax.set_xlabel("Spatial frequency [cycles/mm]")
                analysis_ax.set_ylabel("MTF")
                analysis_ax.set_ylim(0.0, 1.05)
                analysis_ax.set_xlim(0.0, min(max_plot_freq, max(100.0, target_freq * 2.5)))
                analysis_ax.set_aspect("auto")
                analysis_ax.set_box_aspect(0.62)
                analysis_ax.grid(True, alpha=0.2)
                y_top = analysis_ax.get_ylim()[1]
                if label_specs:
                    x_min, x_max = [float(v) for v in analysis_ax.get_xlim()]
                    y_min, _ = [float(v) for v in analysis_ax.get_ylim()]
                    active_x_max = x_max
                    for spec in label_specs:
                        curve_x = np.asarray(spec["curve_x"], dtype=float)
                        curve_y = np.asarray(spec["curve_y"], dtype=float)
                        nonzero = curve_x[curve_y > 0.03]
                        if nonzero.size:
                            active_x_max = min(active_x_max, float(np.max(nonzero)))
                    label_left = max(x_min + 0.06 * (x_max - x_min), min(target_freq + 1.5, active_x_max * 0.25))
                    label_right = max(label_left + 1.0, min(x_max - 0.06 * (x_max - x_min), active_x_max * 0.98))
                    label_x_positions = np.linspace(label_left, label_right, len(label_specs))
                    row_levels = [y_top - 0.02, y_top - 0.07, y_top - 0.12, y_top - 0.17]
                    for index, (spec, label_x) in enumerate(zip(label_specs, label_x_positions)):
                        curve_x = np.asarray(spec["curve_x"], dtype=float)
                        curve_y = np.asarray(spec["curve_y"], dtype=float)
                        marker_value = float(np.interp(label_x, curve_x, curve_y, left=curve_y[0], right=curve_y[-1]))
                        if not np.isfinite(marker_value):
                            continue
                        label_y = max(row_levels[index % len(row_levels)], marker_value + 0.06)
                        analysis_ax.plot(
                            [label_x, label_x],
                            [marker_value, label_y - 0.015],
                            color=str(spec["color"]),
                            linewidth=0.75,
                            linestyle=spec["linestyle"],
                            alpha=0.8,
                            zorder=2.5,
                        )
                        analysis_ax.text(
                            label_x,
                            label_y,
                            str(spec["label"]),
                            rotation=0,
                            ha="center",
                            va="top",
                            fontsize=6.1,
                            color=str(spec["color"]),
                            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.45, "pad": 0.2},
                        )
                if dl_fc is not None and np.isfinite(dl_fc):
                    dl_x = float(np.clip(dl_fc, analysis_ax.get_xlim()[0], analysis_ax.get_xlim()[1]))
                    analysis_ax.axvline(
                        dl_x,
                        color="#475569",
                        linewidth=0.9,
                        linestyle=(0, (1, 2)),
                        alpha=0.7,
                        zorder=1.5,
                    )
                    label_x = min(dl_x, analysis_ax.get_xlim()[1] - 0.5)
                    analysis_ax.text(
                        label_x,
                        0.035,
                        f"DL cutoff {float(dl_fc):.1f} cy/mm",
                        ha="right" if dl_x >= analysis_ax.get_xlim()[1] - 1.0 else "center",
                        va="bottom",
                        fontsize=7,
                        color="#475569",
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.45, "pad": 0.2},
                    )

                self._set_analysis_parallel_status("MTF", max_workers, max_workers > 1)
                if accelerators:
                    accel_summary = "/".join(sorted(accelerators))
                    self._set_analysis_accelerator(accel_summary)
                self._update_analysis_progress("Rendering MTF", total_steps, total_steps)
                self._finish_analysis_progress("MTF analysis", success=True)
            except Exception as exc:
                self._set_analysis_parallel_status("MTF", 1, True)
                self.append_debug(f"MTF analysis error: {exc}")
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

    def _current_spot_view_mode(self) -> str:
        value = getattr(self, "spot_view_mode_var", None)
        if value is None:
            return "Grid"
        mode = value.get().strip()
        if mode in {"Grid", "Absolute", "Centroid"}:
            return mode
        return "Grid"

    def _current_field_value(self) -> float:
        try:
            return float(self.field_value_var.get())
        except ValueError:
            return 0.0

    def _current_field_angle_deg(self) -> float:
        return float(self._field_metrics().get("angle_deg", 0.0))

    def _current_field_height(self) -> float:
        return float(self._field_metrics().get("object_height", 0.0))

    def _field_metrics_for_value(self, field_type: str, raw_value: float) -> dict[str, float]:
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

    def _field_metrics(self) -> dict[str, float]:
        return self._field_metrics_for_value(self._current_field_type(), self._current_field_value())

    def _field_metrics_summary(self) -> dict[str, float]:
        field_type = self._current_field_type()
        sample_values = self._sample_field_values(self._current_field_value())
        if not sample_values:
            sample_values = [self._current_field_value()]
        metrics = [self._field_metrics_for_value(field_type, value) for value in sample_values]
        current_metrics = self._field_metrics()
        max_paraxial = max(abs(float(item.get("paraxial_image_height", 0.0))) for item in metrics) if metrics else 0.0
        max_real = max(abs(float(item.get("real_image_height", 0.0))) for item in metrics) if metrics else 0.0
        return {
            "current_angle_deg": float(current_metrics.get("angle_deg", 0.0)),
            "current_object_height": float(current_metrics.get("object_height", 0.0)),
            "current_paraxial_image_height": float(current_metrics.get("paraxial_image_height", 0.0)),
            "current_real_image_height": float(current_metrics.get("real_image_height", 0.0)),
            "max_paraxial_image_height": float(max_paraxial),
            "max_real_image_height": float(max_real),
            "image_diameter": float(max(2.0 * max_real, 0.0)),
        }

    def _current_effl_estimate(self) -> float:
        try:
            effl, _ppa, _ppp = self._exact_paraxial_cardinals(self._current_wavelength())
            return max(abs(float(effl)), 1e-6)
        except Exception:
            pass
        if self.last_system is not None:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    _a, _b, _c, _d, effl, *_rest = self.last_system.EFL(self._current_wavelength())  # type: ignore[misc]
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
                line.set_zorder(max(float(line.get_zorder()), 20.0))
            else:
                line.set_zorder(min(float(line.get_zorder()), 10.0))

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
        show_clipped_rays = self.show_clipped_rays_var.get()
        final_surface_index = max(0, len(self.rows) - 1)
        ray_count_hint = max(1, self._preview_field_ray_count)
        ray_linewidth = 1.1 if ray_count_hint <= 9 else 0.8
        ray_alpha = 0.92 if ray_count_hint <= 9 else 0.72
        if self._has_off_axis_geometry():
            if show_clipped_rays:
                before = len(self.ax.lines)
                try:
                    Plot2DRays(rays, 0, 0, self.ax, 0)
                    for line in self.ax.lines[before:]:
                        line.set_linewidth(ray_linewidth)
                        line.set_alpha(ray_alpha)
                except Exception:
                    for ray in rays.CC:
                        points = np.asarray(ray, dtype=float)
                        if points.shape[0] < 2:
                            continue
                        self.ax.plot(points[:, 2], points[:, 1], color="#39FF14", linewidth=ray_linewidth, alpha=ray_alpha)
            else:
                for index, ray in enumerate(rays.CC):
                    points = np.asarray(ray, dtype=float)
                    if points.shape[0] < 2:
                        continue
                    try:
                        surf_ids = np.asarray(rays.SURFACE[index], dtype=int).ravel()
                        if surf_ids.size == 0 or int(surf_ids[-1]) != final_surface_index:
                            continue
                    except Exception:
                        continue
                    self.ax.plot(points[:, 2], points[:, 1], color="#39FF14", linewidth=ray_linewidth, alpha=ray_alpha)
            return
        ray_count = max(1, self._preview_field_ray_count)
        field_count = max(1, self._current_field_count())
        colors = self._field_colors(field_count)
        for index, ray in enumerate(rays.CC):
            points = np.asarray(ray)
            if points.shape[0] < 2:
                continue
            if not show_clipped_rays:
                try:
                    surf_ids = np.asarray(rays.SURFACE[index], dtype=int).ravel()
                    if surf_ids.size == 0 or int(surf_ids[-1]) != final_surface_index:
                        continue
                except Exception:
                        continue
            field_index = min(index // ray_count, field_count - 1)
            color = colors[field_index]
            self.ax.plot(points[:, 2], points[:, 1], color=color, linewidth=ray_linewidth, alpha=ray_alpha)

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
    def _snap_display_direction(direction: np.ndarray, tolerance: float = 0.03) -> np.ndarray:
        d = np.asarray(direction, dtype=float)
        norm = np.linalg.norm(d)
        if norm <= 1e-12:
            return np.array([0.0, 1.0], dtype=float)
        d /= norm
        if abs(d[0]) <= tolerance:
            return np.array([0.0, 1.0 if d[1] >= 0.0 else -1.0], dtype=float)
        if abs(d[1]) <= tolerance:
            return np.array([1.0 if d[0] >= 0.0 else -1.0, 0.0], dtype=float)
        return d

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
        ray_paths: list[tuple[list[np.ndarray], bool]] = []
        for field_index, field_value in enumerate(field_values):
            if self._current_object_mode() == "Infinity":
                for pupil_y in pupil_samples:
                    d = direction.copy()
                    origin = point + np.array([float(field_value) + float(pupil_y), 0.0], dtype=float)
                    p = origin.copy()
                    path = [origin.copy()]
                    current_dir = d
                    current_medium = 1.0
                    reached_image = False
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
                            current_dir = self._snap_display_direction(
                                self._reflect_2d(current_dir, float(row.tilt_x))
                            )
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
                            reached_image = True
                            break
                    ray_paths.append((path, reached_image))
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
                    reached_image = False
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
                            current_dir = self._snap_display_direction(
                                self._reflect_2d(current_dir, float(row.tilt_x))
                            )
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
                            reached_image = True
                            break
                    ray_paths.append((path, reached_image))
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
            show_clipped_rays = self.show_clipped_rays_var.get()
            for index, (path, reached_image) in enumerate(ray_paths):
                if not show_clipped_rays and not reached_image:
                    continue
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
                ("Show clipped rays", "Yes" if self.show_clipped_rays_var.get() else "No"),
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
                current_dir = self._snap_display_direction(self._reflect_2d(current_dir, float(row.tilt_x)))
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
        self._refresh_3d_inspector_if_open()

        extent_points: list[np.ndarray] = []
        folded_elements = None
        folded_max_half = max_radius
        native_overlay_counts: dict[str, int] = {}
        native_overlay_curves: dict[int, np.ndarray] = {}
        native_overlay_metrics: list[tuple[str, str]] = []
        if folded_visual_mode:
            _, _, folded_max_half, extent_points, folded_elements = self._draw_folded_layout_geometry()
            if use_native_surfaces and self.show_native_overlays_var.get():
                native_overlay_counts, native_overlay_metrics, native_overlay_curves = self._overlay_native_folded_surfaces(
                    system,
                    folded_elements,
                    extent_points,
                )
                self._overlay_native_folded_lens_bodies(folded_elements, native_overlay_curves, extent_points)
                self._overlay_native_folded_lens_edges(folded_elements, native_overlay_curves, extent_points)
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
        show_clipped_rays = self.show_clipped_rays_var.get()
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
        native_hit_points_by_surface = (
            self._native_folded_hit_points_by_surface(system, rays, folded_elements)
            if folded_elements is not None
            else {}
        )
        for index, ray in enumerate(rays.CC):
            pts = np.asarray(ray, dtype=float)
            if pts.shape[0] < 2:
                continue
            ray_lengths.append(int(pts.shape[0]))
            last_surface: int | None = None
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
            if not show_clipped_rays and last_surface != final_surface_index:
                continue
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
            if self.show_native_hit_labels_var.get():
                self._overlay_native_folded_surface_hits(
                    folded_elements,
                    surface_hit_counts,
                    folded_max_half,
                    extent_points,
                )
            if self.show_native_active_spans_var.get():
                self._overlay_native_folded_active_spans(
                    folded_elements,
                    native_hit_points_by_surface,
                    native_overlay_curves,
                    extent_points,
                )
            if self.show_native_overlays_var.get() or self.show_native_active_spans_var.get() or self.show_native_hit_labels_var.get():
                self._draw_native_folded_legend(
                    native_overlay_counts,
                    show_hit_count=self.show_native_hit_labels_var.get(),
                    show_active_span=self.show_native_active_spans_var.get(),
                )

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
                ("Show clipped rays", "Yes" if self.show_clipped_rays_var.get() else "No"),
                ("Native rays", str(len(rays.CC))),
                ("Rays to image", str(rays_reaching_image)),
                ("Avg ray hits", f"{np.mean(ray_lengths):.3g}" if ray_lengths else "0"),
                ("Max ray hits", str(max(ray_lengths) if ray_lengths else 0)),
                ("Surface hits", self._native_surface_count_summary(surface_hit_counts)),
                ("Last-hit surfaces", self._native_surface_count_summary(last_surface_counts)),
                ("Native mirrors", str(native_overlay_counts.get("mirror", 0))),
                ("Native refractive", str(native_overlay_counts.get("standard", 0))),
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
        if native_overlay_counts:
            overlay_summary = ", ".join(f"{key}={value}" for key, value in sorted(native_overlay_counts.items()))
            self.append_debug(f"Native folded overlays: {overlay_summary}")
        for label, metric in native_overlay_metrics:
            self.append_debug(f"Native surface metric {label}: {metric}")
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

    def _overlay_native_folded_surfaces(
        self, system, elements, extent_points: list[np.ndarray]
    ) -> tuple[dict[str, int], list[tuple[str, str]], dict[int, np.ndarray]]:
        transforms = getattr(system, "TRANS_2A", None)
        surfaces = getattr(system, "AAA", None)
        if transforms is None or surfaces is None or not elements:
            return {}, [], {}
        block_count = min(len(self.rows), getattr(surfaces, "n_blocks", 0), len(transforms))
        counts = {"mirror": 0, "standard": 0}
        metrics: list[tuple[str, str]] = []
        curves: dict[int, np.ndarray] = {}
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
            display_pts, metric = self._map_native_surface_to_folded(surface_type, center2, row, branch_dir, pts)
            if display_pts is None or display_pts.shape[0] < 8:
                continue
            self.ax.plot(display_pts[:, 0], display_pts[:, 1], color="white", linewidth=4.6, alpha=0.96)
            self.ax.plot(display_pts[:, 0], display_pts[:, 1], color="#111111", linewidth=2.2, alpha=0.98)
            extent_points.extend(display_pts.tolist())
            counts["standard" if surface_type == "Standard" else "mirror"] += 1
            curves[surface_index] = np.asarray(display_pts, dtype=float)
            if metric:
                metrics.append((f"{surface_index}:{row.name or row.surface}", metric))
        return counts, metrics, curves

    def _draw_native_folded_legend(
        self,
        native_overlay_counts: dict[str, int],
        show_hit_count: bool,
        show_active_span: bool,
    ) -> None:
        handles = [
            Line2D([0], [0], color="#202020", linewidth=1.2, label="Base object/image"),
            Line2D([0], [0], color="#111111", linewidth=2.2, label="Native optical surface"),
            Line2D([0], [0], color="#39FF14", linewidth=1.8, label="Displayed ray"),
        ]
        if show_hit_count:
            handles.insert(
                2,
                Line2D([0], [0], color="#d97706", linewidth=2.4, linestyle="--", alpha=0.45, label="Native hit count"),
            )
        if show_active_span:
            handles.insert(
                3 if show_hit_count else 2,
                Line2D([0], [0], color="#f59e0b", linewidth=3.2, label="Native active span"),
            )
        title_parts = []
        if native_overlay_counts.get("mirror", 0):
            title_parts.append(f"M{native_overlay_counts['mirror']}")
        if native_overlay_counts.get("standard", 0):
            title_parts.append(f"S{native_overlay_counts['standard']}")
        title = "Native overlays"
        if title_parts:
            title += f" ({', '.join(title_parts)})"
        self.ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9, title=title, title_fontsize=8)

    def _overlay_native_folded_lens_edges(self, elements, native_overlay_curves: dict[int, np.ndarray], extent_points) -> None:
        if not elements or not native_overlay_curves:
            return
        for group in self._native_lens_surface_groups(elements, native_overlay_curves):
            self._draw_native_edge_group(group, native_overlay_curves, extent_points)

    def _overlay_native_folded_lens_bodies(self, elements, native_overlay_curves: dict[int, np.ndarray], extent_points) -> None:
        if not elements or not native_overlay_curves:
            return
        for group in self._native_lens_surface_groups(elements, native_overlay_curves):
            self._fill_native_lens_group(group, native_overlay_curves, extent_points)

    @staticmethod
    def _native_lens_surface_groups(elements, native_overlay_curves: dict[int, np.ndarray]) -> list[list[int]]:
        groups: list[list[int]] = []
        group: list[int] = []
        for surface_index, (surface_type, _center, row, _branch_dir) in enumerate(elements, start=1):
            if surface_type == "Standard" and surface_index in native_overlay_curves:
                group.append(surface_index)
                # `glass` is the post-surface medium. AIR closes the current optical element.
                if str(row.glass).strip().upper() == "AIR":
                    if len(group) >= 2:
                        groups.append(group[:])
                    group = []
            else:
                if len(group) >= 2:
                    groups.append(group[:])
                group = []
        if len(group) >= 2:
            groups.append(group[:])
        return groups

    def _fill_native_lens_group(self, indices: list[int], native_overlay_curves: dict[int, np.ndarray], extent_points) -> None:
        for first, second in zip(indices, indices[1:]):
            curve_a = np.asarray(native_overlay_curves.get(first), dtype=float)
            curve_b = np.asarray(native_overlay_curves.get(second), dtype=float)
            if curve_a.shape[0] < 8 or curve_b.shape[0] < 8:
                continue
            polygon = np.vstack([curve_a, curve_b[::-1]])
            self.ax.fill(
                polygon[:, 0],
                polygon[:, 1],
                facecolor="#f3f4f6",
                edgecolor="none",
                alpha=0.95,
                zorder=0.5,
            )
            extent_points.extend(polygon.tolist())

    def _draw_native_edge_group(self, indices: list[int], native_overlay_curves: dict[int, np.ndarray], extent_points) -> None:
        for first, second in zip(indices, indices[1:]):
            curve_a = np.asarray(native_overlay_curves.get(first), dtype=float)
            curve_b = np.asarray(native_overlay_curves.get(second), dtype=float)
            if curve_a.shape[0] < 2 or curve_b.shape[0] < 2:
                continue
            a_top = curve_a[np.argmin(curve_a[:, 1])]
            a_bottom = curve_a[np.argmax(curve_a[:, 1])]
            b_top = curve_b[np.argmin(curve_b[:, 1])]
            b_bottom = curve_b[np.argmax(curve_b[:, 1])]
            for p0, p1 in ((a_top, b_top), (a_bottom, b_bottom)):
                self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="white", linewidth=4.2, alpha=0.96)
                self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#111111", linewidth=2.0, alpha=0.98)
                extent_points.extend([p0, p1])

    def _overlay_native_folded_active_spans(
        self,
        elements,
        hit_points_by_surface: dict[int, list[np.ndarray]],
        native_overlay_curves: dict[int, np.ndarray],
        extent_points: list[np.ndarray],
    ) -> None:
        for surface_index, points in sorted(hit_points_by_surface.items()):
            if surface_index <= 0 or surface_index > len(elements):
                continue
            if len(points) < 2:
                continue
            surface_type, center, row, branch_dir = elements[surface_index - 1]
            pts = np.asarray(points, dtype=float)
            if surface_type == "Mirror":
                theta = np.deg2rad(float(row.tilt_x))
                tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
                tangent /= max(np.linalg.norm(tangent), 1e-12)
                along = (pts - center[None, :]) @ tangent
                lo = float(np.min(along))
                hi = float(np.max(along))
                p0 = center + tangent * lo
                p1 = center + tangent * hi
                self.ax.plot(
                    [p0[0], p1[0]],
                    [p0[1], p1[1]],
                    color="#f59e0b",
                    linewidth=3.2,
                    alpha=0.95,
                    solid_capstyle="round",
                )
                extent_points.extend([p0, p1])
            elif surface_type == "Standard":
                curve = native_overlay_curves.get(surface_index)
                if curve is not None and curve.shape[0] >= 8:
                    snapped = self._snap_points_to_curve(pts, curve)
                    if snapped.shape[0] >= 2:
                        self.ax.plot(
                            snapped[:, 0],
                            snapped[:, 1],
                            color="#a855f7",
                            linewidth=2.2,
                            alpha=0.9,
                            linestyle="-.",
                        )
                        extent_points.extend(snapped.tolist())
                else:
                    tangent = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
                    tangent /= max(np.linalg.norm(tangent), 1e-12)
                    along = (pts - center[None, :]) @ tangent
                    order = np.argsort(along)
                    pts = pts[order]
                    self.ax.plot(
                        pts[:, 0],
                        pts[:, 1],
                        color="#a855f7",
                        linewidth=2.2,
                        alpha=0.9,
                        linestyle="-.",
                    )
                    extent_points.extend(pts.tolist())

    @staticmethod
    def _snap_points_to_curve(points: np.ndarray, curve: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        crv = np.asarray(curve, dtype=float)
        if pts.size == 0 or crv.size == 0:
            return np.empty((0, 2), dtype=float)
        chosen: list[int] = []
        for point in pts:
            distances = np.sum((crv - point[None, :]) ** 2, axis=1)
            chosen.append(int(np.argmin(distances)))
        lo = min(chosen)
        hi = max(chosen)
        return crv[lo : hi + 1]

    def _map_native_surface_to_folded(
        self, surface_type: str, center2: np.ndarray, row: SurfaceRow, branch_dir: np.ndarray, pts3: np.ndarray
    ) -> tuple[np.ndarray | None, str | None]:
        pts3 = np.asarray(pts3, dtype=float)
        if pts3.size == 0:
            return None, None
        tangent2 = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
        tangent2 /= max(np.linalg.norm(tangent2), 1e-12)
        axis2 = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
        if surface_type == "Mirror":
            theta = np.deg2rad(float(row.tilt_x))
            mirror_tangent = np.array([np.cos(theta), np.sin(theta)], dtype=float)
            mirror_tangent /= max(np.linalg.norm(mirror_tangent), 1e-12)
            aperture_component = pts3[:, 0] - float(np.mean(pts3[:, 0]))
            if aperture_component.size < 8:
                return None, None
            native_half = max(float(np.max(np.abs(aperture_component))), 1e-9)
            target_half = max(float(row.diameter) / 2.0, 0.5)
            scale = target_half / native_half
            aperture_component = aperture_component * scale
            display_pts = center2[None, :] + np.outer(aperture_component, mirror_tangent)
            order = np.argsort(aperture_component)
            return display_pts[order], f"span_scale={scale:.3g}"
        center_x = float(np.median(pts3[:, 0]))
        tolerance = max(0.08, 0.015 * max(np.ptp(pts3[:, 1]), np.ptp(pts3[:, 2]), 1.0))
        mask = np.abs(pts3[:, 0] - center_x) <= tolerance
        sliced = pts3[mask]
        if sliced.shape[0] < 16:
            tolerance = max(0.2, tolerance * 2.0)
            mask = np.abs(pts3[:, 0] - center_x) <= tolerance
            sliced = pts3[mask]
        if sliced.shape[0] < 8:
            return None, None
        center3 = np.mean(sliced, axis=0)
        centered = sliced - center3
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            return None, None
        vectors = [np.asarray(v, dtype=float) for v in vh]
        yz_ranked = [
            (float(np.linalg.norm(vec[1:])), np.asarray(vec, dtype=float) / max(np.linalg.norm(vec), 1e-12))
            for vec in vectors
        ]
        yz_ranked.sort(key=lambda item: item[0], reverse=True)
        if len(yz_ranked) < 2:
            return None, None
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
        center_sag = 0.0
        if np.any(center_mask):
            center_sag = float(np.mean(sag_component[center_mask]))
        if np.any(edge_mask) and np.any(center_mask) and abs(float(row.rc)) > 1e-9:
            edge_sag = float(np.mean(sag_component[edge_mask]))
            expected_sign = 1.0 if float(row.rc) > 0.0 else -1.0
            observed_sign = np.sign(edge_sag - center_sag)
            if observed_sign != 0.0 and observed_sign != expected_sign:
                sag_component = -sag_component
                center_sag = -center_sag
        # Anchor the mapped native surface at the surface vertex, not at the mean mesh sag.
        sag_component = sag_component - center_sag
        native_half = max(float(np.max(np.abs(aperture_component))), 1e-9)
        target_half = max(float(row.diameter) / 2.0, 0.5)
        aperture_scale = target_half / native_half
        aperture_component = aperture_component * aperture_scale
        native_sag = max(float(np.max(np.abs(sag_component))), 1e-9)
        nominal_sag = 0.0
        if abs(float(row.rc)) > target_half + 1e-9:
            rr = abs(float(row.rc))
            sign = 1.0 if float(row.rc) >= 0.0 else -1.0
            nominal_sag = abs(float(row.rc) - sign * np.sqrt(max(rr * rr - target_half * target_half, 0.0)))
        sag_scale = nominal_sag / native_sag if nominal_sag > 0.0 else 1.0
        sag_component = sag_component * sag_scale
        display_pts = center2[None, :] + np.outer(aperture_component, tangent2) + np.outer(sag_component, axis2)
        order = np.argsort(display_pts[:, 1])
        return display_pts[order], f"span_scale={aperture_scale:.3g}, sag_scale={sag_scale:.3g}"

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

    def _native_folded_hit_points_by_surface(self, system, rays, elements) -> dict[int, list[np.ndarray]]:
        if elements is None:
            return {}
        mapping_meta = self._native_surface_mapping_metadata(system)
        if not mapping_meta:
            return {}
        hit_points: dict[int, list[np.ndarray]] = {}
        for ray_index, (surface_ids_raw, cc_raw) in enumerate(zip(rays.SURFACE, rays.CC)):
            cc = np.asarray(cc_raw, dtype=float)
            if cc.shape[0] < 2:
                continue
            surface_ids = [int(v) for v in np.asarray(surface_ids_raw, dtype=int).ravel().tolist()]
            last_id: int | None = None
            for hit_index, surface_index in enumerate(surface_ids, start=1):
                if surface_index == last_id:
                    continue
                last_id = surface_index
                if hit_index >= len(cc):
                    continue
                display_hit = self._map_native_hit_to_folded_surface(surface_index, cc[hit_index], elements, mapping_meta)
                if display_hit is None:
                    continue
                hit_points.setdefault(surface_index, []).append(display_hit)
        return hit_points

    def _native_surface_mapping_metadata(self, system) -> dict[int, dict[str, np.ndarray | float]]:
        transforms = getattr(system, "TRANS_2A", None)
        surfaces = getattr(system, "AAA", None)
        if transforms is None or surfaces is None:
            return {}
        metadata: dict[int, dict[str, np.ndarray | float]] = {}
        block_count = min(len(self.rows), getattr(surfaces, "n_blocks", 0), len(transforms))
        for surface_index in range(1, block_count):
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
            center = np.mean(pts, axis=0)
            centered = pts - center
            try:
                _, _, vh = np.linalg.svd(centered, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            vectors = [np.asarray(v, dtype=float) / max(np.linalg.norm(v), 1e-12) for v in vh]
            yz_ranked = [(float(np.linalg.norm(vec[1:])), vec) for vec in vectors]
            yz_ranked.sort(key=lambda item: item[0], reverse=True)
            if len(yz_ranked) < 2:
                continue
            basis_a = yz_ranked[0][1]
            basis_b = yz_ranked[1][1]
            u_a = centered @ basis_a
            u_b = centered @ basis_b
            if np.ptp(u_a) >= np.ptp(u_b):
                aperture_basis = basis_a
                sag_basis = basis_b
                aperture_component = u_a
                sag_component = u_b
            else:
                aperture_basis = basis_b
                sag_basis = basis_a
                aperture_component = u_b
                sag_component = u_a
            metadata[surface_index] = {
                "center": center,
                "aperture_basis": aperture_basis,
                "sag_basis": sag_basis,
                "native_half": max(float(np.max(np.abs(aperture_component))), 1e-9),
                "native_sag": max(float(np.max(np.abs(sag_component))), 1e-9),
            }
        return metadata

    def _map_native_hit_to_folded_surface(self, surface_index: int, native_hit: np.ndarray, elements, mapping_meta) -> np.ndarray | None:
        if surface_index <= 0 or surface_index > len(elements):
            return None
        meta = mapping_meta.get(surface_index)
        if meta is None:
            return None
        surface_type, center2, row, branch_dir = elements[surface_index - 1]
        center3 = np.asarray(meta["center"], dtype=float)
        aperture_basis = np.asarray(meta["aperture_basis"], dtype=float)
        sag_basis = np.asarray(meta["sag_basis"], dtype=float)
        delta = np.asarray(native_hit, dtype=float) - center3
        aperture_component = float(delta @ aperture_basis)
        sag_component = float(delta @ sag_basis)
        tangent2 = np.array([-branch_dir[1], branch_dir[0]], dtype=float)
        tangent2 /= max(np.linalg.norm(tangent2), 1e-12)
        axis2 = branch_dir / max(np.linalg.norm(branch_dir), 1e-12)
        native_half = float(meta["native_half"])
        target_half = max(float(row.diameter) / 2.0, 0.5)
        aperture_scale = target_half / max(native_half, 1e-9)
        aperture_component *= aperture_scale
        aperture_component = float(np.clip(aperture_component, -target_half, target_half))
        if surface_type == "Mirror":
            return center2 + tangent2 * aperture_component
        native_sag = float(meta["native_sag"])
        nominal_sag = 0.0
        if abs(float(row.rc)) > target_half + 1e-9:
            rr = abs(float(row.rc))
            sign = 1.0 if float(row.rc) >= 0.0 else -1.0
            nominal_sag = abs(float(row.rc) - sign * np.sqrt(max(rr * rr - target_half * target_half, 0.0)))
        sag_scale = nominal_sag / max(native_sag, 1e-9) if nominal_sag > 0.0 else 1.0
        sag_component *= sag_scale
        if nominal_sag > 0.0:
            sag_limit = nominal_sag * 1.15
            sag_component = float(np.clip(sag_component, -sag_limit, sag_limit))
        return center2 + tangent2 * aperture_component + axis2 * sag_component

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

    def _reference_native_ray_index(self, rays, max_half: float) -> int | None:
        starts = self._preview_ray_start_specs(max_half)
        best_index: int | None = None
        best_score: tuple[int, float] | None = None
        for index, surface_ids_raw in enumerate(rays.SURFACE):
            if index >= len(starts):
                break
            seq = [int(v) for v in np.asarray(surface_ids_raw, dtype=int).ravel().tolist()]
            unique_count = 0
            last = None
            for value in seq:
                if value == last:
                    continue
                last = value
                unique_count += 1
            start_origin, _ = starts[index]
            offset_score = abs(float(start_origin[0]))
            score = (unique_count, -offset_score)
            if best_score is None or score > best_score:
                best_score = score
                best_index = index
        return best_index

    def _derive_native_folded_elements(self, rays, folded_elements, max_half: float):
        if not folded_elements:
            return folded_elements
        ref_index = self._reference_native_ray_index(rays, max_half)
        if ref_index is None or ref_index >= len(rays.CC):
            return folded_elements
        cc = np.asarray(rays.CC[ref_index], dtype=float)
        surface_ids = [int(v) for v in np.asarray(rays.SURFACE[ref_index], dtype=int).ravel().tolist()]
        if cc.shape[0] < 2 or not surface_ids:
            return folded_elements

        unique_hits: list[tuple[int, np.ndarray]] = []
        last_id: int | None = None
        for hit_index, surface_index in enumerate(surface_ids, start=1):
            if surface_index == last_id:
                continue
            last_id = surface_index
            if hit_index >= len(cc):
                break
            unique_hits.append((surface_index, np.asarray(cc[hit_index], dtype=float)))
        if not unique_hits:
            return folded_elements

        new_elements = list(folded_elements)
        display_point = np.array([0.0, max(float(self.rows[0].thickness), 0.0)], dtype=float) if self.rows else np.array([0.0, 0.0], dtype=float)
        native_prev = np.asarray(cc[0], dtype=float)
        current_dir = np.array([0.0, 1.0], dtype=float)
        current_medium = 1.0

        for surface_index, native_hit in unique_hits:
            if surface_index <= 0 or surface_index > len(new_elements):
                native_prev = native_hit
                continue
            surface_type, _, row, _ = new_elements[surface_index - 1]
            native_step = float(np.linalg.norm(np.asarray(native_hit[1:3]) - np.asarray(native_prev[1:3])))
            center = display_point + current_dir * native_step
            new_elements[surface_index - 1] = (surface_type, center.copy(), row, current_dir.copy())
            display_point = center.copy()
            native_prev = native_hit

            if surface_type == "Mirror":
                current_dir = self._snap_display_direction(self._reflect_2d(current_dir, float(row.tilt_x)))
            elif surface_type == "Standard":
                axis = current_dir / max(np.linalg.norm(current_dir), 1e-12)
                sphere_center = center + axis * float(row.rc)
                normal = center - sphere_center
                next_medium = self._glass_index_for_preview(row.glass)
                current_dir = self._refract_ray_2d(current_dir, normal, current_medium, next_medium)
                current_medium = next_medium

        return new_elements

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
                    current_dir = self._snap_display_direction(
                        self._reflect_2d(current_dir, float(row.tilt_x))
                    )
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

    def _clear_cardinal_marker_artists(self) -> None:
        for artist in self._cardinal_marker_artists:
            try:
                artist.remove()
            except Exception:
                pass
        self._cardinal_marker_artists.clear()

    def _on_toggle_cardinal_markers(self) -> None:
        self._clear_cardinal_marker_artists()
        if not self.show_cardinals_var.get():
            self.canvas.draw_idle()
            self.status_var.set("PP / EP / XP hidden")
            return

        if self._last_optics_info is None and self.last_system is not None and self.last_rays is not None:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    self._last_optics_info = self._collect_optics_info(
                        self.last_system,
                        self.last_rays,
                        self._current_wavelength(),
                    )
            except Exception:
                self._last_optics_info = None

        if self._last_optics_info is None:
            self.canvas.draw_idle()
            self.status_var.set("PP / EP / XP unavailable for current view")
            return

        self._draw_optics_markers(self._last_optics_info)
        self.canvas.draw_idle()
        self.status_var.set("PP / EP / XP updated")
        self._autosave_plot()

    def _draw_optics_markers(self, optics_info: dict) -> None:
        self._clear_cardinal_marker_artists()
        if not self.show_cardinals_var.get():
            return
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()
        x_min, x_max = min(x0, x1), max(x0, x1)
        y_min, y_max = min(y0, y1), max(y0, y1)
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
                if y_mark < y_min or y_mark > y_max:
                    continue
                line = self.ax.axhline(y_mark, color=color, linewidth=1.0, linestyle=":", alpha=0.9)
                text = self.ax.text(
                    x0 + 0.04 * (x1 - x0),
                    y_mark,
                    label,
                    color=color,
                    fontsize=8,
                    ha="left",
                    va="bottom",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.6},
                )
                self._cardinal_marker_artists.extend((line, text))
            else:
                if z_val < x_min or z_val > x_max:
                    continue
                line = self.ax.axvline(z_val, color=color, linewidth=1.0, linestyle=":", alpha=0.9)
                text = self.ax.text(
                    z_val,
                    y_top,
                    label,
                    color=color,
                    fontsize=8,
                    ha="center",
                    va="top",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.6},
                )
                self._cardinal_marker_artists.extend((line, text))

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
        line = message.rstrip()
        self.debug_text.insert("end", line + "\n")
        self.debug_text.see("end")
        self._append_debug_log(line)
        self.update_idletasks()

    def _bind_text_copy_shortcuts(self, widget: tk.Text) -> None:
        for sequence in ("<Control-c>", "<Control-C>", "<Control-Insert>", "<<Copy>>", "<Control-KeyPress-c>", "<Control-KeyPress-C>"):
            widget.bind(sequence, lambda _e, w=widget: self._copy_selection_from_text_widget(w), add="+")

    def _bind_text_context_menu(self, widget: tk.Text) -> None:
        widget.bind("<Button-3>", lambda e, w=widget: self._show_text_context_menu(e, w), add="+")

    def _bind_global_copy_shortcuts(self) -> None:
        for sequence in ("<Control-c>", "<Control-C>", "<Control-Insert>"):
            self.bind_all(sequence, self._copy_selection_from_focus, add="+")

    def _show_text_context_menu(self, event, widget: tk.Text):
        if self._text_popup_menu is None:
            menu = tk.Menu(self, tearoff=0)
            menu.add_command(label="Copy Selected", command=lambda: self._copy_selection_from_text_widget(widget))
            menu.add_command(label="Copy All", command=lambda: self._copy_all_from_text_widget(widget))
            self._text_popup_menu = menu
        else:
            self._text_popup_menu.entryconfigure(0, command=lambda: self._copy_selection_from_text_widget(widget))
            self._text_popup_menu.entryconfigure(1, command=lambda: self._copy_all_from_text_widget(widget))
        self._text_popup_menu.tk_popup(event.x_root, event.y_root)
        return "break"

    def _copy_selection_from_focus(self, _event=None):
        candidates = []
        focused = self.focus_get()
        if isinstance(focused, tk.Text):
            candidates.append(focused)
        for widget in (getattr(self, "debug_text", None), getattr(self, "progress_text", None)):
            if isinstance(widget, tk.Text) and widget not in candidates:
                candidates.append(widget)
        for widget in candidates:
            try:
                text = widget.get("sel.first", "sel.last")
            except tk.TclError:
                continue
            if not text:
                continue
            try:
                ok, backend = self._copy_text_to_clipboard(text)
                if ok:
                    self.status_var.set(f"Selected text copied to clipboard ({backend})")
                else:
                    self.status_var.set("Copy failed")
                return "break"
            except Exception as exc:
                self.append_debug(f"Copy selected text failed: {exc}")
                return "break"
        return None

    def _copy_selection_from_text_widget(self, widget: tk.Text) -> str:
        try:
            text = widget.get("sel.first", "sel.last")
        except tk.TclError:
            self.status_var.set("No text selected")
            return "break"
        if not text:
            self.status_var.set("No text selected")
            return "break"
        try:
            ok, backend = self._copy_text_to_clipboard(text)
            if ok:
                self.status_var.set(f"Selected text copied to clipboard ({backend})")
            else:
                self.status_var.set("Copy failed")
        except Exception as exc:
            self.append_debug(f"Copy selected text failed: {exc}")
        return "break"

    def _copy_all_from_text_widget(self, widget: tk.Text) -> str:
        text = widget.get("1.0", "end-1c")
        if not text:
            self.status_var.set("No text to copy")
            return "break"
        try:
            ok, backend = self._copy_text_to_clipboard(text)
            if ok:
                self.status_var.set(f"All text copied to clipboard ({backend})")
            else:
                self.status_var.set("Copy failed")
        except Exception as exc:
            self.append_debug(f"Copy all text failed: {exc}")
        return "break"

    def _copy_text_to_clipboard(self, text: str) -> tuple[bool, str]:
        tools = (
            ("wl-copy", ["wl-copy"]),
            ("xclip", ["xclip", "-selection", "clipboard"]),
            ("xsel", ["xsel", "--clipboard", "--input"]),
        )
        encoded = text.encode("utf-8", errors="replace")
        for label, cmd in tools:
            if shutil.which(cmd[0]) is None:
                continue
            try:
                subprocess.run(cmd, input=encoded, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                return True, label
            except Exception:
                continue
        try:
            self.clipboard_clear()
            self.clipboard_append(text)
            self.update()
            return True, "Tk"
        except Exception:
            return False, "none"

    def copy_debug_to_clipboard(self) -> None:
        try:
            self._copy_all_from_text_widget(self.debug_text)
        except Exception as exc:
            self.append_debug(f"Copy debug failed: {exc}")

    def _reset_debug_log(self) -> None:
        try:
            DEBUG_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
            DEBUG_LOG_PATH.write_text("", encoding="utf-8")
        except Exception:
            pass

    def _append_debug_log(self, line: str) -> None:
        try:
            with DEBUG_LOG_PATH.open("a", encoding="utf-8") as fh:
                fh.write(line + "\n")
        except Exception:
            pass

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

    def _set_analysis_parallel_status(self, label: str, workers: int = 1, parallel_capable: bool = False) -> None:
        self._last_analysis_label = str(label)
        self._last_analysis_workers = max(1, int(workers))
        self._last_analysis_parallel_capable = bool(parallel_capable)
        self._last_analysis_accelerator = "CPU"

    def _set_analysis_accelerator(self, label: str) -> None:
        self._last_analysis_accelerator = str(label)

    def _analysis_parallel_summary(self) -> str:
        workers = max(1, int(self._last_analysis_workers))
        if workers <= 1:
            return "workers: 1"
        if self._last_analysis_parallel_capable:
            return f"workers: {workers} (parallel)"
        return f"workers: {workers}"

    def _analysis_compute_summary(self) -> str:
        return f"{self._analysis_parallel_summary()} | {self._last_analysis_accelerator}"

    def _report_compute_backends(self) -> None:
        if self._gpu_backend_reported:
            return
        self._gpu_backend_reported = True
        backend_pref = os.getenv("KRAKEN_POSTPROC_BACKEND", "auto").strip().lower()
        if backend_pref not in {"auto", "torch", "cupy", "cpu"}:
            backend_pref = "auto"
        self.append_debug(f"Post-processing backend preference: {backend_pref}")

        torch = _optional_torch()
        if torch is None:
            self.append_debug("Torch backend: unavailable.")
        else:
            try:
                if bool(torch.cuda.is_available()):
                    self.append_debug(f"Torch backend: CUDA available ({torch.cuda.device_count()} device(s)).")
                else:
                    self.append_debug("Torch backend: installed, CUDA not available.")
            except Exception as exc:
                self.append_debug(f"Torch backend: probe failed: {exc}")

        cp = _optional_cupy()
        if cp is None:
            self.append_debug("GPU backend: CuPy unavailable, PSF/MTF post-processing will use CPU.")
            return
        try:
            device_count = int(cp.cuda.runtime.getDeviceCount())
        except Exception:
            device_count = 0
        if device_count > 0:
            self.append_debug(f"GPU backend: CuPy available, detected {device_count} CUDA device(s).")
        else:
            self.append_debug("GPU backend: CuPy import succeeded, but no CUDA devices were detected.")

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

    def _serializable_row_specs(self) -> list[dict]:
        return [asdict(row) for row in self.rows]

    def _mtf_worker_count(self, ray_count: int) -> int:
        cpu_total = os.cpu_count() or 1
        if cpu_total <= 1 or ray_count < 2048:
            return 1
        return max(1, min(cpu_total - 1, ray_count // 2048))

    def _optimization_worker_count(self) -> int:
        cpu_total = max(1, int(os.cpu_count() or 1))
        if hasattr(self, "optimization_workers_var"):
            selected = self.optimization_workers_var.get().strip()
            if selected and selected.lower() != "auto":
                try:
                    parsed = int(selected)
                    if parsed > 0:
                        return max(1, min(parsed, cpu_total))
                except ValueError:
                    pass
        configured = os.getenv("KRAKEN_OPT_WORKERS", "").strip()
        if configured:
            try:
                parsed = int(configured)
                if parsed > 0:
                    return max(1, min(parsed, cpu_total))
            except ValueError:
                pass
        return 1 if cpu_total <= 1 else max(2, cpu_total - 1)

    def _ensure_analysis_executor(self, worker_count: int) -> ProcessPoolExecutor | None:
        worker_count = max(1, int(worker_count))
        if worker_count <= 1:
            return None
        if self._analysis_executor is not None and self._analysis_executor_workers == worker_count:
            return self._analysis_executor
        self._shutdown_analysis_executor()
        self._analysis_executor = ProcessPoolExecutor(max_workers=worker_count)
        self._analysis_executor_workers = worker_count
        return self._analysis_executor

    def _shutdown_analysis_executor(self) -> None:
        if self._analysis_executor is not None:
            self._analysis_executor.shutdown(wait=False, cancel_futures=True)
            self._analysis_executor = None
            self._analysis_executor_workers = 0

    @staticmethod
    def _terminate_process_group(pid: int, *, force: bool = False) -> None:
        if pid <= 0 or os.name != "posix":
            return
        signals = (signal.SIGTERM, signal.SIGKILL) if force else (signal.SIGTERM,)
        for sig in signals:
            try:
                os.killpg(pid, sig)
            except ProcessLookupError:
                return
            except Exception:
                break
            if sig == signal.SIGTERM:
                time.sleep(0.05)

    def _shutdown_optimization_worker(self, force: bool = False) -> None:
        if self._optimization_stop_event is not None:
            try:
                self._optimization_stop_event.set()
            except Exception:
                pass
        process = self._optimization_process
        if process is not None:
            try:
                process.join(timeout=0.15)
            except Exception:
                pass
            try:
                if force and process.is_alive():
                    self._terminate_process_group(int(process.pid or 0), force=True)
            except Exception:
                pass
            try:
                if process.is_alive():
                    process.join(timeout=0.15)
            except Exception:
                pass
            try:
                if not process.is_alive():
                    process.close()
            except Exception:
                pass
        queue = self._optimization_queue
        if queue is not None:
            try:
                queue.close()
            except Exception:
                pass
        self._optimization_process = None
        self._optimization_queue = None
        self._optimization_stop_event = None

    def _trace_pattern_chunks_parallel(
        self,
        wavelength: float,
        bundles: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    ) -> tuple[np.ndarray, np.ndarray, int]:
        total_rays = int(sum(len(np.asarray(bundle[0])) for bundle in bundles))
        worker_count = self._mtf_worker_count(total_rays)
        if worker_count <= 1:
            x_parts: list[np.ndarray] = []
            y_parts: list[np.ndarray] = []
            row_specs = self._serializable_row_specs()
            for bundle in bundles:
                x_local, y_local = _trace_analysis_chunk(row_specs, wavelength, *bundle)
                if x_local.size:
                    x_parts.append(x_local)
                    y_parts.append(y_local)
            if not x_parts:
                return np.asarray([], dtype=float), np.asarray([], dtype=float), 1
            return np.concatenate(x_parts), np.concatenate(y_parts), 1

        row_specs = self._serializable_row_specs()
        futures = []
        x_parts = []
        y_parts = []
        executor = self._ensure_analysis_executor(worker_count)
        if executor is None:
            return np.asarray([], dtype=float), np.asarray([], dtype=float), 1
        for bundle in bundles:
            ray_total = len(np.asarray(bundle[0]))
            if ray_total == 0:
                continue
            indices = np.array_split(np.arange(ray_total), worker_count)
            for chunk in indices:
                if chunk.size == 0:
                    continue
                chunk_bundle = tuple(np.asarray(values)[chunk] for values in bundle)
                futures.append(executor.submit(_trace_analysis_chunk, row_specs, wavelength, *chunk_bundle))
        for future in futures:
            x_local, y_local = future.result()
            if x_local.size:
                x_parts.append(x_local)
                y_parts.append(y_local)
        if not x_parts:
            return np.asarray([], dtype=float), np.asarray([], dtype=float), worker_count
        return np.concatenate(x_parts), np.concatenate(y_parts), worker_count

    def _trace_pattern_chunks_parallel_full(
        self,
        wavelength: float,
        bundles: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
        total_rays = int(sum(len(np.asarray(bundle[0])) for bundle in bundles))
        worker_count = self._mtf_worker_count(total_rays)
        if worker_count <= 1:
            parts = [[] for _ in range(6)]
            row_specs = self._serializable_row_specs()
            for bundle in bundles:
                outputs = _trace_analysis_chunk_full(row_specs, wavelength, *bundle)
                if outputs[0].size:
                    for idx, arr in enumerate(outputs):
                        parts[idx].append(arr)
            if not parts[0]:
                empty = np.asarray([], dtype=float)
                return empty, empty, empty, empty, empty, empty, 1
            merged = [np.concatenate(group) for group in parts]
            return (*merged, 1)

        row_specs = self._serializable_row_specs()
        futures = []
        parts = [[] for _ in range(6)]
        executor = self._ensure_analysis_executor(worker_count)
        if executor is None:
            empty = np.asarray([], dtype=float)
            return empty, empty, empty, empty, empty, empty, 1
        for bundle in bundles:
            ray_total = len(np.asarray(bundle[0]))
            if ray_total == 0:
                continue
            indices = np.array_split(np.arange(ray_total), worker_count)
            for chunk in indices:
                if chunk.size == 0:
                    continue
                chunk_bundle = tuple(np.asarray(values)[chunk] for values in bundle)
                futures.append(executor.submit(_trace_analysis_chunk_full, row_specs, wavelength, *chunk_bundle))
        for future in futures:
            outputs = future.result()
            if outputs[0].size:
                for idx, arr in enumerate(outputs):
                    parts[idx].append(arr)
        if not parts[0]:
            empty = np.asarray([], dtype=float)
            return empty, empty, empty, empty, empty, empty, worker_count
        merged = [np.concatenate(group) for group in parts]
        return (*merged, worker_count)

    def _compute_diffraction_mtf_sample(
        self,
        system,
        *,
        wavelength: float,
        surface_index: int,
        aperture_type: str,
        aperture_value: float,
        field_type: str,
        field_x: float,
        field_y: float,
    ) -> dict[str, object]:
        def _sanitize_phase_arrays(px, py, phase):
            px = np.asarray(px, dtype=float)
            py = np.asarray(py, dtype=float)
            phase = np.asarray(phase, dtype=float)
            finite = np.isfinite(px) & np.isfinite(py) & np.isfinite(phase)
            return px[finite], py[finite], phase[finite]

        def _phase_is_degenerate(px: np.ndarray, py: np.ndarray) -> bool:
            if px.size < 6:
                return True
            unique_x = np.unique(np.round(px, 10)).size
            unique_y = np.unique(np.round(py, 10)).size
            return unique_x <= 1 or unique_y <= 1

        pupil = Kos.PupilCalc(
            system,
            int(surface_index),
            float(wavelength),
            str(aperture_type),
            float(aperture_value),
        )
        pupil.Samp = max(11, min(21, self._current_ray_count()))
        pupil.Ptype = "hexapolar"
        pupil.FieldType = str(field_type)
        pupil.FieldX = float(field_x)
        pupil.FieldY = float(field_y)
        phase_method = "Phase"
        px, py, phase, _p2v = Kos.Phase(pupil)
        px, py, phase = _sanitize_phase_arrays(px, py, phase)
        if _phase_is_degenerate(px, py):
            capture = io.StringIO()
            with redirect_stdout(capture), redirect_stderr(capture):
                px, py, phase, _p2v = Kos.Phase2(pupil)
            phase2_log = capture.getvalue().strip()
            if phase2_log:
                self.append_debug(phase2_log)
            px, py, phase = _sanitize_phase_arrays(px, py, phase)
            phase_method = "Phase2"
        if _phase_is_degenerate(px, py):
            raise RuntimeError("Degenerate pupil sample from Phase/Phase2")
        if px.size < 6:
            raise RuntimeError("Not enough finite pupil samples for MTF fitting")
        phase_ptv = float(np.ptp(phase)) if phase.size else 0.0
        phase_peak = float(np.max(np.abs(phase))) if phase.size else 0.0
        if phase_ptv <= 1e-14 and phase_peak <= 1e-14:
            raise RuntimeError(f"Flat phase map from {phase_method}")

        zcoef = None
        used_terms = None
        last_error = None
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

        focal = max(0.01, abs(float(getattr(pupil, "EFFL", 0.0))))
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
        # Convert the FFT axis to image-plane cycles/mm using fc = D / (lambda * f).
        freq_max = diameter / max(wavelength * 1e-3 * focal, 1e-12)
        freq = np.linspace(0.0, freq_max, samples // 2)
        tangential = np.asarray(mtf[samples // 2, samples // 2:], dtype=float)
        sagittal = np.asarray(mtf[samples // 2 :, samples // 2], dtype=float)
        count = min(len(freq), len(tangential), len(sagittal))
        return {
            "plot_freq": np.asarray(freq[:count], dtype=float),
            "plot_tan": np.asarray(tangential[:count], dtype=float),
            "plot_sag": np.asarray(sagittal[:count], dtype=float),
            "method": "Diffraction",
            "worker_count": 1,
            "accelerator": "CPU",
            "sample_count": int(px.size),
            "used_terms": int(used_terms) if used_terms is not None else 0,
            "phase_method": phase_method,
        }

    def _compute_geometric_mtf_sample(
        self,
        system,
        *,
        wavelength: float,
        surface_index: int,
        aperture_type: str,
        aperture_value: float,
        field_type: str,
        field_x: float,
        field_y: float,
        algorithm: str = "psf_fft",
    ) -> dict[str, object]:
        dense_count = max(24, self._current_ray_count() * 6)
        x_local, y_local, worker_count = self._build_geometric_image_samples(
            system,
            wavelength,
            sample_count=dense_count,
            pattern="hexapolar",
            surface_index=int(surface_index),
            aperture_type=str(aperture_type),
            aperture_value=float(aperture_value),
            field_type=str(field_type),
            field_x=float(field_x),
            field_y=float(field_y),
        )
        if x_local.size < 4:
            raise RuntimeError("Not enough image-plane ray samples for geometric MTF")

        span_x = max(float(np.ptp(x_local)), 1e-3)
        span_y = max(float(np.ptp(y_local)), 1e-3)
        span = max(span_x, span_y) * 1.25
        if span <= 0:
            span = 1.0
        bins = 128
        if str(algorithm).strip().lower() == "lsf_fft":
            positive, tangential, sagittal, accelerator = self._compute_lsf_mtf_arrays(
                x_local,
                y_local,
                bins,
                span,
            )
            method_name = "Geometric-LSF"
        else:
            mtf, freq, _xedges, _unused, accelerator = self._compute_geometric_mtf_arrays(
                x_local,
                y_local,
                bins,
                span,
            )
            center = bins // 2
            positive = np.asarray(freq[center:], dtype=float)
            tangential = np.asarray(mtf[center, center:], dtype=float)
            sagittal = np.asarray(mtf[center:, center], dtype=float)
            method_name = "Geometric-PSF"
        count = min(len(positive), len(tangential), len(sagittal))
        return {
            "plot_freq": np.asarray(positive[:count], dtype=float),
            "plot_tan": np.asarray(tangential[:count], dtype=float),
            "plot_sag": np.asarray(sagittal[:count], dtype=float),
            "method": method_name,
            "worker_count": int(worker_count),
            "accelerator": str(accelerator),
            "sample_count": int(x_local.size),
            "pupil_samp": int(dense_count),
        }

    @staticmethod
    def _compute_lsf_mtf_arrays(
        x_local: np.ndarray,
        y_local: np.ndarray,
        bins: int,
        span: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
        lower = -span / 2.0
        upper = span / 2.0
        hist_x, xedges = np.histogram(x_local, bins=bins, range=(lower, upper))
        hist_y, _yedges = np.histogram(y_local, bins=bins, range=(lower, upper))
        lsf_tan = np.asarray(hist_x, dtype=float)
        lsf_sag = np.asarray(hist_y, dtype=float)
        if float(np.sum(lsf_tan)) <= 0.0 or float(np.sum(lsf_sag)) <= 0.0:
            raise RuntimeError("Not enough line-spread samples for LSF MTF")

        # Apply a smooth window before FFT to stabilize high-frequency ringing.
        window = np.hanning(bins)
        lsf_tan *= window
        lsf_sag *= window
        lsf_tan /= max(float(np.sum(lsf_tan)), 1e-12)
        lsf_sag /= max(float(np.sum(lsf_sag)), 1e-12)
        mtf_tan = np.abs(np.fft.fft(lsf_tan))
        mtf_sag = np.abs(np.fft.fft(lsf_sag))
        mtf_tan /= max(float(mtf_tan[0]), 1e-12)
        mtf_sag /= max(float(mtf_sag[0]), 1e-12)
        dx = float(xedges[1] - xedges[0])
        freq = np.fft.fftfreq(bins, d=dx)
        positive = freq >= 0.0
        return (
            np.asarray(freq[positive], dtype=float),
            np.asarray(mtf_tan[positive], dtype=float),
            np.asarray(mtf_sag[positive], dtype=float),
            "CPU",
        )

    def _build_geometric_image_samples(
        self,
        system,
        wavelength: float,
        sample_count: int,
        pattern: str = "hexapolar",
        *,
        surface_index: int,
        aperture_type: str,
        aperture_value: float,
        field_type: str,
        field_x: float,
        field_y: float,
    ) -> tuple[np.ndarray, np.ndarray, int]:
        pupil = Kos.PupilCalc(
            system,
            int(surface_index),
            wavelength,
            str(aperture_type),
            float(aperture_value),
        )
        pupil.Samp = max(2, int(sample_count))
        pupil.Ptype = str(pattern)
        pupil.FieldType = str(field_type)
        field_pairs = [(float(field_x), float(field_y))]
        bundles = []
        for fx, fy in field_pairs:
            pupil.FieldX = fx
            pupil.FieldY = fy
            bundles.append(tuple(np.asarray(values, dtype=float) for values in pupil.Pattern2Field()))
        x_local, y_local, worker_count = self._trace_pattern_chunks_parallel(wavelength, bundles)
        finite = np.isfinite(x_local) & np.isfinite(y_local)
        return x_local[finite], y_local[finite], worker_count

    def _build_geometric_image_samples_full(
        self,
        system,
        wavelength: float,
        sample_count: int,
        pattern: str = "hexapolar",
        *,
        surface_index: int,
        aperture_type: str,
        aperture_value: float,
        field_type: str,
        field_x: float,
        field_y: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
        pupil = Kos.PupilCalc(
            system,
            int(surface_index),
            wavelength,
            str(aperture_type),
            float(aperture_value),
        )
        pupil.Samp = max(2, int(sample_count))
        pupil.Ptype = str(pattern)
        pupil.FieldType = str(field_type)
        pupil.FieldX = float(field_x)
        pupil.FieldY = float(field_y)
        bundles = [tuple(np.asarray(values, dtype=float) for values in pupil.Pattern2Field())]
        x_local, y_local, z_local, l_local, m_local, n_local, worker_count = self._trace_pattern_chunks_parallel_full(
            wavelength,
            bundles,
        )
        finite = (
            np.isfinite(x_local)
            & np.isfinite(y_local)
            & np.isfinite(z_local)
            & np.isfinite(l_local)
            & np.isfinite(m_local)
            & np.isfinite(n_local)
        )
        return (
            x_local[finite],
            y_local[finite],
            z_local[finite],
            l_local[finite],
            m_local[finite],
            n_local[finite],
            worker_count,
        )

    def _compute_psf_histogram(
        self,
        x_local: np.ndarray,
        y_local: np.ndarray,
        bins: int,
        span: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
        backend_pref = os.getenv("KRAKEN_POSTPROC_BACKEND", "auto").strip().lower()
        gpu_min_samples = max(1, int(os.getenv("KRAKEN_GPU_MIN_SAMPLES", "1000000")))
        allow_auto_gpu = x_local.size >= gpu_min_samples
        if backend_pref == "auto" and not allow_auto_gpu:
            backend_pref = "cpu"
        if backend_pref in {"auto", "torch"}:
            torch = _optional_torch()
            if torch is not None:
                try:
                    if bool(torch.cuda.is_available()):
                        device = torch.device("cuda")
                        lower = -span / 2.0
                        upper = span / 2.0
                        step = (upper - lower) / float(bins)
                        x_t = torch.as_tensor(x_local, dtype=torch.float64, device=device)
                        y_t = torch.as_tensor(y_local, dtype=torch.float64, device=device)
                        ix = torch.floor((x_t - lower) / step).to(torch.int64)
                        iy = torch.floor((y_t - lower) / step).to(torch.int64)
                        valid = (ix >= 0) & (ix < bins) & (iy >= 0) & (iy < bins)
                        ix = ix[valid]
                        iy = iy[valid]
                        lin = ix * bins + iy
                        hist_t = torch.zeros(bins * bins, dtype=torch.float64, device=device)
                        hist_t.scatter_add_(0, lin, torch.ones_like(lin, dtype=torch.float64))
                        hist_t = hist_t.view(bins, bins)
                        xedges_t = torch.linspace(lower, upper, bins + 1, dtype=torch.float64, device=device)
                        yedges_t = torch.linspace(lower, upper, bins + 1, dtype=torch.float64, device=device)
                        return (
                            hist_t.detach().cpu().numpy(),
                            xedges_t.detach().cpu().numpy(),
                            yedges_t.detach().cpu().numpy(),
                            "GPU-Torch",
                        )
                except Exception:
                    pass
        if backend_pref in {"auto", "cupy"}:
            cp = _optional_cupy()
            if cp is not None:
                try:
                    x_gpu = cp.asarray(x_local, dtype=cp.float64)
                    y_gpu = cp.asarray(y_local, dtype=cp.float64)
                    hist_gpu, xedges_gpu, yedges_gpu = cp.histogram2d(
                        x_gpu,
                        y_gpu,
                        bins=bins,
                        range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
                    )
                    return (
                        cp.asnumpy(hist_gpu),
                        cp.asnumpy(xedges_gpu),
                        cp.asnumpy(yedges_gpu),
                        "GPU-CuPy",
                    )
                except Exception:
                    pass
        hist, xedges, yedges = np.histogram2d(
            x_local,
            y_local,
            bins=bins,
            range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
        )
        return hist, xedges, yedges, "CPU"

    def _compute_geometric_mtf_arrays(
        self,
        x_local: np.ndarray,
        y_local: np.ndarray,
        bins: int,
        span: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
        backend_pref = os.getenv("KRAKEN_POSTPROC_BACKEND", "auto").strip().lower()
        gpu_min_samples = max(1, int(os.getenv("KRAKEN_GPU_MIN_SAMPLES", "1000000")))
        allow_auto_gpu = x_local.size >= gpu_min_samples
        if backend_pref == "auto" and not allow_auto_gpu:
            backend_pref = "cpu"
        if backend_pref in {"auto", "torch"}:
            torch = _optional_torch()
            if torch is not None:
                try:
                    if bool(torch.cuda.is_available()):
                        device = torch.device("cuda")
                        lower = -span / 2.0
                        upper = span / 2.0
                        step = (upper - lower) / float(bins)
                        x_t = torch.as_tensor(x_local, dtype=torch.float64, device=device)
                        y_t = torch.as_tensor(y_local, dtype=torch.float64, device=device)
                        ix = torch.floor((x_t - lower) / step).to(torch.int64)
                        iy = torch.floor((y_t - lower) / step).to(torch.int64)
                        valid = (ix >= 0) & (ix < bins) & (iy >= 0) & (iy < bins)
                        ix = ix[valid]
                        iy = iy[valid]
                        lin = ix * bins + iy
                        hist_t = torch.zeros(bins * bins, dtype=torch.float64, device=device)
                        hist_t.scatter_add_(0, lin, torch.ones_like(lin, dtype=torch.float64))
                        hist_t = hist_t.view(bins, bins)
                        psf_t = hist_t / torch.clamp(torch.sum(hist_t), min=1.0)
                        otf_t = torch.fft.fftshift(torch.fft.fft2(psf_t))
                        mtf_t = torch.abs(otf_t)
                        mtf_t = mtf_t / torch.clamp(torch.max(mtf_t), min=1e-12)
                        xedges_t = torch.linspace(lower, upper, bins + 1, dtype=torch.float64, device=device)
                        dx = float((xedges_t[1] - xedges_t[0]).detach().cpu().item())
                        freq_t = torch.fft.fftshift(torch.fft.fftfreq(bins, d=dx, device=device))
                        return (
                            mtf_t.detach().cpu().numpy(),
                            freq_t.detach().cpu().numpy(),
                            xedges_t.detach().cpu().numpy(),
                            np.asarray([], dtype=float),
                            "GPU-Torch",
                        )
                except Exception:
                    pass
        if backend_pref in {"auto", "cupy"}:
            cp = _optional_cupy()
            if cp is not None:
                try:
                    x_gpu = cp.asarray(x_local, dtype=cp.float64)
                    y_gpu = cp.asarray(y_local, dtype=cp.float64)
                    hist_gpu, xedges_gpu, _yedges_gpu = cp.histogram2d(
                        x_gpu,
                        y_gpu,
                        bins=bins,
                        range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
                    )
                    psf_gpu = hist_gpu / cp.maximum(cp.sum(hist_gpu), 1.0)
                    otf_gpu = cp.fft.fftshift(cp.fft.fft2(psf_gpu))
                    mtf_gpu = cp.abs(otf_gpu)
                    mtf_gpu /= cp.maximum(cp.max(mtf_gpu), 1e-12)
                    dx = float(cp.asnumpy(xedges_gpu[1] - xedges_gpu[0]))
                    freq_gpu = cp.fft.fftshift(cp.fft.fftfreq(bins, d=dx))
                    return (
                        cp.asnumpy(mtf_gpu),
                        cp.asnumpy(freq_gpu),
                        cp.asnumpy(xedges_gpu),
                        np.asarray([], dtype=float),
                        "GPU-CuPy",
                    )
                except Exception:
                    pass
        hist, xedges, _yedges = np.histogram2d(
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
        return mtf, freq, xedges, np.asarray([], dtype=float), "CPU"

    def _collect_optics_info(self, system, rays, wavelength: float) -> dict:
        info: dict[str, float | None | str] = {
            "effl": None,
            "magnification": None,
            "ppa": None,
            "ppp": None,
            "paraxial_image_size": None,
            "sensor_fill": None,
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
            effl, ppa, ppp = self._exact_paraxial_cardinals(wavelength)
            info.update(
                {
                    "effl": float(effl),
                    "ppa": float(ppa),
                    "ppp": float(ppp),
                }
            )
            if self._current_object_mode() == "Finite" and len(self.rows) >= 2:
                object_gap = max(float(self.rows[0].thickness), 1e-9)
                object_size = max(float(self.rows[0].diameter), 0.0)
                sensor_size = max(float(self.rows[-1].diameter), 0.0)
                object_principal = object_gap + float(ppa)
                if np.isfinite(object_principal) and abs(object_principal) > 1e-12:
                    power_balance = (1.0 / float(effl)) - (1.0 / object_principal)
                    if abs(power_balance) > 1e-12:
                        image_principal = 1.0 / power_balance
                        magnification = image_principal / object_principal
                        image_size = abs(magnification) * object_size
                        info["paraxial_image_size"] = float(image_size)
                        info["magnification"] = float(magnification)
                        if sensor_size > 1e-12:
                            info["sensor_fill"] = float(image_size / sensor_size)
        except Exception:
            pass
        try:
            _, _, _, a, b, c, d, _effl, _ppa, _ppp, _, _, _ = system.Parax(wavelength)
            info["parax_a"] = float(a)
            info["parax_b"] = float(b)
            info["parax_c"] = float(c)
            info["parax_d"] = float(d)
            if info.get("magnification") is None:
                info["magnification"] = float(a)
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
        items.append(("Analysis mode", self._last_analysis_label))
        items.append(("Analysis workers", self._analysis_compute_summary()))
        items.append(("Analysis surface", str(self._analysis_surface_index())))
        items.append(("Aperture type", self._current_aperture_type()))
        items.append(("Aperture value", f"{self._current_aperture_value():.4g}"))
        field_metrics = self._field_metrics_summary()
        items.append(("Field type", self._current_field_type()))
        items.append(("Field angle [deg]", f"{field_metrics['current_angle_deg']:.4g}"))
        items.append(("Object height [mm]", f"{field_metrics['current_object_height']:.4g}"))
        items.append(("Paraxial img semi-ht [mm]", f"{field_metrics['current_paraxial_image_height']:.4g}"))
        items.append(("Real img semi-ht [mm]", f"{field_metrics['current_real_image_height']:.4g}"))
        if self._current_field_count() > 1:
            items.append(("Max parax img semi-ht [mm]", f"{field_metrics['max_paraxial_image_height']:.4g}"))
            items.append(("Max real img semi-ht [mm]", f"{field_metrics['max_real_image_height']:.4g}"))
            items.append(("Required image dia [mm]", f"{field_metrics['image_diameter']:.4g}"))
        if self.rows and self.rows[-1].surface == "Image":
            items.append(("Image row diameter [mm]", f"{float(self.rows[-1].diameter):.4g}"))

        total_length = sum(max(float(row.thickness), 0.0) for row in self.rows)
        items.append(("Total length [mm]", f"{total_length:.4g}"))

        if optics_info.get("effl") is not None:
            items.append(("Imaging", ""))
            items.append(("EFFL [mm]", f"{float(optics_info['effl']):.4g}"))
            items.append(("Magnification", f"{float(optics_info['magnification']):.4g}"))
            if optics_info.get("paraxial_image_size") is not None:
                items.append(("Paraxial image size [mm]", f"{float(optics_info['paraxial_image_size']):.4g}"))
            if optics_info.get("sensor_fill") is not None:
                items.append(("Sensor fill", f"{100.0 * float(optics_info['sensor_fill']):.3g}%"))
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
        if self._current_object_mode() == "Infinity":
            return "angle"
        return "height"

    def _resolved_field_coordinate(self, field_basis: str, raw_value: float, resolved_field_type: str) -> float:
        metrics = self._field_metrics_for_value(field_basis, raw_value)
        if resolved_field_type == "angle":
            return float(metrics["angle_deg"])
        return float(metrics["object_height"])

    def _resolved_mtf_field_samples(self, label: str) -> list[dict[str, float | str]]:
        field_basis = self._current_field_type()
        resolved_field_type = "angle" if self._current_object_mode() == "Infinity" else "height"
        basis_label = field_basis
        if field_basis == "Angle" and resolved_field_type == "height":
            basis_label = "Angle->Object Height"
        unit = self._field_type_unit(field_basis)
        raw_x = 0.0
        resolved_x = 0.0
        raw_values = self._sample_field_values(self._current_field_value())
        if not raw_values:
            raw_values = [self._current_field_value()]

        samples: list[dict[str, float | str]] = []
        for raw_value in raw_values:
            resolved_y = self._resolved_field_coordinate(field_basis, raw_value, resolved_field_type)
            samples.append(
                {
                    "basis": basis_label,
                    "unit": unit,
                    "display_x": float(raw_x),
                    "display_y": float(raw_value),
                    "field_type": resolved_field_type,
                    "field_x": float(resolved_x),
                    "field_y": float(resolved_y),
                    "legend": f"{self._format_field_sample_value(raw_value)} {unit}".strip(),
                }
            )
        return samples

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
            if row.surface == "Image":
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

        self.append_progress(
            "Optimization start | operands: "
            + ", ".join(spec.label for spec in merit_specs)
        )
        self.append_progress(f"Variables: {', '.join(v.normalized_name() for v in variables)}")
        population_size = 12
        optimization_workers = 1
        parallel_pref = os.getenv("KRAKEN_OPT_PARALLEL", "1").strip().lower()
        parallel_enabled = parallel_pref not in {"0", "false", "off", "no"}
        if parallel_enabled:
            optimization_workers = self._optimization_worker_count()
        self.status_var.set("Optimization starting...")
        self.append_progress("Preparing optimization worker...")
        self.optimization_running = True
        self.optimization_cancel_requested = False
        self.optimization_context = {
            "variables": variables,
            "champion_x": list(x0),
            "generations_total": 12,
            "generation_done": 0,
            "verbosity_every": 1,
            "workers": optimization_workers,
            "compute_backend": "pending",
            "initial_total": None,
            "last_best": None,
        }
        self._update_progress_indicators()
        self._start_progress_spinner()
        ctx = mp.get_context("spawn")
        self._optimization_queue = ctx.Queue()
        self._optimization_stop_event = ctx.Event()
        self._optimization_process = ctx.Process(
            target=_run_optimization_job,
            args=(
                self._optimization_queue,
                self._optimization_stop_event,
                self._serializable_row_specs(),
                merit,
                variables,
                list(map(float, x0)),
                int(self.optimization_context["generations_total"]),
                int(self.optimization_context["verbosity_every"]),
                int(population_size),
                int(optimization_workers),
                bool(parallel_enabled),
            ),
        )
        self._optimization_process.start()
        self.after(75, self._poll_optimization_worker)

    def stop_optimization(self) -> None:
        if not self.optimization_running:
            self.append_progress("Stop ignored: no optimization is running.")
            return
        self.optimization_cancel_requested = True
        if self._optimization_stop_event is not None:
            try:
                self._optimization_stop_event.set()
            except Exception:
                pass
        self.append_progress("Stop requested. Applying the latest completed generation and terminating the worker.")
        partial_result = None
        if self.optimization_context is not None:
            partial_result = {
                "champion_x": list(self.optimization_context.get("champion_x", [])),
                "initial_total": self.optimization_context.get("initial_total"),
                "final_total": self.optimization_context.get("last_best", self.optimization_context.get("initial_total")),
                "compute_backend": self.optimization_context.get("compute_backend", "pending"),
                "workers": self.optimization_context.get("workers", 1),
                "operands": [],
            }
        self._finish_optimization(cancelled=True, result=partial_result)

    def _poll_optimization_worker(self) -> None:
        if not self.optimization_running or self.optimization_context is None:
            return

        ctx = self.optimization_context
        queue = self._optimization_queue
        process = self._optimization_process
        if queue is None or process is None:
            self.append_progress("Optimization failed: worker process not available.")
            self._finish_optimization(cancelled=True)
            return

        completed = False
        while True:
            try:
                payload = queue.get_nowait()
            except Empty:
                break

            message_type = str(payload.get("type", ""))
            if message_type == "bootstrap":
                for line in payload.get("debug_messages", []) or []:
                    self.append_debug(str(line))
                initial_total = float(payload.get("initial_total", 0.0))
                ctx["initial_total"] = initial_total
                ctx["compute_backend"] = str(payload.get("compute_backend", "sequential"))
                ctx["workers"] = max(1, int(payload.get("workers", 1)))
                self.status_var.set(f"Optimization running: initial merit = {initial_total:.6g}")
                self.append_progress(f"Initial merit: {initial_total:.6g}")
                self.append_progress(f"Optimization compute: {ctx['compute_backend']}")
            elif message_type == "generation":
                capture_text = str(payload.get("debug", ""))
                if capture_text:
                    self.append_debug(capture_text)
                ctx["generation_done"] = max(ctx["generation_done"], int(payload.get("generation_done", 0)))
                champion_x = list(payload.get("champion_x", []))
                if champion_x:
                    ctx["champion_x"] = champion_x
                self._update_progress_indicators()
                if "log_best" in payload:
                    ctx["last_best"] = float(payload.get("log_best", 0.0))
                    if (
                        ctx["generation_done"] == 1
                        or ctx["generation_done"] == ctx["generations_total"]
                        or ctx["generation_done"] % ctx["verbosity_every"] == 0
                    ):
                        self.append_progress(
                            "Gen {gen:>3} | fevals {fevals:>4} | best {best:.6g} | dx {dx:.6g} | df {df:.6g}".format(
                                gen=int(ctx["generation_done"]),
                                fevals=int(payload.get("log_fevals", 0)),
                                best=float(payload.get("log_best", 0.0)),
                                dx=float(payload.get("log_dx", 0.0)),
                                df=float(payload.get("log_df", 0.0)),
                            )
                        )
                    if champion_x:
                        for variable, value in zip(ctx["variables"], champion_x):
                            row = self.rows[variable.surface_index]
                            if variable.parameter == "Rc":
                                row.rc = float(value)
                            elif variable.parameter == "Thickness":
                                row.thickness = float(value)
                        self._sync_table()
                        self.status_var.set(
                            f"Optimization running: generation {ctx['generation_done']}/{ctx['generations_total']} | best merit = {ctx['last_best']:.6g}"
                        )
            elif message_type == "complete":
                self._finish_optimization(cancelled=bool(payload.get("cancelled", False)), result=payload)
                completed = True
                break
            elif message_type == "error":
                tb = str(payload.get("traceback", ""))
                if tb:
                    self.append_debug(tb)
                self.append_progress(f"Optimization failed: {payload.get('message', 'unknown error')}")
                self._finish_optimization(cancelled=True)
                completed = True
                break

        if completed:
            return
        if process.is_alive():
            self.after(75, self._poll_optimization_worker)
            return
        self.append_progress("Optimization worker exited unexpectedly.")
        self._finish_optimization(cancelled=True)

    def _finish_optimization(self, cancelled: bool, result: dict | None = None) -> None:
        self._shutdown_optimization_worker(force=bool(cancelled))
        if self.optimization_context is None:
            self.optimization_running = False
            self.optimization_cancel_requested = False
            self._stop_progress_spinner()
            return

        ctx = self.optimization_context
        if result is not None:
            champion_x = list(result.get("champion_x", []))
            ctx["champion_x"] = champion_x
            ctx["compute_backend"] = str(result.get("compute_backend", ctx.get("compute_backend", "sequential")))
            ctx["workers"] = max(1, int(result.get("workers", ctx.get("workers", 1))))
            if result.get("initial_total") is not None:
                ctx["initial_total"] = float(result["initial_total"])

        champion_x = list(ctx.get("champion_x", []))
        for variable, value in zip(ctx["variables"], champion_x):
            row = self.rows[variable.surface_index]
            if variable.parameter == "Rc":
                row.rc = float(value)
            elif variable.parameter == "Thickness":
                row.thickness = float(value)

        self._sync_table()
        self.refresh_plot()
        initial_total = float(ctx.get("initial_total") or 0.0)
        final_total_value = initial_total
        if result is not None:
            candidate_final = result.get("final_total", initial_total)
            if candidate_final is not None:
                final_total_value = float(candidate_final)
        final_total = float(final_total_value)
        compute_backend = str(ctx.get("compute_backend", "sequential"))
        compute_workers = max(1, int(ctx.get("workers", 1)))
        if cancelled:
            self.status_var.set(
                f"Optimization stopped: {initial_total:.6g} -> {final_total:.6g}"
            )
            self.append_progress(
                f"Optimization stopped | merit {initial_total:.6g} -> {final_total:.6g}"
            )
        else:
            self.status_var.set(
                f"Optimization finished: {initial_total:.6g} -> {final_total:.6g}"
            )
            self.append_progress(
                f"Optimization finished | merit {initial_total:.6g} -> {final_total:.6g}"
            )
        self.append_progress(f"Optimization compute: {compute_backend} | workers={compute_workers}")
        for operand in (result.get("operands", []) if result is not None else []):
            self.append_progress(
                f"  {operand['name']}: value={float(operand['value']):.6g} weighted={float(operand['weighted']):.6g}"
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
        if count == 1:
            return [float(maximum)]
        span = abs(float(maximum))
        if span <= 1e-9:
            if self._current_object_mode() == "Finite" and self.rows:
                # Finite mode fallback: half object diameter.
                span = max(float(self.rows[0].diameter) * 0.5, 0.0)
            elif self._current_object_mode() == "Infinity":
                # Infinity mode fallback: default half-angle fan.
                span = 5.0
        if span <= 1e-9:
            return [0.0]
        return list(np.linspace(-span, span, count))

    def _current_image_diameter_mode(self) -> str:
        if not hasattr(self, "image_diameter_mode_var"):
            return "Auto"
        value = self.image_diameter_mode_var.get().strip()
        return value if value in {"Auto", "Manual"} else "Auto"

    def _auto_image_diameter_value(self) -> float:
        if not self.rows:
            return 3.0
        current_diameter = max(float(self.rows[-1].diameter), 1.0)
        sample_values = self._sample_field_values(self._current_field_value())
        if not sample_values:
            sample_values = [self._current_field_value()]
        image_heights = [
            abs(float(self._field_metrics_for_value(self._current_field_type(), value).get("real_image_height", 0.0)))
            for value in sample_values
        ]
        if not image_heights:
            return current_diameter
        max_height = max(image_heights)
        if max_height <= 1e-9:
            return current_diameter
        diameter = 2.0 * max_height
        return max(float(diameter), 1.0)

    def _apply_image_diameter_mode(self) -> None:
        if not self.rows or self.rows[-1].surface != "Image":
            return
        if self._current_image_diameter_mode() != "Auto":
            return
        self.rows[-1].diameter = self._auto_image_diameter_value()

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
        info: dict[str, object] = {"surfaces": [], "settings": {}}
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
        self._apply_layout_settings(info.get("settings", {}))
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
        settings = self._collect_layout_settings()
        settings_text = pformat(settings, sort_dicts=False, width=100)
        lines = [
            "#!/usr/bin/env python3",
            f'TITLE = "{title}"',
            "",
            f"SETTINGS = {settings_text}",
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
        self._apply_image_diameter_mode()

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
