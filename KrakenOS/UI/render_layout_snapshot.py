from __future__ import annotations

import argparse
from pathlib import Path

from KrakenOS.UI.layout_editor import AUTO_PLOT_PATH, KrakenLayoutEditor


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render a Kraken layout snapshot without opening the UI.")
    parser.add_argument("--mode", choices=["2d", "native", "mtf"], default="native", help="Render mode")
    parser.add_argument("--layout", default=None, help="Common layout title to load")
    parser.add_argument("--output", type=Path, default=AUTO_PLOT_PATH, help="Output image path")
    parser.add_argument("--dpi", type=int, default=180, help="Output DPI")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    app = KrakenLayoutEditor(headless=True)
    try:
        if args.layout:
            app.load_layout_by_name(args.layout)
        if args.mode == "native":
            app.analysis_mode = "native_off_axis"
        elif args.mode == "mtf":
            app.analysis_mode = "mtf"
        else:
            app.analysis_mode = "none"
        app.auto_save_plot_var.set(False)
        try:
            app.attributes("-alpha", 0.0)
        except Exception:
            pass
        app.geometry("1800x1100")
        app.update()
        app.refresh_plot()
        app.update()
        app.figure.set_size_inches(16, 9, forward=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        app.figure.savefig(args.output, dpi=args.dpi)
        print(args.output)
    finally:
        app.destroy()


if __name__ == "__main__":
    main()
