{ pkgs, lib, ... }:

let
  runtimeLibs = with pkgs; [
    stdenv.cc.cc.lib
    zlib
    glib
    libGL
    fontconfig
    freetype
    dbus
    wayland
    libxkbcommon
    xorg.libX11
    xorg.libXcursor
    xorg.libXext
    xorg.libXfixes
    xorg.libXi
    xorg.libXrandr
    xorg.libXtst
    xorg.libXrender
    xorg.libxcb
  ];
in
{
  env = {
    GREET = "KrakenOS devenv";
    LD_LIBRARY_PATH = lib.makeLibraryPath runtimeLibs;
  };

  languages.python = {
    enable = true;
    package = pkgs.python312;
    uv.enable = true;
  };

  packages = with pkgs; [
    git
    pkg-config
  ] ++ runtimeLibs;

  enterShell = ''
    VENV_DIR="$PWD/.devenv/state/venv"
    REQ_HASH_FILE="$PWD/.devenv/state/kraken-requirements.hash"
    REQ_HASH="krakenos-v4"

    if [ ! -x "$VENV_DIR/bin/python" ]; then
      python -m venv "$VENV_DIR"
    fi

    if [ ! -f "$REQ_HASH_FILE" ] || [ "$(cat "$REQ_HASH_FILE" 2>/dev/null)" != "$REQ_HASH" ]; then
      "$VENV_DIR/bin/python" -m pip install --upgrade pip setuptools wheel
      "$VENV_DIR/bin/python" -m pip install \
        -e . \
        numpy scipy matplotlib pandas pyvista vtk \
        PyVTK csv342 ipython ipykernel jupyter jupyterlab pyzmq \
        packaging setuptools basedpyright ruff PyQt5 sip
      printf '%s\n' "$REQ_HASH" > "$REQ_HASH_FILE"
    fi

    if [ -n "''${WAYLAND_DISPLAY:-}" ] || [ -n "''${DISPLAY:-}" ]; then
      export MPLBACKEND=qtagg
    fi

    echo "$GREET"
    "$VENV_DIR/bin/python" --version
  '';
}
