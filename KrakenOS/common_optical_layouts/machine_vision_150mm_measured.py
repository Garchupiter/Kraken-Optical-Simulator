TITLE = "Machine Vision 150 mm (Measured)"

SETTINGS = {
    "object_mode": "Finite",
    "display_orientation": "Vertical",
    "field_type": "Object Height",
    "field_value": 9.75,
    "field_count": 3,
    "analysis_surface": 2,
    "aperture_type": "STOP",
    "aperture_value": 26.8,
    "wavelength": 0.546,
}

SURFACES = [
    {
        "surface": "Object",
        "name": "Object",
        "rc": 0.0,
        "thickness": 268.0,
        "diameter": 19.5,
        "glass": "AIR",
    },
    {
        "surface": "Standard",
        "name": "Lens Front",
        "rc": 114.60971480633905,
        "thickness": 24.405,
        "diameter": 35.0,
        "glass": "BK7",
    },
    {
        "surface": "Standard",
        "name": "Stop Plane",
        "rc": 0.0,
        "thickness": 24.405,
        "diameter": 26.8,
        "glass": "BK7",
    },
    {
        "surface": "Standard",
        "name": "Lens Back",
        "rc": -211.2053946277084,
        "thickness": 308.19,
        "diameter": 35.0,
        "glass": "AIR",
    },
    {
        "surface": "Image",
        "name": "Image",
        "rc": 0.0,
        "thickness": 0.0,
        "diameter": 23.0,
        "glass": "AIR",
    },
]
