TITLE = "Machine Vision 150 mm (Datasheet 1X)"

SETTINGS = {
    "object_mode": "Finite",
    "display_orientation": "Vertical",
    "field_type": "Object Height",
    "field_value": 16.5,
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
        "thickness": 275.0,
        "diameter": 33.0,
        "glass": "AIR",
    },
    {
        "surface": "Standard",
        "name": "Lens Front",
        "rc": 289.54325068870634,
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
        "rc": -99.78246835443024,
        "thickness": 291.3844936708861,
        "diameter": 35.0,
        "glass": "AIR",
    },
    {
        "surface": "Image",
        "name": "Image",
        "rc": 0.0,
        "thickness": 0.0,
        "diameter": 33.0,
        "glass": "AIR",
    },
]
