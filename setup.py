from setuptools import setup

setup(
    name='KrakenOS',
    version='1.0.0.24',
    packages=[
        'KrakenOS',
        'KrakenOS.AstroAtmosphere',
        'KrakenOS.Cat',
        'KrakenOS.Examples',
        'KrakenOS.Docs',
        'KrakenOS.LensCat',
        'KrakenOS.Optimization',
        'KrakenOS.Optimization.adapters',
        'KrakenOS.Optimization.examples',
        'KrakenOS.UI',
        'KrakenOS.common_optical_layouts',
    ],

    install_requires=['pyvista','PyVTK','vtk','numpy','scipy','matplotlib', 'csv342', 'pandas' ],
    package_data={'': ['LICENSE.txt', '../*.AGF', '../*.agf', '../*.stl', '../*.scad', '../*.pdf', '../*.zmx', '../*.ZMX', '../*.csv']},
    include_package_data=True,
    url='https://github.com/Garchupiter/Kraken-Optical-Simulator',
    download_url='https://github.com/Garchupiter/Kraken-Optical-Simulator',
    keywords=['Optical', 'Simulator', 'lens'],
    license='GNU General Public License v3.0',
    author='joel Herrera et al.',
    author_email='joel@astro.unam.mx',
    description='Optical Simulation and ray tracing '
)
