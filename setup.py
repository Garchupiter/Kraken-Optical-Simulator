from setuptools import setup

setup(
    name='KrakenOS',
    version='0.0.32',
    packages=['KrakenOS', 'KrakenOS.AstroAtmosphere', 'KrakenOS.Cat', 'KrakenOS.Examples', 'KrakenOS.Docs'],

    install_requires=['pyvista','PyVTK','vtk','numpy','scipy','matplotlib', 'csv342'],
    package_data={'': ['LICENSE.txt', '../*.AGF', '../*.agf', '../*.stl', '../*.scad', '../*.pdf']},
    include_package_data=True,
    url='https://github.com/Garchupiter/Kraken-Optical-Simulator',
    license='GNU General Public License v3.0',
    author='joel Herrera et al.',
    author_email='joel@astro.unam.mx',
    description='Optical Simulation and ray tracing '
)
