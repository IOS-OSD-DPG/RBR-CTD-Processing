import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ios_rbr_processing",
    version="0.0.1",
    author="Lu Guan, Samantha Huntington, Hana Hourston",
    author_email="DFO.PAC.SCI.IOSData-DonneesISO.SCI.PAC.MPO@dfo-mpo.gc.ca",
    description="A Python package for processing RBR CTD data and exporting in IOS header format.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IOS-OSD-DPG/RBR-CTD-Processing",
    install_requires=['cartopy', 'pyrsktools', 'numpy', 'pandas', 'scipy', 'datetime', 'gsw',
                      'matplotlib', 'scipy', 'ocean_data_parser', 'seawater', 'mpl_toolkits'],
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 3"],
    python_requires='>=3.8',
)