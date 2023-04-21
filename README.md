# RBR-CTD-Processing
To process raw RBR CTD data in Python. 

RBR_Processing.py is deprecated and has been replaced by RBR_CTD_IOS.py.

RBR_CTD_IOS.py contains processing steps:

- Export raw data from .rsk to .csv files

- Create Metadata dictionary

- Prepare CTD_DATA.csv with 6 line headers for IOSShell

- Plot and check profile locations

- Plot and Check for zero-order hold

- CALIB: Pressure/Depth correction

- CLIP: remove measurements near sea surface and bottom

- FILTER: apply a low pass filter

- SHIFT: shift conductivity and recalculate salinity

- SHIFT: shift oxygen

- DELETE: remove the pressure reversal

- BINAVE: calculate bin averages

- EDIT: apply final editing

- Prepare .ctd files with IOS Header File

##### Requirements
- Python >= 3.8
- pyrsktools >= 1.1.1
- Hakai ocean-data-parser package for converting oxygen saturation:\
`pip install git+https://github.com/HakaiInstitute/ocean-data-parser.git`

##### References
Halverson, M., Jackson, J., Richards, C., Melling, H., Brunsting, R., 
Dempsey, M., Gatien, G., Hamilton, A., Hunt, B., Jacob, W., and 
Zimmerman, S. 2017. Guidelines for processing RBR CTD profiles. Can. 
Tech. Rep. Hydrogr. Ocean Sci. 314: iv + 38 p. 

Bittig Henry, Körtzinger Arne, Johnson Ken, Claustre Herve, Emerson 
Steve, Fennel Katja, Garcia Hernan, Gilbert Denis, Gruber Nicolas, 
Kang Dong-Jin, Naqvi Wajih, Prakash Satya, Riser Steven, Thierry 
Virginie, Tilbrook Bronte, Uchida Hiroshi, Ulloa Osvaldo, Xing 
Xiaogang (2018). **SCOR WG 142: Quality Control Procedures for Oxygen 
and Other Biogeochemical Sensors on Floats and Gliders. Recommendations 
on the conversion between oxygen quantities for Bio-Argo floats and 
other autonomous sensor platforms.** https://doi.org/10.13155/45915