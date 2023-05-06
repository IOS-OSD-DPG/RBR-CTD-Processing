# RBR-CTD-Processing
To process raw RBR CTD data in Python. The main reference for this routine is Halverson et al. (2017). 

RBR_Processing.py is deprecated and has been replaced by RBR_CTD_IOS.py.

RBR_CTD_IOS.py contains the following processing steps, which are accompanied by before/after plots.

- Export raw data from .rsk or .xlsx to .csv files

- Create metadata dictionary

- Prepare CTD_DATA.csv with 6 line headers for IOSShell

- Plot and check profile locations

- Plot first-order pressure differences to check need for zero-order hold

- (If needed) Correct for zero-order holds by replacing holds with NaNs

- (If needed) CALIB: Pressure/Depth correction if there are negative pressure values with corresponding conductivity values above about 30 mS/cm

- CLIP: Remove measurements near sea surface and bottom

- FILTER: Apply a low pass filter to pressure, temperature, conductivity, and fluorescence (if available)

- SHIFT: Shift conductivity and recalculate salinity

- SHIFT: Shift oxygen saturation

- DERIVE: Derive oxygen concentration in mL/L and umol/kg from oxygen saturation following SCOR WG 142

- DELETE: Remove pressure reversals

- BINAVE: Average the data into 1-dbar bins

- EDIT: Apply final editing, including converting conductivity units from mS/cm to S/m

- Prepare .ctd files following the IOS Header file ASCII format

##### Requirements
- Python >= 3.8
- pyrsktools >= 1.1.1
- Hakai ocean-data-parser package for converting oxygen saturation: https://github.com/HakaiInstitute/ocean-data-parser

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