# GEE DisALEXI Input Data

The following is a brief summary of the data sets that are needed to compute Landsat DisALEXI ET.
The GEE DisALEXI class is is intended to be applied in a function that is mapped over a "prepped" image collection (see EXAMPLE).  The prepped image collection must at minimum have the following 5 bands ('albedo', 'cfmask', 'lai', 'lst', 'ndvi') and the following 2 properties ('SUN_ELEVATION' 'system:time_start').  The intended use of this code is for computing ET from Landsat images, but it should be possible to use it with other image collections with the appropriate data/bands.

## Landsat

ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA')
ee.ImageCollection('LANDSAT/LE07/C01/T1_RT_TOA')
ee.ImageCollection('LANDSAT/LC08/C01/T1_RT_TOA')

A custom Landsat class (disalexi/landsat.py) was built in order to "prep" the Landsat top-of-atmosphere reflectance images for use in the GEE DisALEXI code.  T

## ALEXI Evapotranspiration (ET)

Daily ALEXI ET data is not currently available in GEE as a native asset.  The CONUS ALEXI ET for 2014 were ingested into GEE.  The global ALEXI ET images are still being actively developed and were not available for ingestion into GEE.  A single global image for the validation date (2014-05-21) was also ingested, but this asset will need to be replaced once the global ALEXI ET data is regenerated.

If an ET collection is not set when the DisALEXI class is instantiated, the code will default to the CONUS ALEXI et collection.  The default will likely be changed to "global" once the global ALEXI ET images are

#### Global
Collection ID: "projects/climate-engine/alexi/global/daily/et"
Units: mm
Band: "ET"
Example: https://code.earthengine.google.com/16b4c9387497663fb60c920d23f9c09b

#### CONUS
Collection ID: "projects/climate-engine/alexi/conus/daily/et"
Units: mm
Band: "ET"
Example: https://code.earthengine.google.com/c3e3623202cd257325ddb814f2eccf0b

## Land Cover

The code currently supports two land cover datasets, GlobeLand30 and the National Land Cover Database (NLCD).  Both of these datasets have a nominal resolution of 30m.  If a land cover image and type are not defined when the DisALEXI class is instantiated, the code will default to the CONUS NLCD 2011 image.

In the python DisALEXI code, the land cover values were remapped to various properties (aleafv, adeadv, leaf width, etc) based on the values "landcover.xlsx" file.  For the GEE DisALEXI code, this data is now being imported directly into the code from the file "disalexi/landcover.py", where the values are stored in nested dictionaries by land cover type, property, and class.

To change the land cover from the CONUS default, the user will need to set both the lc_img and lc_type properties on the DisALEXI class object.  The "lc_type" property is needed for the lookup into the dictionaries in the landcover.py file.
```
disalexi = disalexi.DisALEXI(prep_img)
disalexi.lc_img = ee.Image('USGS/NLCD/NLCD2011')
disalexi.lc_type = 'NLCD'
```

#### GlobeLand30

[GlobeLand30](http://www.globallandcover.com/GLC30Download/index.aspx) is a global land cover dataset that is not currently available in Earth Engine as a native asset.  Images for 2010 for the globe were ingested into GEE as personal assets.

The GlobeLand30 collection should be mosaiced to a single image before use because the images in the collection overlap and have different projections.  The GlobeLand30 values (in the landcover.xlsx file) were only defined for the 10 main landuse classes (10, 20, 30, etc.) whereas the GlobeLand30 images have additional stratification.  In the GEE DisALEXI code, the images are rounded down to the closest even multiple of ten (e.g. 33 -> 30 and 11 -> 10) (see the following example).

```
disalexi = disalexi.Image(prep_img)
disalexi.lc_img = ee.Image(
    ee.ImageCollection('users/nbearson/GlobeLand30').mosaic()) \
        .divide(10).floor().multiply(10)
disalexi.lc_type = 'GLOBELAND30'
```

Collection ID: "users/cgmorton/GlobeLand30"
lc_type keyword: "GlobeLand30"
Band: "b1"
Example: https://code.earthengine.google.com/8d65a42248d813eaf076fcbbb928d57e

#### National Land Cover Database (NLCD)

[NLCD](https://www.mrlc.gov/nlcd2011.php) is avaiable in Earth Engine as a native asset.

Image ID: "USGS/NLCD/NLCD2011"
lc_type keyword: "NLCD"

## MERRA2

MERRA2 daily and hourly solar insolation values are not currently available in GEE as native assets.  Images for 2014 were ingested into GEE and it would be relatively easy to ingest images for any date in the MERRA2 archive (1980-01-01 to 2-6 weeks prior to present).

For simplicity, the MERRA2 image with the same UTC date and hour as the Landsat is used.  The drawback to this approach is that the image may not be the one closest to the Landsat image collection time.  This approach should be revisited before the code is used operationally.

Collection ID: "projects/climate-engine/merra2/rs_daily"
Collection ID: "projects/climate-engine/merra2/rs_hourly"
Units: W m-2 (total per hour and per day respectively)
Band: "rs"
Example: https://code.earthengine.google.com/5ac23430894de5487ddf2496005ab0f0

## CFSR

Wind speed is being computed from the CFSR 6 hour wind speed vector images.  For simplicity, the code is computing a daily average wind speed from the 6 hour data by selecting all images in the same UTC day as the Landsat image (even if the CFSR images are well after the Landsat image).  For example, the Landsat image LC08_028031_20140521 was captured at ~17:05 UTC.  The CFSR image times that would be used in the daily average would be 2014-05-21 00:00, 06:00, 12:00, & 18:00. (https://code.earthengine.google.com/043ed92d80fd641fabc8a5d1ac2fa5cb).  This approach should be revisited before the code is used operationally.

Collection ID: "NOAA/CFSV2/FOR6H"
Bands: "u-component_of_wind_height_above_ground" and "v-component_of_wind_height_above_ground"

## Elevation

Elevation is currently only being used to compute air pressure.  There are numerous elevation images available in GEE that could be used, but the code currently defaults to using the void filled SRTM 30m dataset (https://www2.jpl.nasa.gov/srtm/).

Image ID: "USGS/SRTMGL1_003"
Units: m
