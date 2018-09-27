# Google Earth Engine DisALEXI

The GEE DisALEXI code was derived from the projectMAS DisALEXI version 0.4.0.

## Approach

The general approach is to map the GEE DisALEXI function over a filtered/merged/prepped Landsat collection (or any collection with thermal data).  Because the Landsat bands are different between Landsat 5, 7, and 8, the bands need to be renamed to common band names before merging.

## Usage

Before computing ET images, the GEE DisALEXI object needs to be instantiated.  The purpose of this step is to set the Alexi ET image collection and landcover image that will be used for all images.

```
disalexi = DisALEXI(
    et_coll=ee.ImageCollection('projects/climate-engine/alexi/conus/daily/et'),
    lc_img=ee.Image('USGS/NLCD/NLCD2011'), lc_type='NLCD')
```

Example call using the prep and ET functions directly on a single Landsat 8 collection 1 TOA image:
```
input_image = ee.Image(Landsat(ee.Image(
    'LANDSAT/LC08/C01/T1_TOA/LC08_043033_20150805')).prep_landsat8())
et_image = ee.Image(disalexi.compute_et(input_image))
```

Example call mapping the prep and ET functions over a Landsat 8 TOA collection (with only one image):
```
landsat_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA') \
    .filterDate('2016-08-05', '2016-08-06') \
    .filterMetadata('WRS_PATH', 'equals', 43) \
    .filterMetadata('WRS_ROW', 'equals', 33)
input_coll = ee.ImageCollection(landsat_coll.map(
    lambda x: ee.Image(landsat.Landsat(ee.Image(x)).prep_landsat8())))
et_coll = ee.ImageCollection(input_coll.map(disalexi.compute_et))
```

or
```
def apply_prep(img):
    return disalexi.Landsat(img).prep_landsat8()
input_coll = ee.ImageCollection(landsat_coll.map(apply_prep))
```

## Tests

The initial test values were generated using a slightly modified version of the [DisALEXI python code](https://gitlab.com/EE_pydisalexi/projectMAS).  See the [validation documentation](validation/VALIDATION.md) for additional details on how the test values were generated and the modifications that were made to the Python code.

The full test suite can be run using Pytest:
```
> python -m pytest
```

To get additional logging info, use the "-v" or "-s" tags:
```
py.test -v
```

To run tests for a single module:
```
python -m pytest tests\test_disalexi.py
```

## Python

See the [Python documentation](PYTHON.md) for additional details on installing Python and the external modules that are needed.
