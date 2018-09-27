# GEE DisALEXI Conversion Status

## Global Default

Once the global ALEXI ET images are available in GEE, switch the code to default to the global ALEXI ET and land cover collections.

## LAI

Currently the code is computing LAI based on the METRIC NDVI/LAI empirical relationship.  Nick is looking into implementing a GEE approach for computing Landsat scale LAI similar to what is being done in the projectMASlai repository.  The main .

## LST Sharpening

Currently the code is using the LST band directly from Landsat with no additional sharpening.

## Global ALEXI ET

Follow up with Chris about getting global ALEXI ET images or having them operationally ingest the images into GEE.

## CONUS ALEXI ET

Ingest additional CONUS ALEXI ET images or follow up with Chris about having them operationally ingest the images into GEE

## MERRA2

Ingest additional MERRA2 daily and hourly images into GEE.
Look into interpolating the hourly Rs images to the Landsat image collection time and/or selecting a different daily image.  Currently the code is selecting the image with the same UTC hour and date as the Landsat image.

## CFSR

Look into interpolating CFSR wind speed or aggregating over a different time period (instead of the UTC day)

## Testing

Additional test cases should be written to validate the corner cases and extremes of the intermediate functions.
