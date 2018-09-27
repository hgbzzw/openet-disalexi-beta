import ee

# Use US-NE1 as default sample point
default_xy = (-96.47672812080845, 41.16506126041818)
default_crs = 'EPSG:32614'
default_geo = [30, 0, 632685, 0, -30, 4742715]
default_scale = 0.1


def constant_image_value(image, crs=default_crs):
    return image.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=ee.Geometry.Rectangle([0, 0, 10, 10], crs, False),
        scale=1).getInfo()


def image_value(image, xy=default_xy, scale=None,
                crs=default_crs, crsTransform=default_geo):
    # Default to using Landsat crsTransform unless scale is set
    if scale is not None:
        return image.reduceRegion(
            reducer=ee.Reducer.first(),
            geometry=ee.Geometry.Point(xy),
            scale=default_scale).getInfo()
    else:
        return image.reduceRegion(
            reducer=ee.Reducer.first(),
            geometry=ee.Geometry.Point(xy),
            crs=crs, crsTransform=crsTransform).getInfo()



def coll_value(coll, xy=default_xy, crs=default_crs,
               crsTransform=default_geo):
    # values = coll.getRegion(
    #     geometry=ee.Geometry.Point(xy),
    #     crs=crs, crsTransform=crsTransform).getInfo()
    values = coll.getRegion(
        geometry=ee.Geometry.Point(xy),
        scale=default_scale).getInfo()
    return [item for item in values]
