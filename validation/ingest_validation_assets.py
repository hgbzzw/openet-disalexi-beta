from builtins import input
import datetime
import logging
import os
import shutil
import subprocess


def main():
    """Process DisAlexi test datasets

    LC80280312014189LGN00 -> LC08_028031_20140708
    """
    input_id = 'LC08_028031_20140708'
    input_ws = os.path.join(
        os.getcwd(), 'IDL_Dataset', 'landsat_LC8_2014189', 'INPUT')
    output_ws = os.path.join(
        os.getcwd(), 'IDL_Dataset', 'landsat_LC8_2014189', 'OUTPUT')

    # Hardcoding image time for now
    image_time = '1705'
    image_dt = datetime.datetime.strptime(
        '{} {}{}'.format(input_id[12:20], image_time, '50'), '%Y%m%d %H%M%S')
    image_date = image_dt.strftime('%Y%j')
    image_date_start = image_dt.strftime('%Y-%m-%d')
    image_time_start = image_dt.strftime('%Y-%m-%dT%H:%M:%S')

    image_id = 'LC08_028031_20140708'
    # image_id = '{}0{}_{}_{}'.format(
    #     input_id[:2], input_id[2], input_id[3:9],
    #     image_dt.strftime('%Y%m%d'))

    adjust_ws = os.path.join(os.getcwd(), '{}_64bit'.format(image_id))
    if not os.path.isdir(adjust_ws):
        os.makedirs(adjust_ws)
    bucket_ws = 'gs://{}/{}'.format('disalexi', image_id)
    asset_folder = '{}/{}'.format('disalexi', image_id)
    asset_ws = 'users/cgmorton/{}'.format(asset_folder, image_id)

    # Set as strings for subprocess calls
    nodata_value = '-9999'
    # input_nodata = '-9999'
    # output_nodata = '-340282306073709653000000000000000000000'

    alexi_ws = os.path.join(input_ws, 'ALEXI')
    landsat_ws = os.path.join(input_ws, 'Surface')
    met_ws = os.path.join(input_ws, 'Meteo')
    aux_ws = os.path.join(input_ws, 'auxiliar')
    etd_ws = os.path.join(output_ws, 'Daily', 'Evapotranspiration')
    ta_ws = os.path.join(output_ws, 'TA')

    images = {
        'albedo': {
            'input_path': os.path.join(
                landsat_ws, 'Albedo', 'Albedo_{}'.format(image_date)),
            'adjust_path': os.path.join(adjust_ws, 'albedo.tif'),
            'bucket_path': '{}/albedo.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'albedo'),
            'time_start': image_time_start},
        'lai': {
            'input_path': os.path.join(
                landsat_ws, 'LAI', 'LAI_{}'.format(image_date)),
            'adjust_path': os.path.join(adjust_ws, 'lai.tif'),
            'bucket_path': '{}/lai.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'lai'),
            'time_start': image_time_start},
        'lst': {
            'input_path': os.path.join(
                landsat_ws, 'LST', 'LST_{}{}'.format(image_date, image_time)),
            'adjust_path': os.path.join(adjust_ws, 'lst.tif'),
            'bucket_path': '{}/lst.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'lst'),
            'time_start': image_time_start,
            # To match IDL codes exactly, this should be 273.16
            'calc': 'A+273.15',
            'nodata': '-9999'},
        'ndvi': {
            'input_path': os.path.join(
                landsat_ws, 'NDVI', 'NDVI_{}'.format(image_date)),
            'adjust_path': os.path.join(adjust_ws, 'ndvi.tif'),
            'bucket_path': '{}/ndvi.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'ndvi'),
            'time_start': image_time_start},
        # Compute CFMask based on nodata pixels in albedo
        'cfmask': {
            'input_path': os.path.join(
                landsat_ws, 'Albedo', 'Albedo_{}'.format(image_date)),
            'adjust_path': os.path.join(adjust_ws, 'cfmask.tif'),
            'bucket_path': '{}/cfmask.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'cfmask'),
            'type': "--type=UInt16",
            'nodata': '65535',
            'calc': 'A/-9999.0',
            'time_start': image_time_start},
        'landcover': {
            'input_path': os.path.join(aux_ws, 'landcover'),
            'adjust_path': os.path.join(adjust_ws, 'landcover.tif'),
            'bucket_path': '{}/landcover.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'landcover'),
            'type': "--type=UInt16",
            'calc': "A*1",
            'nodata': '255',
            'time_start': image_date_start},
        'Insol1': {
            'input_path': os.path.join(
                met_ws, 'Solar_Radiation',
                'Rs_{}{}'.format(image_date, image_time)),
            'adjust_path': os.path.join(adjust_ws, 'Insol1.tif'),
            'bucket_path': '{}/Insol1.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'Insol1'),
            'calc': "A*1",
            'time_start': image_time_start},
        # For now, just convert hourly Rs to TIF
        'Insol24': {
            'timestep': 'daily',
            'input_path': os.path.join(
                met_ws, 'Daily', 'Solar_Radiation',
                image_dt.strftime('%j'), 'Rs_{}'.format(image_date)),
            'adjust_path': os.path.join(
                adjust_ws, 'Rs', 'Rs_{}.tif'.format(image_date)),
            'bucket_path': '{}/Insol24.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'Insol24'),
            'time_start': image_date_start},
        'windspeed': {
            'input_path': os.path.join(
                met_ws, 'Wind_Speed', 'U_{}{}'.format(image_date, image_time)),
            'adjust_path': os.path.join(adjust_ws, 'u.tif'),
            'bucket_path': '{}/u.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'u'),
            'time_start': image_date_start},
        'Ta': {
            'input_path': os.path.join(
                ta_ws, 'Ta_{}{}'.format(image_date, image_time)),
            'adjust_path': os.path.join(adjust_ws, 'Ta.tif'),
            'bucket_path': '{}/Ta.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'Ta'),
            'calc': 'A+273.15',
            'time_start': image_date_start},
        'ETd': {
            'input_path': os.path.join(etd_ws, 'ETd_{}'.format(image_date)),
            'adjust_path': os.path.join(adjust_ws, 'ETd.tif'),
            'bucket_path': '{}/ETd.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'ETd'),
            'time_start': image_time_start},
        'ETalexi': {
            'input_path': os.path.join(alexi_ws, 'ETD_{}'.format(image_date)),
            'adjust_path': os.path.join(adjust_ws, 'ETalexi.tif'),
            'bucket_path': '{}/alexiET.tif'.format(bucket_ws),
            'asset_id': '{}/{}'.format(asset_ws, 'alexiET'),
            'time_start': image_date_start},
        # # 'mask': {
        # #     # Note - bad doy in filename
        # #     'input_path': os.path.join(aux_ws, 'mask'),
        # #     'adjust_path': os.path.join(adjust_ws, 'mask.tif'),
        # #     'bucket_path': '{}/mask.tif'.format(bucket_ws),
        # #     'asset_id': '{}/{}'.format(asset_ws, 'mask'),
        # #     'time_start': image_time_start,
        # #     'type': "--type=UInt16"},
        # # # Converting pressure from X to mb
        # # 'pressure': {
        # #     'input_path': os.path.join(met_ws, '{}_pSub.tiff'.format(input_id)),
        # #     'adjust_path': os.path.join(adjust_ws, 'p.tif'),
        # #     'bucket_path': '{}/p.tif'.format(bucket_ws),
        # #     'asset_id': '{}/{}'.format(asset_ws, 'p'),
        # #     'calc': 'A*0.1',
        # #     'time_start': image_date_start},
        # # 'q2': {
        # #     'input_path': os.path.join(met_ws, '{}_q2Sub.tiff'.format(input_id)),
        # #     'adjust_path': os.path.join(adjust_ws, 'q2.tif'),
        # #     'bucket_path': '{}/q2.tif'.format(bucket_ws),
        # #     'asset_id': '{}/{}'.format(asset_ws, 'q2'),
        # #     'calc': "A*1",
        # #     'time_start': image_date_start},
    }


    # Get a list of assets in the folder
    # Assume disalexi folder was created
    test_image_id_list = subprocess.check_output(
        ['earthengine', 'ls', os.path.dirname(asset_ws)],
        universal_newlines=True, shell=True)
    test_image_id_list = [
        x.strip() for x in test_image_id_list.split('\n') if x]
    if image_id not in test_image_id_list:
        subprocess.check_output([
            'earthengine', 'create', 'folder', asset_ws])


    # Get existing asset ID list
    asset_id_list = subprocess.check_output(
        ['earthengine', 'ls', asset_ws],
        universal_newlines=True, shell=True)
    asset_id_list = [
        x.strip() for x in asset_id_list.split('\n') if x]
    logging.info('Existing assets: {}'.format(', '.join(asset_id_list)))


    # Process each band
    for band_name, image in images.items():
        logging.info('Band: {}'.format(band_name))

        # Set defaults
        if 'input_path' not in image.keys():
            image['input_path'] = None
        if 'nodata' not in image.keys():
            image['nodata'] = nodata_value
        if 'type' not in image.keys():
            image['type'] = "--type=Float64"
        if 'calc' not in image.keys():
            image['calc'] = None

        # Make a copy of the file in a separate directory
        if 'timestep' in image.keys() and image['timestep'] == 'daily':
            # Process each hour
            for i in range(24):
                input_path = image['input_path'] + '{:02d}00'.format(i)
                adjust_path = image['adjust_path'].replace(
                    '.tif', '{:02d}00.tif'.format(i))
                if not os.path.isdir(os.path.dirname(adjust_path)):
                    os.makedirs(os.path.dirname(adjust_path))

                # Use gdal_translate to copy in case inputs are not TIF
                logging.info('  Copying: {} -> {}'.format(
                    os.path.basename(input_path),
                    os.path.basename(adjust_path)))
                if not os.path.isfile(input_path):
                    logging.info(' Image {} doesn\'t exist, skipping'.format(
                        input_path))
                    continue
                subprocess.check_output(
                    ['gdalwarp', '-overwrite', '-ot', 'Float64', input_path, adjust_path],
                    shell=True)
            # For now, don't upload to GCS or ingest to GEE
            continue
        elif image['input_path'] and not os.path.isfile(image['input_path']):
            logging.info('  The input file doesn\'t exist, skipping')
            continue
        elif 'timestep' not in image.keys() or image['timestep'] != 'daily':
            # GDAL translate won't overwrite so remove first
            if os.path.isfile(image['adjust_path']):
                logging.info('  Removing: {}'.format(
                    os.path.basename(image['adjust_path'])))
                os.remove(image['adjust_path'])

            # Use gdal_translate to copy in case inputs are not TIF
            logging.info('  Copying: {} -> {}'.format(
                os.path.basename(image['input_path']),
                os.path.basename(image['adjust_path'])))
            subprocess.check_output(
                ['gdal_translate', image['input_path'], image['adjust_path']],
                shell=True)


        # Set nodata value
        if band_name in ['albedo', 'ndvi', 'lai', 'lst', 'ETd', 'Ta']:
            logging.info('  Setting nodata')
            temp_path = image['adjust_path'].replace('.tif', '_temp.tif')
            logging.info('    {}'.format(image['adjust_path']))
            logging.info('    {}'.format(temp_path))
            # GDAL translate won't overwrite so remove first
            if os.path.isfile(temp_path):
                logging.info('  Removing: {}'.format(
                    os.path.basename(temp_path)))
                os.remove(temp_path)

            subprocess.check_output([
                'gdal_translate', '-a_nodata', image['nodata'],
                '-ot', 'Float64', image['adjust_path'], temp_path])
            # subprocess.check_output([
            #     'gdalwarp', '-srcnodata', '-9999', '-dstnodata', nodata_value,
            #     image['adjust_path'], temp_path])
            os.remove(image['adjust_path'])
            shutil.move(temp_path, image['adjust_path'])


        # Scale values
        if image['calc']:
            logging.info('  Adjust values')
            temp_path = image['adjust_path'].replace('.tif', '_temp.tif')
            logging.info('    {}'.format(image['adjust_path']))
            logging.info('    {}'.format(temp_path))
            # GDAL translate won't overwrite so remove first
            if os.path.isfile(temp_path):
                logging.info('  Removing: {}'.format(
                    os.path.basename(temp_path)))
                os.remove(temp_path)

            subprocess.check_output(
                [
                    'gdal_calc',
                    '--calc={}'.format(image['calc']),
                    '-A', image['adjust_path'],
                    '--outfile={}'.format(temp_path),
                    '{}'.format(image['type']),
                    #'--co="COMPRESS=LZW"'
                    '--NoDataValue={}'.format(image['nodata']),
                ],
                shell=True)
            os.remove(image['adjust_path'])
            shutil.move(temp_path, image['adjust_path'])


        # # Remove existing images from bucket?
        # image_list = subprocess.check_output(
        #     ['gsutil', 'ls', 'gs://' + bucket_name],
        #     universal_newlines=True, shell=True)
        # image_list = [x.strip() for x in image_list.split('\n') if x]
        # image_list = sorted([
        #     image_path for image_path in image_list
        #     if image_re.match(os.path.basename(image_path))])


        # Copy to the storage bucket
        logging.info('  Uploading: {}'.format(band_name))
        logging.info('    {}'.format(image['adjust_path']))
        logging.info('    {}'.format(image['bucket_path']))
        subprocess.check_output(
            ['gsutil', 'cp', image['adjust_path'], image['bucket_path']],
            shell=True)


        # Remove existing assets
        if image['asset_id'] in asset_id_list:
            logging.info('  Removing: {}'.format(image['asset_id']))
            subprocess.call(['earthengine', 'rm', image['asset_id']])


        # Ingest into Earth Engine
        logging.info('  Ingesting: {}'.format(image['asset_id']))
        subprocess.call(
            [
                'earthengine', 'upload', 'image',
                '--bands', band_name,
                '--asset_id', image['asset_id'],
                '--time_start', image['time_start'],
                '--nodata_value', image['nodata'],
                '--property', '(string)SCENE_ID={}'.format(image_id),
                # '--property', '(string)DOY={}'.format(start_dt.strftime('%j')),
                image['bucket_path']
            ])


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    main()
