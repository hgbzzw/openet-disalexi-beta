import base64
import os

import ee

# On Travis-CI authenticate using private key environment variable
if 'EE_PRIVATE_KEY_B64' in os.environ:
    print('Writing privatekey.json from environmental variable ...')
    content = base64.b64decode(os.environ['EE_PRIVATE_KEY_B64']).decode('ascii')
    EE_PRIVATE_KEY_FILE = 'privatekey.json'
    with open(EE_PRIVATE_KEY_FILE, 'w') as f:
        f.write(content)
    EE_CREDENTIALS = ee.ServiceAccountCredentials(
        '', key_file=EE_PRIVATE_KEY_FILE)
    ee.Initialize(EE_CREDENTIALS)
else:
    ee.Initialize()
