# Python

## Dependencies
The EE DisALEXI functions have been tested using Python 2.7 (using the  "future" module, see below).

The following module must be present to run the EE DisALEXI function:
* [earthengine-api](https://github.com/google/earthengine-api)

The following module is used to run the test suite
* [pytest](http://doc.pytest.org/en/latest/)

The following modules must be present if using Python 2.7
* [future](http://python-future.org/)

## Anaconda/Miniconda

The easiest way of obtaining Python and all of the necessary external modules, is to install [Miniconda](https://www.continuum.io/downloads).

It is important to double check that you are calling the Anaconda version, especially if you have two or more version of Python installed (e.g. Anaconda and ArcGIS).

+ Windows: "where python"
+ Linux/Mac: "which python"

#### Installing/Updating Python Modules

Additional python modules will need to be installed (and/or updated) using "conda".  For example to install the pytest and future modules, enter the following in a command prompt or terminal window:

```
conda install pytest future
```

To update the pytest and future modules to the latest version, enter the following in a command prompt or terminal window:

```
conda update pytest future
```

External modules can also be installed/updated one at at a time:
```
> conda install pytest
> conda install future
```

#### EarthEngine-API / PIP

The EarthEngine API must be installed through pip (not conda):
```
> pip install earthengine-api
```

After installing the EarthEngine API module, you will need to authenticate the Earth Engine API (see [setting-up-authentication-credentials](https://developers.google.com/earth-engine/python_install#setting-up-authentication-credentials)):
```
> python -c "import ee; ee.Initialize()"
```
