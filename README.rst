=================
OpenET - DisALEXI
=================

|version| |build|

This repository provides an Earth Engine Python API based implementation of the DisALEXI ET model.

Examples
========

Jupyter notebooks are provided in the "examples" folder that show various approaches for calling the OpenET SSEBop model.

Installation
============

The OpenET DisALEXI python module can be installed via pip:

.. code-block:: console

    pip install openet-disalexi

Dependencies
============

Modules needed to run the model:

 * `earthengine-api <https://github.com/google/earthengine-api>`__
 * `openet <https://github.com/Open-ET/openet-core-beta>`__

Modules needed to run the test suite:

 * `pytest <https://docs.pytest.org/en/latest/>`__

Running Tests
=============

.. code-block:: console

    python -m pytest

OpenET Namespace Package
========================

Each OpenET model should be stored in the "openet" folder (namespace).  The benefit of the namespace package is that each ET model can be tracked in separate repositories but called as a "dot" submodule of the main openet module,

.. code-block:: console

    import openet.disalexi as disalexi


References
==========

 * `Senay et al., 2013 <http://onlinelibrary.wiley.com/doi/10.1111/jawr.12057/abstract>`__
 * `Senay et al., 2016 <http://www.sciencedirect.com/science/article/pii/S0034425715302650>`__

.. |build| image:: https://travis-ci.org/Open-ET/openet-disalexi-beta.svg?branch=master
   :alt: Build status
   :target: https://travis-ci.org/Open-ET/openet-disalexi-beta
.. |version| image:: https://badge.fury.io/py/openet-disalexi.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/openet-disalexi