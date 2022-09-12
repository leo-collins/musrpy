Usage
=====

.. _installation:

Installation
------------

To use musrpy, first install using pip:

.. code:: console

    $ pip install musrpy

Then import

.. code:: python

    from musrpy.instruments import MuonInstrument

Creating a Muon instrument
--------------------------

First we have to define our muon instrument as follows

.. code:: python

    CHRONUS = MuonInstrument("CHRONUS", num_detectors=606, detector_groups=MuonInstrument.chronus_standard_groups)

Then we can load our data in using the .load_data method:

.. code:: python

    CHRONUS.load_data(1, 1, path="directory", bins=200)

Note that the directory must be structured as follows:

::

    directory
    ├── data
    │   ├── his_1_1.v1190.root
    │   └── musr_1.root
    ├── 1.mac
    ├── 1.v1190
    ├── g4_1.wrl

By default musrSim/musrSimAna place the root trees inside the "data" folder.