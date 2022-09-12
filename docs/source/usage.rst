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

.. _creatinginstrument:

Creating a Muon instrument
--------------------------

First we have to define our muon instrument as follows

.. code:: python

    CHRONUS = MuonInstrument("CHRONUS", num_detectors=606, detector_groups=group_dict)

where the "group_dict" object is a dictionary with the following structures:

::

    {"roup_one": (start_slice, end_slice),
     "group_two": ...  }

or

::

    {"group_one": [detector_ids],
     "group_two": ...  }

The name of the group is the key, and the value can either be a list of detector ids, or a start slice and end slice of detector ids. Note that detectors are indexed from 0 in the order they are defined inside the .v1190 file.
The MuonInstrument object comes with some standard groupings for the CHRONUS and EMU instruments at ISIS. See the API for details.

We can then load our data in using the .load_data method:

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

By default musrSim/musrSimAna place the root trees inside the "data" folder. After the data has been loaded we can use the :any:`group_data` method,
which will group all the detector histograms as defined in the detector_groups dictionary. For convenience we can use the :any:`load_and_group_data` method which
will combine these two methods into one.

Creating Pairs
--------------

We can create pairs of groups and find the asymmetry between them. To do this we use the :any:`create_pair` method as follows:

.. code:: python

    CHRONUS.create_pair("one", "group_one", "group_two")

This will create a new attribute of the CHRONUS object which can be accessed by

.. code:: python

    CHRONUS.pair_one

From this we can calculate the asymmetry using the :any:`get_asymmetry` method.

Creating Plots
--------------

Curve Fitting
-------------


