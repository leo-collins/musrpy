API
===

.. py:class:: musrpy.instruments.MuonInstrument(name, num_detectors, detector_groups, model_dict={}, mac=None, v1190=None, data=None, num_events=None, path=None, simana_histograms=None, sim_histograms=None)

    The MuonInstrument object represents a muon spin spectroscopy instrument. It can be loaded with data from
    musrSim musrSimAna. The object has several methods to load this data, group the detectors, and to visualise
    and analyse the data.

    :param name: Name of muon instrument
    :type name: str
    :param num_detectors: Total number of detectors of the instrument. Used when loading and grouping detectors.
    :type num_detectors: int
    :param detector_groups: How the detectors will be grouped.
    :type detector_groups: dict[str, list[int]] | dict[str, tuple[int, int]]
    :param model_dict: Optional. Where models fitted to data are stored.
    :type model_dict: dict[str, Fitfunction]
    :param mac: Optional. Mac file number. Used when loading the data.
    :type mac: int
    :param v1190: Optional. v1190 file number. Used when loading the data.
    :type v1190: int
    :param data: Optional. The raw or grouped detector data.
    :type data: pd.DataFrame
    :param num_events: Optional. Total number of events simulated. Defined in mac file via /run/beamOn command. Automatically obtained when loading data in.
    :type num_events: int
    :param path: Optional. Path to directory. See section :ref:`creatinginstrument` for details.
    :type path: str
    :param simana_histograms: Optional. List of histograms inside the .v1190.root tree. Obtained using get_histograms method.
    :type simana_histograms: list[str]
    :param sim_histograms: Optional. List of histograms inside the .root tree. Obtained using get_histograms method.
    :type sim_histograms: list[str]

    .. py:method:: load_data(mac, v1190, path, bins=2000)

        Loads data from ROOT histogram and parses into pandas dataframe with shape:
        ::

            time detector1 detector2 ...
              .       .         .
              .       .         .
              .       .         .

        Path must be to the folder containing the mac file and the data folder. See section :ref:`creatinginstrument` for details.
        Data can be rebinned using the bins parameter. Note that this must be *less* than the number of bins defined inside the v1190 file.

        :param mac: Mac file number of simulation to be loaded.
        :type mac: int
        :param v1190: v1190 file number of simulation to be loaded.
        :type v1190: int
        :param path: Path to directory containing mac file and data folder
        :type path: str
        :param bins: Number of bins to put data into.
        :type bins: int


    .. py:method:: group_data()

        Groups the detector histograms according to the detector_groups parameter. There must be data already loaded in using the :any:`load_data` method.
        This methods transforms the data parameter dataframe into the following structure:
        ::

            time    group1    group2 ...
              .       .         .
              .       .         .
              .       .         .


    .. py:method:: load_and_group_data(mac, v1190, path, bins=2000)

        A convenience method that combines the :any:`load_data` and :any:`group_data` methods into a single method. It first calls :any:`load_data`, then :any:`group_data`.

        :param mac: Mac file number of simulation to be loaded.
        :type mac: int
        :param v1190: v1190 file number of simulation to be loaded.
        :type v1190: int
        :param path: Path to directory containing mac file and data folder
        :type path: str
        :param bins: Optional. Number of bins to put data into.
        :type bins: int

    .. py:method:: fit(name, function, start_time=0, end_time=15, initial_guess=None, bounds=None)

        Fits curve to detector group histograms. Uses `scipy's curve_fit method <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_
        Used by xx and xx methods. The resulting model will then be saved inside the model_dict attribute with key equal to the name parameter.

        :param name: Name of detector group/pair to fit.
        :type name: str
        :param function: Function to fit to data. Can be one of xxx found in xxx.
        :type function: str
        :param start_time: Optional. Start time to regress from.
        :type start_time: float
        :param end_time: Optional. End time to regress to.
        :type end_time: float
        :param initial_guess: Optional. List of initial parameters for regression. See curve_fit documentation for details. Default is obtained from xxx
        :type initial_guess: list[float]
        :param bounds: Optional. Bounds can be placed on the parameters. See curve_fit documentation for details.
        :type bounds: tuple[list[float], list[float]]

    .. py:method:: plot_counts(group, plot_fit, start_time=0, end_time=15, save_path=None, show_plot=None, initial_guess=None, bounds=None)

        Plots detector counts against time for a group. Option to plot fitted curve to data. Can plot just one group or multiple groups on one plot.
        Plots can be saved into a chosen directory.

        :param group: Detector grouping or list of groups to plot.
        :type group: str | list[str]
        :param plot_fit: If true then the :any:`fit` method is called for the group(s) and is shown on the plot.
        :type plot_fit: bool
        :param start_time: Optional. Start time to regress from.
        :type start_time: float
        :param end_time: Optional. End time to regress to.
        :type end_time: float
        :param save_path: Optional. Will create a folder inside the specified directory and saves plot.
        :type save_path: str
        :param show_plot: Optional. If true then the plot will be shown when the method is called.
        :type show_plot: bool
        :param initial_guess: Optional. List of initial parameters for regression. See curve_fit documentation for details. Default is obtained from xxx
        :type initial_guess: list[float]
        :param bounds: Optional. Bounds can be placed on the parameters. See curve_fit documentation for details.
        :type bounds: tuple[list[float], list[float]]

    .. py:method:: plot_asymmetry(pair, plot_fit, start_time=0, end_time=15, save_path=None, show_plot=None, initial_guess=None, bounds=None)

        Plots asymmetry for a pair of groups. Option to plot fitted curve to data. Pairs can be created using the create_pair method.

        :param pair: Group pair to plot.
        :type pair: str
        :param plot_fit: If true then the :any:`fit` method is called for the pair and is shown on the plot.
        :type plot_fit: bool
        :param start_time: Optional. Start time to regress from.
        :type start_time: float
        :param end_time: Optional. End time to regress to.
        :type end_time: float
        :param save_path: Optional. Will create a folder inside the specified directory and saves plot.
        :type save_path: str
        :param show_plot: Optional. If true then the plot will be shown when the method is called.
        :type show_plot: bool
        :param initial_guess: Optional. List of initial parameters for regression. See curve_fit documentation for details. Default is obtained from xxx
        :type initial_guess: list[float]
        :param bounds: Optional. Bounds can be placed on the parameters. See curve_fit documentation for details.
        :type bounds: tuple[list[float], list[float]]

    .. py:method:: get_histograms()

        Returns lists of histograms contained in the .root and .v1190.root trees. These can then be plotted using the :any:`plot_histogram` method.
        These lists are stored inside the simana_histograms and sim_histograms attributes.

        :return: Lists of histograms
        :rtype: tuple[list[str], list[str]]

    .. py:method:: plot_histogram(hist_name, bins=None, data_range=None, show_plot=None, save_path=None)

        Method to plot other histograms stored inside the .root and .v1190.root trees. Can plot 2D histograms, 1D histograms and bar charts.

        :param hist_name: Name of histogram to plot. Names of all histograms can be found using the :any:`get_histograms` method.
        :type hist_name: str
        :param bins: Optional. Can rebin data. Default number of bins is determined from the data.
        :type bins: int
        :param data_range: Optional. Range of data to plot. Format is (start, end) for 1D histograms, and ((start_x, end_x), (start_y, end_y)) for 2D histograms. Default is the full range of data in the tree.
        :type data_range: tuple[float, float] | tuple[tuple[float, float], tuple[float, float]]
        :param show_plot: Optional. If true then the plot will be shown when the method is called.
        :type show_plot: bool
        :param save_path: Optional. Will create a folder inside the specified directory and saves plot.
        :type save_path: str

    .. py:method:: create_pair(pair_name, group_1, group_2)

        Creates a pair object from the two given groups. This is saved as a new attribute of the :any:`MuonInstrument`