import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import uproot as up
from models import FitFunction


# noinspection DuplicatedCode
class MuonInstrument:
    """Class representing a muon instrument. Instance attributes are name of instrument, number of detectors, and how
    detectors are grouped.
    """

    emu_standard_groups = {
        "forward_outer": (0, 16),
        "forward_middle": (16, 32),
        "forward_inner": (32, 48),
        "backward_outer": (48, 64),
        "backward_middle": (64, 80),
        "backward_inner": (80, 96),
        "forward_total": (0, 48),
        "backward_total": (48, 96)
    }

    chronus_standard_groups = {
        "forward": (0, 303),
        "backward": (303, 606)
    }

    def __init__(self, name: str, num_detectors: int, detector_groups: dict):
        """Dictionary of detector groups has following format:
        {"detector_group_one": (start_slice, end_slice),
        "detector_group_two": ...  }

        For example, the first 16 detectors would be (0, 16).

        :param name: Name of muon instrument
        :param num_detectors: Total number of detectors for the instrument
        :param detector_groups: Dictionary of detector groups. Names are keys, values are slices.
        """
        self.num_detectors = num_detectors
        self.groups = detector_groups
        self.name = name
        self.model_dict = {}
        self.mac = None
        self.v1190 = None
        self.data = None
        self.num_events = None
        self.data_path = None

    def __repr__(self):
        return f"MuonInstrument({self.name})"

    def group_data(self):
        """Groups the detector histograms according to the instrument grouping.

        """
        dataframe = self.data.copy()
        for grouping in self.groups:
            start, end = self.groups[grouping]
            dataframe[grouping] = self.data.iloc[:, (start + 1):(end + 1)].sum(axis=1)
        dataframe.drop(list(range(1, self.num_detectors + 1)), axis=1, inplace=True)
        self.data = dataframe

    def load_data(self, mac: int, v1190: int, num_events: int, path: str = "../data"):
        """Loads data from ROOT histogram and parses into pandas dataframe with shape
        time detector1 detector2 ...
          .       .         .
          .       .         .
          .       .         .

        :param mac: mac file of simulation
        :param v1190: v1190 file of simulation
        :param num_events: Number of events simulated
        :param path: Path to ROOT file
        :return: Dataframe of detector counts
        """
        detectors = []
        with up.open(fr"{path}/his_{mac}_{v1190}.v1190.root") as f:  # opens ROOT histogram
            for detector_num in range(self.num_detectors):  # loops through all detector histograms
                detector_hist = f[f.keys()[2 + detector_num]]
                if detector_num == 0:
                    # gets the times from the first detector histogram
                    detectors.append(detector_hist.axis().centers())
                    detectors.append(detector_hist.values())
                else:
                    detectors.append(detector_hist.values())
        data = pd.DataFrame(np.stack(detectors, axis=1), columns=["time"] + list(range(1, self.num_detectors + 1)))
        self.data = data
        self.mac = mac
        self.v1190 = v1190
        self.num_events = num_events
        self.data_path = path

    def load_and_group_data(self, mac: int, v1190: int, num_events: int, path: str = "../data"):
        """Combines load_data and group_data into a single method.

        :param mac: mac file of simulation
        :param v1190: v1190 file of simulation
        :param num_events: Number of events simulated
        :param path: path to ROOT file
        :return: Dataframe of detector group counts
        """
        self.load_data(mac, v1190, num_events, path=path)
        self.group_data()

    def fit(self, name: str, function: str,
            start_time: float = 0, end_time: float = 15, initial_guess: list = None, bounds: tuple = None):
        """Fits curve to detector group histograms. Uses scipy's curve_fit method.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html

        :param name: Name of detector group/pair to fit.
        :param function: Function to fit to data
        :param start_time: Start time to regress from
        :param end_time: End time to regress from
        :param initial_guess: List of initial parameters for regression.
        :param bounds: Option to place bounds on parameters during fit
        """
        if hasattr(self, f"pair_{name}"):
            group_model = FitFunction(function)
            group_model.fit(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                            getattr(self, f"pair_{name}").asymmetry[0][(self.data["time"] >= start_time)
                                                                       & (self.data["time"] <= end_time)],
                            start_time=start_time, end_time=end_time, initial_guess=initial_guess, bounds=bounds)
        elif name in self.groups:
            group_model = FitFunction(function)
            group_model.fit(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                            self.data[name][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                            start_time=start_time, end_time=end_time, initial_guess=initial_guess, bounds=bounds)
        else:
            print("Enter valid group or pair")
            return
        self.model_dict[name] = group_model

    def plot_counts(self, group: str, plot_fit: bool,
                    start_time: float = 0, end_time: float = 15, save_path: str = None, show_plot: bool = False,
                    initial_guess: list = None, bounds: tuple = None):
        """Plots detector counts against time for a group. Option to plot fitted curve to data.

        :param group: Detector grouping to plot
        :param plot_fit: Fits and plots curve
        :param start_time: time to start plotting from
        :param end_time: time to end plotting at
        :param save_path: Optional directory to save figures
        :param show_plot: Option to output plot to terminal
        :param initial_guess: Option to change starting point for iteration. Defaults defined in models.py
        :param bounds: Option to place bounds on parameters.
        """
        if self.model_dict.get(group) is None:
            self.fit(group, "ExpDecayOsc", start_time=start_time, end_time=end_time, initial_guess=initial_guess,
                     bounds=bounds)
        plt.scatter(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                    self.data[group][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                    s=1, label=f"rate = {self.data[group].sum() / self.num_events}")
        if plot_fit:
            plt.plot(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                     self.model_dict[group].curve(
                         self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                         *self.model_dict[group].model[0]),
                     "r-", label=self.model_dict[group].graph_label)
        plt.title(f"{self.name}, {self.mac}.mac, {group}")
        plt.xlabel(r"time ($\mu$s)")
        plt.ylabel("n")
        plt.legend(loc="upper center", fontsize="x-small")
        if save_path is None:
            pass
        else:
            os.makedirs(f"{save_path}/{self.mac}", exist_ok=True)
            plt.savefig(f"{save_path}/{self.mac}/{self.mac}_{group}.png", dpi=200)
        if show_plot:
            plt.show()

    def plot_asymmetry(self, pair: str, plot_fit: bool,
                       start_time: float = 0, end_time: float = 15, save_path: str = None, show_plot: bool = False,
                       initial_guess: list = None, bounds: tuple = None):
        """Plots asymmetry for a pair of groups. Option to plot fitted curve to data.

        :param pair: Group pair to plot
        :param plot_fit: Fits and plots curve
        :param start_time: time to start plotting from
        :param end_time: time to end plotting at
        :param save_path: Optional directory to save figures
        :param show_plot: Option to output plot to terminal
        :param initial_guess: Option to change starting point for iteration. Defaults defined in models.py
        :param bounds: Option to place bounds on parameters.
        """
        if self.model_dict.get(pair) is None:
            self.fit(pair, "Sinusoid", start_time=start_time, end_time=end_time,
                     initial_guess=initial_guess, bounds=bounds)
        plt.scatter(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                    getattr(self, f"pair_{pair}").asymmetry[0][(self.data["time"] >= start_time)
                                                               & (self.data["time"] <= end_time)],
                    s=1, label=f"rate = {getattr(self, f'pair_{pair}').asymmetry.sum() / self.num_events}")
        if plot_fit:
            plt.plot(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                     self.model_dict[pair].curve(
                         self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                         *self.model_dict[pair].model[0]),
                     "r-", label=self.model_dict[pair].graph_label)
        plt.title(f"{self.name}, {self.mac}, pair_{pair} asymmetry")
        plt.xlabel(r"time ($\mu$s)")
        plt.ylabel("asymmetry")
        plt.legend(loc="upper center", fontsize="x-small")
        if save_path is None:
            pass
        else:
            os.makedirs(f"{save_path}/{self.mac}", exist_ok=True)
            plt.savefig(f"{save_path}/{self.mac}/{self.mac}_pair_{pair}_asymmetry.png", dpi=200)
        if show_plot:
            plt.show()

    def plot_histogram(self, hist_name: str, condition: int, show_plot: bool = True) -> tuple:
        with up.open(fr"{self.data_path}/his_{self.mac}_{self.v1190}.v1190.root")[f"{hist_name}_{condition};1"] as f:
            if f.axes.__len__() == 2:
                fig, ax = plt.subplots(dpi=200)
                w, x, y = f.to_numpy()
                mesh = ax.pcolormesh(x, y, w.T)
                fig.colorbar(mesh)
                ax.set_xlabel(f.axes[0].all_members["fTitle"])
                ax.set_ylabel(f.axes[1].all_members["fTitle"])
                ax.set_title(f"{f.title}, {f.name}, {self.mac}")
                if show_plot:
                    plt.show()
                return w, x, y
            else:
                if f.axis().labels() is None:
                    fig, ax = plt.subplots(dpi=200)
                    time, values = f.axis().centers(), f.values()
                    ax.errorbar(time, values, [np.sqrt(x) for x in values],
                                fmt="o", elinewidth=0.6, capsize=1, capthick=0.6, markersize=0.5)
                    ax.set_xlabel(f.axis().all_members["fTitle"])
                    ax.set_ylabel("N")
                    ax.set_title(f"{f.title}, {f.name}, {self.mac}")
                    if show_plot:
                        plt.show()
                    return time, values
                else:
                    fig, ax = plt.subplots(dpi=200)
                    labels, heights = f.axis().labels(), f.values()
                    ax.bar(labels, heights, yerr=[np.sqrt(x) for x in heights], capsize=4)
                    ax.set_ylabel("N")
                    ax.set_title(f'{f.title}, {f.name}, {self.mac}')
                    plt.xticks(rotation=-25, ha="left", size=7)
                    plt.tight_layout()
                    if show_plot:
                        plt.show()
                    return labels, heights

    def create_pair(self, pair_name: str, group_1: str, group_2: str):
        """Creates a pairing of detector groups. Asymmetry between the groups can then be calculated.

        :param pair_name: Name of pair
        :param group_1: First detector group in pair
        :param group_2: Second detector group in pair
        """
        if hasattr(self, f"pair_{pair_name}"):
            print("Pair already exists with this name!")
        else:
            setattr(self, f"pair_{pair_name}", Pair(self.data[group_1], self.data[group_2]))


class Pair:
    """Object representing a pairing of detector groups."""

    def __init__(self, group_1: pd.DataFrame, group_2: pd.DataFrame):
        """A pairing of detector groups.

        :param group_1: First detector group in the pair, known as "forward"
        :param group_2: Second detector group in the pair, known as "backward
        """
        self.group_1 = group_1
        self.group_2 = group_2
        self.alpha = None
        self.asymmetry = None

    def __repr__(self):
        return f"Pair({self.group_1}, {self.group_2}"

    def get_alpha(self) -> float:
        """Estimates alpha (balance parameter) between the two groups. Used in calculation for asymmetry.

        Defined as sum(forward) / sum(backward)
        """
        self.alpha = self.group_1.sum() / self.group_2.sum()
        return self.alpha

    def get_asymmetry(self, alpha: float = None) -> pd.DataFrame:
        """Calculates the asymmetry between the two groups.

        Defined as (forward - alpha * backward) / (forward + alpha * backward)

        :param alpha: The balance parameter used in the calculation. By default, it is estimated using get_alpha method.
        """
        if alpha is None:
            alpha = self.get_alpha()
        self.asymmetry = pd.DataFrame([0 if x == y == 0 else (x - alpha * y) / (x + alpha * y) for x, y in
                                       zip(self.group_1, self.group_2)])
        return self.asymmetry
