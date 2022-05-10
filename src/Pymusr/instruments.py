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

    chronus_tf_groups = {
        "trans1": [23, 24, 25, 26, 55, 56, 57, 58, 93, 94, 95, 96, 97, 136, 137, 138, 139, 140, 141, 184, 185, 186, 187,
                   188, 189, 237, 238, 239, 240, 241, 242, 293, 294, 295, 296, 297, 298, 299, 326, 327, 328, 329, 358,
                   359, 360, 361, 396, 397, 398, 399, 400, 439, 440, 441, 442, 443, 444, 487, 488, 489, 490, 491, 492,
                   540, 541, 542, 543, 544, 545, 596, 597, 598, 599, 600, 601, 602],
        "trans2": [27, 28, 29, 30, 59, 60, 61, 62, 63, 98, 99, 100, 101, 102, 142, 143, 144, 145, 146, 147, 190, 191,
                   192, 193, 194, 195, 196, 243, 244, 245, 246, 247, 248, 249, 300, 301, 302, 0, 1, 330, 331, 332, 333,
                   362, 363, 364, 365, 366, 401, 402, 403, 404, 405, 445, 446, 447, 448, 449, 450, 493, 494, 495, 496,
                   497, 498, 499, 546, 547, 548, 549, 550, 551, 552, 603, 604, 605, 303, 304],
        "trans3": [2, 3, 4, 5, 31, 32, 33, 34, 64, 65, 66, 67, 68, 103, 104, 105, 106, 107, 108, 148, 149, 150, 151,
                   152, 153, 197, 198, 199, 200, 201, 202, 250, 251, 252, 253, 254, 255, 256, 305, 306, 307, 308, 334,
                   335, 336, 337, 367, 368, 369, 370, 371, 406, 407, 408, 409, 410, 411, 451, 452, 453, 454, 455, 456,
                   500, 501, 502, 503, 504, 505, 553, 554, 555, 556, 557, 558, 559],
        "trans4": [6, 7, 8, 35, 36, 37, 38, 69, 70, 71, 72, 73, 109, 110, 111, 112, 113, 154, 155, 156, 157, 158, 159,
                   203, 204, 205, 206, 207, 208, 209, 257, 258, 259, 260, 261, 262, 263, 309, 310, 311, 338, 339, 340,
                   341, 372, 373, 374, 375, 376, 412, 413, 414, 415, 416, 457, 458, 459, 460, 461, 462, 506, 507, 508,
                   509, 510, 511, 512, 560, 561, 562, 563, 564, 565, 566],
        "trans5": [9, 10, 11, 12, 39, 40, 41, 42, 74, 75, 76, 77, 78, 114, 115, 116, 117, 118, 119, 160, 161, 162, 163,
                   164, 165, 210, 211, 212, 213, 214, 215, 216, 264, 265, 266, 267, 268, 269, 270, 312, 313, 314, 315,
                   342, 343, 344, 345, 377, 378, 379, 380, 381, 417, 418, 419, 420, 421, 422, 463, 464, 465, 466, 467,
                   468, 513, 514, 515, 516, 517, 518, 519, 567, 568, 569, 570, 571, 572, 573],
        "trans6": [13, 14, 15, 43, 44, 45, 46, 79, 80, 81, 82, 120, 121, 122, 123, 124, 166, 167, 168, 169, 170, 171,
                   217, 218, 219, 220, 221, 222, 271, 272, 273, 274, 275, 276, 277, 278, 316, 317, 318, 346, 347, 348,
                   349, 382, 383, 384, 385, 423, 424, 425, 426, 427, 469, 470, 471, 472, 473, 474, 520, 521, 522, 523,
                   524, 525, 574, 575, 576, 577, 578, 579, 580, 581],
        "trans7": [16, 17, 18, 19, 47, 48, 49, 50, 83, 84, 85, 86, 87, 125, 126, 127, 128, 129, 130, 172, 173, 174, 175,
                   176, 177, 223, 224, 225, 226, 227, 228, 229, 279, 280, 281, 282, 283, 284, 285, 319, 320, 321, 322,
                   350, 351, 352, 353, 386, 387, 388, 389, 390, 428, 429, 430, 431, 432, 433, 475, 476, 477, 478, 479,
                   480, 526, 527, 528, 529, 530, 531, 532, 582, 583, 584, 585, 586, 587, 588],
        "trans8": [20, 21, 22, 51, 52, 53, 54, 88, 89, 90, 91, 92, 131, 132, 133, 134, 135, 178, 179, 180, 181, 182,
                   183, 230, 231, 232, 233, 234, 235, 236, 286, 287, 288, 289, 290, 291, 292, 323, 324, 325, 354, 355,
                   356, 357, 391, 392, 393, 394, 395, 434, 435, 436, 437, 438, 481, 482, 483, 484, 485, 486, 533, 534,
                   535, 536, 537, 538, 539, 589, 590, 591, 592, 593, 594, 595]
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
            if isinstance(self.groups[grouping], tuple):
                start, end = self.groups[grouping]
                dataframe[grouping] = self.data.iloc[:, 1:(self.num_detectors + 1)].iloc[:, start:end].sum(axis=1)
            elif isinstance(self.groups[grouping], list):
                grouping_list = self.groups[grouping]
                dataframe[grouping] = self.data.iloc[:, 1:(self.num_detectors + 1)].iloc[:, grouping_list].sum(axis=1)
            else:
                print("Group must be defined using list or tuple")
                return
        dataframe.drop(list(range(self.num_detectors)), axis=1, inplace=True)
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
        data = pd.DataFrame(np.stack(detectors, axis=1), columns=["time"] + list(range(self.num_detectors)))
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
                            getattr(self, f"pair_{name}").asymmetry[0][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                            start_time=start_time, end_time=end_time, initial_guess=initial_guess, bounds=bounds)
        elif name in self.groups:
            group_model = FitFunction(function)
            group_model.fit(self.data["time"][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                            self.data[name][(self.data["time"] >= start_time) & (self.data["time"] <= end_time)],
                            start_time=start_time, end_time=end_time, initial_guess=initial_guess, bounds=bounds)
        else:
            print("Enter valid group or pair!")
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
            os.makedirs(f"{save_path}/{self.mac}_plots", exist_ok=True)
            plt.savefig(f"{save_path}/{self.mac}_plots/{self.mac}_{group}.png", dpi=200)
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
        if getattr(self, f"pair_{pair}").asymmetry is None:
            getattr(self, f"pair_{pair}").get_asymmetry()
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
            os.makedirs(f"{save_path}/{self.mac}_plots", exist_ok=True)
            plt.savefig(f"{save_path}/{self.mac}_plots/{self.mac}_pair_{pair}_asymmetry.png", dpi=200)
        if show_plot:
            plt.show()

    def get_histograms(self) -> list:
        with up.open(fr"{self.data_path}/his_{self.mac}_{self.v1190}.v1190.root") as f:
            hist_names = [hist[:-2] for hist in f.keys()]
        return hist_names

    def plot_histogram(self, hist_name: str, show_plot: bool = True, save_path: str = None) -> tuple:
        with up.open(fr"{self.data_path}/his_{self.mac}_{self.v1190}.v1190.root")[f"{hist_name};1"] as f:
            if f.axes.__len__() == 2:
                fig, ax = plt.subplots(dpi=200)
                w, x, y = f.to_numpy()
                mesh = ax.pcolormesh(x, y, w.T)
                fig.colorbar(mesh)
                ax.set_xlabel(f.axes[0].all_members["fTitle"])
                ax.set_ylabel(f.axes[1].all_members["fTitle"])
                ax.set_title(f"{f.title}, {f.name}, {self.mac}.mac", fontsize="small")
                plt.tight_layout()
                hist_data = w, x, y
            else:
                if f.axis().labels() is None:
                    fig, ax = plt.subplots(dpi=200)
                    time, values = f.axis().centers(), f.values()
                    ax.errorbar(time, values, [np.sqrt(x) for x in values],
                                fmt="o", elinewidth=0.6, capsize=1, capthick=0.6, markersize=0.5)
                    ax.set_xlabel(f.axis().all_members["fTitle"])
                    ax.set_ylabel("N")
                    ax.set_title(f"{f.title}, {f.name}, {self.mac}.mac", fontsize="small")
                    plt.tight_layout()
                    hist_data = time, values
                else:
                    fig, ax = plt.subplots(dpi=200)
                    labels, heights = f.axis().labels(), f.values()
                    ax.bar(labels, heights, yerr=[np.sqrt(x) for x in heights], capsize=4)
                    ax.set_ylabel("N")
                    ax.set_title(f'{f.title}, {f.name}, {self.mac}.mac', fontsize="small")
                    plt.xticks(rotation=-25, ha="left", size=7)
                    plt.tight_layout()
                    hist_data = labels, heights
            if save_path is None:
                pass
            else:
                os.makedirs(f"{save_path}/{self.mac}_plots", exist_ok=True)
                fig.savefig(f"{save_path}/{self.mac}_plots/{self.mac}_{f.name}.png", dpi=200)
            if show_plot:
                plt.show()
            plt.clf()
        return hist_data

    def create_pair(self, pair_name: str, group_1: str, group_2: str):
        """Creates a pairing of detector groups. Asymmetry between the groups can then be calculated.

        :param pair_name: Name of pair
        :param group_1: First detector group in pair
        :param group_2: Second detector group in pair
        """
        if hasattr(self, f"pair_{pair_name}"):
            print("Pair already exists with this name!")
        else:
            setattr(self, f"pair_{pair_name}", Pair(pair_name, self, group_1, group_2))


class Pair:
    """Object representing a pairing of detector groups."""

    def __init__(self, name: str, instrument: MuonInstrument, group_1: str, group_2: str):
        """A pairing of detector groups.

        :param group_1: First detector group in the pair, known as "forward"
        :param group_2: Second detector group in the pair, known as "backward
        """
        self.name = name
        self.instrument = instrument
        self.group_1 = group_1
        self.group_2 = group_2
        self.alpha = None
        self.asymmetry = None

    def __repr__(self):
        return f"Pair_{self.name}({self.group_1}, {self.group_2})"

    def get_alpha(self) -> float:
        """Estimates alpha (balance parameter) between the two groups. Used in calculation for asymmetry.

        Defined as sum(forward) / sum(backward)
        """
        self.alpha = self.instrument.data[self.group_1].sum() / self.instrument.data[self.group_2].sum()
        return self.alpha

    def get_asymmetry(self, alpha: float = None) -> pd.DataFrame:
        """Calculates the asymmetry between the two groups.

        Defined as (forward - alpha * backward) / (forward + alpha * backward)

        :param alpha: The balance parameter used in the calculation. By default, it is estimated using get_alpha method.
        """
        if alpha is None:
            alpha = self.get_alpha()
        self.asymmetry = pd.DataFrame([0 if x == y == 0 else (x - alpha * y) / (x + alpha * y) for x, y in
                                       zip(self.instrument.data[self.group_1], self.instrument.data[self.group_2])])
        return self.asymmetry
