from src.musrpy.instruments import MuonInstrument

EMU = MuonInstrument("EMU", num_detectors=96, detector_groups=MuonInstrument.emu_standard_groups)
EMU.load_data(1, 1, path="test_files", bins=200)
EMU.group_data()
