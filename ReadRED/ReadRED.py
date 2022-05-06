import uproot
import format_plot
import os
import sndisplay as sn
from ROOT import TFile


c_tdc2sec = 0.625E-08
t_tdc2sec = 1.25E-08
tdc2ns = 0.390625
adc2mv = 0.610352
MAX_GG_TIMES = 10


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-d', required=False, type=str, help='Output data path for saving plots and files')
    args = parser.parse_args()
    return args


def process_timestamps(timestamps, n_tracker_hits):
    if len(timestamps) != n_tracker_hits*MAX_GG_TIMES:
        raise ValueError("Length of the timestamps: {}".format(len(timestamps)))
    output = [[] for i in range(int(n_tracker_hits))]
    for index, timestamp in enumerate(timestamps):
        j = index // MAX_GG_TIMES
        output[j].append(timestamp)
    return output


def parse_root_file(file_name):
    root_file = TFile(file_name, "READ")
    tree = root_file.RED

    # Build a dictionary full of the data you want to keep
    data = {}

    n_events = tree.GetEntries()
    counter = 0
    print(f">>> Events in tree: {n_events}")
    print(">>> Parsing ...")
    print(">")
    for event in tree:
        ################################################################################################################
        #
        #  Access the branches of the root file
        #  To know what is in the TTree either use:
        #               tree.Print()
        #               OR
        #               branches = tree.GetListOfBranches()
        #               print([branch.GetName() for branch in branches])
        #
        ################################################################################################################
        run_id = event.run_id
        event_id = event.event_id
        # event_clock = event.event_clock                                       # int The clock of the events seconds
        # event_ticks = event.event_ticks                                       # Event times - TODO: not yet implemented
        n_calo_hits = int(event.nb_calo_hits)                                   # int Number of calorimeter hits
        # calo_hit_id = list(event.calo_hit_id)                                 # list(int) List of the calo hit IDs
        calo_om_side_id = list(event.calo_om_side_id)                           # list(int) List of OM sides (1: French, 0: Italian)
        calo_om_wall_id = list(event.calo_om_wall_id)                           # list(int) List of OM walls (0,1) see url:https://nemo.lpc-caen.in2p3.fr/wiki/NEMO/SuperNEMO/Calorimeter
        calo_om_column_id = list(event.calo_om_column_id)                       # list(int) List of OM columns
        calo_om_row_id = list(event.calo_om_row_id)                             # list(int) List of OM rows
        # calo_om_ref_id = list(event.calo_om_ref_id)                           # Don't know?
        # calo_clock = list(event.calo_clock)                                   # Don't know?
        calo_ticks = list(event.calo_ticks)                                     # list(int) Calo TDC (6.25ns)
        calo_lto = list(event.calo_lto)                                         # list(bool) Calo low threshold flag
        calo_ht = list(event.calo_ht)                                           # list(bool) Calo high threshold flag
        # calo_fcr = list(event.calo_fcr)                                       # list(int) Dont' know
        # calo_lt_trigger_counter = list(event.calo_lt_trigger_counter)         # Don't know
        # calo_lt_time_counter = list(event.calo_lt_time_counter)               # Don't know
        calo_fwmeas_baseline = list(event.calo_fwmeas_baseline)                 # list(int) Calo baseline LSB: ADC unit/16
        calo_fwmeas_peak_amplitude = list(event.calo_fwmeas_peak_amplitude)     # list(int) Calo amplitude LSB: ADC unit/8
        calo_fwmeas_peak_cell = list(event.calo_fwmeas_peak_cell)               # list(int) Calo peak cell TDC: 0-1023
        calo_fwmeas_charge = list(event.calo_fwmeas_charge)                     # list(int) Calo charge NO IDEA OF UNITS - do not trust this
        calo_fwmeas_rising_cell = list(event.calo_fwmeas_rising_cell)           # list(int) Calo rising cell (right edge) TDC unit/256
        calo_fwmeas_falling_cell = list(event.calo_fwmeas_falling_cell)         # list(int) Calo falling cell (left edge) TDC unit/256

        calo_waveform = [list(event.calo_waveform)[1024 * i:1024 * (i + 1)] for i in range(n_calo_hits)]  # list(list(int)) Calo waveform

        nb_tracker_hits = event.nb_tracker_hits                                 # int Number of tracker cells that have at least one hit
        # tracker_hit_id = list(event.tracker_hit_id)                           # list(int) TR hit ids
        tracker_cell_side_id = list(event.tracker_cell_side_id)                 # list(int) TR side (1: French, 0: Italian)
        tracker_cell_row_id = list(event.tracker_cell_row_id)                   # list(int) TR row (0-112)
        tracker_cell_layer_id = list(event.tracker_cell_layer_id)               # list(int) TR layer (0-8)
        # tracker_clock = list(event.tracker_clock)                             # list(int) Don't know
        tracker_anode_R0_ticks = process_timestamps(list(event.tracker_anode_R0_ticks), nb_tracker_hits)
        tracker_anode_R1_ticks = process_timestamps(list(event.tracker_anode_R1_ticks), nb_tracker_hits)
        tracker_anode_R2_ticks = process_timestamps(list(event.tracker_anode_R2_ticks), nb_tracker_hits)
        tracker_anode_R3_ticks = process_timestamps(list(event.tracker_anode_R3_ticks), nb_tracker_hits)
        tracker_anode_R4_ticks = process_timestamps(list(event.tracker_anode_R4_ticks), nb_tracker_hits)
        tracker_bottom_cathode_R5_ticks = process_timestamps(list(event.tracker_bottom_cathode_R5_ticks), nb_tracker_hits)
        tracker_top_cathode_R6_ticks = process_timestamps(list(event.tracker_top_cathode_R6_ticks), nb_tracker_hits)

        '''tracker_anode_R0_ticks = [list(event.tracker_anode_R0_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode time LSB: 12.5 ns
        tracker_anode_R1_ticks = [list(event.tracker_anode_R1_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 1st low threshold LSB: 12.5 ns
        tracker_anode_R2_ticks = [list(event.tracker_anode_R2_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 2nd low threshold LSB: 12.5 ns
        tracker_anode_R3_ticks = [list(event.tracker_anode_R3_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 1st high threshold LSB: 12.5 ns
        tracker_anode_R4_ticks = [list(event.tracker_anode_R4_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 2st high threshold LSB: 12.5 ns
        tracker_bottom_cathode_R5_ticks = [list(event.tracker_bottom_cathode_R5_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                           range(nb_tracker_hits)]  # list(int) Cathode bottom LSB: 12.5 ns
        tracker_top_cathode_R6_ticks = [list(event.tracker_top_cathode_R6_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                        range(nb_tracker_hits)]  # list(int) Cathode top LSB: 12.5 ns'''

        ################################################################################################################
        #  Process data here either store in a container like in the data dictionary
        #  OR process directly here
        ################################################################################################################

        ################
        # Quality Cuts #
        ################
        if n_calo_hits == 0 or nb_tracker_hits == 0 or nb_tracker_hits > 50:
            # Cut out noisy data
            continue

        calo_om_nums = [None for i in range(n_calo_hits)]
        for i_om in range(n_calo_hits):
            if calo_om_row_id[i_om] == -1:
                # GVETO OM
                calo_om_num = 520 + 128 + calo_om_side_id[i_om] * 32 + calo_om_wall_id[i_om] * 16 + calo_om_column_id[i_om]
            elif calo_om_wall_id[i_om] == -1:
                # MAIN WALL OM
                calo_om_num = calo_om_side_id[i_om] * 20 * 13 + calo_om_column_id[i_om] * 13 + calo_om_row_id[i_om]
            else:
                # XWALL OM
                calo_om_num = 520 + calo_om_side_id[i_om] * 64 + calo_om_wall_id[i_om] * 32 + calo_om_column_id[i_om] * 16 + calo_om_row_id[i_om]
            calo_om_nums[i_om] = calo_om_num

        try:
            # Take the first calorimeter time
            calo_time = min(calo_ticks) * c_tdc2sec
        except ValueError:
            # If there are no calo hits (its likely noise) set the calo time to 0
            calo_time = 0

        tracker_cell_nums = [None for i in range(nb_tracker_hits)]
        tracker_ppts = [None for i in range(nb_tracker_hits)]
        tracker_cell_times = [[] for i in range(nb_tracker_hits)]
        tracker_cell_timestamps = [[] for i in range(nb_tracker_hits)]

        for i_cell in range(nb_tracker_hits):
            tracker_cell_num = tracker_cell_side_id[i_cell] * 113 * 9 + tracker_cell_row_id[i_cell] * 9 + tracker_cell_layer_id[i_cell]
            # Note: There are 10 entries for each of the cell R0-6. I have chosen to use the 0th entry
            # >>> If the value is -9223372036854775808 this is the default value
            # >>> It is when there is a value for one RX but not another RY for a specific cell - ignore these values

            # TDC RX
            rs = [tracker_anode_R0_ticks[i_cell][0],
                  tracker_anode_R1_ticks[i_cell][0], tracker_anode_R2_ticks[i_cell][0],
                  tracker_anode_R3_ticks[i_cell][0], tracker_anode_R4_ticks[i_cell][0],
                  tracker_bottom_cathode_R5_ticks[i_cell][0], tracker_top_cathode_R6_ticks[i_cell][0]]
            tracker_cell_timestamps[i_cell] = rs
            # Time (RX - calo hit time) in µs
            ts = []
            for r in rs:
                if r == -9223372036854775808:
                    ts.append(None)
                elif r == 0:
                    ts.append(None)
                else:
                    ts.append((r*t_tdc2sec - calo_time)*1E6)
            tracker_cell_times[i_cell] = ts
            tracker_cell_nums[i_cell] = tracker_cell_num
            try:
                tracker_ppts[i_cell] = ts[5] + ts[6]
            except TypeError:
                pass

        # Fill the dictionary
        data[event_id] = {
            "cells": {cell: {
                "rs": tracker_cell_timestamps[index],
                "ts": tracker_cell_times[index],
                "ppts": tracker_ppts[index]
            } for index, cell in enumerate(tracker_cell_nums)},
            "oms": {om: {
                "amplitudes": -1 * calo_fwmeas_peak_amplitude[index] * adc2mv/8
            } for index, om in enumerate(calo_om_nums)}
        }

        '''if event_id == 0:
            for index in range(len(tracker_cell_timestamps)):
                print(tracker_cell_nums[index], tracker_cell_timestamps[index])'''

        # Progress
        new_counter = (event_id + 1) / n_events * 100 // 10
        if new_counter != counter:
            counter = new_counter
            print("> Progress: {:.1f}".format(event_id / n_events * 100) + f"%. Events {len(data.keys())}")
    return data


def plot_event_map(data: dict, output_path):
    sntracker = sn.tracker(new_name='tracker_event_map', with_palette=True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{}')

    for event in data:
        cell_info = data[event]["cells"]
        for cell in cell_info:
            sntracker.fill(cell, 1)

    sntracker.draw()
    sntracker.save(output_path)
    del sntracker

    sncalo = sn.calorimeter(new_name='calo_event_map', with_palette=True)
    sncalo.draw_content_label('{}')
    sncalo.draw_omid_label()
    for event in data:
        om_info = data[event]["oms"]
        for om in om_info:
            sncalo.fill(om)
    sncalo.draw()
    sncalo.save(output_path)


def main():
    args = parse_arguments()
    input_file = args.i
    output_path = args.d
    if output_path is None:
        output_path = os.getcwd()
    print(f">>> File to process {input_file}")

    data = parse_root_file(input_file)
    plot_event_map(data=data, output_path=output_path)


if __name__ == "__main__":
    main()
