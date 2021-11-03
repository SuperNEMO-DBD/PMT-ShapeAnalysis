import ROOT
import numpy as np


class Calo_Data:
    def __init__(self):
        self.fee_crate = None
        self.fee_board = None
        # self.fee_chip = -1
        self.fee_channel = None

        self.om_side = None
        self.om_wall = None
        self.om_column = None
        self.om_row = None

        self.om_num = None
        self.flag = None

        self.tdc = None
        self.time = None  # // v2

        self.fee_baseline = None
        self.fee_amplitude = None
        self.fee_charge = None

        self.baseline = None
        self.amplitude = None
        self.charge = None

    def set(self, index: int, event):
        self.fee_crate = event.calo_fee_crate[index]
        self.fee_board = event.calo_fee_board[index]
        # self.fee_chip = event.calo_fee_chip[index]
        self.fee_channel = event.calo_fee_channel[index]

        self.om_side = event.calo_om_side[index]
        self.om_wall = event.calo_om_wall[index]
        self.om_column = event.calo_om_column[index]
        self.om_row = event.calo_om_row[index]

        self.om_num = event.calo_om_num[index]
        # self.flag = event.calo_flag[index]

        self.tdc = event.calo_tdc[index]
        self.time = event.calo_time[index]

        self.fee_baseline = event.calo_fee_baseline[index]
        self.fee_amplitude = event.calo_fee_amplitude[index]
        self.fee_charge = event.calo_fee_charge[index]

        self.baseline = event.calo_baseline[index]
        self.amplitude = event.calo_amplitude[index]
        self.charge = event.calo_charge[index]


class Tracker_Data:
    def __init__(self):
        self.fee_crate = None
        self.fee_board = None
        # self.fee_chip = 0
        self.fee_channel = None

        self.cell_side = None
        self.cell_row = None
        self.cell_layer = None

        self.cell_num = None

        self.timestamp_r0 = None
        self.timestamp_r1 = None
        self.timestamp_r2 = None
        self.timestamp_r3 = None
        self.timestamp_r4 = None
        self.timestamp_r5 = None
        self.timestamp_r6 = None

        self.time_anode = None
        self.time_top_cathode = None
        self.time_bottom_cathode = None

    def set(self, index: int, event):
        self.fee_crate = event.tracker_fee_crate[index]
        self.fee_board = event.tracker_fee_board[index]
        # self.fee_chip = event.tracker_fee_chip[index]
        self.fee_channel = event.tracker_fee_channel[index]

        self.cell_side = event.tracker_cell_side[index]
        self.cell_row = event.tracker_cell_row[index]
        self.cell_layer = event.tracker_cell_layer[index]

        self.cell_num = event.tracker_cell_num[index]

        self.timestamp_r0 = event.tracker_timestamp_r0[index]
        self.timestamp_r1 = event.tracker_timestamp_r1[index]
        self.timestamp_r2 = event.tracker_timestamp_r2[index]
        self.timestamp_r3 = event.tracker_timestamp_r3[index]
        self.timestamp_r4 = event.tracker_timestamp_r4[index]
        self.timestamp_r5 = event.tracker_timestamp_r5[index]
        self.timestamp_r6 = event.tracker_timestamp_r6[index]

        self.time_anode = event.tracker_time_anode[index]
        self.time_top_cathode = event.tracker_time_top_cathode[index]
        self.time_bottom_cathode = event.tracker_time_bottom_cathode[index]


def event_builder(input_filename: str):

    input_file = ROOT.TFile(input_filename, "READ")
    input_tree = input_file.event_tree
    print(">>> Number of events in input file:", input_tree.GetEntries())

    first_event = True
    previous_last_time = -86400

    events = []
    event_time = 0
    new_calo_events = []
    new_tracker_events = []

    n_failures = 0

    for event in input_tree:
        try:
            event_time = event.time
            calo_oms = list(event.calo_om_num)
            tracker_cells = list(event.tracker_cell_num)
            tracker_timestamp_r0s = list(event.tracker_timestamp_r0)
            tracker_timestamp_r5s = list(event.tracker_timestamp_r5)
            tracker_timestamp_r6s = list(event.tracker_timestamp_r6)
        except SystemError:
            n_failures += 1
            print("SystemError", n_failures)
            continue

        '''if len(calo_oms) == 0:
            continue'''
        for i_cal in range(len(calo_oms)):
            new_calo_data = Calo_Data()
            new_calo_data.set(i_cal, event)
            new_calo_events.append(new_calo_data)

        for i_tr in range(len(tracker_cells)):
            if tracker_timestamp_r6s[i_tr] != 0 and tracker_timestamp_r5s[i_tr] != 0 and tracker_timestamp_r0s[i_tr] != 0:
                new_tracker_data = Tracker_Data()
                new_tracker_data.set(i_tr, event)
                new_tracker_events.append(new_tracker_data)

        # if len(new_tracker_events) > 0:
        events.append([new_calo_events, new_tracker_events, event_time])

        new_calo_events = []
        new_tracker_events = []

    '''for event in input_tree:
        first_time = 86400  # sec
        last_time = -86400  # sec
        calo_oms = list(event.calo_om_num)
        calo_times = list(event.calo_time)

        for i_cal in range(len(calo_oms)):
            if calo_times[i_cal] < first_time:
                first_time = calo_times[i_cal]

            if calo_times[i_cal] > last_time:
                last_time = calo_times[i_cal]

        try:
            tracker_cells = list(event.tracker_cell_num)
            tracker_sides = list(event.tracker_cell_side)
            tracker_rows = list(event.tracker_cell_row)
            tracker_layers = list(event.tracker_cell_layer)
            tracker_timestamp_r0s = list(event.tracker_timestamp_r0)
            tracker_timestamp_r1s = list(event.tracker_timestamp_r1)
            tracker_timestamp_r2s = list(event.tracker_timestamp_r2)
            tracker_timestamp_r3s = list(event.tracker_timestamp_r3)
            tracker_timestamp_r4s = list(event.tracker_timestamp_r4)
            tracker_timestamp_r5s = list(event.tracker_timestamp_r5)
            tracker_timestamp_r6s = list(event.tracker_timestamp_r6)

            tracker_time_anodes = list(event.tracker_time_anode)
            tracker_time_top_cathodes = list(event.tracker_time_top_cathode)
            tracker_time_bottom_cathodes = list(event.tracker_time_top_cathode)
        except SystemError:
            n_failures += 1
            print("SystemError", n_failures)
            continue

        for i_tr in range(len(tracker_cells)):
            if tracker_timestamp_r0s[i_tr] != 0:
                if tracker_time_anodes[i_tr] < first_time:
                    first_time = tracker_time_anodes[i_tr]
                if tracker_time_anodes[i_tr] > last_time:
                    last_time = tracker_time_anodes[i_tr]
            if tracker_timestamp_r5s[i_tr] != 0:
                if tracker_time_bottom_cathodes[i_tr] < first_time:
                    first_time = tracker_time_bottom_cathodes[i_tr]
                if tracker_time_bottom_cathodes[i_tr] > last_time:
                    last_time = tracker_time_bottom_cathodes[i_tr]
            if tracker_timestamp_r6s[i_tr] != 0:
                if tracker_time_top_cathodes[i_tr] < first_time:
                    first_time = tracker_time_top_cathodes[i_tr]
                if tracker_time_top_cathodes[i_tr] > last_time:
                    last_time = tracker_time_top_cathodes[i_tr]

        print(first_time - previous_last_time, 's')

        if (first_time - previous_last_time) > 60e-6 and not first_event:
            events.append([new_calo_events, new_tracker_events])
            new_calo_events = []
            new_tracker_events = []

        previous_last_time = last_time
        if first_event:
            first_event = False

        # output_event_trigger_v.push_back(input_event_trigger);

        for i_cal in range(len(calo_oms)):
            new_calo_data = Calo_Data()
            new_calo_data.set(i_cal, event)
            new_calo_events.append(new_calo_data)

        for i_tr in range(len(tracker_cells)):
            output_tracker_index = -1

            for index in range(len(new_tracker_events)):
                if tracker_cells[i_tr] != new_tracker_events[index].cell_num:
                    continue
                output_tracker_index = index
                break

            if output_tracker_index == -1:
                new_tracker_data = Tracker_Data()
                new_tracker_data.set(i_tr, event)
                new_tracker_events.append(new_tracker_data)
            else:
                if tracker_timestamp_r0s[i_tr] != 0:
                    if new_tracker_events[output_tracker_index].timestamp_r0 != 0:
                        print(">>> R0 already set for GG:{}.{}.{} in output entry {}\n".format(tracker_sides[i_tr],
                                                                                               tracker_rows[i_tr],
                                                                                               tracker_layers[i_tr],
                                                                                               len(events)))
                        pass
                    else:
                        new_tracker_events[output_tracker_index].timestamp_r0 = tracker_timestamp_r0s[i_tr]
                        new_tracker_events[output_tracker_index].time_anode = tracker_time_anodes[i_tr]
                if tracker_timestamp_r5s[i_tr] != 0:
                    if new_tracker_events[output_tracker_index].timestamp_r5 != 0:
                        print(">>> R5 already set for GG:{}.{}.{} in output entry {}\n".format(tracker_sides[i_tr],
                                                                                               tracker_rows[i_tr],
                                                                                               tracker_layers[i_tr],
                                                                                               len(events)))
                        pass
                    else:
                        new_tracker_events[output_tracker_index].timestamp_r5 = tracker_timestamp_r5s[i_tr]
                        new_tracker_events[output_tracker_index].time_bottom_cathode = tracker_time_bottom_cathodes[i_tr]
                if tracker_timestamp_r6s[i_tr] != 0:
                    if new_tracker_events[output_tracker_index].timestamp_r6 != 0:
                        print(">>> R6 already set for GG:{}.{}.{} in output entry {}\n".format(tracker_sides[i_tr],
                                                                                               tracker_rows[i_tr],
                                                                                               tracker_layers[i_tr],
                                                                                               len(events)))
                        pass
                    else:
                        new_tracker_events[output_tracker_index].timestamp_r6 = tracker_timestamp_r6s[i_tr]
                        new_tracker_events[output_tracker_index].time_top_cathode = tracker_time_top_cathodes[i_tr]'''

    input_file.Close()

    return events


if __name__ == '__main__':
    file = ""
    events = event_builder(file)
