import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sndisplay as sn

is_none = -9999


def plot_event_map(name: str, hits: list, output_directory):
    sntracker = sn.tracker(name, with_palette=True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(hits)):
        if hits[i_cell] == 0:
            continue
        sntracker.setcontent(i_cell, hits[i_cell])

    # sntracker.setrange(0, 35)

    sntracker.draw()
    sntracker.save(output_directory)


def main(cut: bool):
    root_file = ROOT.TFile("/path/to/file.root", "READ")
    event_tree = root_file.event_tree

    top_cathode_events = [0 for i in range(2034)]
    bottom_cathode_events = [0 for i in range(2034)]
    anode_events = [0 for i in range(2034)]

    for event in event_tree:
        if cut:
            if len(event.tracker_time_anode) > 50:
                continue
            if len(event.calo_time) == 0:
                continue
            ht = False
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    ht = True
            if not ht:
                continue

        for index, i_cell in enumerate(event.tracker_cell_num):
            if event.tracker_time_top_cathode[index] == is_none:
                pass
            else:
                top_cathode_events[event.tracker_cell_num[index]] += 1
            if event.tracker_time_bottom_cathode[index] == is_none:
                pass
            else:
                bottom_cathode_events[event.tracker_cell_num[index]] += 1
            if event.tracker_time_anode[index] == is_none:
                pass
            else:
                anode_events[event.tracker_cell_num[index]] += 1

    plot_event_map("top_name", top_cathode_events, ".")
    plot_event_map("bot_name", bottom_cathode_events, ".")
    plot_event_map("anode_name", anode_events, ".")


if __name__ == "__main__":
    main(True)
