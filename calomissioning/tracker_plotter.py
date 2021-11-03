import sys
import ROOT
sys.path.insert(1, '../')

import matplotlib.pyplot as plt
import numpy as np
import sndisplay
from functions.other_functions import *
from event_builder import *


def draw_anode_events(run: int, events):
    sn_tracker = sndisplay.tracker('tr_anode_' + str(run))
    sn_tracker.draw_cellid_label()
    sn_tracker.draw_content_label('{}')

    n_events = [0]*2034

    for event in events:
        tracker_events = event[1]
        for i_tr in tracker_events:
            if i_tr.time_anode != 0:
                n_events[i_tr.cell_num] += 1
    for i in range(2034):
        if n_events[i] > 0:
            sn_tracker.setcontent(i, n_events[i])

    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/" + str(run))


def draw_top_cathode_events(run: int, events):
    sn_tracker = sndisplay.tracker('tr_top_cathode_' + str(run))
    sn_tracker.draw_cellid_label()
    sn_tracker.draw_content_label('{}')

    n_events = [0] * 2034

    for event in events:
        tracker_events = event[1]
        for i_tr in tracker_events:
            # print(i_tr.time_anode, i_tr.time_top_cathode)
            if i_tr.time_top_cathode != 0:
                n_events[i_tr.cell_num] += 1
    for i in range(2034):
        if n_events[i] > 0:
            sn_tracker.setcontent(i, n_events[i])

    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/" + str(run))


def draw_bottom_cathode_events(run: int, events):
    sn_tracker = sndisplay.tracker('tr_bottom_cathode_' + str(run))
    sn_tracker.draw_cellid_label()
    sn_tracker.draw_content_label('{}')

    n_events = [0] * 2034

    for event in events:
        tracker_events = event[1]
        for i_tr in tracker_events:
            if i_tr.time_bottom_cathode != 0:
                n_events[i_tr.cell_num] += 1
    for i in range(2034):
        if n_events[i] > 0:
            sn_tracker.setcontent(i, n_events[i])

    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/" + str(run))


def draw_events(run: int, events):
    for i in range(10):
        sn_demo = sndisplay.demonstrator(str(run) + '_demo_' + str(i))

        for calo_event in events[i][0]:
            sn_demo.setomcontent(calo_event.om_num, 1)
        for tracker_event in events[i][1]:
            if tracker_event.time_anode != 0 and tracker_event.time_top_cathode != 0 and tracker_event.time_bottom_cathode != 0:
                sn_demo.setggcontent(tracker_event.cell_num, 1)
        sn_demo.draw_top()
        sn_demo.save("/Users/williamquinn/Desktop/" + str(run))


def draw_calo_events(run: int, events):
    sn_calo = sndisplay.calorimeter('calo_om_' + str(run))
    sn_calo.draw_omid_label()
    sn_calo.draw_content_label('{}')

    n_events = [0] * 712

    for event in events:
        calo_events = event[0]
        for i_calo in calo_events:
            n_events[i_calo.om_num] += 1
    for i in range(712):
        if n_events[i] > 0:
            sn_calo.setcontent(i, n_events[i])

    sn_calo.draw()
    sn_calo.save("/Users/williamquinn/Desktop/" + str(run))


def draw_full_cell_events(run: int, events):
    sn_tracker = sndisplay.tracker('tr_full_cell_' + str(run))
    sn_tracker.draw_cellid_label()
    sn_tracker.draw_content_label('{}')

    n_events = [0] * 2034

    for event in events:
        tracker_events = event[1]
        for i_tr in tracker_events:
            if i_tr.time_anode != 0 and i_tr.time_top_cathode != 0 and i_tr.time_bottom_cathode != 0:
                n_events[i_tr.cell_num] += 1
    for i in range(2034):
        if n_events[i] > 0:
            sn_tracker.setcontent(i, n_events[i])

    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/" + str(run))


def draw_t56_all_cells_2D(run: int, events):
    lower = 0
    higher = 100
    ROOT.gStyle.SetOptStat(0)
    hist = ROOT.TH2F("", "", 40, lower, higher,
                     40, lower, higher)
    can = ROOT.TCanvas()
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time)*1e6
            t6 = (tracker_event.time_top_cathode - event_time)*1e6
            if t5 < 10 and t6 < 10:
                continue
            if lower < t5 < higher and lower < t6 < higher:
                hist.Fill(t5, t6)
    can.cd()
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("t5 /us")
    hist.GetYaxis().SetTitle("t6 /us")
    can.SetGrid()
    can.Draw()
    can.SaveAs("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t56_2D.pdf")

    del hist
    del can


def draw_t56_all_cells(run: int, events):
    lower = -50
    higher = 100

    f = []
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time)*1e6
            t6 = (tracker_event.time_top_cathode - event_time)*1e6
            if t5 < 10 and t6 < 10:
                continue
            if lower < t5 < higher and lower < t6 < higher:
                f.append(t5+t6)

    fig = plt.figure(figsize=(9, 6), facecolor='white')
    freq, bin_edges = np.histogram(f, 100, range=(0, 100))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, color="blue", width=width, label='')
    plt.xlabel('t5 + t6 /µs')
    plt.ylabel('counts')
    plt.grid()
    # plt.axvline(65, ls='--', color='red')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t56.pdf")
    plt.close()


def draw_t56_by_layer(run: int, events):
    lower = -50
    higher = 100

    f = [[] for i in range(9)]
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            cell_num = tracker_event.cell_num
            cell_side = cell_num // (9 * 113)
            cell_row = cell_num % (9 * 113) // 9
            cell_layer = cell_num % (9 * 113) % 9
            t5 = (tracker_event.time_bottom_cathode - event_time)*1e6
            t6 = (tracker_event.time_top_cathode - event_time)*1e6
            if t5 < 10 and t6 < 10:
                continue
            if lower < t5 < higher and lower < t6 < higher:
                f[cell_layer].append(t5+t6)

    for i in range(len(f)):
        fig = plt.figure(figsize=(9, 6), facecolor='white')
        freq, bin_edges = np.histogram(f[i], 100, range=(0, 100))
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, color="blue", width=width, label='')
        plt.xlabel('t5 + t6 /µs')
        plt.ylabel('counts')
        plt.grid()
        # plt.axvline(65, ls='--', color='red')
        plt.tight_layout()
        plt.savefig("/Users/williamquinn/Desktop/" + str(run) + "/t56_layer_{}.pdf".format(i))
        plt.close()


def draw_av_drift(run: int, events):
    lower = -50
    higher = 100

    sn_tracker = sndisplay.tracker('tr_av_drift_' + str(run))
    sn_tracker.draw_cellid_label()
    sn_tracker.draw_content_label('{:.2f}')

    drift_times = [[] for i in range(2034)]

    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time)*1e6
            t6 = (tracker_event.time_top_cathode - event_time)*1e6
            if t5 < 10 and t6 < 10:
                continue
            if lower < t5 < higher and lower < t6 < higher:
                drift_times[int(tracker_event.cell_num)].append(t5 + t6)

    for i_cell in range(len(drift_times)):
        average = 0
        if len(drift_times[i_cell]) == 0:
            continue
        for i_drift in range(len(drift_times[i_cell])):
            average += drift_times[i_cell][i_drift]
        average = average/len(drift_times[i_cell])
        sn_tracker.setcontent(i_cell, average)

    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/" + str(run))


def draw_std_drift(run: int, events):
    lower = -50
    higher = 100

    sn_tracker = sndisplay.tracker('tr_std_drift_' + str(run))
    sn_tracker.draw_cellid_label()
    sn_tracker.draw_content_label('{:.2f}')

    drift_times = [[] for i in range(2034)]

    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time)*1e6
            t6 = (tracker_event.time_top_cathode - event_time)*1e6
            if t5 < 10 and t6 < 10:
                continue
            if lower < t5 < higher and lower < t6 < higher:
                drift_times[tracker_event.cell_num].append(t5+t6)

    for i_cell in range(2034):
        temp_drift = np.array(drift_times[i_cell], dtype='float')
        if temp_drift.size == 0:
            continue
        std = np.std(temp_drift)
        sn_tracker.setcontent(i_cell, std)

    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/" + str(run))


def draw_t056_comp(run: int, events):
    lower = -50
    higher = 100
    c = []
    d = []
    e = []
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t0 = (tracker_event.time_anode - event_time)*1e6
            t5 = (tracker_event.time_bottom_cathode - event_time)*1e6
            t6 = (tracker_event.time_top_cathode - event_time)*1e6
            if lower < t0 < higher and lower < t5 < higher and lower < t6 < higher:
                c.append(t6)
                d.append(t5)
                e.append(t0)

    fig = plt.figure(figsize=(9, 6), facecolor='white')
    freq, bin_edges = np.histogram(c, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "b-", label='t6')

    freq, bin_edges = np.histogram(d, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "g-", label='t5')

    freq, bin_edges = np.histogram(e, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "r-", label='t0')

    plt.grid()
    plt.xlabel("drift time /µs")
    plt.ylabel("counts")
    plt.xlim(lower, higher)
    plt.title("All Cells")
    plt.yscale('log')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t056_comp.pdf")
    plt.close()


def draw_t012_comp(run: int, events):
    lower = -50
    higher = 100
    c = []
    d = []
    e = []
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                c.append(t1)
                d.append(t2)
                e.append(t0)

    fig = plt.figure(figsize=(9, 6), facecolor='white')
    freq, bin_edges = np.histogram(c, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "b-", label='t1')

    freq, bin_edges = np.histogram(d, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "g-", label='t2')

    freq, bin_edges = np.histogram(e, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "r-", label='t0')

    plt.grid()
    plt.xlabel("drift time /µs")
    plt.ylabel("counts")
    plt.xlim(lower, higher)
    plt.title("All Cells")
    plt.yscale('log')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t012_comp.pdf")
    plt.close()


def draw_top_bottom_comp(run: int, events):
    lower = -50
    higher = 100
    c = []
    d = []
    e = []
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                if t5 < t6:
                    if t1 < t2:
                        t_bottom = t1
                        t_top = t2
                    else:
                        t_bottom = t2
                        t_top = t1
                else:
                    if t1 < t2:
                        t_bottom = t2
                        t_top = t1
                    else:
                        t_bottom = t1
                        t_top = t2

                c.append(t_bottom)
                d.append(t_top)
                e.append(t0)

    fig = plt.figure(figsize=(9, 6), facecolor='white')
    freq, bin_edges = np.histogram(c, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "b-", label='t_bottom')

    freq, bin_edges = np.histogram(d, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "g-", label='t_top')

    freq, bin_edges = np.histogram(e, 250, range=(lower, higher))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq, "r-", label='t0')

    plt.grid()
    plt.xlabel("drift time /µs")
    plt.ylabel("counts")
    plt.xlim(lower, higher)
    plt.title("All Cells")
    plt.yscale('log')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_top_bottom_comp.pdf")
    plt.close()


def draw_bottom_signal_comp(run: int, events):
    lower = 0
    higher = 100
    ROOT.gStyle.SetOptStat(0)
    hist = ROOT.TH2F("", "", 40, lower, higher,
                     40, lower, higher)

    bottom = [[] for i in range(2034)]

    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                if t5 < t6:
                    if t1 < t2:
                        t_bottom = t1
                        t_top = t2
                    else:
                        t_bottom = t2
                        t_top = t1
                else:
                    if t1 < t2:
                        t_bottom = t2
                        t_top = t1
                    else:
                        t_bottom = t1
                        t_top = t3
                bottom[tracker_event.cell_num].append([t5, t_bottom])

    for i_cellnum in range(len(bottom)):
        for i in range(len(bottom[i_cellnum])):
            hist.Fill(bottom[i_cellnum][i][0], bottom[i_cellnum][i][1])

    can = ROOT.TCanvas()
    can.cd()
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("t5 /us")
    hist.GetYaxis().SetTitle("t1 or t2 /us")
    can.SetGrid()
    can.Draw()
    can.SaveAs("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_bottom_signal_comp.pdf")
    del hist
    del can


def draw_top_signal_comp(run: int, events):
    lower = 0
    higher = 100
    ROOT.gStyle.SetOptStat(0)
    hist = ROOT.TH2F("", "", 40, lower, higher,
                     40, lower, higher)

    top = [[] for i in range(2034)]

    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                if t5 < t6:
                    if t1 < t2:
                        t_bottom = t1
                        t_top = t2
                    else:
                        t_bottom = t2
                        t_top = t1
                else:
                    if t1 < t2:
                        t_bottom = t2
                        t_top = t1
                    else:
                        t_bottom = t1
                        t_top = t2
                top[tracker_event.cell_num].append([t6, t_top])

    for i_cellnum in range(len(top)):
        for i in range(len(top[i_cellnum])):
            hist.Fill(top[i_cellnum][i][0], top[i_cellnum][i][1])

    can = ROOT.TCanvas()
    can.cd()
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("t6 /us")
    hist.GetYaxis().SetTitle("t1 or t2 /us")
    can.SetGrid()
    can.Draw()
    can.SaveAs("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_top_signal_comp.pdf")
    del hist
    del can


def draw_t12_all_cells(run: int, events):
    lower = -50
    higher = 100

    f = []
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                f.append(t1+t2)

    fig = plt.figure(figsize=(9, 6), facecolor='white')
    freq, bin_edges = np.histogram(f, 100, range=(0, 100))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, color="blue", width=width, label='')
    plt.xlabel('t1 + t2 /µs')
    plt.ylabel('counts')
    plt.grid()
    # plt.axvline(65, ls='--', color='red')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t12.pdf")
    plt.close()


def draw_t12_all_cells_2D(run: int, events):
    lower = 0
    higher = 100
    ROOT.gStyle.SetOptStat(0)
    hist = ROOT.TH2F("", "", 40, lower, higher,
                     40, lower, higher)
    can = ROOT.TCanvas()
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                hist.Fill(t1, t2)
    can.cd()
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("t1 /us")
    hist.GetYaxis().SetTitle("t2 /us")
    can.SetGrid()
    can.Draw()
    can.SaveAs("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t12_2D.pdf")

    del hist
    del can


def draw_top_bottom_all_cells_2D(run: int, events):
    lower = 0
    higher = 100
    ROOT.gStyle.SetOptStat(0)
    hist = ROOT.TH2F("", "", 40, lower, higher,
                     40, lower, higher)
    can = ROOT.TCanvas()
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                if t5 < t6:
                    if t1 < t2:
                        t_bottom = t1
                        t_top = t2
                    else:
                        t_bottom = t2
                        t_top = t1
                else:
                    if t1 < t2:
                        t_bottom = t2
                        t_top = t1
                    else:
                        t_bottom = t1
                        t_top = t2

                hist.Fill(t_bottom, t_top)
    can.cd()
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("t_bottom /us")
    hist.GetYaxis().SetTitle("t_top /us")
    can.SetGrid()
    can.Draw()
    can.SaveAs("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_top_bottom_2D.pdf")

    del hist
    del can


def draw_t1256_comp(run: int, events):
    lower = 0
    higher = 100
    ROOT.gStyle.SetOptStat(0)
    hist = ROOT.TH2F("", "", 40, lower, higher,
                     40, lower, higher)
    can = ROOT.TCanvas()
    for event in events:
        event_time = []
        for calo_event in event[0]:
            event_time.append(calo_event.time)
        event_time = np.array(event_time)
        if event_time.size == 0:
            continue
        event_time = np.min(event_time)
        tracker_events = event[1]
        for tracker_event in tracker_events:
            t5 = (tracker_event.time_bottom_cathode - event_time) * 1e6
            t6 = (tracker_event.time_top_cathode - event_time) * 1e6

            t0 = (tracker_event.time_anode - event_time) * 1e6
            t1 = (tracker_event.timestamp_r1 * 1.25e-8 - event_time) * 1e6
            t2 = (tracker_event.timestamp_r2 * 1.25e-8 - event_time) * 1e6
            t3 = (tracker_event.timestamp_r3 * 1.25e-8 - event_time) * 1e6
            t4 = (tracker_event.timestamp_r4 * 1.25e-8 - event_time) * 1e6

            if t5 < 10 and t6 < 10:
                continue

            if lower < t1 < higher and lower < t2 < higher and lower < t5 < higher and lower < t6 < higher:
                hist.Fill(t5 + t6, t1 + t2)
    can.cd()
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("t5 + t6 /us")
    hist.GetYaxis().SetTitle("t1 + t2 /us")
    can.SetGrid()
    can.Draw()
    can.SaveAs("/Users/williamquinn/Desktop/" + str(run) + "/all_cell_t1256_2D.pdf")

    del hist
    del can


def main():
    args = io_parse_arguments()
    input_file = args.i
    run_num = args.o

    if input_file is None:
        input_file = "/Users/williamquinn/Desktop/test.root"

    print(">>> input file: ", input_file)
    events = event_builder(input_file)
    print(">>> Number of merged events:", len(events))

    draw_anode_events(run_num, events)
    draw_top_cathode_events(run_num, events)
    draw_bottom_cathode_events(run_num, events)
    '''draw_calo_events(run_num, events)
    draw_full_cell_events(run_num, events)
    draw_events(run_num, events)
    draw_t56_all_cells(run_num, events)
    draw_av_drift(run_num, events)
    draw_std_drift(run_num, events)
    draw_t56_all_cells_2D(run_num, events)
    draw_t056_comp(run_num, events)
    draw_bottom_signal_comp(run_num, events)
    draw_top_signal_comp(run_num, events)
    draw_t012_comp(run_num, events)
    draw_t1256_comp(run_num, events)
    draw_t12_all_cells(run_num, events)
    draw_t12_all_cells_2D(run_num, events)
    draw_top_bottom_comp(run_num, events)
    draw_top_bottom_all_cells_2D(run_num, events)
    draw_t56_by_layer(run_num, events)'''

    print(">>> Finished")


if __name__ == '__main__':
    main()
