import sys

import json
from multiprocessing import Process

import numpy as np

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


adc2mv = 0.610352


def merge_files():

    chain = ROOT.TChain("T")
    for i in [521, 522]:
        file = "/Users/williamquinn/Desktop/commissioning/run_{}_waveform.root".format(i)
        chain.Add(file)

    chain.Merge("/Users/williamquinn/Desktop/commissioning/demon_ap.root")


def plot_event_map(events):
    sncalo = sn.calorimeter("demo_event_map", with_palette=True)
    sncalo.draw_omid = False
    # sncalo.draw_omnum_label()
    sncalo.draw_content = False
    for om, val in enumerate(events):
        if val > 0:
            sncalo.setcontent(om, val)
    sncalo.draw()
    sncalo.save("/Users/williamquinn/Desktop/commissioning")
    del sncalo


def store_waveform(waveform, name, freq):
    with open("/Users/williamquinn/Desktop/commissioning/" + name + ".csv", "w") as file:
        for i, val in enumerate(waveform):
            file.write('{},{}\n'.format(i/freq, val))


def get_bad_temp():
    # The peak of this pulse is about the 10th bin
    template = []
    root_file = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/templates_521.root")
    temp_hist = root_file.Get("Template_Ch1")
    for i_bin in range(temp_hist.GetNbinsX()):
        template.append(temp_hist.GetBinContent(i_bin + 1))
    template = np.array(template)
    template = template/np.sqrt(np.dot(template, template))
    return template


def mf_sweep(waveform, template, shape_cut, amp_cut):
    # the template needs to be normalised
    # The waveform needs to be baseline subtracted
    sweep_start = 250
    shape = []
    amplitude = []
    for i in range(250, len(waveform) - len(template)):
        test = waveform[i:i+len(template)]
        test_norm = np.sqrt(np.dot(test, test))
        amplitude_index = np.dot(template, test)
        shape_index = amplitude_index/test_norm
        shape.append(shape_index)
        amplitude.append(amplitude_index)
    ap_shapes, _ = find_peaks(shape, height=shape_cut, distance=len(template))
    ap_amps, _ = find_peaks(amplitude, height=amp_cut, distance=len(template))

    ap_info = {}
    for i in ap_shapes:
        if i in ap_amps:
            ap_info[i + sweep_start] = {"shape": shape[i], "amplitude": amplitude[i]}

    return ap_info


def create_templates():
    template = get_bad_temp()
    filename = "/Users/williamquinn/Desktop/commissioning/demon_ap.root"
    file = ROOT.TFile(filename, "READ")
    tree = file.T
    events = [0 for i in range(4)]
    n_average = 1000
    i_e = 0
    n_events = tree.GetEntries()
    ap_templates = [[0 for i in range(len(template))] for j in range(4)]
    for event in tree:
        '''if sum(events) >= 4*n_average:
            break'''
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events, "aps:", events)
        for index, om in enumerate(list(event.OM_ID)):
            om_t, om_s = om_type(om)
            if events[om_t] >= n_average:
                continue
            waveform = list(event.waveform)[index * 1024:1024 * (1 + index)]
            baseline = get_baseline(waveform, 20)
            amplitude = get_amplitude(waveform, baseline) * adc2mv * -1
            if amplitude < 300:
                continue
            else:
                # store_waveform((np.array(waveform) - baseline)*adc2mv, "demon_waveform", 0.64)
                ap_info = mf_sweep([(val - baseline) * adc2mv for val in waveform], template, shape_cut=0.85, amp_cut=20)
                if len(ap_info.keys()) == 0:
                    # No after-pulses
                    continue
                elif len(ap_info.keys()) > 1:
                    # Too many after-pulses
                    continue
                for ap_time in ap_info.keys():
                    ap_pulse = [val - baseline for index, val in enumerate(waveform[ap_time:ap_time+len(template)])]
                    '''plt.plot([i for i in range(len(waveform))], waveform, 'C0.-')
                    plt.plot([i for i in range(len(waveform))][ap_time:ap_time + len(template)], ap_pulse, 'C1.-')
                    plt.plot([i for i in range(len(waveform))][ap_time + 10], waveform[ap_time + 10], 'C2x')
                    plt.axhline(baseline, ls='--')
                    plt.show()'''
                    # amp = ap_time[ap_time]["amplitude"]
                    amp = -1*np.min(ap_pulse)
                    ap_templates[om_t] = [val + ap_pulse[index]/amp for index, val in enumerate(ap_templates[om_t])]
                    events[om_t] += 1

    out_file = open("/Users/williamquinn/Desktop/commissioning/ap_templates.csv", "w")
    for i in range(4):
        for j in range(len(ap_templates[i])):
            out_file.write('{}\n'.format(np.array(ap_templates[i])[j]/events[i]))
        plt.plot(np.array(ap_templates[i])/events[i])
    plt.show()

    file.Close()


def alter_template():
    templates = [[], [], [], []]
    with open("/Users/williamquinn/Desktop/commissioning/ap_templates.csv", "r") as file:
        fl = file.readlines()
        i = 0
        for index, line in enumerate(fl):
            if index < 40:
                i = 0
            elif index < 80:
                i = 1
            elif index < 120:
                i = 2
            elif index < 160:
                i = 3
            else:
                continue
            templates[i].append(float(line.strip()))
    for i in range(4):
        baseline = np.average(templates[i][:5])
        plt.plot(np.array(templates[i]) - baseline, ".-")
        templates[i] = np.array(templates[i]) - baseline
        templates[i] = templates[i]/np.sqrt(np.dot(templates[i], templates[i]))
    plt.show()
    with open("/Users/williamquinn/Desktop/commissioning/ap_templates_altered.csv", "w") as file:
        for i in range(len(templates[0])):
            file.write('{},{},{},{}\n'.format(templates[0][i], templates[1][i], templates[2][i], templates[3][i]))


def get_templates():
    templates = [[], [], [], []]
    with open("/Users/williamquinn/Desktop/commissioning/ap_templates_altered.csv", "r") as file:
        fl = file.readlines()
        for line in fl:
            line_list = line.split(",")
            for index, val in enumerate(line_list):
                templates[index].append(float(val.strip()))
    templates = [np.array(templates[i]) for i in range(4)]
    return templates


def plot_aan(ap_info):
    sncalo = sn.calorimeter("aan", with_palette=True)
    sncalo.draw_content_label('{:.2f}')
    sncalo.draw_omid = False
    avs = [[], [], [], []]
    for om in range(712):
        if om not in ap_info.keys():
            continue
        om_t, om_s = om_type(om)
        nums = []
        for i in range(len(ap_info[om]["ap_num"])):
            nums.append(ap_info[om]["ap_num"][i])
        av = np.average(nums)
        avs[om_t].append(av)
        sncalo.setcontent(om, av)
    sncalo.draw()
    sncalo.save("/Users/williamquinn/Desktop/commissioning")
    del sncalo

    plt.figure(figsize=figsize)
    names = ['5inchMW', '8inchMW', '5inchXW', '5inchGV']
    for i in range(4):
        freq, bin_edges = np.histogram(avs[i], bins=10, range=(0, 2.5))
        width = bin_edges[2] - bin_edges[1]
        bin_centres = bin_edges[:-1] + width/2
        plt.bar(bin_centres, freq/np.sum(freq), color='C{}'.format(i), alpha=0.5, width=width,
                label=names[i] + r' $\mu$={:.2f}'.format(np.average(avs[i])))
        plt.axvline(np.average(avs[i]), ls='--', color='C{}'.format(i))
    plt.xlabel('Average After-pulse Number')
    plt.ylabel('OM Count')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/aan_demon_all_oms.pdf")


def store_test_bench_waveform():
    file = open("/Users/williamquinn/Desktop/PMT_Project/200117_A1400_B1400_t1003.xml", "r")
    fl = file.readlines()
    for index, line in enumerate(fl):
        line_list = line.split(" ")
        if len(line_list) > 1000:
            waveform = []
            for i in range(3, len(line_list)):
                waveform.append(int(line_list[i].strip()))
            baseline = get_baseline(waveform, 200)
            amplitude = -1 * get_amplitude(waveform, baseline)

            if amplitude > 300 and 0 not in waveform:
                with open("/Users/williamquinn/Desktop/commissioning/test_bench_waveform.csv", "w") as outfile:
                    for i, val in enumerate(waveform):
                        outfile.write('{},{}\n'.format(i, (val - baseline)/amplitude))
                break
    file.close()


def plot_waveform_comparisons():
    # Get the afterpulse average waveform for the 8" PMTs
    after_pulse = [np.array([i/0.64 for i in range(len(get_templates()[1]))]), np.array(get_templates()[1])]
    ap_peak = np.argmin(after_pulse[1])

    demon_pulse = [[], []]
    with open("/Users/williamquinn/Desktop/commissioning/demon_waveform.csv", "r") as temp_file:
        fl = temp_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")
            demon_pulse[0].append(float(line_list[0].strip()))
            demon_pulse[1].append(float(line_list[1].strip()))
    demon_pulse = [np.array(demon_pulse[0]), np.array(demon_pulse[1])]
    d_peak = np.argmin(demon_pulse[1])

    test_bench_pulse = [[], []]
    with open("/Users/williamquinn/Desktop/commissioning/test_bench_waveform.csv", "r") as temp_file:
        fl = temp_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")
            test_bench_pulse[0].append(float(line_list[0].strip()))
            test_bench_pulse[1].append(float(line_list[1].strip()))
    test_bench_pulse = [np.array(test_bench_pulse[0]), np.array(test_bench_pulse[1])]
    t_peak = np.argmin(test_bench_pulse[1])

    plt.figure(figsize=figsize)
    plt.plot(demon_pulse[0][d_peak-10:d_peak+30] - demon_pulse[0][d_peak],
             (-1*demon_pulse[1]/np.min(demon_pulse[1]))[d_peak-10:d_peak+30], "C0.-",
             label='Demonstrator')
    plt.plot(test_bench_pulse[0][t_peak-10:t_peak+50] - test_bench_pulse[0][t_peak],
             -1*test_bench_pulse[1][t_peak-10:t_peak+50]/np.min(test_bench_pulse[1][t_peak-10:t_peak+50]), "C1.--",
             label='Test Bench')
    plt.plot(after_pulse[0] - after_pulse[0][ap_peak],
             -1*after_pulse[1]/np.min(after_pulse[1]), "C2.-", label='After-Pulse')
    plt.legend()
    plt.xlim(-15, 44)
    plt.xlabel("Timestamp /ns")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/demon_test_bench_comp.pdf")


def reduce_data(n_record):
    template = get_bad_temp()
    filename = "/Users/williamquinn/Desktop/commissioning/demon_ap.root"
    file = ROOT.TFile(filename, "READ")
    tree = file.T

    out_file = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/demon_ap_reduced.root", "RECREATE")
    out_file.cd()
    new_tree = ROOT.TTree('new_tree', 'Tree containing reduced RTD data')
    b_om = array('i', [0])
    b_event_num = array('i', [0])
    b_waveform = array('i', 1024 * [0])
    new_tree.Branch('OM_ID', b_om, 'OM_ID/I')
    new_tree.Branch('event_num', b_om, 'event_num/I')
    new_tree.Branch('waveform', b_waveform, 'waveform/I')

    events = [0 for i in range(712)]
    i_e = 0
    n_events = tree.GetEntries()
    for event in tree:
        '''if sum(events) >= 4*n_average:
            break'''
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events)
        for index, om in enumerate(list(event.OM_ID)):
            if events[om] == n_record:
                continue
            waveform = list(event.waveform)[index * 1024:1024 * (1 + index)]
            baseline = get_baseline(waveform, 20)
            amplitude = get_amplitude(waveform, baseline) * adc2mv * -1
            if amplitude < 300:
                continue
            else:
                b_om[0] = om
                b_event_num[0] = event.event_num
                for i, val in enumerate(waveform):
                    b_waveform[i] = val
                new_tree.Fill()
                events[om] += 1

    file.Close()

    out_file.cd()
    new_tree.Write()
    out_file.Close()


def get_pulse_time(waveform, amplitude, baseline, peak):
    cell = 0
    for i in range(len(waveform[peak-20:peak])):
        j = i + peak - 20
        if (waveform[j] - baseline)*adc2mv < -1*amplitude*0.4:
            cell = j
            break
    return cell


def store_ap_info(n_record):
    templates = get_templates()
    output_file = open("/Users/williamquinn/Desktop/commissioning/ap_nums.txt", "w")

    # filename = "/Users/williamquinn/Desktop/commissioning/demon_ap_reduced.root"
    filename = "/Users/williamquinn/Desktop/commissioning/demon_ap.root"
    file = ROOT.TFile(filename, "READ")
    tree = file.T

    # out_dicts = [{"om": i, "events": {}} for i in range(712)]

    n_events = tree.GetEntries()
    events = [0 for i in range(712)]
    i_e = 0
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events, '{:.2f}%'.format(np.sum(events)/(712*n_record)*100))
        for index, om in enumerate(list(event.OM_ID)):
            if events[om] == n_record:
                continue
            om_t, om_s = om_type(om)
            eventn = event.event_num
            '''if eventn not in out_dicts[om]["events"].keys():
                out_dicts[om]["events"][eventn] = {}'''
            waveform = list(event.waveform)[index*1024:1024*(1 + index)]
            baseline = get_baseline(waveform, 20)
            amplitude = get_amplitude(waveform, baseline) * adc2mv * -1
            if amplitude < 300:
                continue
            peak = np.argmin(waveform)
            pulse_time = get_pulse_time(waveform, amplitude, baseline, peak)
            ap_info = mf_sweep([(val - baseline) * adc2mv for val in waveform], templates[om_t],
                               shape_cut=0.95, amp_cut=1)
            '''out_dicts[om]["events"][eventn] = {"amplitude": float(amplitude),
                                               "baseline": float(baseline),
                                               "pulse_time": pulse_time,
                                               "ap_times": [int(i) for i in ap_info.keys()],
                                               "ap_num": int(len(ap_info.keys()))}'''
            events[om] += 1
            # ap_nums[om].append(len(ap_info.keys()))
            output_file.write('{}:{}:{}:{}:{}:{}:{}\n'.format(om, eventn, amplitude, baseline, pulse_time,
                                                              len(ap_info.keys()), [val for val in ap_info.keys()]))

    output_file.close()


def read_ap_info():
    # file = open("/Users/williamquinn/Desktop/commissioning/ap_nums.csv", "r")
    file = open("/Users/williamquinn/Desktop/commissioning/ap_nums.txt", "r")
    fl = file.readlines()
    ap_info = {}
    for line in fl:
        # OM, Event Number, amplitude, baseline, pulse time, number of apulses, apulse_times
        line_list = line.split(":")
        if len(line_list) < 3:
            continue
        om = int(line_list[0])
        event_num = int(line_list[1])
        amplitude = float(line_list[2])
        baseline = float(line_list[3])
        pulse_time = int(line_list[4])
        ap_num = int(line_list[5].strip())
        if ap_num == 0:
            ap_times = None
        else:
            new_list = line_list[6].split("[")[-1].split("]")[0].split(",")
            ap_times = [int(i) for i in new_list]

        if om not in ap_info.keys():
            ap_info[om] = {"ap_num": [], "event_num": [], "amplitude": [], "baseline": [], "pulse_time": [],
                           "ap_times": []}
        ap_info[om]["event_num"].append(event_num)
        ap_info[om]["amplitude"].append(amplitude)
        ap_info[om]["ap_num"].append(ap_num)
        ap_info[om]["pulse_time"].append(pulse_time)
        ap_info[om]["ap_times"].append(ap_times)
        ap_info[om]["baseline"].append(baseline)

    file.close()
    return ap_info


def plot_pile_up_example():
    templates = get_templates()
    filename = "/Users/williamquinn/Desktop/commissioning/demon_ap.root"
    file = ROOT.TFile(filename, "READ")
    tree = file.T
    events = [0 for i in range(712)]
    i_e = 0
    n_events = tree.GetEntries()
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events, "aps:", events)
        for index, om in enumerate(list(event.OM_ID)):
            om_t, om_s = om_type(om)
            waveform = list(event.waveform)[index * 1024:1024 * (1 + index)]
            baseline = get_baseline(waveform, 20)
            amplitude = get_amplitude(waveform, baseline) * adc2mv * -1
            if amplitude < 300:
                continue
            else:
                # store_waveform((np.array(waveform) - baseline)*adc2mv, "demon_waveform", 0.64)
                ap_info = mf_sweep([(val - baseline) * adc2mv for val in waveform], templates[om_t], shape_cut=0.85,
                                   amp_cut=1)
                times = []
                for i in ap_info.keys():
                    times.append(i + 10)
                if len(ap_info.keys()) < 1:
                    # Too many after-pulses
                    continue
                for ap_time in ap_info.keys():
                    x = np.array([i for i in range(len(waveform))])
                    waveform = (np.array(waveform) - baseline) * adc2mv
                    ap_pulse = [val - baseline for index, val in enumerate(waveform[ap_time:ap_time + len(templates[om_t])])]
                    plt.figure(figsize=figsize)
                    plt.plot(x, waveform, 'C0.-', markersize=1.5)
                    plt.plot(x[times], waveform[times] , 'C2x', label='MF times')
                    plt.axhline(0, ls='--', color='C1', label='Baseline')
                    plt.xlabel("Timestamp /ns")
                    plt.ylabel("Voltage /mV")
                    plt.legend(loc='best')
                    plt.tight_layout()
                    plt.show()
                    events[om_t] += 1
                    return

    file.Close()


def plot_HV_map(plot):
    file = open("/Users/williamquinn/Desktop/commissioning/calorimeter_hv_211207.txt", "r")
    fl = file.readlines()
    HV_map = {}
    for line in fl:
        if '#' in line:
            continue
        line_list = line.split("\t")
        om = int(line_list[0].strip())
        '''om_id_str = line_list[0]
        om_typ = om_id_str.split(":")[0]
        if om_typ == 'M':
            side = int(om_id_str.split(":")[1].split(".")[0])
            col = int(om_id_str.split(":")[1].split(".")[1])
            row = int(om_id_str.split(":")[1].split(".")[2])
            om = get_om_num(om_typ, side=side, col=col, row=row, wall=None)
        elif om_typ == 'X':
            side = int(om_id_str.split(":")[1].split(".")[0])
            wall = int(om_id_str.split(":")[1].split(".")[1])
            col = int(om_id_str.split(":")[1].split(".")[2])
            row = int(om_id_str.split(":")[1].split(".")[3])
            om = get_om_num(om_typ, side=side, col=col, row=row, wall=wall)
        elif om_typ == 'X':
            side = int(om_id_str.split(":")[1].split(".")[0])
            wall = int(om_id_str.split(":")[1].split(".")[1])
            col = int(om_id_str.split(":")[1].split(".")[2])
            om = get_om_num(om_typ, side=side, col=col, row=None, wall=wall)
        else:
            continue'''
        HV = float(line_list[-1].strip())
        HV_map[om] = HV
    file.close()
    if plot:
        sncalo = sn.calorimeter("HV_map_march2021", with_palette=True)
        sncalo.draw_content = False
        sncalo.draw_omnum = False
        sncalo.draw_omid = False
        for om in range(712):
            if om not in HV_map.keys():
                continue
            if HV_map[om] == 0:
                continue
            sncalo.setcontent(om, HV_map[om])
        sncalo.draw()
        sncalo.save("/Users/williamquinn/Desktop/commissioning")
        del sncalo
    return HV_map


def plot_ap_times(ap_info, HV_map):
    ap_tot = [[], [], [], []]
    hvs = []
    for om in ap_info.keys():
        if om not in HV_map.keys():
            continue
        om_t, om_s = om_type(om)
        for index, aps in enumerate(ap_info[om]["ap_times"]):
            if ap_info[om]["ap_times"][index] is None:
                continue
            for ap in aps:
                new_ap = (ap - ap_info[om]["pulse_time"][index])/0.64 * np.sqrt(HV_map[om]/1000)
                ap_tot[om_t].append(new_ap)
                hvs.append(HV_map[om])

    labels = ['5" MW', '8" MW', '5" XW', '5" GV']
    colours = ['C0', 'k', 'C1', 'C2']
    alphas = [0.5, 1, 0.5, 0.5]
    z = [1, 10, 1, 1]
    plt.figure(figsize=figsize)
    for index in range(4):
        if len(ap_tot[index]) == 0:
            continue
        freq, bin_edges = np.histogram(ap_tot[index], bins=100, range=(0, 1024/0.64*np.sqrt(np.max(hvs)/1000)))
        width = bin_edges[2] - bin_edges[1]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='{}'.format(colours[index]), alpha=alphas[index],
                label=labels[index], zorder=z[index])
    # plt.ylabel("")
    plt.xlabel("After-pulse times /ns")
    plt.xlim(0, 1024/0.64*np.sqrt(np.max(hvs)/1000))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/om_type_time_dist.pdf")


def main():
    # The data I have used 521 and 522 need to be merged
    # merge_files()
    #
    # Get some templates using a bad template
    # create_templates()
    #
    # These new templates are still bad so need to zero them etc
    # alter_template()
    #
    # This saves some time but will reduce the file size mainly - 1.4M events (712*2000)
    # reduce_data(2000)
    #
    # This function will take a while to run so store the relevant data in a ascii file
    # store_ap_info(1000)
    #
    ap_info = read_ap_info()
    # print(ap_info)
    # plot_aan(ap_info)
    #
    # store_test_bench_waveform()
    # plot_waveform_comparisons()
    # plot_pile_up_example()
    HV_map = plot_HV_map(False)
    plot_ap_times(ap_info, HV_map)


if __name__ == "__main__":
    main()
    '''jobs = []
    for i in range(4):
        p = Process(target=main)
        jobs.append(p)
        p.start()'''
