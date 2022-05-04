from format_plot import *
import sndisplay as sn


def read_short_map(filename):
    file = open(filename, "r")
    fl = file.readlines()
    shorts = {"top": {}, "bottom": {}}
    for index in range(len(fl)):
        if index == 0:
            continue
        line_list = fl[index].split(",")
        side = int(line_list[0])
        end = int(line_list[1])
        pp = line_list[2].strip()
        pp_pos = int(line_list[3])
        pins = []
        for i in line_list[4:]:
            if i == '' or i == '\n':
                pins.append(0)
            else:
                pins.append(int(i))
        data = {"side": side, "end": end, "pp": pp, "pp_pos": pp_pos, "pins": pins}
        info = get_info(data)
        if info is None:
            continue
        for i in range(len(info["cells"])):
            shorts[info["tob"][i]][info["cells"][i]] = info["short"][i]
        # print(side, end, pp, pp_pos, pins)
    file.close()
    return shorts


def get_row(side, end, pp, pp_pos):
    spare = False
    special_cassette = False
    smallest_row = None
    if side == 1:
        if end == 0:
            if pp == 'C':
                if pp_pos < 11:
                    if pp_pos in [1]:
                        spare = True
                    elif pp_pos == 2:
                        special_cassette = True
                        smallest_row = 56
                    else:
                        i = pp_pos - 3
                        smallest_row = 54 - i * 2
                else:
                    if pp_pos in [21, 20, 19, 18]:
                        spare = True
                    elif pp_pos == 17:
                        smallest_row = 52
                    elif pp_pos == 16:
                        special_cassette = True
                        smallest_row = 56
                    else:
                        i = pp_pos - 11
                        smallest_row = 46 + (i + 1) % 2 * (2 * i) + i % 2 * (-2 + 2 * (i - 1))
            elif pp == 'D':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 38 - i * 2
                elif pp_pos == 11:
                    smallest_row = 20
                else:
                    i = pp_pos - 12
                    smallest_row = 26 + (i+1) % 2 * (2*i) + i % 2 * (-2 + 2*(i - 1))
            elif pp == 'E':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 18 - i * 2
                else:
                    i = pp_pos - 11
                    smallest_row = 2 + (i+1) % 2 * (2*i) + i % 2 * (-2 + 2*(i - 1))
            else:
                raise ValueError('pp = {}'.format(pp))
        elif end == 1:
            if pp == 'C':
                if pp_pos < 11:
                    if pp_pos in [1, 2]:
                        spare = True
                    else:
                        i = pp_pos - 3
                        smallest_row = 57 + 2*i
                else:
                    if pp_pos in [20, 19, 18, 17]:
                        spare = True
                    else:
                        i = pp_pos - 16
                        smallest_row = 57 - 2*i
            elif pp == 'D':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 73 + 2 * i
                else:
                    i = pp_pos - 21
                    smallest_row = 69 - 2 * i
            elif pp == 'E':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 93 + 2 * i
                else:
                    i = pp_pos - 21
                    smallest_row = 91 - 2 * i
            else:
                raise ValueError('pp = {}'.format(pp))
        else:
            raise ValueError('end = {}'.format(end))
    elif side == 0:
        if end == 0:
            if pp == 'C':
                if pp_pos < 11:
                    if pp_pos in [1, 2]:
                        spare = True
                    else:
                        i = pp_pos - 3
                        smallest_row = 54 - 2*i
                else:
                    if pp_pos in [20, 19, 18, 17]:
                        spare = True
                    else:
                        i = pp_pos - 16
                        smallest_row = 54 + 2*i
            elif pp == 'D':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 38 - 2 * i
                else:
                    i = pp_pos - 21
                    smallest_row = 42 + 2 * i
            elif pp == 'E':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 18 - 2 * i
                else:
                    i = pp_pos - 21
                    smallest_row = 20 + 2 * i
            else:
                raise ValueError('pp = {}'.format(pp))
        elif end == 1:
            if pp == 'C':
                if pp_pos < 11:
                    if pp_pos in [1]:
                        spare = True
                    elif pp_pos == 2:
                        special_cassette = True
                        smallest_row = 56
                    else:
                        i = pp_pos - 3
                        smallest_row = 57 + 2*i
                else:
                    if pp_pos in [21, 20, 19, 18]:
                        spare = True
                    else:
                        i = pp_pos - 11
                        smallest_row = 65 - (i + 1) % 2 * (2 * i) - i % 2 * (-2 + 2 * (i - 1))
            elif pp == 'D':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 73 + 2 * i
                elif pp_pos == 11:
                    smallest_row = 91
                else:
                    i = pp_pos - 12
                    smallest_row = 85 - (i + 1) % 2 * (2 * i) - i % 2 * (-2 + 2 * (i - 1))
            elif pp == 'E':
                if pp_pos < 11:
                    i = pp_pos - 1
                    smallest_row = 93 + 2 * i
                else:
                    i = pp_pos - 11
                    smallest_row = 109 - (i + 1) % 2 * (2 * i) - i % 2 * (-2 + 2 * (i - 1))
            else:
                print(pp, type(pp), pp=='C', len(pp))
                raise ValueError('pp = {}'.format(pp))
        else:
            raise ValueError('end = {}'.format(end))
    else:
        raise ValueError('side = {}'.format(side))

    if pp_pos < 11:
        tob = 'bottom'
    else:
        tob = 'top'

    return smallest_row, special_cassette, spare, tob


def get_info(data):
    side = data["side"]
    end = data["end"]
    pp = data["pp"]
    pp_pos = data["pp_pos"]
    pins = data["pins"]
    small_row, special, spare, tob = get_row(side, end, pp, pp_pos)
    output = {"cells": [], "tob": [], "short": []}
    if special:
        '''output["cells"].append(cell_num)
        output["tob"].append(tob)
        output["short"].append(pin)'''
        return None
    if spare:
        return None
    if small_row is None:
        return None
    for index, pin in enumerate(pins):
        if index == 18:
            continue
        if side == 1:
            if index < 9:
                layer = 8 - index
                row = small_row
            elif 9 <= index < 18:
                layer = 8 - (index - 9)
                row = small_row + 1
        elif side == 0:
            if index < 9:
                layer = 8 - index
                row = small_row + 1
            elif 9 <= index < 18:
                layer = 8 - (index - 9)
                row = small_row
        else:
            pass
        cell_num = side * 113 * 9 + row * 9 + layer
        output["cells"].append(cell_num)
        output["tob"].append(tob)
        output["short"].append(pin)
    return output


def plot_short_map(short_map, name):
    data_top = short_map["top"]
    data_bot = short_map["bottom"]
    new_data = {}

    for cell_num in data_top.keys():
        if cell_num not in data_bot.keys():
            # This shouldn't happen
            raise ValueError
        top = bool(data_top[cell_num])
        bot = bool(data_bot[cell_num])
        if top and bot:
            new_data[cell_num] = 3
        elif top and not bot:
            new_data[cell_num] = 2
        elif not top and bot:
            new_data[cell_num] = 1
        else:
            new_data[cell_num] = 0

    sn_tracker = sn.tracker("short_" + name, with_palette=True)
    sn_tracker.draw_content = False
    sn_tracker.draw_cellnum = False
    sn_tracker.draw_cellid = False

    for cell_num in new_data.keys():
        # cell_num = side * 113 * 9 + row * 9 + layer
        if cell_num == 0 * 113 * 9 + 108 * 9 + 4:
            continue
        elif cell_num == 0 * 113 * 9 + 89 * 9 + 1:
            continue
        sn_tracker.setcontent(cell_num, new_data[cell_num])
    sn_tracker.draw()
    sn_tracker.save("/Users/williamquinn/Desktop/tracksioning")
    del sn_tracker


def main():
    short_map_before = read_short_map("/Users/williamquinn/Desktop/tracksioning/shorts_before.csv")
    short_map_after = read_short_map("/Users/williamquinn/Desktop/tracksioning/shorts_after.csv")
    plot_short_map(short_map_before, "200908_before")
    plot_short_map(short_map_after, "201124_after")


if __name__ == "__main__":
    main()
