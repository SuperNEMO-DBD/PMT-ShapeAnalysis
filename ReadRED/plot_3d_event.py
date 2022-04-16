import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def om_id(omnum: int):
    id_ = {"side": None, "wall": None, "col": None, 'row': None}
    if 0 <= omnum < 260:
        id_['side'] = 0
        id_['col'] = omnum // 13
        id_['row'] = omnum % 13
    elif omnum < 520:
        omnum = omnum - 260
        id_['side'] = 1
        id_['col'] = omnum // 13
        id_['row'] = omnum % 13
    elif omnum < 584:
        omnum = omnum - 520
        id_['side'] = 0
        id_['wall'] = omnum // 32
        id_['col'] = omnum % 32 // 16
        id_['row'] = omnum % 32 % 16
    elif omnum < 648:
        omnum = omnum - 520 - 64
        id_['side'] = 1
        id_['wall'] = omnum // 32
        id_['col'] = omnum % 32 // 16
        id_['row'] = omnum % 32 % 16
    elif omnum < 680:
        omnum = omnum - 520 - 128
        id_['side'] = 0
        id_['wall'] = omnum % 16
        id_['col'] = omnum // 16
    elif omnum < 712:
        omnum = omnum - 520 - 128 - 32
        id_['side'] = 1
        id_['wall'] = omnum % 16
        id_['col'] = omnum // 16
    return id_


def cell_id(cellnum):
    cell_side = cellnum // (9 * 113)
    cell_row = cellnum % (9 * 113) // 9
    cell_layer = cellnum % (9 * 113) % 9

    return cell_side, cell_row, cell_layer


class coord:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        if isinstance(other, coord):
            X = [self.x, other.x]
            Y = [self.y, other.y]
            Z = [self.z, other.z]

            return coord(X, Y, Z)
        else:
            raise TypeError


def main():
    oms = {}
    trs = {}
    file = open("/Users/williamquinn/Desktop/read_red/event_88.csv")
    fl = file.readlines()
    for index, line in enumerate(fl):
        line_list = line.split(",")
        if line_list[0] == "OM":
            oms[int(line_list[1].strip())] = None
        elif line_list[0] == "GG":
            trs[int(line_list[1].strip())] = float(line_list[2].strip())

    for om in oms.keys():
        oms[om] = om_id(om)

    for cell in trs.keys():
        pos = trs[cell]
        side, row, layer = cell_id(cell)
        trs[cell] = {"side": side, "row": row, "layer": layer, "pos": pos}

    demo_height = 200 # cm
    demo_length = 113 * 5
    cell_diameter = 5
    om_dim = {"height": demo_height/13, "width": 10, "length":demo_length/20}
    demo_width = (18*cell_diameter) + om_dim["width"]*2

    # Plot figure
    fig = plt.figure(figsize=(12,6), facecolor='grey')
    ax = fig.add_subplot(111, projection='3d')

    for om in range(520):
        blank_oms = om_id(om)
        xl = blank_oms["col"] * om_dim["length"]
        xh = xl + om_dim["length"]
        yl = blank_oms["side"] * (18*cell_diameter + om_dim["width"])
        yh = yl + om_dim["width"]
        zl = blank_oms["row"] * om_dim["height"]
        zh = zl + om_dim["height"]

        A = coord(xl, yl, zl)
        B = coord(xh, yl, zl)
        C = coord(xh, yh, zl)
        D = coord(xl, yh, zl)
        E = coord(xl, yl, zh)
        F = coord(xh, yl, zh)
        G = coord(xh, yh, zh)
        H = coord(xl, yh, zh)
        '''X = np.array([xl, xh, xh, xl, xl, xl, xl, xl, xl, xh, xh, xh, xh, xh])
        Y = np.array([yl, yl, yh, yh, yl, yl, yh, yh, yh, yh, yh, yh, yl, yl])
        Z = np.array([zl, zl, zl, zl, zl, zh, zh, zl, zh, zh, zl, zh, zh, zl])'''
        ax.plot((A + B).x, (A + B).y, (A + B).z, color='grey', alpha=0.2)
        ax.plot((B + C).x, (B + C).y, (B + C).z, color='grey', alpha=0.2)
        ax.plot((C + D).x, (C + D).y, (C + D).z, color='grey', alpha=0.2)
        ax.plot((D + A).x, (D + A).y, (D + A).z, color='grey', alpha=0.2)
        ax.plot((E + H).x, (E + H).y, (E + H).z, color='grey', alpha=0.2)
        ax.plot((H + G).x, (H + G).y, (H + G).z, color='grey', alpha=0.2)
        ax.plot((G + F).x, (G + F).y, (G + F).z, color='grey', alpha=0.2)
        ax.plot((F + E).x, (F + E).y, (F + E).z, color='grey', alpha=0.2)
        ax.plot((A + E).x, (A + E).y, (A + E).z, color='grey', alpha=0.2)
        ax.plot((B + F).x, (B + F).y, (B + F).z, color='grey', alpha=0.2)
        ax.plot((G + C).x, (G + C).y, (G + C).z, color='grey', alpha=0.2)
        ax.plot((H + D).x, (H + D).y, (H + D).z, color='grey', alpha=0.2)

    for om in oms.keys():
        full_oms = oms[om]
        xl = full_oms["col"] * om_dim["length"]
        xh = xl + om_dim["length"]
        yl = full_oms["side"] * (18*cell_diameter + om_dim["width"])
        yh = yl + om_dim["width"]
        zl = full_oms["row"] * om_dim["height"]
        zh = zl + om_dim["height"]
        aplha = 0.9
        # Plot the surfaces #6
        # x
        ax.add_collection3d(plt.fill_between(x=[xl, xh], y1=yl, y2=yh, color='red', alpha=aplha), zl, zdir='z')
        ax.add_collection3d(plt.fill_between(x=[xl, xh], y1=yl, y2=yh, color='red', alpha=aplha), zh, zdir='z')
        # y
        ax.add_collection3d(plt.fill_between(x=[xl, xh], y1=zl, y2=zh, color='red', alpha=aplha), yl, zdir='y')
        ax.add_collection3d(plt.fill_between(x=[xl, xh], y1=zl, y2=zh, color='red', alpha=aplha), yh, zdir='y')
        # z
        ax.add_collection3d(plt.fill_between(x=[yl, yh], y1=zl, y2=zh, color='red', alpha=aplha), xl, zdir='x')
        ax.add_collection3d(plt.fill_between(x=[yl, yh], y1=zl, y2=zh, color='red', alpha=aplha), xh, zdir='x')

    xs = []
    ys = []
    zs = []
    for index, cell in enumerate(trs.keys()):
        if index == 1:
            continue
        x = trs[cell]["row"]*5 + 2.5
        if trs[cell]["side"] == 0:
            trs[cell]["layer"] = 8 - trs[cell]["layer"]
        y = om_dim["width"] + (trs[cell]["side"] * 9 * cell_diameter) + trs[cell]["layer"]*cell_diameter + 2.5
        z = trs[cell]["pos"]*100
        xs.append(x)
        ys.append(y)
        zs.append(z)
        print(trs[cell], x, y, z)
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    ax.scatter3D(xs, ys, zs, marker='o', color='b')
    ax.set(xlabel="x", ylabel='y', zlabel='z')
    ax.add_collection3d(plt.fill_between(x=[0, demo_length], y1=0, y2=demo_height, color='grey', alpha=0.1), demo_width/2, zdir='y')

    # Voxels is used to customizations of the
    # sizes, positions and colors.
    plt.tight_layout()
    ax.view_init(20, 140)
    # rotate the axes and update
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('grey')
    ax.yaxis.pane.set_edgecolor('grey')
    ax.zaxis.pane.set_edgecolor('grey')
    plt.savefig("event_88_3d.pdf")


if __name__ == "__main__":
    main()
