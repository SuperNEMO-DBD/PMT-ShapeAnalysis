import ROOT
import numpy as np
from random import gauss


class palette:
    def __init__(self):
        nRGBs = 6
        stops = np.array([0.00, 0.20, 0.40, 0.60, 0.80, 1.00], dtype='float64')
        red = np.array([0.25, 0.00, 0.20, 1.00, 1.00, 0.90], dtype='float64')
        green = np.array([0.25, 0.80, 1.00, 1.00, 0.80, 0.00], dtype='float64')
        blue = np.array([1.00, 1.00, 0.20, 0.00, 0.00, 0.00], dtype='float64')
        self.palette_index = ROOT.TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, 100)

    def get_index(self):
        return self.palette_index


'''// // // // // // // // // // // // // //
// sndisplay::calorimeter //
// // // // // // // // // // // // // //'''


class calorimeter:
    def __init__(self, new_name: str, with_palette: bool):
        self.palette_index = palette()
        self.name = new_name
        # self.canvas_it = ROOT.TCanvas(f'C_it_{self.name}', self.name, 900, 600)
        # self.canvas_fr = ROOT.TCanvas(f'C_fr_{self.name}', self.name, 900, 600)
        self.canvas_it = None
        self.canvas_fr = None
        # TODO: check this
        self.nmwall = 520
        self.nxwall = 128
        self.ngveto = 64
        self.nb_om = 712
        self.content_err = []
        self.content = []
        self.ombox = []
        self.omid_text_v = []
        self.omnum_text_v = []
        self.content_err_text_v = []
        self.content_text_v = []

        self.draw_omid = True
        self.draw_omnum = False
        self.draw_content_err = False
        self.draw_content = False
        self.draw_content_err_format = '± {:.0f}'
        self.draw_content_format = '{:.0f}'

        self.has_italy_data = False
        self.has_french_data = False

        for omnum in range(self.nb_om):
            self.content.append(None)
            self.content_err.append(None)

        self.range_min = -1
        self.range_max = -1

        spacery = 0.0125
        if with_palette:
            spacerx = 0.0093458
            mw_sizex = (1 - 5 * spacerx) / (20 + 4 + 1.5)
        else:
            spacerx = 0.0100
            mw_sizex = (1-4*spacerx)/(20+4)

        mw_sizey = (1 - 4 * spacery) / (13 + 2)
        gv_sizey = mw_sizey
        xw_sizey = mw_sizey * 13. / 16.

        gv_sizex = mw_sizex * 20. / 16.
        xw_sizex = mw_sizex


        '''// // // // // // // // // // // // //
        // MWALL
        initialisation //
        // // // // // // // // // // // // //'''

        for mw_side in range(2):

            for mw_column in range(20):

                for mw_row in range(13):
                    omnum = mw_side * 20 * 13 + mw_column * 13 + mw_row

                    x1 = spacerx + 2 * xw_sizex + spacerx

                    if mw_side == 0:
                        x1 += mw_sizex * (19 - mw_column)
                    else:
                        x1 += mw_sizex * mw_column

                    y1 = spacery + gv_sizey + spacery + mw_sizey * mw_row
                    x2 = x1 + mw_sizex
                    y2 = y1 + mw_sizey

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.ombox.append(box)

                    omid_string = f'M:{mw_side}.{mw_column}.{mw_row}'
                    omid_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.7 * mw_sizey, omid_string)
                    omid_text.SetTextSize(0.014)
                    omid_text.SetTextAlign(22)
                    self.omid_text_v.append(omid_text)

                    omnum_string = f'{omnum}'
                    omnum_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.7 * mw_sizey, omnum_string)
                    omnum_text.SetTextFont(42)
                    omnum_text.SetTextSize(0.013)
                    omnum_text.SetTextAlign(22)
                    self.omnum_text_v.append(omnum_text)

                    content_err_text_string = '±'
                    content_err_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.6 * mw_sizey, content_err_text_string)
                    content_err_text.SetTextSize(0.014)
                    content_err_text.SetTextAlign(22)
                    self.content_err_text_v.append(content_err_text)

                    content_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.3 * mw_sizey, "")
                    content_text.SetTextSize(0.02)
                    content_text.SetTextAlign(22)
                    self.content_text_v.append(content_text)

        '''// // // // // // // // // // // // //
        // XWALL initialisation //
        // // // // // // // // // // // // //'''

        for xw_side in range(2):

            for xw_wall in range(2):

                for xw_column in range(2):

                    for xw_row in range(16):

                        omnum = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row

                        if xw_side == 0:
                            if xw_wall == 0:
                                x1 = spacerx + 2 * xw_sizex + spacerx + 20 * mw_sizex + spacerx + (
                                            1 - xw_column) * xw_sizex
                            else:
                                x1 = spacerx + xw_sizex * xw_column
                        else:
                            if xw_wall == 0:
                                x1 = spacerx + xw_sizex * xw_column
                            else:
                                x1 = spacerx + 2 * xw_sizex + spacerx + 20 * mw_sizex + spacerx + (
                                            1 - xw_column) * xw_sizex

                        x2 = x1 + xw_sizex
                        y1 = spacery + gv_sizey + spacery + xw_sizey * xw_row
                        y2 = spacery + gv_sizey + spacery + xw_sizey * (xw_row + 1)

                        box = ROOT.TBox(x1, y1, x2, y2)
                        box.SetFillColor(0)
                        box.SetLineWidth(1)
                        self.ombox.append(box)

                        omid_string = f'X:{xw_side}.{xw_wall}.{xw_column}.{xw_row}'
                        omid_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.7 * mw_sizey, omid_string)
                        omid_text.SetTextFont(42)
                        omid_text.SetTextSize(0.013)
                        omid_text.SetTextAlign(22)
                        self.omid_text_v.append(omid_text)

                        omnum_string = f'{omnum}'
                        omnum_text = ROOT.TText(x1 + 0.5 * xw_sizex, y1 + 0.333 * xw_sizey, omnum_string)
                        omnum_text.SetTextFont(42)
                        omnum_text.SetTextSize(0.013)
                        omnum_text.SetTextAlign(22)
                        self.omnum_text_v.append(omnum_text)

                        content_err_text_string = '±'
                        content_err_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.6 * mw_sizey, content_err_text_string)
                        content_err_text.SetTextSize(0.014)
                        content_err_text.SetTextAlign(22)
                        self.content_err_text_v.append(content_err_text)

                        content_text = ROOT.TText(x1 + 0.5 * xw_sizex, y1 + 0.3 * xw_sizey, "")
                        content_text.SetTextSize(0.02)
                        content_text.SetTextAlign(22)
                        self.content_text_v.append(content_text)

        '''// // // // // // // // // // // // //
        // GVETO initialisation //
        // // // // // // // // // // // // //'''

        for gv_side in range(2):

            for gv_wall in range(2):

                for gv_column in range(16):

                    omnum = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column

                    if gv_side == 0:
                        x1 = spacerx + 2 * xw_sizex + spacerx + gv_sizex * (16 - 1 - gv_column)
                    else:

                        x1 = spacerx + 2 * xw_sizex + spacerx + gv_sizex * gv_column

                    x2 = x1 + gv_sizex
                    y1 = spacery + gv_wall * (gv_sizey + spacery + 13 * mw_sizey + spacery)
                    y2 = y1 + gv_sizey

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.ombox.append(box)

                    omid_string = f'G:{gv_side}.{gv_wall}.{gv_column}'
                    omid_text = ROOT.TText(x1 + 0.5 * gv_sizex, y1 + 0.7 * gv_sizey, omid_string)
                    omid_text.SetTextFont(42)
                    omid_text.SetTextSize(0.013)
                    omid_text.SetTextAlign(22)
                    self.omid_text_v.append(omid_text)

                    omnum_string = f'{omnum}'
                    omnum_text = ROOT.TText(x1 + 0.5 * gv_sizex, y1 + 0.7 * gv_sizey, omnum_string)
                    omnum_text.SetTextFont(42)
                    omnum_text.SetTextSize(0.02)
                    omnum_text.SetTextAlign(22)
                    self.omnum_text_v.append(omnum_text)

                    content_err_text_string = '±'
                    content_err_text = ROOT.TText(x1 + 0.5 * mw_sizex, y1 + 0.6 * mw_sizey, content_err_text_string)
                    content_err_text.SetTextSize(0.014)
                    content_err_text.SetTextAlign(22)
                    self.content_err_text_v.append(content_err_text)

                    content_text = ROOT.TText(x1 + 0.5 * gv_sizex, y1 + 0.3 * gv_sizey, "")
                    content_text.SetTextSize(0.02)
                    content_text.SetTextAlign(22)
                    self.content_text_v.append(content_text)

        self.it_label = ROOT.TText(spacerx+xw_sizex, spacery + gv_sizey + spacery + 13 * mw_sizey + spacery + 0.5 * gv_sizey,
                                   "  ITALY")
        self.it_label.SetTextSize(0.028)
        self.it_label.SetTextAlign(22)

        self.fr_label = ROOT.TText(spacerx+xw_sizex, spacery + gv_sizey + spacery + 13 * mw_sizey + spacery + 0.5 * gv_sizey,
                                   "FRANCE")
        self.fr_label.SetTextSize(0.028)
        self.fr_label.SetTextAlign(22)

        self.palette_histo = None
        self.palette_axis = None
        if with_palette:
            palette_sizey = mw_sizey * 13
            palette_sizex = mw_sizex * 1.5

            palette_histo = ROOT.TH2D("{}_palette_histo".format(self.name), "", 1, 0, 1, 1, 0, 1)
            palette_histo.GetZaxis().SetNdivisions(509)
            palette_histo.GetZaxis().SetLabelSize(0.024)
            palette_histo.GetZaxis().SetLabelFont(62)
            palette_histo.SetMinimum(self.range_min)
            palette_histo.SetMaximum(self.range_max)
            palette_histo.SetContour(100)
            self.palette_histo = palette_histo

            palette_axis = ROOT.TPaletteAxis()
            palette_axis.SetHistogram(palette_histo)
            palette_axis.SetX1NDC(1-spacerx-palette_sizex * (5. / 6))
            palette_axis.SetY1NDC(spacery+gv_sizey+spacery+palette_sizey / 8.)
            palette_axis.SetX2NDC(1-spacerx-palette_sizex * (3. / 6))
            palette_axis.SetY2NDC(1-palette_axis.GetY1NDC())
            palette_axis.SetY2NDC(1-palette_axis.GetY1NDC())
            self.palette_axis = palette_axis

    def setrange(self, xmin: float, xmax: float):
        self.range_min = xmin
        self.range_max = xmax

    def draw_omid_label(self):
        self.draw_omid = True
        self.draw_omnum = False

    def draw_omnum_label(self):
        self.draw_omnum = True
        self.draw_omid = False

    def draw_content_err_label(self, string: str):
        self.draw_content_err_format = string
        self.draw_content_err = True

    def draw_content_label(self, string: str):
        self.draw_content_format = string
        self.draw_content = True

    def draw(self):

        if self.palette_axis is not None:
            canvas_width = 1284
        else:
            canvas_width = 1200
        canvas_height = 780

        if self.canvas_it is None:
            self.canvas_it = ROOT.TCanvas("{}_canvas_it".format(self.name), "{} (IT side)".format(self.name),
                                          canvas_width, canvas_height)
            decoration_width = canvas_width - self.canvas_it.GetWw()
            decoration_height = canvas_height - self.canvas_it.GetWh()
            self.canvas_it.SetWindowSize(canvas_width + decoration_width, canvas_height + decoration_height)
            self.canvas_it.SetFixedAspectRatio()

        if self.canvas_fr is None:
            self.canvas_fr = ROOT.TCanvas("{}_canvas_fr".format(self.name), "{} (FR side)".format(self.name),
                                          canvas_width, canvas_height)
            decoration_width = canvas_width - self.canvas_fr.GetWw()
            decoration_height = canvas_height - self.canvas_fr.GetWh()
            self.canvas_fr.SetWindowSize(canvas_width + decoration_width, canvas_height + decoration_height)
            self.canvas_fr.SetFixedAspectRatio()

        if self.draw_content:
            for omnum in range(self.nb_om):
                if self.content[omnum] is not None:
                    ttext = self.content_text_v[omnum]
                    ttext.SetText(ttext.GetX(), ttext.GetY(), self.draw_content_format.format(self.content[omnum]))

        if self.draw_content_err:
            for omnum in range(self.nb_om):
                if self.content[omnum] is not None:
                    ttext = self.content_err_text_v[omnum]
                    ttext.SetText(ttext.GetX(), ttext.GetY(),
                                  self.draw_content_err_format.format(self.content_err[omnum]))

        '''// // // // // // /
        // Draw IT //
        // // // // // // /'''

        self.canvas_it.cd()
        self.canvas_it.SetEditable(True)

        mw_side = 0
        for mw_column in range(20):
            for mw_row in range(13):
                id_ = mw_side * 20 * 13 + mw_column * 13 + mw_row
                self.ombox[id_].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id_].Draw()
                elif self.draw_omnum:
                    self.omnum_text_v[id_].Draw()
                if self.draw_content_err and self.content_err[id_] is not None:
                    self.content_err_text_v[id_].Draw()
                if self.draw_content and self.content[id_] is not None:
                    self.content_text_v[id_].Draw()

        xw_side = 0
        for xw_wall in range(2):
            for xw_column in range(2):
                for xw_row in range(16):
                    id_ = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row
                    self.ombox[id_].Draw("l")
                    if self.draw_omid:
                        self.omid_text_v[id_].Draw()
                    if self.draw_omnum:
                        self.omnum_text_v[id_].Draw()
                    if self.draw_content_err and self.content_err[id_] is not None:
                        self.content_err_text_v[id_].Draw()
                    if self.draw_content and self.content[id_] is not None:
                        self.content_text_v[id_].Draw()

        gv_side = 0
        for gv_wall in range(2):
            for gv_column in range(16):
                id_ = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column
                self.ombox[id_].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id_].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id_].Draw()
                if self.draw_content_err and self.content_err[id_] is not None:
                    self.content_err_text_v[id_].Draw()
                if self.draw_content and self.content[id_] is not None:
                    self.content_text_v[id_].Draw()

        self.it_label.Draw()
        if self.palette_axis is not None:
            self.palette_axis.Draw()
        self.canvas_it.SetEditable(False)

        '''// // // // // // /
        // Draw FR //
        // // // // // // /'''

        self.canvas_fr.cd()
        self.canvas_fr.SetEditable(True)

        mw_side = 1
        for mw_column in range(20):
            for mw_row in range(13):
                id_ = mw_side * 20 * 13 + mw_column * 13 + mw_row
                self.ombox[id_].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id_].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id_].Draw()
                if self.draw_content_err and self.content_err[id_] is not None:
                    self.content_err_text_v[id_].Draw()
                if self.draw_content and self.content[id_] is not None:
                    self.content_text_v[id_].Draw()

        xw_side = 1
        for xw_wall in range(2):
            for xw_column in range(2):
                for xw_row in range(16):
                    id_ = 520 + xw_side * 2 * 2 * 16 + xw_wall * 2 * 16 + xw_column * 16 + xw_row
                    self.ombox[id_].Draw("l")
                    if self.draw_omid:
                        self.omid_text_v[id_].Draw()
                    if self.draw_omnum:
                        self.omnum_text_v[id_].Draw()
                    if self.draw_content_err and self.content_err[id_] is not None:
                        self.content_err_text_v[id_].Draw()
                    if self.draw_content and self.content[id_] is not None:
                        self.content_text_v[id_].Draw()

        gv_side = 1
        for gv_wall in range(2):
            for gv_column in range(16):
                id_ = 520 + 128 + gv_side * 2 * 16 + gv_wall * 16 + gv_column
                self.ombox[id_].Draw("l")
                if self.draw_omid:
                    self.omid_text_v[id_].Draw()
                if self.draw_omnum:
                    self.omnum_text_v[id_].Draw()
                if self.draw_content_err and self.content_err[id_] is not None:
                    self.content_err_text_v[id_].Draw()
                if self.draw_content and self.content[id_] is not None:
                    self.content_text_v[id_].Draw()

        self.fr_label.Draw()
        if self.palette_axis is not None:
            self.palette_axis.Draw()
        self.canvas_fr.SetEditable(False)

        self.update()

    def reset(self):
        for omnum in range(self.nb_om):
            self.content[omnum] = None
            self.content_text_v[omnum].Clear()
            self.ombox[omnum].SetFillColor(14)

        self.canvas_it.Modified()
        self.canvas_it.Update()
        self.canvas_fr.Modified()
        self.canvas_fr.Update()
        ROOT.gSystem.ProcessEvents()

    def getcontent(self, omnum: int):
        return self.content[omnum]

    def setcontent(self, omnum: int, value: float):
        if 0 <= omnum < 260:
            self.has_italy_data = True
        elif omnum < 520:
            self.has_french_data = True
        elif omnum < 584:
            self.has_italy_data = True
        elif omnum < 648:
            self.has_french_data = True
        elif omnum < 680:
            self.has_italy_data = True
        elif omnum < 712:
            self.has_french_data = True

        self.content[omnum] = value

    def set_err_content(self, omnum: int, value: float):
        if 0 <= omnum < 260:
            self.has_italy_data = True
        elif omnum < 520:
            self.has_french_data = True
        elif omnum < 584:
            self.has_italy_data = True
        elif omnum < 648:
            self.has_french_data = True
        elif omnum < 680:
            self.has_italy_data = True
        elif omnum < 712:
            self.has_french_data = True

        self.content_err[omnum] = value

    def setcontent_(self, om_side: int, om_wall: int, om_column: int, om_row: int, value: float):
        omnum = -1

        # auto detect MW
        if om_side != -1 and om_wall == -1 and om_column != -1 and om_row != -1:
            omnum = 260 * om_side + 13 * om_column + om_row
        # auto detect XW
        elif om_side != -1 and om_wall != -1 and om_column != -1 and om_row != -1:
            omnum = 520 + 64 * om_side + 32 * om_wall + 16 * om_column + om_row
        # auto detect GV
        elif om_side != -1 and om_wall != -1 and om_column != -1 and om_row == -1:
            omnum = 520 + 128 + 32 * om_side + 16 * om_wall + om_column
        else:
            print(f'+++ sndisplay: skipping OM ({om_side}.{om_wall}.{om_column}.{om_row})\n')
            return

        if om_side == 0:
            self.has_italy_data = True
        elif om_side == 1:
            self.has_french_data = True

        self.content[omnum] = value

    def set_colour(self, omnum: int, colour):
        self.ombox[omnum].SetFillColor(colour)

    def fill(self, omnum: int):
        if self.content[omnum] is None:
            val = 1
        else:
            val = self.content[omnum] + 1
        self.setcontent(omnum, val)

    def update(self):

        for i_om in range(self.nb_om):
            if self.content[i_om] is not None:
                content_min = self.content[i_om]
                content_max = self.content[i_om]
                break

        if self.range_min == -1:
            for omnum in range(1, self.nb_om):
                if self.content[omnum] is None:
                    continue
                if self.content[omnum] < content_min:
                    content_min = self.content[omnum]
                if self.content[omnum] > content_max:
                    content_max = self.content[omnum]
        else:
            content_min = self.range_min
            content_max = self.range_max

        for omnum in range(self.nb_om):
            if self.content[omnum] is not None:
                colour_index = int(99*(self.content[omnum] - content_min)/(content_max - content_min))
                self.ombox[omnum].SetFillColor(self.palette_index.get_index() + int(colour_index))
            else:
                self.ombox[omnum].SetFillColor(14)

        if self.palette_axis is not None:
            self.palette_histo.SetMinimum(content_min)
            self.palette_histo.SetMaximum(content_max)
            self.palette_histo.SetContour(100)

        self.canvas_it.Modified()
        self.canvas_it.Update()

        self.canvas_fr.Modified()
        self.canvas_fr.Update()

        ROOT.gSystem.ProcessEvents()

    def save(self, location: str):
        self.canvas_it.SaveAs(location + "/" + self.name + '_it.pdf')
        self.canvas_fr.SaveAs(location + "/" + self.name + '_fr.pdf')


'''// // // // // // // // // // // //
// sndisplay::tracker //
// // // // // // // // // // // //'''


class tracker:
    def __init__(self, new_name: str, with_palette):
        self.palette_index = palette()
        self.name = new_name
        self.canvas = None
        self.draw_cellid = False
        self.draw_cellnum = False
        self.draw_content = False
        self.draw_content_format = '{:.0f}'  # "%.0f"
        self.range_min = -1
        self.range_max = -1
        self.nb_cell = 2034
        self.content = [None for i in range(self.nb_cell)]
        self.spacerx = 0.0082
        self.spacery = 0.0500
        if with_palette:
            self.cell_sizex = (1 - 3 * self.spacerx) / (113. + 2)
        else:
            self.cell_sizex = (1 - 2 * self.spacerx) / 113.0
        self.cell_sizey = (1 - 3 * self.spacery) / (9 * 2)
        self.cellbox = []
        self.cellid_text_v = []
        self.cellnum_text_v = []
        self.content_text_v = []
        self.row_text_v = []

        '''// // // // // // // // // // // // //
        // CELLS initialisation //
        // // // // // // // // // // // // //'''

        for cell_side in range(2):
            for cell_row in range(113):
                x1 = self.spacerx + cell_row * self.cell_sizex
                for cell_layer in range(9):
                    cellnum = cell_side * 113 * 9 + cell_row * 9 + cell_layer

                    y1 = self.spacery

                    if cell_side == 0:
                        y1 += 9 * self.cell_sizey + self.spacery + cell_layer * self.cell_sizey
                    else:
                        y1 += (8 - cell_layer) * self.cell_sizey

                    x2 = x1 + self.cell_sizex
                    y2 = y1 + self.cell_sizey
                    width = x2 - x1

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.cellbox.append(box)

                    cellid_string = "{}.{}.{}".format(cell_side, cell_row, cell_layer)
                    cellid_text = ROOT.TText(x1 + 0.5 * self.cell_sizex, y1 + 0.667 * self.cell_sizey, cellid_string)
                    cellid_text.SetTextSize(0.01)
                    cellid_text.SetTextAlign(22)
                    self.cellid_text_v.append(cellid_text)

                    cellnum_string = "{}".format(cellnum)
                    cellnum_text = ROOT.TText(x1 + 0.5 * self.cell_sizex, y1 + 0.333 * self.cell_sizey, cellnum_string)
                    cellnum_text.SetTextSize(0.01)
                    cellnum_text.SetTextAlign(22)
                    self.cellnum_text_v.append(cellnum_text)

                    content_text = ROOT.TText(x1 + 0.5 * self.cell_sizex, y1 + 0.333 * self.cell_sizey, "")
                    content_text.SetTextSize(0.01)
                    content_text.SetTextAlign(22)
                    self.content_text_v.append(content_text)

                if (cell_row % 5) == 0:
                    row_text = ROOT.TText(x1 + 0.5 * self.cell_sizex, 0.5, "{}".format(cell_row))
                    row_text.SetTextSize(0.03)
                    row_text.SetTextAngle(90)
                    row_text.SetTextAlign(22)
                    self.row_text_v.append(row_text)

        self.it_label = ROOT.TText(self.spacerx, 2.5*self.spacery+2*9*self.cell_sizey, "ITALY")
        self.it_label.SetTextSize(0.036)
        self.it_label.SetTextAlign(12)

        self.fr_label = ROOT.TText(self.spacerx, 0.5*self.spacery, "FRANCE")
        self.fr_label.SetTextSize(0.036)
        self.fr_label.SetTextAlign(12)

        self.palette_histo = None
        self.palette_axis = None

        if with_palette:
            palette_sizex = 2 * self.cell_sizex
            palette_sizey = self.cell_sizey * 9 * 2 + self.spacery

            palette_histo = ROOT.TH2D("{}_palette_histo".format(self.name), "", 1, 0, 1, 1, 0, 1)
            palette_histo.GetZaxis().SetNdivisions(509)
            palette_histo.GetZaxis().SetLabelSize(0.032)
            palette_histo.GetZaxis().SetLabelFont(62)
            palette_histo.GetZaxis().SetTickLength(0.009)
            palette_histo.SetMinimum(self.range_min)
            palette_histo.SetMaximum(self.range_max)
            palette_histo.SetContour(100)
            self.palette_histo = palette_histo

            palette_axis = ROOT.TPaletteAxis()
            palette_axis.SetHistogram(palette_histo)
            palette_axis.SetX1NDC(1 - self.spacerx - palette_sizex)
            palette_axis.SetY1NDC(self.spacery)
            palette_axis.SetX2NDC(1 - self.spacerx - palette_sizex / 2)
            palette_axis.SetY2NDC(self.spacery + palette_sizey)
            palette_axis.SetY2NDC(1 - palette_axis.GetY1NDC())
            self.palette_axis = palette_axis

    def setrange(self, zmin: float, zmax: float):
        self.range_min = zmin
        self.range_max = zmax

    def draw_cellid_label(self):
        self.draw_cellid = True

    def draw_cellnum_label(self):
        self.draw_cellnum = True

    def draw_content_label(self, string: str):
        self.draw_content_format = string
        self.draw_content = True

    def draw(self):
        if self.palette_axis is not None:
            canvas_width = 1231 * 2
        else:
            canvas_width = 1200 * 2
        canvas_height = 221 * 2

        if self.canvas is None:
            self.canvas = ROOT.TCanvas("{}_canvas".format(self.name), "{}".format(self.name), int(canvas_width / 2), int(canvas_height / 2))

            decoration_width = canvas_width / 2 - self.canvas.GetWw()
            decoration_height = canvas_height / 2 - self.canvas.GetWh()
            scroll_height = 16
            self.canvas.SetWindowSize(int(canvas_width * 3 / 4 + decoration_width),
                                      int(canvas_height + decoration_height + scroll_height))
            self.canvas.SetCanvasSize(canvas_width, canvas_height)

        if self.draw_content:
            for cellnum in range(self.nb_cell):
                if self.content[cellnum] is not None:
                    ttext = self.content_text_v[cellnum]
                    ttext.SetText(ttext.GetX(), ttext.GetY(), self.draw_content_format.format(self.content[cellnum]))

        self.canvas.cd()
        self.canvas.SetEditable(True)

        for cell_side in range(2):
            for cell_row in range(113):
                for cell_layer in range(9):
                    cellnum = cell_side * 113 * 9 + cell_row * 9 + cell_layer
                    self.cellbox[cellnum].Draw("l")
                    if self.draw_cellid:
                        self.cellid_text_v[cellnum].Draw()

                    if self.draw_cellnum:
                        self.cellnum_text_v[cellnum].Draw()

                    if self.draw_content and self.content[cellnum] is not None:
                        self.content_text_v[cellnum].Draw()

        if self.palette_axis is not None:
            self.palette_axis.Draw()

        for row in range(len(self.row_text_v)):
            self.row_text_v[row].Draw()

        self.it_label.Draw()
        self.fr_label.Draw()

        self.canvas.SetEditable(False)

        self.update()

    def reset(self):
        for cellnum in range(self.nb_cell):
            self.content[cellnum] = None
            self.cellbox[cellnum].SetFillColor(14)
            self.content_text_v[cellnum].Clear()

        self.canvas.Modified()
        self.canvas.Update()

        ROOT.gSystem.ProcessEvents()

    def getcontent(self, cellnum: int):
        return self.content[cellnum]

    def setcontent(self, cellnum: int, value: float):
        if cellnum < self.nb_cell:
            self.content[cellnum] = value
        else:
            print("*** wrong cell ID\n")

    def setcontent_(self, cell_side: int, cell_row: int, cell_layer: int, value: float):
        cellnum = cell_side * 9 * 113 + cell_row * 9 + cell_layer
        self.setcontent(cellnum, value)

    def fill(self, cellnum: int, value =1):
        if self.content[cellnum] is None:
            self.setcontent(cellnum, value)
        else:
            self.setcontent(cellnum, self.content[cellnum] + value)

    def update(self):
        # content_min = self.content[0]
        # content_max = self.content[0]
        content_min = 0
        content_max = 0

        for cellnum in range(1, self.nb_cell):
            if self.content[cellnum] is None:
                continue
            if self.content[cellnum] < content_min:
                content_min = self.content[cellnum]
            if self.content[cellnum] > content_max:
                content_max = self.content[cellnum]

        content_min = 0
        if self.range_min != -1:
            content_min = self.range_min
        if self.range_max != -1:
            content_max = self.range_max

        for cellnum in range(self.nb_cell):
            if self.content[cellnum] is not None:
                color_index = np.floor(99 * (self.content[cellnum] - content_min) / (content_max - content_min))
                if color_index < 0:
                    color_index = 0
                elif color_index >= 100:
                    color_index = 99
                self.cellbox[cellnum].SetFillColor(int(self.palette_index.get_index() + color_index))
            else:
                self.cellbox[cellnum].SetFillColor(14)

        if self.palette_axis is not None:
            self.palette_histo.SetMinimum(content_min)
            self.palette_histo.SetMaximum(content_max)
            self.palette_histo.SetContour(100)

        self.canvas.Modified()
        self.canvas.Update()

        ROOT.gSystem.ProcessEvents()

    def save(self, location: str):
        self.canvas.SaveAs(location + "/" + self.name + '_tr.pdf')


'''// // // // // // // // // // // // // // /
// sndisplay::demonstrator //
// // // // // // // // // // // // // // /'''


class demonstrator:

    def __init__(self, new_name: str):
        self.palette_index = palette()
        self.name = new_name
        self.demonstrator_canvas = None
        self.range_min = -1
        self.range_max = -1

        '''// TOP_VIEW //'''

        self.spacerx = 0.005
        self.spacery = 0.025
        self.title_sizey = 0.0615

        self.mw_sizey = (1-2*self.spacery-self.title_sizey)/(2.0 + 4*1.035 + 0.125)
        self.xw_sizey = 1.035 * self.mw_sizey
        self.se_sizey = 0.125 * self.mw_sizey
        self.gg_sizey = (1-2*self.spacery-self.title_sizey-2*self.mw_sizey-self.se_sizey)/18.0

        self.mw_sizex = (1 - 2 * self.spacerx) / (20 + 2 * 0.5 * 0.720)
        self.xw_sizex = (1 - 2 * self.spacerx - 20 * self.mw_sizex)
        self.se_sizex = (1 - 2 * self.spacerx - 2 * self.xw_sizex)
        self.gg_sizex = self.se_sizex / 113.0

        self.top_om_content = []
        self.top_om_box = []
        self.top_om_text = []
        self.top_gg_content = []
        self.top_gg_box = []
        self.top_gg_ellipse = []

        '''// MW(columnonly)'''

        for mw_side in range(2):
            for mw_column in range(20):
                x1 = self.spacerx + 0.5 * self.xw_sizex + mw_column * self.mw_sizex
                y1 = self.spacery + (1-mw_side) * (self.mw_sizey+4 * self.xw_sizey+self.se_sizey)
                x2 = x1 + self.mw_sizex
                y2 = y1 + self.mw_sizey

                self.top_om_content.append(0)

                box = ROOT.TBox(x1, y1, x2, y2)
                box.SetFillColor(0)
                box.SetLineWidth(1)
                self.top_om_box.append(box)

                omid_string = "M:{}.{}.*".format(mw_side, mw_column)
                omid_text = ROOT.TText(x1+0.5 * self.mw_sizex, y1+0.667 * self.mw_sizey, omid_string)
                omid_text.SetTextSize(0.032)
                omid_text.SetTextAlign(22)
                self.top_om_text.append(omid_text)

        '''// XW (column only)'''

        for xw_side in range(2):
            for xw_wall in range(2):
                for xw_column in range(2):
                    # omnum = 40 + xw_side * 2 * 2 + xw_wall * 2 + xw_column

                    x1 = self.spacerx + xw_wall * (self.xw_sizex+113 * self.gg_sizex)
                    x2 = x1 + self.xw_sizex
                    y1 = self.spacery + self.mw_sizey

                    if xw_side == 0:
                        y1 += 2 * self.xw_sizey + self.se_sizey + xw_column * self.xw_sizey
                    else:
                        y1 += (1-xw_column) * self.xw_sizey

                    y2 = y1 + self.xw_sizey

                    self.top_om_content.append(0)

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.top_om_box.append(box)

                    omid_string = "X:{}.{}.{}.*".format(xw_side, xw_wall, xw_column)
                    omid_text = ROOT.TText(x1+0.5 * self.xw_sizex, y1+0.6 * self.xw_sizey, omid_string)
                    omid_text.SetTextSize(0.032)
                    omid_text.SetTextAlign(22)
                    self.top_om_text.append(omid_text)

        for gg_side in range(2):
            for gg_row in range(113):
                for gg_layer in range(9):
                    # ggnum = gg_side * 113 * 9 + gg_row * 9 + gg_layer
                    x1 = self.spacerx + self.xw_sizex + gg_row * self.gg_sizex

                    y1 = self.spacery + self.mw_sizey
                    if gg_side == 0:
                        y1 += 9 * self.gg_sizey + self.se_sizey + gg_layer * self.gg_sizey
                    else:
                        y1 += (8-gg_layer) * self.gg_sizey

                    x2 = x1 + self.gg_sizex
                    y2 = y1 + self.gg_sizey

                    self.top_gg_content.append(0)

                    box = ROOT.TBox(x1, y1, x2, y2)
                    box.SetFillColor(0)
                    box.SetLineWidth(1)
                    self.top_gg_box.append(box)

                    ellipse = ROOT.TEllipse((x1+x2) / 2, (y1+y2) / 2, self.gg_sizex / 2, self.gg_sizey / 2)
                    ellipse.SetFillColor(0)
                    ellipse.SetLineWidth(1)
                    self.top_gg_ellipse.append(ellipse)
        self.title = ROOT.TText(self.spacerx, 1 - self.title_sizey * 3 / 4, "")
        self.title.SetTextSize(0.056)
        self.title.SetTextAlign(12)

    def setrange(self, zmin: float, zmax: float):
        self.range_min = zmin
        self.range_max = zmax

    def draw_top(self):
        if self.demonstrator_canvas is None:
            canvas_width = 1600
            canvas_height = 400

            self.demonstrator_canvas = ROOT.TCanvas("C_demonstrator_{}".format(self.name), "{}".format(self.name),
                                                    canvas_width, canvas_height)
            decoration_width = canvas_width - self.demonstrator_canvas.GetWw()
            decoration_height = canvas_height - self.demonstrator_canvas.GetWh()
            self.demonstrator_canvas.SetWindowSize(canvas_width + decoration_width, canvas_height + decoration_height)
            self.demonstrator_canvas.SetFixedAspectRatio()

        for mw_side in range(2):
            for mw_column in range(20):
                top_om_num = mw_side * 20 + mw_column
                self.top_om_box[top_om_num].Draw("l")
                self.top_om_text[top_om_num].Draw()

        for xw_side in range(2):
            for xw_wall in range(2):
                for xw_column in range(2):
                    top_om_num = 40 + xw_side * 2 * 2 + xw_wall * 2 + xw_column
                    self.top_om_box[top_om_num].Draw("l")
                    self.top_om_text[top_om_num].Draw()

        for gg_side in range(2):
            for gg_row in range(113):
                for gg_layer in range(9):
                    top_gg_num = gg_side * 113 * 9 + gg_row * 9 + gg_layer
                    self.top_gg_box[top_gg_num].Draw("l")
                    self.top_gg_ellipse[top_gg_num].Draw("l")
        self.title.Draw()
        self.update()

    def setomcontent(self, om_num: int, value: float):
        top_om_num = -1
        if om_num < 260:  # // MW IT
            om_side = 0
            om_column = (om_num / 13)
            top_om_num = om_side * 20 + om_column
        elif om_num < 520:  # // MW IT
            om_side = 1
            om_column = (om_num-260) / 13
            top_om_num = om_side * 20 + om_column
        elif om_num < 648:  # // XW
            if om_num < 584:
                om_side = 0
            else:
                om_side = 1
            # om_side = (om_num < 584) ? 0: 1
            om_wall = (om_num - 520 - om_side * 64) / 32
            om_column = (om_num - 520 - om_side * 64 - om_wall * 32) / 16
            top_om_num = 40 + om_side * 2 * 2 + om_wall * 2 + om_column

        self.top_om_content[int(top_om_num)] = value

    def setggcontent(self, cell_num: int, value: float):
        if cell_num < 2034:
            self.top_gg_content[cell_num] = value
        else:
            print("*** wrong cell ID\n")

    def setggcontent_(self, cell_side: int, cell_row: int, cell_layer: int, value: float):
        cell_num = cell_side * 9 * 113 + cell_row * 9 + cell_layer
        self.setggcontent(cell_num, value)

    def setggcolor(self, cell_num: int, colour: ROOT.TColor):
        if cell_num < 2034:
            self.top_gg_ellipse[cell_num].SetFillColor(colour)
        else:
            print("*** wrong cell ID\n")

    def setggcolor_(self, cell_side: int, cell_row: int, cell_layer: int, color: ROOT.TColor):
        cell_num = cell_side * 9 * 113 + cell_row * 9 + cell_layer
        self.setggcolor(cell_num, color)

    def reset(self):
        for om in range(len(self.top_om_content)):
            self.top_om_content[om] = None
            self.top_om_box[om].SetFillColor(14)

        for gg in range(len(self.top_gg_content)):
            self.top_gg_content[gg] = None
            self.top_gg_ellipse[gg].SetFillColor(14)

        self.demonstrator_canvas.Modified()
        self.demonstrator_canvas.Update()

        ROOT.gSystem.ProcessEvents()

    def settitle(self, text: str):
        self.title.SetText(self.title.GetX(), self.title.GetY(), text)

    def update(self):
        top_content_min = self.top_om_content[0]
        top_content_max = self.top_om_content[0]

        for om in range(len(self.top_om_content)):
            if self.top_om_content[om] < top_content_min:
                top_content_min = self.top_om_content[om]
            if self.top_om_content[om] > top_content_max:
                top_content_max = self.top_om_content[om]

        for gg in range(len(self.top_gg_content)):
            if self.top_gg_content[gg] < top_content_min:
                top_content_min = self.top_gg_content[gg]
            if self.top_gg_content[gg] > top_content_max:
                top_content_max = self.top_gg_content[gg]

        top_content_min = 0
        if self.range_min != -1:
            top_content_min = self.range_min
        if self.range_max != -1:
            top_content_max = self.range_max

        for om in range(len(self.top_om_content)):
            if self.top_om_content[om] != 0:
                color_index = np.floor(99 * (self.top_om_content[om]-top_content_min) / (top_content_max-top_content_min))
                if color_index < 0:
                    color_index = 0
                elif color_index >= 100:
                    color_index = 99
                self.top_om_box[om].SetFillColor(int(self.palette_index.get_index() + color_index))
            else:
                self.top_om_box[om].SetFillColor(0)

        for gg in range(len(self.top_gg_content)):
            if self.top_gg_content[gg] != 0:
                color_index = np.floor(99 * (self.top_gg_content[gg]-top_content_min) / (top_content_max-top_content_min))
                if color_index < 0:
                    color_index = 0
                elif color_index >= 100:
                    color_index = 99
                self.top_gg_ellipse[gg].SetFillColor(int(self.palette_index.get_index() + color_index))
            else:
                self.top_gg_ellipse[gg].SetFillColor(0)

        self.demonstrator_canvas.Modified()
        self.demonstrator_canvas.Update()

        ROOT.gSystem.ProcessEvents()

    def save(self, location: str):
        self.demonstrator_canvas.SaveAs(location + "/" + self.name + '_d.pdf')
        self.reset()


def sndisplay_calo_test():
    sncalo = calorimeter('test', with_palette=True)

    sncalo.draw_omid_label()
    sncalo.draw_content_label('{:.1f}')

    for omnum in range(520):
        sncalo.setcontent(omnum, gauss(100.0, 10.0))

    '''for omnum in range(712):
        sncalo.setcontent(omnum, gauss(50.0, 10.0))'''

    sncalo.draw()
    sncalo.save(".")


def sndisplay_tracker_test():
    sntracker = tracker('test', with_palette=True)

    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.1f}')

    for cell in range(2034):
        sntracker.setcontent(cell, gauss(100.0, 10.0))

    '''for omnum in range(712):
        sncalo.setcontent(omnum, gauss(50.0, 10.0))'''

    sntracker.draw()
    sntracker.save(".")


def sndisplay_demon_test():
    sndemonstrator = demonstrator("demonstrator_test")
    anode_and_two_cathodes = 1
    anode_and_one_cathode = 0.85
    anode_and_no_cathode = 0.7
    two_cathodes_only = 0.5
    one_cathode_only = 0.2

    sndemonstrator.setggcontent_(0,  7, 7, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(0,  7, 8, anode_and_one_cathode)
    sndemonstrator.setggcontent_(0,  8, 6, one_cathode_only)
    sndemonstrator.setggcontent_(0,  8, 7, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(0,  9, 2, anode_and_one_cathode)
    sndemonstrator.setggcontent_(0,  9, 3, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(0,  9, 4, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(0, 10, 0, anode_and_no_cathode)
    sndemonstrator.setggcontent_(0, 10, 1, anode_and_one_cathode)

    sndemonstrator.setggcontent_(1,  8, 6, one_cathode_only)
    sndemonstrator.setggcontent_(1,  8, 7, anode_and_one_cathode)
    sndemonstrator.setggcontent_(1,  9, 3, one_cathode_only)
    sndemonstrator.setggcontent_(1,  9, 4, anode_and_one_cathode)
    sndemonstrator.setggcontent_(1,  9, 5, anode_and_one_cathode)
    sndemonstrator.setggcontent_(1,  9, 6, anode_and_one_cathode)
    sndemonstrator.setggcontent_(1,  9, 7, anode_and_one_cathode)
    sndemonstrator.setggcontent_(1, 10, 0, anode_and_no_cathode)
    sndemonstrator.setggcontent_(1, 10, 1, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(1, 10, 2, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(1, 10, 3, anode_and_two_cathodes)
    sndemonstrator.setggcontent_(1, 11, 0, anode_and_one_cathode)

    sndemonstrator.settitle("RUN 609 // TRIGGER 9")
    sndemonstrator.draw_top()

    for side in range(2):
        for row in range(15, 113):
            for layer in range(9):
                sndemonstrator.setggcolor_(side, row, layer, 14)
    sndemonstrator.save(".")


if __name__ == '__main__':
    sndisplay_calo_test()
    sndisplay_tracker_test()
    sndisplay_demon_test()
