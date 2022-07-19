#!/usr/bin/env python3
"""Make mapping svg file.

usage:
  map1.py [-h|--help] -m efl_per_cellsize -o out.svg

options:
  --help                show this help message and exit
  -m efl_per_cellsize   parameter efl_per_cellsize (efl / Cell scale)
  -o out.svg            output file
"""
import numpy as np
import math
import svgwrite
from docopt import docopt

efl_per_cellsize = 6320
outfile = 'map.svg'

if __name__ == '__main__':
    args = docopt(__doc__)
    efl_per_cellsize = int(args['-m'])
    outfile = args['-o']

# Drawing parameters
mm = 3.543307
deg = 14*mm  # 1 deg = 40mm
paper_width = 297*mm
paper_hight = 210*mm
xcenter = paper_width/2
ycenter = paper_hight/2

# Model parameters
cells_per_pix = 13
pixel_size = 10
efl = efl_per_cellsize*pixel_size/cells_per_pix
pfov = 10e-3/efl * 180/math.pi
num_pixels = 1952
chip_spacing_a = 22.4/efl * 180/math.pi
chip_spacing_b = 22.4/efl * 180/math.pi

text = 'EFL={:.3f}'.format(efl)

# Field parameters
widex = 2
widey = 2
gc = 0.7
gdw = 1.5
gdh = 0.3


def degxys(x, y):
    return (x*deg*mm, y*deg*mm)


def degxy(x, y):
    return ((x*deg*mm+xcenter, y*deg*mm+ycenter))


# Prepare a container for all elements
dwg = svgwrite.Drawing(outfile, size=(297*mm, 210*mm))

# Add an axis group
axis = dwg.add(dwg.g(id='axis'))
axis.add(dwg.line(start=degxy(-widex, -widey),
         end=degxy(widex, -widey), stroke='black', stroke_width=2))
axis.add(dwg.line(start=degxy(-widex,     0),
         end=degxy(widex,     0), stroke='black', stroke_width=1))
axis.add(dwg.line(start=degxy(-widex, widey),
         end=degxy(widex, widey), stroke='black', stroke_width=2))
axis.add(dwg.line(start=degxy(-widex, -widey),
         end=degxy(-widex, widey), stroke='black', stroke_width=2))
axis.add(dwg.line(start=degxy(0, -widey),
         end=degxy(0, widey), stroke='black', stroke_width=1))
axis.add(dwg.line(start=degxy(widex, -widey),
         end=degxy(widex, widey), stroke='black', stroke_width=2))

# Add Text
title = dwg.add(dwg.g(id='axis'))
title.add(dwg.text(text, insert=degxy(-gc, 0.8)))

# Add the central region
cent = dwg.add(dwg.g(id='cent'))
cent.add(dwg.rect(insert=degxy(-gc, -gc), size=degxys(2*gc, 2*gc)))
cent.fill('green', opacity=0.1)

# Add the disk region
disk = dwg.add(dwg.g(id='disk'))
disk.add(dwg.rect(insert=degxy(-gc, -gdh), size=degxys(gc+gdw, 2*gdh)))
disk.fill('blue', opacity=0.1)

# Add a field


def addfield(center_x, center_y, position_angle, dwg, id):
    field = dwg.add(dwg.g(id='field{}'.format(id)))
    for i in range(4):
        x = np.empty(4)
        y = np.empty(4)
        if i == 0 or i == 1:
            x[0] = -chip_spacing_a/2-num_pixels*pfov/2
        else:
            x[0] = chip_spacing_a/2-num_pixels*pfov/2
        if i == 0 or i == 2:
            y[0] = -chip_spacing_b/2-num_pixels*pfov/2
        else:
            y[0] = chip_spacing_b/2-num_pixels*pfov/2
        x[1] = x[0]+num_pixels*pfov
        x[2] = x[1]
        x[3] = x[0]
        y[1] = y[0]
        y[2] = y[1]+num_pixels*pfov
        y[3] = y[2]
        xd = x+center_x
        yd = y+center_y
        points = (degxy(xd[0], yd[0]), degxy(xd[1], yd[1]),
                  degxy(xd[2], yd[2]), degxy(xd[3], yd[3]))
        field.add(dwg.polygon(points=points))
        field.fill('red', opacity=0.2)


addfield(-gc+chip_spacing_a/2+num_pixels*pfov/2, -gc +
         chip_spacing_b/2+num_pixels*pfov/2, 0, dwg, 0)
addfield(-gc+chip_spacing_a/2+num_pixels*pfov/2, -gc +
         chip_spacing_b/2*3+num_pixels*pfov/2, 0, dwg, 1)
addfield(-gc+chip_spacing_a/2+num_pixels*pfov/2, -gc +
         chip_spacing_b/2*7+num_pixels*pfov/2, 0, dwg, 1)
addfield(-gc+chip_spacing_a*3/2+num_pixels*pfov/2, -gc +
         chip_spacing_b/2*7+num_pixels*pfov/2, 0, dwg, 1)

dwg.save()
