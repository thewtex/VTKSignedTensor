#!/usr/bin/env python

import numpy as np



def write_uniaxial_strain(filename, axial_strain):
    poissons_ratio = 0.3
    lateral_strain = -poissons_ratio * axial_strain
    side_points = 4
    with open(filename, 'w') as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Strain tensors for uniaxial loading of a cube.\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_POINTS\n')
        f.write('DIMENSIONS ' + str(side_points) + ' ' + str(side_points) + ' ' + str(side_points) + '\n')
        f.write('ORIGIN 0.0 0.0 0.0\n')
        f.write('SPACING 1.0 1.0 1.0\n')
        f.write('\nPOINT_DATA ' + str(side_points**3) + '\n')
        f.write('TENSORS strain double\n')

        # we will take y as the axial direction
        for idx in range(side_points**3):
            f.write('{0:.20g} 0.0 0.0\n'.format(lateral_strain))
            f.write('0.0 {0:.20g} 0.0\n'.format(axial_strain))
            f.write('0.0 0.0 {0:.20g}\n'.format(lateral_strain))
            f.write('\n')

axial_strain = -0.05
write_uniaxial_strain('uniaxial_compression.vtk', axial_strain)
axial_strain = 0.05
write_uniaxial_strain('uniaxial_tension.vtk', axial_strain)
