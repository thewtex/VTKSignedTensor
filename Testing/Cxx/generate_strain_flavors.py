#!/usr/bin/env python

import numpy as np

# Eigenvalues a, b, c
e_a = 1.0
e_b = 0.6
e_c = 0.4
eigenvalues = [[e_a, e_b, e_c],
        [-e_a, e_b, e_c],
        [-e_a, -e_b, e_c],
        [-e_a, e_b, -e_c],
        [e_a, -e_b, e_c],
        [e_a, -e_b, -e_c],
        [e_a, e_b, -e_c],
        [-e_a, -e_b, -e_c]]

with open('strain_flavors_eigenvalues.txt', 'w') as f:
    for eigenval in eigenvalues:
        for ii in range(2):
            f.write('{0:g}\t'.format(eigenval[ii]))
        f.write('{0:g}\n'.format(eigenval[2]))

angles = [0.0]
while(angles[-1] < 300.0):
    angles.append(angles[-1] + 60.0)

with open('strain_flavors_angles.txt', 'w') as f:
    for a in angles:
        f.write('{0}\n'.format(a))

with open('strain_flavors_strain.vtk', 'w') as f:
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Different strain flavors.\n')
    f.write('ASCII\n')
    f.write('DATASET STRUCTURED_POINTS\n')
    f.write('DIMENSIONS ' + str(len(angles)) + ' ' + str(len(eigenvalues)) + ' 1\n')
    f.write('ORIGIN 0.0 0.0 0.0\n')
    f.write('SPACING 1.0 1.0 1.0\n')
    f.write('\nPOINT_DATA ' + str(len(angles) * len(eigenvalues)) + '\n')
    f.write('TENSORS strain double\n')

    for eigenval in eigenvalues:
        for theta in angles:
            angle = np.pi/180*theta
            eigenvector_a = np.array([np.cos(angle), np.sin(angle), 0.0])
            eigenvector_a.shape = (3,1)
            eigenvector_b = np.array([np.cos(angle+np.pi/2), np.sin(angle+np.pi/2), 0.0])
            eigenvector_b.shape = (3,1)
            eigenvector_c = np.array([0.0, 0.0, 1.0])
            eigenvector_c.shape = (3,1)

            lam = np.concatenate((eigenvector_a, eigenvector_b, eigenvector_c), axis=1)
            tensor = np.dot(np.dot(lam,np.identity(3)*eigenval),np.linalg.inv(lam))
            for row in tensor:
                f.write('{0:.20g} {1:.20g} {2:.20g}\n'.format(row[0], row[1], row[2]))
            f.write('\n')
