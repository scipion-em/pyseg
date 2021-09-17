# **************************************************************************
# *
# * Authors:  Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from os.path import dirname

import numpy as np
from pwem.emlib.image import ImageHandler
from pwem.objects.data import Transform, String
import pwem.convert.transformations as tfs
from os.path import join

from pyworkflow.object import List, Float
from pyworkflow.utils import removeBaseExt
from relion.convert import Table
from reliontomo.convert.convert30_tomo import TOMO_NAME, SUBTOMO_NAME, COORD_X, COORD_Y, COORD_Z, ROT, TILT, PSI, \
    RELION_TOMO_LABELS, TILT_PRIOR, PSI_PRIOR, CTF_MISSING_WEDGE, SHIFTX, SHIFTY, SHIFTZ
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SubTomogram, Coordinate3D, TomoAcquisition, Tomogram

FILE_NOT_FOUND = 'file_not_found'
PS_SEG_IMAGE = '_psSegImage'

PYSEG_PICKING_LABELS = [TOMO_NAME,
                        SUBTOMO_NAME,
                        PS_SEG_IMAGE,
                        COORD_X,
                        COORD_Y,
                        COORD_Z,
                        ROT,
                        TILT,
                        PSI,
                        ]

# Star files coding
RELION_SUBTOMO_STAR = 0
PYSEG_PICKING_STAR = 1


def readStarFile(prot, outputSetObject, fileType, starFile=None, invert=True, returnTable=False):
    warningMsg = None
    # # Star file can be provided by the user or not, depending on the protocol invoking this method
    # if not starFile:
    #     starFile = prot.starFile.get()
    #
    # # If the star file is currently in the extra folder of another protocol execution, the paths
    # # generated with method _getExtraPath() will become wrong, but it doesn't have to be located there
    # if 'extra' in starFile:
    #     starPath = ''
    # else:
    #     starPath = dirname(starFile) + '/'

    tomoTable = Table()
    tomoTable.read(starFile)

    if fileType == RELION_SUBTOMO_STAR:
        labels = RELION_TOMO_LABELS
        _relionTomoStar2Subtomograms(prot,outputSetObject, tomoTable, invert)
    else:  # fileType == PYSEG_PICKING_STAR:
        labels = PYSEG_PICKING_LABELS
        _pysegStar2Coords3D(prot, outputSetObject, tomoTable, invert)

    if not tomoTable.hasAllColumns(labels):
        missingCols = [name for name in labels if name not in tomoTable.getColumnNames()]
        warningMsg = 'Columns %s\nwere not found in the star file provided.\nThe corresponding numerical ' \
                     'values will be considered as 0.' \
                     % '  '.join(['*' + colName + '*' for colName in missingCols])

    if returnTable:
        return warningMsg, tomoTable
    else:
        return warningMsg


def manageIhDims(fileName, z, n):
    if fileName.endswith('.mrc') or fileName.endswith('.map'):
        # fileName += ':mrc'
        if z == 1 and n != 1:
            zDim = n
        else:
            zDim = z
    else:
        zDim = z

    return zDim, fileName


def _relionTomoStar2Subtomograms(prot, outputSubTomogramsSet, tomoTable, invert):
    ih = ImageHandler()
    samplingRate = outputSubTomogramsSet.getSamplingRate()
    for row, inSubtomo in zip(tomoTable, prot.inputSubtomos.get()):
        subtomo = SubTomogram()
        coordinate3d = Coordinate3D()
        transform = Transform()
        origin = Transform()

        volname = row.get(TOMO_NAME, FILE_NOT_FOUND)
        subtomoFn = row.get(SUBTOMO_NAME, FILE_NOT_FOUND)

        subtomo.setVolName(managePath4Sqlite(volname))
        subtomo.setTransform(transform)
        subtomo.setAcquisition(TomoAcquisition())
        subtomo.setClassId(row.get('rlnClassNumber', 0))
        subtomo.setSamplingRate(samplingRate)
        subtomo.setCoordinate3D(coordinate3d)  # Needed to get later the tomogram pointer via getCoordinate3D()

        x = row.get(COORD_X, 0)
        y = row.get(COORD_Y, 0)
        z = row.get(COORD_Z, 0)
        tiltPrior = row.get(TILT_PRIOR, 0)
        psiPrior = row.get(PSI_PRIOR, 0)
        # ctf3d = row.get(CTF_MISSING_WEDGE, FILE_NOT_FOUND)
        coordinate3d = subtomo.getCoordinate3D()
        coordinate3d.setX(float(x), BOTTOM_LEFT_CORNER)
        coordinate3d.setY(float(y), BOTTOM_LEFT_CORNER)
        coordinate3d.setZ(float(z), BOTTOM_LEFT_CORNER)
        coordinate3d._3dcftMrcFile = inSubtomo.getCoordinate3D()._3dcftMrcFile  # Used for the ctf3d in Relion 3.0 (tomo)
        M = _getTransformMatrix(row, invert)
        transform.setMatrix(M)
        subtomo.setCoordinate3D(coordinate3d)
        subtomo._tiltPriorAngle = Float(tiltPrior)
        subtomo._psiPriorAngle = Float(psiPrior)

        # Set the origin and the dimensions of the current subtomogram
        x, y, z, n = ih.getDimensions(subtomoFn)
        zDim, _ = manageIhDims(subtomoFn, z, n)
        origin.setShifts(x / -2. * samplingRate, y / -2. * samplingRate, zDim / -2. * samplingRate)
        subtomo.setOrigin(origin)

        subtomo.setFileName(managePath4Sqlite(subtomoFn))
        # if subtomo is in a vesicle
        if 'tid_' in subtomoFn:
            vesicleId = subtomoFn.split('tid_')[1]
            vesicleId = vesicleId[0]
            scoor = subtomo.getCoordinate3D()
            scoor.setGroupId(vesicleId)
            subtomo.setCoordinate3D(scoor)

        # Add current subtomogram to the output set
        outputSubTomogramsSet.append(subtomo)


def _getTransformMatrix(row, invert):
    shiftx = row.get(SHIFTX, 0)
    shifty = row.get(SHIFTY, 0)
    shiftz = row.get(SHIFTZ, 0)
    tilt = row.get(TILT, 0)
    psi = row.get(PSI, 0)
    rot = row.get(ROT, 0)
    shifts = (float(shiftx), float(shifty), float(shiftz))
    angles = (float(rot), float(tilt), float(psi))
    radAngles = -np.deg2rad(angles)
    M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if invert:
        M[0, 3] = -shifts[0]
        M[1, 3] = -shifts[1]
        M[2, 3] = -shifts[2]
        M = np.linalg.inv(M)
    else:
        M[0, 3] = shifts[0]
        M[1, 3] = shifts[1]
        M[2, 3] = shifts[2]

    return M


def managePath4Sqlite(fpath):
    return fpath if fpath != FILE_NOT_FOUND else fpath


def getTomoSetFromStar(prot, starFile):
    samplingRate = prot.pixelSize.get()
    imgh = ImageHandler()
    tomoTable = Table()
    tomoTable.read(starFile)
    tomoList = [row.get(TOMO_NAME, FILE_NOT_FOUND) for row in tomoTable]
    prot.tomoList = List(tomoList)
    tomoNamesUnique = list(set(tomoList))

    # Create a Volume template object
    tomo = Tomogram()
    tomo.setSamplingRate(samplingRate)
    for fileName in tomoNamesUnique:
        x, y, z, n = imgh.getDimensions(fileName)
        if fileName.endswith('.mrc') or fileName.endswith('.map'):
            if z == 1 and n != 1:
                zDim = n
                n = 1
            else:
                zDim = z
        else:
            zDim = z
        origin = Transform()
        origin.setShifts(x / -2. * samplingRate,
                         y / -2. * samplingRate,
                         zDim / -2. * samplingRate)

        tomo.setOrigin(origin)  # read origin from form

        for index in range(1, n + 1):
            tomo.cleanObjId()
            tomo.setLocation(index, fileName)
            tomo.setAcquisition(TomoAcquisition(**prot.acquisitionParams))
            prot.tomoSet.append(tomo)


def _pysegStar2Coords3D(prot, output3DCoordSet, tomoTable, invert):
    for tomoNum, tomo in enumerate(prot.tomoSet.iterItems()):
        tomoName = tomo.getFileName().replace(':mrc', '')
        for row in tomoTable:
            # Create the set of coordinates referring each of them to its corresponding tomogram (ancestor)
            if row.get(TOMO_NAME) == tomoName:
                coordinate3d = Coordinate3D()
                coordinate3d.setVolId(tomoNum)
                coordinate3d.setVolume(tomo)
                x = row.get(COORD_X, 0)
                y = row.get(COORD_Y, 0)
                z = row.get(COORD_Z, 0)
                M = _getTransformMatrix(row, invert)
                coordinate3d.setX(float(x), BOTTOM_LEFT_CORNER)
                coordinate3d.setY(float(y), BOTTOM_LEFT_CORNER)
                coordinate3d.setZ(float(z), BOTTOM_LEFT_CORNER)
                coordinate3d.setMatrix(M)
                coordinate3d.setGroupId(_getVesicleIdFromSubtomoName(row.get(SUBTOMO_NAME, FILE_NOT_FOUND)))

                # Add current subtomogram to the output set
                output3DCoordSet.append(coordinate3d)


def _getVesicleIdFromSubtomoName(subtomoName):
    """PySeg adds the vesicle index to the name of the subtomogram, with a suffix of type
    _tid_[VesicleNumber].mrc. Example: Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_0.mrc.
    This function returns that vesicle number for a given """
    pattern = '_tid_'
    baseName = removeBaseExt(subtomoName)
    pos = baseName.find(pattern)
    # regexPattern = re.compile('^_tid_[0-9]{1, 3}')
    # res = regexPattern.match(baseName)
    return baseName[pos + len(pattern)::]
