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
from os.path import join
from emtable import Table
from pwem.emlib.image import ImageHandler
from pwem.objects.data import Transform
from pyseg.constants import NOT_FOUND, GRAPHS_OUT, VESICLE
from pyseg.utils import manageDims
from pyworkflow.object import Float
from pyworkflow.utils import removeBaseExt, createLink
from reliontomo.constants import TILT_PRIOR, PSI_PRIOR, SUBTOMO_NAME, TOMO_NAME_30
from reliontomo.convert import RELION_30_TOMO_LABELS
from reliontomo.convert.convert30_tomo import Reader
from reliontomo.convert.convertBase import getTransformMatrixFromRow
from tomo.objects import SubTomogram, TomoAcquisition


class PysegStarReader(Reader):

    def __init__(self, starFile, dataTable, **kwargs):
        super().__init__(starFile, dataTable)

    def starFile2Coords3D(self, coordsSet, precedentsSet, scaleFactor=1):
        precedentDict = {removeBaseExt(tomo.getFileName()): tomo.clone() for tomo in precedentsSet}
        for row in self.dataTable:
            coord3d = self.gen3dCoordFromStarRow(row, precedentDict, scaleFactor)
            # GroupId stuff
            vsicleName = row.get(VESICLE, None)
            if 'tid_' in vsicleName:
                vesicleId = vsicleName.split('tid_')[1]
                vesicleId = vesicleId[0]
                coord3d.setGroupId(vesicleId)

            coordsSet.append(coord3d)

    def starFile2Subtomograms(self, inSubtomos, outputSubtomos):
        warningMsg = None
        labels = RELION_30_TOMO_LABELS
        self.genOutputSubtomograms(inSubtomos, outputSubtomos)
        if not self.dataTable.hasAllColumns(labels):
            missingCols = [name for name in labels if name not in self.dataTable.getColumnNames()]
            warningMsg = 'Columns %s\nwere not found in the star file provided.\nThe corresponding numerical ' \
                         'values will be considered as 0.' \
                         % '  '.join(['*' + colName + '*' for colName in missingCols])
        return warningMsg, self.dataTable

    def genOutputSubtomograms(self, inSubtomos, outputSubtomos):
        ih = ImageHandler()
        samplingRate = outputSubtomos.getSamplingRate()
        for row, inSubtomo in zip(self.dataTable, inSubtomos):
            subtomo = SubTomogram()
            transform = Transform()
            origin = Transform()

            volname = row.get(TOMO_NAME_30, NOT_FOUND)
            subtomoFn = row.get(SUBTOMO_NAME, NOT_FOUND)
            transform.setMatrix(getTransformMatrixFromRow(row))

            subtomo.setVolName(managePath4Sqlite(volname))
            subtomo.setTransform(transform)
            subtomo.setAcquisition(TomoAcquisition())
            subtomo.setClassId(row.get('rlnClassNumber', 0)+1)
            subtomo.setSamplingRate(samplingRate)

            tiltPrior = row.get(TILT_PRIOR, 0)
            psiPrior = row.get(PSI_PRIOR, 0)
            subtomo.setCoordinate3D(inSubtomo.getCoordinate3D())
            subtomo._tiltPriorAngle = Float(tiltPrior)
            subtomo._psiPriorAngle = Float(psiPrior)

            # Set the origin and the dimensions of the current subtomogram
            x, y, z, n = ih.getDimensions(subtomoFn)
            zDim = manageDims(subtomoFn, z, n)
            origin.setShifts(x / -2. * samplingRate,
                             y / -2. * samplingRate,
                             zDim / -2. * samplingRate)
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
            outputSubtomos.append(subtomo)


def splitPysegStarFile(inStar, outDir, j=1, prefix=GRAPHS_OUT + '_', fileCounter=1):
    """Split a star file which one line for each membrane into n files of one membrane, in order to make the
    filament protocol runs faster. If the input star only has one row, it will be linked. Attribute fileCOunter is used
    manage the generated files enumeration, being possible to consider that previous files have been generated by
    using a number higher than 1."""
    outStarFiles = []
    tomoTable = Table()

    def _addRow():
        values = [vesicleRow.get(label, NOT_FOUND) for label in labels]
        outTable.addRow(*values)

    def _writeStarFile():
        outStarFile = join(outDir, '%s%03d.star' % (prefix, fileCounter))
        outStarFiles.append(outStarFile)
        outTable.write(outStarFile)
        outTable.clearRows()

    tomoTable.read(inStar)
    if len(tomoTable) == 1:
        outStarFile = join(outDir, '%s%03d.star' % (prefix, fileCounter))
        createLink(inStar, outStarFile)
        outStarFiles.append(outStarFile)
    else:
        nVesicles = tomoTable.size()
        labels = tomoTable.getColumnNames()
        outTable = Table(columns=labels)
        counter = 1
        for vesicleRow in tomoTable:
            _addRow()
            if j == 1 or (counter > 1 and counter % j == 0):
                _writeStarFile()
                fileCounter += 1
            counter += 1

        rem = nVesicles % j
        if rem > 0:
            for vesicleRow in tomoTable[-rem:-1]:
                _addRow()
            _writeStarFile()

    return outStarFiles


def managePath4Sqlite(fpath):
    return fpath if fpath != NOT_FOUND else fpath


def getVesicleIdFromSubtomoName(subtomoName):
    """PySeg adds the vesicle index to the name of the subtomogram, with a suffix of type
    _tid_[VesicleNumber].mrc. Example: Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_tid_0.mrc. In case of splitting
    into slices, the name is slightly different: Pertuzumab_1_defocus_25um_tomo_7_aliSIRT_EED_id_2_split_2.mrc.
    This function returns that vesicle number for a given subtomogram name."""
    splitPattern = '_split_'
    tidPatten = '_tid_'
    idPattern = '_id_'
    baseName = removeBaseExt(subtomoName)
    if splitPattern in baseName:
        posIni = baseName.find(idPattern) + len(idPattern)
        posEnd = baseName.find(splitPattern)
        return baseName[posIni:posEnd]
    else:
        posIni = baseName.find(tidPatten) + len(tidPatten)
        return baseName[posIni::]