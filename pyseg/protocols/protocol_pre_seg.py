# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import glob
import logging
import os
from enum import Enum
from os.path import abspath, exists

import mrcfile

from pwem.convert.headers import fixVolume
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyseg.convert.convert import getVesicleIdFromSubtomoName

from pyworkflow.object import String, Integer
from pyworkflow.protocol import IntParam, FloatParam, GT, LEVEL_ADVANCED, PointerParam
from pyworkflow.utils import Message, removeBaseExt, removeExt
from scipion.constants import PYTHON
from tomo.objects import SetOfTomoMasks, TomoMask, SetOfSubTomograms, SubTomogram, SetOfTomograms
from pyseg import Plugin
from pyseg.constants import PRESEG_SCRIPT, TOMOGRAM, PYSEG_LABEL, VESICLE, NOT_FOUND, \
    PYSEG_OFFSET_X, PYSEG_OFFSET_Y, PYSEG_OFFSET_Z, SEGMENTATION, RLN_ORIGIN_X, RLN_ORIGIN_Y, \
    RLN_ORIGIN_Z, PYSEG_PSI, PYSEG_TILT, PYSEG_ROT
from relion.convert import Table
import numpy as np
from tomo.utils import getObjFromRelation


logger = logging.getLogger(__name__)
PRE_STAR_SUFFIX = '_pre.star'
N_SEGS_REMOVED = '_nSegmentationsRemoved'


class outputObjects(Enum):
    vesicles = SetOfSubTomograms
    segmentations = SetOfTomoMasks


class ProtPySegPreSegParticles(EMProtocol):
    """Segment membranes into membranes, inner surroundings and outer surroundings"""

    _label = 'preseg membranes'
    _starFile = None

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inTomoMasks', PointerParam,
                      pointerClass='SetOfTomoMasks',
                      label='Tomomasks (segmentations)',
                      important=True,
                      allowsNull=False,
                      help='Pointer to segmented and annotated tomograms within Scipion.')
        group = form.addGroup('Membrane segmentation')
        group.addParam('spOffVoxels', IntParam,
                       label='Offset voxels',
                       default=1,
                       validators=[GT(0)],
                       help='Margin to ensure that the desired entities, e. g. membranes, proteins, are included.')
        group.addParam('sgThreshold', IntParam,
                       default=-1,
                       label='Density threshold',
                       expertLevel=LEVEL_ADVANCED,
                       help='All the voxels with density equal to or higher than the threshold are set to 1. '
                            'The remaining voxels are set to 0.')
        group.addParam('sgSizeThreshold', IntParam,
                       default=-1,
                       label='Size threshold (voxels)',
                       condition='sgThreshold > 0',
                       help='It sets the minimal size for a component to be considered as membrane.')
        group.addParam('sgMembThk', FloatParam,
                       label='Segmented membrane thickness (Å)',
                       default=40,
                       allowsNull=False,
                       validators=[GT(0)],
                       help='Value introduced will be divided by 2 internally, because it is expected like that '
                            'by PySeg.')
        group.addParam('sgMembNeigh', FloatParam,
                       label='Segmented membrane neighbours (Å)',
                       allowsNull=False,
                       validators=[GT(0)],
                       help='Thickness around the membrane to represent the in-membrane and out-membrane surroundings '
                            'desired to be included in the analysis.')
        line = group.addLine('Artifact removal', help='When the segmentations do not come from the a manual annotation '
                                                      'tool, like the Membrane Annotator from scipion-em-tomosegmemtv, '
                                                      'there may be spurious structures that should be removed.')
        line.addParam('minXY', IntParam,
                      label='Min size (px)',
                      default=-1,
                      help='Minimum length, in pixels, allowed in height and width. Value -1 ignore the size '
                           'filtering.')
        line.addParam('minDepth', IntParam,
                      label='Min depth (px)',
                      default=-1,
                      help='Minimum depth, in pixels, allowed. Value -1 ignore the depth filtering.')

    def _insertAllSteps(self):
        self._starFile = self._getExtraPath('inStar.star')
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.pysegPreSegStep)
        self._insertFunctionStep(self.getMembraneCenterStep)
        self._insertFunctionStep(self.pysegPreSegCenteredStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertInputStep(self):
        from pwem import Domain
        xmipp3 = Domain.importFromPlugin('xmipp3')
        # Generate the star file with the vesicles centered for the second pre_seg execution
        outputTable = self._createTable(isConvertingInput=True)
        for tomoMask in self.inTomoMasks.get():
            # Convert to MRC the tomograms to which the tomomasks are referred to if they are not, because it's
            # the format searched by pyseg
            tomoFile = tomoMask.getVolName()
            if not tomoFile.endswith('.mrc'):
                mrcTomoFile = self._getExtraPath(removeBaseExt(tomoFile) + '.mrc')
                args = '-i %s -o %s -t vol ' % (tomoFile, mrcTomoFile)
                xmipp3.Plugin.runXmippProgram('xmipp_image_convert', args)
                # Set the volName to the generated volume
                tomoMask.setVolName(mrcTomoFile)

            vesicle = tomoMask.getFileName()
            materials = self._getMaterialsList(vesicle)
            if materials:
                for materialIndex in materials:  # Get annotated materials from txt file and add one
                    # Add row to output table
                    outputTable.addRow(tomoMask.getVolName(),
                                       vesicle,
                                       int(materialIndex),
                                       vesicle)
        outputTable.write(self._starFile)

    def pysegPreSegStep(self):
        outDir = self._getExtraPath()
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPreSegCmd(self._starFile, outDir))

    def getMembraneCenterStep(self):
        inStar = self.getPresegOutputFile(self._starFile)
        self._findVesicleCenter(self._starFile, inStar)

    def pysegPreSegCenteredStep(self):
        inStar = abspath(self.getVesiclesCenteredStarFile())
        outDir = self._getExtraPath()
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPreSegCmd(inStar, outDir))

    def createOutputStep(self):
        segVesSet, vesSet = self._genOutputSetOfTomoMasks()
        self._defineOutputs(**{outputObjects.vesicles.name: vesSet,
                               outputObjects.segmentations.name: segVesSet})

        self._defineSourceRelation(self.inTomoMasks.get(), vesSet)
        self._defineSourceRelation(self.inTomoMasks.get(), segVesSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        valMsg = []
        if not self._getTomoFromRelations():
            valMsg.append('Unable to get trough the relations the tomograms corresponding to the introduced tomomasks')
        return valMsg

    def _summary(self):
        msgs = []
        nSegmentationsRemoved = getattr(self, N_SEGS_REMOVED, None)
        if nSegmentationsRemoved:
            msgs.append(f'Number of segmentations removed by size filtering was {nSegmentationsRemoved}')
    # --------------------------- UTIL functions -----------------------------------

    def _getPreSegCmd(self, inStar, outDir):
        preSegCmd = ' '
        preSegCmd += '%s ' % Plugin.getHome(PRESEG_SCRIPT)
        preSegCmd += '--inStar %s ' % inStar
        preSegCmd += '--outDir %s ' % outDir
        preSegCmd += '--spOffVoxels %s ' % self.spOffVoxels.get()
        preSegCmd += '--sgVoxelSize %s ' % (float(self._getSamplingRate())/10)  # required in nm
        preSegCmd += '--sgThreshold %s ' % self.sgThreshold.get()
        preSegCmd += '--sgSizeThreshold %s ' % self.sgSizeThreshold.get()
        preSegCmd += '--sgMembThk %s ' % self._checkValue4PySeg(self.sgMembThk.get()/2)  # half of the thickness in nm
        preSegCmd += '--sgMembNeigh %s ' % self._checkValue4PySeg(self.sgMembNeigh.get())  # required in nm

        return preSegCmd

    @staticmethod
    def _getMaterialsList(vesicle):
        # Get annotated materials from txt file and add one line for each one
        materialsFile = removeExt(vesicle) + '.txt'
        if exists(materialsFile):  # Segmentations come from the Membrane Annotator app
            logger.info(f'==> Found material list file for current vesicle: {materialsFile}')
            with open(materialsFile) as matFile:
                materialsList = matFile.read()
                # Expected format is a string like 'ind1,ind2,...,indn\n, so it's necessary to transform it into
                # a list of material indices, which may have been annotated more than once (incomplete membranes
                # that were annotated part by part)
                return set(materialsList.replace('\n', '').split(','))
        else:
            with mrcfile.mmap(vesicle) as mrc:
                materialIndices = sorted(np.unique(mrc.data))
                logger.info(f'==> Material list file not found. Simulating it. Indices found are {materialIndices}')
                return materialIndices

    def _findVesicleCenter(self, starFileInit, starFilePreseg1):
        ih = ImageHandler()
        outputTable = self._createTable()
        # Read preseg (vesicles not centered) table
        presegTable = Table()
        presegTable.read(starFilePreseg1)
        # Read initial data
        initTable = Table()
        initTable.read(starFileInit)

        # Generate the star file with the vesicles centered for the second pre_seg execution
        for row, rowp in zip(initTable, presegTable):
            tomo = row.get(TOMOGRAM, NOT_FOUND)
            vesicle = row.get(VESICLE, NOT_FOUND)
            vesicle_seg = rowp.get(SEGMENTATION, NOT_FOUND)
            materialIndex = row.get(PYSEG_LABEL, NOT_FOUND)
            # Get the box dimensions, be sure that the Dimensions format is Dimensions: 239 x 1 x 298 x 298
            # ((N)Objects x (Z)Slices x (Y)Rows x (X)Columns)
            tomoX, tomoY, tomoZ, _ = ih.getDimensions(tomo)
            x, y, z, _ = ih.getDimensions(vesicle_seg)

            # Get the upper left corner of the vesicle in the original tomogram: _psSegOffX #7 _psSegOffY #8 _psSegOffZ
            # #9 from output star file of pre_tomos_seg.py that do not distinguish inner and outer densities
            xdimCorner = rowp.get(PYSEG_OFFSET_X, 0)
            ydimCorner = rowp.get(PYSEG_OFFSET_Y, 0)
            zdimCorner = rowp.get(PYSEG_OFFSET_Z, 0)

            # Add row to output table
            outputTable.addRow(tomo,
                               vesicle,
                               materialIndex,
                               vesicle,
                               xdimCorner + x / 2,
                               ydimCorner + y / 2,
                               zdimCorner + z / 2)

        outputTable.write(self.getVesiclesCenteredStarFile())

    def getPresegOutputFile(self, inStar):
        return self._getExtraPath(removeBaseExt(inStar) + PRE_STAR_SUFFIX)

    def getVesiclesCenteredStarFile(self):
        return self._getExtraPath('presegVesiclesCentered.star')

    @ staticmethod
    def _createTable(isConvertingInput=False):
        if isConvertingInput:
            # Headers for the input star file to the first pre seg when the input data is from Scipion
            return Table(columns=[TOMOGRAM,
                                  VESICLE,
                                  PYSEG_LABEL,
                                  SEGMENTATION])
        else:
            # Headers for pySeg pre_seg centered star file
            return Table(columns=[TOMOGRAM,
                                  VESICLE,
                                  PYSEG_LABEL,
                                  SEGMENTATION,
                                  RLN_ORIGIN_X,
                                  RLN_ORIGIN_Y,
                                  RLN_ORIGIN_Z
                                  ])

    @staticmethod
    def _createOutputTable():
        # Headers for pySeg second preseg filtered by size
        return Table(columns=[TOMOGRAM,
                              VESICLE,
                              SEGMENTATION,
                              PYSEG_ROT,
                              PYSEG_TILT,
                              PYSEG_PSI,
                              PYSEG_OFFSET_X,
                              PYSEG_OFFSET_Y,
                              PYSEG_OFFSET_Z
                              ])

    @staticmethod
    def _checkValue4PySeg(value):
        return value if value == -1 else value/10

    def _genOutputSetOfTomoMasks(self):
        ih = ImageHandler()
        minSize = self.minXY.get()
        minDepth = self.minDepth.get()
        sRate = self._getSamplingRate()

        finalPreSegTable = Table()
        finalPreSegTable.read(self.getPresegOutputFile(self.getVesiclesCenteredStarFile()))
        finalStarFile = self._getExtraPath('presegResults.star')
        finalResultTable = self._createOutputTable()

        tomoMaskSet = SetOfTomoMasks.create(self._getPath(), template='tomomasks%s.sqlite', suffix='segVesicles')
        subTomoSet = SetOfSubTomograms.create(self._getPath(), template='subtomograms%s.sqlite', suffix='vesicles')
        inTomoMaskSet = self.inTomoMasks.get()
        tomoMaskSet.copyInfo(inTomoMaskSet)
        tomoMaskSet._pysegPresegVesicleCStar = String(finalStarFile)
        subTomoSet.copyInfo(inTomoMaskSet)

        nSegmentationsRemoved = 0
        for row in finalPreSegTable:
            tomogram = row.get(TOMOGRAM)
            vesicle = row.get(VESICLE)
            vesicleSeg = row.get(SEGMENTATION)
            tomoX, tomoY, tomoZ, _ = ih.getDimensions(tomogram)
            x, y, z, _ = ih.getDimensions(vesicle)
            if (x, y, z) == (tomoX, tomoY, tomoZ):  # Sometimes one of the segmented structures is the whole tomogram
                # Remove the files generated in this intermediate step
                self._removeVesicleFiles(vesicleSeg)
                logger.info(f'==> Removed vesicle files {vesicleSeg} of size {x} x {y} x {z} px')
                nSegmentationsRemoved += 1
                continue
            if minSize > 0 or minDepth > 0:  # Apply the size filtering
                # Only the regions with sizes larger than the minimun specified will be considered from this point
                if x < minSize or y < minSize or z < minDepth:
                    # Remove the files generated in this intermediate step
                    self._removeVesicleFiles(vesicleSeg)
                    logger.info(f'==> Removed vesicle files {vesicleSeg} of size {x} x {y} x {z} px')
                    nSegmentationsRemoved += 1
                    continue

            classId = int(getVesicleIdFromSubtomoName(vesicle))
            # Populate the set of tomomasks
            tomoMask = TomoMask()
            tomoMask.setSamplingRate(sRate)
            fixVolume([vesicle, vesicleSeg])
            tomoMask.setFileName(vesicleSeg)
            tomoMask.setVolName(vesicle)
            tomoMask.setClassId(classId)
            tomoMask.setSamplingRate(sRate)
            tomoMaskSet.append(tomoMask)

            # Populate the set of subtomograms
            subtomo = SubTomogram()
            subtomo.setFileName(vesicle)
            subtomo.setSamplingRate(sRate)
            subtomo.setClassId(classId)
            subtomo.setVolName(tomogram)
            subTomoSet.append(subtomo)

            # Add the corresponding row to the final star file
            finalResultTable._rows.append(row)

        finalResultTable.write(finalStarFile)
        if nSegmentationsRemoved > 0:
            setattr(self, N_SEGS_REMOVED, Integer(nSegmentationsRemoved))
        return tomoMaskSet, subTomoSet


    def _getSamplingRate(self):
        return self.inTomoMasks.get().getFirstItem().getSamplingRate()

    def _getTomoFromRelations(self):
        return getObjFromRelation(self.inTomoMasks.get(), self, SetOfTomograms)

    @staticmethod
    def _getMatchingTomo2Vesicle(tomoBNameDict, vesicleBaseName):
        for tomoBName in list(tomoBNameDict.keys()):
            if tomoBName in vesicleBaseName:
                return tomoBNameDict[tomoBName]
        return None

    @staticmethod
    def _removeVesicleFiles(vesicle_seg):
        os.remove(vesicle_seg)
        # Also remove the non-segmented version of the cropped area generated by the preseg, which has the
        # same name and path of the vesicle_seg, but without the suffix '_seg'
        os.remove(vesicle_seg.replace('_seg.mrc', '.mrc'))
