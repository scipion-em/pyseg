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
import logging
from collections import OrderedDict
from enum import Enum
from os.path import basename, join
import xml.etree.ElementTree as ET

from emtable import Table

from pwem.protocols import EMProtocol
from pyseg.convert import readPysegCoordinates
from pyworkflow.object import String
from pyworkflow.protocol import FloatParam, EnumParam, PointerParam, IntParam, LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.utils import Message, removeBaseExt, copyFile, moveFile
from scipion.constants import PYTHON
from tomo.objects import SetOfCoordinates3D, SetOfTomograms, SetOfTomoMasks
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition
from pyseg import Plugin
from pyseg.constants import FILS_SOURCES, FILS_TARGETS, PICKING_SCRIPT, PICKING_SLICES, PRESEG_AREAS_LIST, MEMBRANE, \
    OUT_STARS_DIR, IN_STARS_DIR, FILS_OUT, PICKING_OUT, VESICLE
from pyseg.utils import encodePresegArea, getPrevPysegProtOutStarFiles, createStarDirectories
from tomo.utils import getObjFromRelation

logger = logging.getLogger(__name__)

# Fils slices xml fields
SIDE = 'side'
CONT = 'cont'

# CONT param codification
CUTTING_POINT = 0
PROJECTIONS = 1


class outputObjects(Enum):
    coordinates = SetOfCoordinates3D


class ProtPySegPicking(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """extract particles from a filament network of a oriented single membrane graph"""

    _label = 'picking'
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        self._tomoSet = None
        self.tomoSet = None
        self.acquisitionParams = {
                'angleMin': 90,
                'angleMax': -90,
                'step': None
        }
        self._inStarDir = None
        self._outStarDir = None
        self._xmlSlices = None
        self._outStarFilesList = []
        self.failedVesicles = String('')

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inFilsProt', PointerParam,
                      pointerClass='ProtPySegFils',
                      label='Fils protocol',
                      important=True,
                      allowsNull=False,
                      help='Pointer to fils protocol.')
        form.addParam('inTomoSet', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms to refer the coordinates',
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='Tomograms to which the coordinates should be referred to. If empty, it is assumed that '
                           'the tomograms are the same from were the vesicles were segmented (pre-seg).')
        form.addParam('boxSize', IntParam,
                      label='Box size (pixels)',
                      default=20,
                      important=True,
                      allowsNull=False,
                      expertLevel=LEVEL_ADVANCED)

        form.addSection(label='Picking')
        self._defineFilsXMLParams(form, self._getSlicesXMLDefaultVals())

        form.addSection(label='Refinement')
        group = form.addGroup('Peaks cleaning')
        group.addParam('peakTh', FloatParam,
                       label='Percentile of points to discard by their density level.',
                       default=0,
                       allowsNull=False)
        group.addParam('peakNs', FloatParam,
                       label='Min distance between selected points (nm).',
                       default=0.5,
                       allowsNull=False,
                       help='Scale suppression in nm, two selected points cannot be closer than this distance.')

        form.addParallelSection(threads=3, mpi=1)

    @staticmethod
    def _defineFilsXMLParams(form, d):
        """d is a disctionary with the default values"""
        paramList = list(d.keys())
        valList = list(d.values())
        form.addParam(paramList[0], EnumParam,
                      label='Segmentation area for picking',
                      choices=PRESEG_AREAS_LIST,
                      default=valList[0],
                      allowsNull=False,
                      help='Area in which the cutting point or cutting point + projections of the '
                           'filament will be considered for the picking coordinates.')
        form.addParam(paramList[1], EnumParam,
                      choices=['Cutting point', 'Projected local minima'],
                      default=0,
                      label='Find on two surfaces',
                      display=EnumParam.DISPLAY_HLIST,
                      help="Track fiducials differentiating in which side of the sample are located.")

    @staticmethod
    def _getSlicesXMLDefaultVals():
        d = OrderedDict()
        d[SIDE] = MEMBRANE
        d[CONT] = CUTTING_POINT
        return d

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        allOutputId = []
        starFileList = self._convertInputStep()
        for starFile in starFileList:
            pId = self._insertFunctionStep(self.pysegPicking, starFile,
                                           prerequisites=[],
                                           needsGPU=False)
            allOutputId.append(pId)
        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=allOutputId,
                                 needsGPU=False)

    def _convertInputStep(self):
        outDir = self._getExtraPath()
        # Generate directories for input and output star files
        self._outStarDir, self._inStarDir = createStarDirectories(outDir)
        # Generate slices xml
        self._createPickingXmlFile(Plugin.getHome(PICKING_SLICES), outDir)
        # Get the star files generated in the fils protocol
        return getPrevPysegProtOutStarFiles(self.inFilsProt.get()._getExtraPath(OUT_STARS_DIR),
                                            self._getExtraPath(IN_STARS_DIR))

    def pysegPicking(self, starFile):
        try:
            # Script called
            Plugin.runPySeg(self, PYTHON, self._getPickingCommand(starFile))
            # Move output files to the corresponding directory
            outFile = self._getExtraPath(removeBaseExt(starFile) + '_parts.star')
            newFileName = join(self._outStarDir, basename(outFile).replace(FILS_OUT, PICKING_OUT))
            self._outStarFilesList.append(newFileName)
            moveFile(outFile, newFileName)
        except Exception as e:
            inStarTable = Table()
            inStarTable.read(starFile)
            vesicleName = removeBaseExt(inStarTable[0].get(VESICLE, None))  # Files of only one line
            prevMsg = self.failedVesicles.get() if self.failedVesicles.get() else ''
            self.failedVesicles.set(prevMsg + f', {vesicleName}')
            with self._lock:
                self._store(self.failedVesicles)
            logger.error(f'{e}')

    def createOutputStep(self):
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(self._getTomoSet(), suffix)
        coordsSet.setSamplingRate(self._getSamplingRate())
        coordsSet.setBoxSize(self.boxSize.get())

        # Read the data from all the out star files
        tomoSet = self._getTomoSet()
        for outStar in self._outStarFilesList:
            readPysegCoordinates(outStar, coordsSet, tomoSet)

        if not coordsSet:
            raise Exception('ERROR! No coordinates were picked.')
        self._defineOutputs(**{outputObjects.coordinates.name: coordsSet})
        self._defineSourceRelation(self._getTomoSet(), coordsSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        valMsg = []
        if not self._getTomoFromRelations():
            valMsg.append("Unable to find the corresponding tomograms using the relations.")
        return valMsg

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if self.isFinished():
            outputCoords = getattr(self, outputObjects.coordinates.name, None)
            summary.append('*Picking*:\n\t- Picking area = %s\n\t- Particles picked = %i\n' %
                           (PRESEG_AREAS_LIST[int(self.side.get())], outputCoords.getSize()))

        return summary

    # --------------------------- UTIL functions -----------------------------------
    def _getFilsStarFileName(self):
        prot = self.inFilsProt.get()
        return prot._getExtraPath('fil_' + removeBaseExt(FILS_SOURCES)
                                  + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')

    def _getPickingStarFileName(self):
        filsStar = self._getExtraPath('fil_' + removeBaseExt(FILS_SOURCES)
                                      + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')
        return self._getExtraPath(removeBaseExt(filsStar) + '_parts.star')

    def _getPickingCommand(self, starFile):
        pickingCmd = ' '
        pickingCmd += '%s ' % Plugin.getHome(PICKING_SCRIPT)
        pickingCmd += '--inStar %s ' % starFile
        pickingCmd += '--outDir %s ' % self._getExtraPath()
        pickingCmd += '--slicesFile %s ' % self._xmlSlices
        pickingCmd += '--peakTh %s ' % self.peakTh.get()
        pickingCmd += '--peakNs %s ' % self.peakNs.get()
        return pickingCmd

    def _createPickingXmlFile(self, templateFile, outDir):
        # Copy the template xml to extra
        filename = join(outDir, basename(templateFile))
        self._xmlSlices = filename
        copyFile(templateFile, filename)

        # Edit the corresponding fields
        xmlTree = ET.parse(filename)
        rootElement = xmlTree.getroot()
        mb_slice = rootElement.findall("mb_slice")
        mb_slice = mb_slice[0]
        mb_slice.find(SIDE).text = str(encodePresegArea(getattr(self, SIDE).get()))
        mb_slice.find(CONT).text = self._decodeContValue(getattr(self, CONT).get())

        # Write the modified xml file.
        xmlTree.write(filename, encoding='UTF-8', xml_declaration=True)

    @staticmethod
    def _decodeContValue(val):
        """Decode the cont values and represent them as expected by pySeg"""
        return '+' if val == CUTTING_POINT else '-'

    def _getSamplingRate(self):
        return self._getTomoSet().getFirstItem().getSamplingRate()

    def _getTomoSet(self):
        if self._tomoSet is None:
            if self.inTomoSet.get():
                self._tomoSet = self.inTomoSet.get()
            else:
                self._tomoSet = self._getTomoFromRelations()
        return self._tomoSet

    def _getTomoFromRelations(self):
        # Get the tomograms climbing from this point of the workflow until the pre-seg and if there aren't tomograms
        # at that point, use the relations to go to the corresponding tomograms
        inObj = self.inFilsProt.get().inGraphsProt.get().inSegProt.get()
        tomoMaks = inObj if type(inObj) is SetOfTomoMasks else inObj.inTomoMasks.get()
        tomoSet = getObjFromRelation(tomoMaks, self, SetOfTomograms)
        return tomoSet
