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

from collections import OrderedDict
from os.path import basename, join
import xml.etree.ElementTree as ET

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FloatParam, EnumParam, PointerParam, IntParam, FileParam, LEVEL_ADVANCED
from pyworkflow.utils import Message, makePath, removeBaseExt, copyFile
from scipion.constants import PYTHON
from tomo.objects import SetOfCoordinates3D, SetOfTomograms
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition

from pyseg import Plugin
from pyseg.constants import FILS_SOURCES, FILS_TARGETS, PICKING_SCRIPT, PICKING_SLICES, PRESEG_AREAS_LIST, MEMBRANE
from pyseg.convert import readStarFile, PYSEG_PICKING_STAR

# Fils slices xml fields
from pyseg.utils import encodePresegArea

SIDE = 'side'
CONT = 'cont'

# CONT param codification
CUTTING_POINT = 0
PROJECTIONS = 1


class ProtPySegPicking(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """extract particles from a filament network of a oriented single membrane graph"""

    _label = 'picking'
    _devStatus = BETA
    tomoSet = None
    acquisitionParams = {
            'angleMin': 90,
            'angleMax': -90,
            'step': None,
            'angleAxis1': None,
            'angleAxis2': None
        }

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
                      label='Tomograms',
                      important=True,
                      allowsNull=False,
                      help='Tomograms used for the graphs and filaments calculation.')
        form.addParam('boxSize', IntParam,
                      label='Box size (pixels)',
                      default=20,
                      important=True,
                      allowsNull=False,
                      expertLevel=LEVEL_ADVANCED,)

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
                      choices=['Cutting point', 'Cutting point + projections'],
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

    def _insertAllSteps(self):
        self._insertFunctionStep(self.pysegPicking)
        self._insertFunctionStep(self.createOutputStep)

    def pysegPicking(self):
        # Generate output dir
        outDir = self._getExtraPath()
        makePath(outDir)

        # Generate slices xml
        self._createPickingXmlFile(Plugin.getHome(PICKING_SLICES), outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPickingCommand(outDir))

    def createOutputStep(self):
        self.tomoSet = self.inTomoSet.get()  # Required for the convert
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(self.inTomoSet.get(), suffix)
        coordsSet.setSamplingRate(self._getSamplingRate())
        coordsSet.setBoxSize(self.boxSize.get())
        readStarFile(self, coordsSet, PYSEG_PICKING_STAR,
                     starFile=self._getPickingStarFileName(), invert=True)

        self._defineOutputs(outputCoordinates=coordsSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if self.isFinished():
            outputCoords = getattr(self, 'outputCoordinates', None)
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

    def _getPickingCommand(self, outDir):
        pickingCmd = ' '
        pickingCmd += '%s ' % Plugin.getHome(PICKING_SCRIPT)
        pickingCmd += '--inStar %s ' % self._getFilsStarFileName()
        pickingCmd += '--outDir %s ' % outDir
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
        return self.inTomoSet.get().getFirstItem().getSamplingRate()
