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
import shutil
from glob import glob
from os.path import basename
from emtable import Table
from pwem.protocols import EMProtocol
from pyseg.convert.convert import splitPysegStarFile
from pyseg.protocols.protocol_pre_seg import outputObjects as presegOutputObjects, ProtPySegPreSegParticles
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyseg.utils import createStarDirectories, genOutSplitStarFileName
from pyworkflow.protocol import FloatParam, PointerParam, LEVEL_ADVANCED, BooleanParam, IntParam
from pyworkflow.utils import Message, moveFile
from scipion.constants import PYTHON
from tomo.objects import SetOfTomoMasks
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition
from pyseg import Plugin
from pyseg.constants import GRAPHS_SCRIPT, SEGMENTATION

logger = logging.getLogger(__name__)


class ProtPySegGraphs(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """analyze a GraphMCF (Mean Cumulative Function) from a segmented membrane"""

    _label = 'graphs'
    stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._outStarDir = None
        self._inStarDir = None
        self.starFileList = None

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSegProt', PointerParam,
                      pointerClass='ProtPySegPreSegParticles, SetOfTomoMasks',
                      label='Pre segmentation',
                      important=True,
                      allowsNull=False,
                      help='Pointer to preseg protocol.')
        form.addParam('vesiclePkgSize', IntParam,
                      label='Vesicles packaging size',
                      allowsNull=False,
                      help='The input set of particles will be split into packages of N vesicles. Each package will '
                           'be processed as a different step, allowing to continue the execution from the last step in '
                           'case of the protocol fails. On the other hand, more packages implies more calls to PySeg, '
                           'which can affect to performance.')
        form.addParam('keepOnlyReqFiles', BooleanParam,
                      label='Keep only required files?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to No, all the intermediate Disperse program resulting directories '
                           'will be kept in the extra folder.')

        group = form.addGroup('Graphs parameters')
        group.addParam('sSig', FloatParam,
                       label='Sigma for gaussian filtering',
                       default=1,
                       allowsNull=False,
                       help='Sigma for Gaussian foltering input tomograms. It allows to smooth '
                            'small and irrelevant features and increases teh signal noise ratio (SNR). '
                            'Higher values will provide less dense graphs (lower execution time), so they should be '
                            'used when picking large particles, like ribosomes.')
        group.addParam('vDen', FloatParam,
                       label='Vertex density within membranes (nm³)',
                       default=0.0035,
                       allowsNull=False,
                       expertLevel=LEVEL_ADVANCED,
                       help='Vertex density within membranes. It allows to adjust simplification '
                            'adaptively for every tomogram.')
        group.addParam('vRatio', FloatParam,
                       label='Avg ratio vertex/edge of graph within membrane',
                       default=4,
                       allowsNull=False,
                       expertLevel=LEVEL_ADVANCED,
                       help='Averaged ratio vertex/edge in the graph within membrane.')
        group.addParam('maxLen', FloatParam,
                       label='Maximum distance to membrane (Å)',
                       allowsNull=False,
                       help='Maximum euclidean distance to membrane in Å.')
        form.addParam('binThreads', IntParam,
                      label='PySeg threads',
                      default=2,
                      help='Number of threads used by EMAN each time it is called in the protocol execution. For '
                           'example, if 2 Scipion threads and 3 PySeg threads are set, the tomograms will be '
                           'processed in groups of 2 at the same time with a call of PySeg with 3 threads each, so '
                           '6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of tomograms specified by Scipion threads.')
        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        self._initialize()
        stepIdList = []
        for starFile in self.starFileList:
            stepId = self._insertFunctionStep(self.pysegGraphs, starFile,
                                              prerequisites=[],
                                              needsGPU=False)
            stepIdList.append(stepId)
        self._insertFunctionStep(self.removeUnusedFilesStep,
                                 prerequisites=stepIdList,
                                 needsGPU=False)

    def _initialize(self):
        # Generate directories for input and output star files
        # Split the input file into n (threads) files
        inStar = self._getPreSegStarFile()
        self._outStarDir, self._inStarDir = createStarDirectories(self._getExtraPath())
        subsetStarFile = self._manageInStarFile(inStar)
        inStar = subsetStarFile if subsetStarFile else inStar
        self.starFileList = splitPysegStarFile(inStar, self._inStarDir, j=self.vesiclePkgSize.get())

    def pysegGraphs(self, starFile):
        # Script called
        Plugin.runPySeg(self, PYTHON, self._getGraphsCommand(starFile))
        # Fils returns the same star file name, so it will be renamed to avoid overwriting
        moveFile(self._getExtraPath(basename(starFile)).replace('.star', '_mb_graph.star'),
                 genOutSplitStarFileName(self._outStarDir, starFile))

    def removeUnusedFilesStep(self):
        # Remove Disperse program intermediate result directories if requested
        if self.keepOnlyReqFiles.get():
            disperseDirs = glob(self._getExtraPath('disperse_*'))
            [shutil.rmtree(disperseDir) for disperseDir in disperseDirs]

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summaryMsg = []
        if self.isFinished():
            summaryMsg.append('Graphs were correctly generated.')

    # --------------------------- UTIL functions -----------------------------------
    def _getGraphsCommand(self, starFile):
        graphsCmd = ' '
        graphsCmd += '%s ' % Plugin.getHome(GRAPHS_SCRIPT)
        graphsCmd += '--inStar %s ' % starFile
        graphsCmd += '--outDir %s ' % self._getExtraPath()
        graphsCmd += '--pixelSize %s ' % (self._getSamplingRate() / 10)  # PySeg requires it in nm
        graphsCmd += '--sSig %s ' % self.sSig.get()
        graphsCmd += '--vDen %s ' % self.vDen.get()
        graphsCmd += '--veRatio %s ' % self.vRatio.get()
        graphsCmd += '--maxLen %s ' % (self.maxLen.get() / 10)  # PySeg requires it in nm
        graphsCmd += '-j %s ' % self.binThreads.get()
        return graphsCmd

    def _getPreSegStarFile(self):
        inputObj = self.inSegProt.get()
        if type(inputObj) is ProtPySegPreSegParticles:
            return inputObj.getPresegOutputFile(inputObj.getVesiclesCenteredStarFile())
        else:
            # Set of TomoMasks
            return inputObj._pysegPresegVesicleCStar.get()

    def _getSamplingRate(self):
        inputObj = self.inSegProt.get()
        if type(inputObj) is ProtPySegPreSegParticles:
            inVesicles = getattr(inputObj, presegOutputObjects.vesicles.name)
            return inVesicles.getSamplingRate()
        else:
            return inputObj.getSamplingRate()

    def _manageInStarFile(self, inStar):
        """Since the graphs protocol admits a set of tomo masks, there may have been subset actions in between. In
        that case, a new star file containing only the present masked vesicles (tomo masks) must be generated.
        """
        outStarFile = ''
        inObj = self.inSegProt.get()
        if type(inObj) is SetOfTomoMasks:
            tomoTable = Table()
            tomoTable.read(inStar)
            if len(tomoTable) != len(inObj):
                outStarFile = self._getExtraPath('inSegVesicles.star')
                labels = tomoTable.getColumnNames()
                outTable = Table(columns=labels)
                tableSegVesicleFileDict = {row.get(SEGMENTATION): i for i, row in enumerate(tomoTable)}
                logger.info('==> Subset detected. Writing a new star file with the selected sefmented vesicles.')
                for segVesicle in inObj:
                    matchRowInd = tableSegVesicleFileDict[segVesicle.getFileName()]
                    outTable._rows.append(tomoTable[matchRowInd])
                    outTable.write(outStarFile)
        return outStarFile