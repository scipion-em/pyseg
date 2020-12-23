from os import environ
from os.path import join

from pwem.protocols import EMProtocol, FileParam
from pyworkflow.protocol import FloatParam, NumericListParam, EnumParam
from pyworkflow.utils import Message, makePath, removeBaseExt, moveFile
from scipion.constants import PYTHON
from tomo.objects import SetOfCoordinates3D
from tomo.protocols import ProtTomoBase
from tomo.protocols.protocol_base import ProtTomoImportAcquisition

from pyseg import Plugin
from pyseg.constants import GRAPHS_OUT, GRAPHS_SCRIPT, FILS_OUT, FILS_SCRIPT, FILS_SOURCES, FILS_TARGETS, \
    PICKING_OUT, PICKING_SCRIPT, PICKING_SLICES
from pyseg.convert import readStarFile, PYSEG_PICKING_STAR, getTomoSetFromStar

TH_MODE_IN = 0
TH_MODE_OUT = 1


class ProtPySegGFP(EMProtocol, ProtTomoBase, ProtTomoImportAcquisition):
    """"""

    _label = 'Graphs-Fils-Picking'
    tomoSet = None
    warningMsg = None
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
        form.addParam('inStar', FileParam,
                      label='Seg particles star file',
                      important=True,
                      allowsNull=False,
                      help='Star file obtained in PySeg seg step.')
        form.addParam('pixelSize', FloatParam,
                      label='Pixel size (Å/voxel)',
                      default=1,
                      important=True,
                      allowsNull=False,
                      help='Input tomograms voxel size (Å/voxel)')

        form.addSection(label='Graphs')
        group = form.addGroup('GraphMCF')
        group.addParam('sSig', FloatParam,
                       label='Sigma for gaussian filtering',
                       default=1,
                       allowsNull=False,
                       help='Sigma for Gaussian fltering input tomograms. It allows to smooth '
                            'small and irrelevant features and increases SNR.')
        group.addParam('vDen', FloatParam,
                       label='Vertex density within membranes (nm³)',
                       default=0.0035,
                       allowsNull=False,
                       help='Vertex density within membranes. It allows to adjust simplifcation '
                            'adaptively for every tomogram.')
        group.addParam('vRatio', FloatParam,
                       label='Avg ratio vertex/edge of graph within membrane',
                       default=4,
                       allowsNull=False,
                       help='Averaged ratio vertex/edge in the graph within membrane.')
        group.addParam('maxLen', FloatParam,
                       label='Shortest distance to membrane (nm)',
                       default=10,
                       allowsNull=False,
                       help='Maximum euclidean shortest distance to membrane in nm.')

        form.addSection(label='Fils')
        group = form.addGroup('Graph thresholding')
        group.addParam('thMode', EnumParam,
                       default=0,
                       choices=['in', 'out'],
                       label='Orientation with respect to the membrane/filament',
                       display=EnumParam.DISPLAY_HLIST)
        group = form.addGroup('Filament geometry')
        group.addParam('gRgEud', NumericListParam,
                       label='Euclidean distance of vertices source-target (nm)',
                       default='1 25',
                       allowsNull=False,
                       help='Euclidean distance between source and target vertices in nm.')
        group.addParam('gRgLen', NumericListParam,
                       label='Geodesic distance of vertices source-target (nm)',
                       default='1 60',
                       allowsNull=False,
                       help='Geodesic distance trough the graph between source and target vertices in nm.')
        group.addParam('gRgSin', NumericListParam,
                       label='Filament sinuosity',
                       default='0 3',
                       allowsNull=False,
                       help='Filament sinuosity, geodesic/euclidean distances ratio.')

        form.addSection(label='Picking')
        group = form.addGroup('Peaks cleaning')
        group.addParam('peakTh', FloatParam,
                       label='Percentile of points to discard by their density level.',
                       default=0,
                       allowsNull=False)
        group.addParam('peakNs', FloatParam,
                       label='Min distance beyween selected points (nm).',
                       default=0.5,
                       allowsNull=False,
                       help='Scale suppression in nm, two selected points cannot be closer than this distance.')

        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        self._insertFunctionStep('pysegGraphs')
        self._insertFunctionStep('pysegFils')
        self._insertFunctionStep('pysegPicking')
        self._insertFunctionStep('createOutputStep')

    def pysegGraphs(self):
        # Generate output dir
        outDir = self._getExtraPath(GRAPHS_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getGraphsCommand(outDir))

    def pysegFils(self):
        # Generate output dir
        outDir = self._getExtraPath(FILS_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getFilsCommand(outDir))

    def pysegPicking(self):
        # Generate output dir
        outDir = self._getExtraPath(PICKING_OUT)
        makePath(outDir)

        # Script called
        Plugin.runPySeg(self, PYTHON, self._getPickingCommand(outDir))

    def createOutputStep(self):
        pickingStarFile = self.getPickingStarFileName()
        samplingRate = self.pixelSize.get()

        tomoSet = self._createSetOfTomograms()
        tomoSet.setSamplingRate(samplingRate)
        self.tomoSet = tomoSet
        getTomoSetFromStar(self, pickingStarFile)
        self._defineOutputs(outputTomograms=self.tomoSet)

        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coordsSet = self._createSetOfCoordinates3D(tomoSet, suffix)
        coordsSet.setSamplingRate(samplingRate)
        # coordsSet.setPrecedents(tomoSet)
        readStarFile(self, coordsSet, PYSEG_PICKING_STAR, starFile=pickingStarFile, invert=True)

        self._defineOutputs(outputCoordinates=coordsSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        pass
        # summary = []
        # if self.isFinished():
        #     summary.append('Generated files location:\n'
        #                    'Subtomograms files: *%s*\n'
        #                    'Star file: *%s*' %
        #                    (self._getExtraPath(POST_REC_OUT),
        #                     self._getExtraPath(POST_REC_OUT + '.star')))
        # return summary

    # --------------------------- UTIL functions -----------------------------------
    # @staticmethod
    # def _getPysegScript(scriptName):
    #     return Plugin.getHome(scriptName)
    #
    # @staticmethod
    # def _ListAsStr2ListOfNum(inList):
    #     """Convert string of integers separated by spaces to a list of integers"""
    #     return [int(i) for i in inList.split()]
    #
    def getPickingStarFileName(self):
        filsStar = self._getExtraPath(FILS_OUT, 'fil_' + removeBaseExt(FILS_SOURCES)
                                      + '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')
        return self._getExtraPath(PICKING_OUT, removeBaseExt(filsStar) + '_parts.star')

    def _getGraphsCommand(self, outDir):
        pixSize = self.pixelSize.get()/10
        graphsCmd = ' '
        graphsCmd += '%s ' % Plugin.getHome(GRAPHS_SCRIPT)
        graphsCmd += '--inStar %s ' % self.inStar.get()
        graphsCmd += '--outDir %s ' % outDir
        graphsCmd += '--pixelSize %s ' % pixSize  # PySeg requires it in nm
        graphsCmd += '--sSig %s ' % self.sSig.get()
        graphsCmd += '--vDen %s ' % self.vDen.get()
        graphsCmd += '--veRatio %s ' % self.vRatio.get()
        graphsCmd += '--maxLen %s ' % self.maxLen.get()
        graphsCmd += '-j %s ' % self.numberOfThreads.get()
        return graphsCmd

    def _getFilsCommand(self, outDir):
        filsCmd = ' '
        filsCmd += '%s ' % Plugin.getHome(FILS_SCRIPT)
        filsCmd += '--inStar %s ' % self._getExtraPath(GRAPHS_OUT, removeBaseExt(self.inStar.get()) + '_mb_graph.star')
        filsCmd += '--outDir %s ' % outDir
        filsCmd += '--inSources %s ' % Plugin.getHome(FILS_SOURCES)
        filsCmd += '--inTargets %s ' % Plugin.getHome(FILS_TARGETS)
        filsCmd += '--thMode %s ' % self._parseThModeSelection()
        filsCmd += '--gRgLen %s ' % self.gRgLen.get()
        filsCmd += '--gRgSin %s ' % self.gRgSin.get()
        filsCmd += '--gRgEud %s ' % self.gRgEud.get()
        return filsCmd

    def _getPickingCommand(self, outDir):
        pickingCmd = ' '
        pickingCmd += '%s ' % Plugin.getHome(PICKING_SCRIPT)
        pickingCmd += '--inStar %s ' % \
                      self._getExtraPath(FILS_OUT, 'fil_' + removeBaseExt(FILS_SOURCES) +
                                         '_to_' + removeBaseExt(FILS_TARGETS) + '_net.star')

        pickingCmd += '--outDir %s ' % outDir
        pickingCmd += '--slicesFile %s ' % Plugin.getHome(PICKING_SLICES)
        pickingCmd += '--peakTh %s ' % self.peakTh.get()
        pickingCmd += '--peakNs %s ' % self.peakNs.get()
        return pickingCmd

    def _parseThModeSelection(self):
        if self.thMode.get() == TH_MODE_IN:
            return 'in'
        else:
            return 'out'
