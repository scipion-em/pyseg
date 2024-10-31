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
from glob import glob
from os.path import exists
from imod.protocols import ProtImodTomoNormalization
from pyseg.protocols import ProtPySegGraphs, ProtPySegFils
from pyseg.protocols.protocol_picking import PROJECTIONS, ProtPySegPicking
from pyseg.protocols.protocol_picking import outputObjects as pickingOutputs
from pyseg.protocols.protocol_pre_seg import outputObjects as presegOutputs, ProtPySegPreSegParticles
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from pyseg.constants import FROM_SCIPION, MEMBRANE_OUTER_SURROUNDINGS, MEMBRANE, OUT_STARS_DIR, FILS_FILES
from tomo.protocols import ProtImportTomograms, ProtImportTomomasks
from tomo.tests import EMD_10439, DataSetEmd10439


class TestFromPresegToPicking(BaseTest):

    outputPath = None
    ds = None
    boxSize = 44
    samplingRate = 13.68
    nVesicles = 3
    vesiclesPackagesSize = 3
    inTomoSet = None
    inTomoSetBinned = None
    inTomomaskSetBinned = None
    protPreseg = None
    ProtGraphs = None
    ProtFils = None
    ProtPicking = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        ds = DataSet.getDataSet(EMD_10439)
        cls.ds = ds
        cls.inTomoSet = cls._importTomograms()
        cls.inTomoSetBinned = cls._normalizeTomo()
        cls.inTomomaskSetBinned = cls._ImportTomoMasks()
        cls.protPreseg = cls._runPreseg()
        cls.ProtGraphs = cls._runGraphs()
        cls.ProtFils = cls._runFils()
        cls.ProtPicking = cls._runPicking()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.value),
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, 'Tomograms', None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated.')

        return outputTomos

    @classmethod
    def _normalizeTomo(cls):
        print(magentaStr("\n==> Tomogram normalization:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=cls.inTomoSet,
                                                binning=2)

        cls.launchProtocol(protTomoNormalization)
        outputTomos = getattr(protTomoNormalization, 'Tomograms', None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated in tomo normalization.')

        return outputTomos

    @classmethod
    def _ImportTomoMasks(cls):
        print(magentaStr("\n==> Importing data - tomoMasks"
                         ":"))
        protImportTomomasks = cls.newProtocol(ProtImportTomomasks,
                                              filesPath=cls.ds.getFile(DataSetEmd10439.tomoMaskAnnotated.value),
                                              inputTomos=cls.inTomoSetBinned)

        cls.launchProtocol(protImportTomomasks)
        tomoMaskSet = getattr(protImportTomomasks, protImportTomomasks._possibleOutputs.tomomasks.name, None)
        cls.assertIsNotNone(tomoMaskSet, 'No tomograms were genetated.')

        return tomoMaskSet

    @classmethod
    def _runPreseg(cls):
        print(magentaStr("\n==> Running preSeg:"))
        protPreseg = cls.newProtocol(
            ProtPySegPreSegParticles,
            inTomoMasks=cls.inTomomaskSetBinned,
            spOffVoxels=22,
            sgMembThk=60,
            sgMembNeigh=330
        )
        protPreseg.setObjLabel('Preseg')
        protPreseg = cls.launchProtocol(protPreseg)

        return protPreseg

    def testPreseg(self):
        # Check that resulting files are created as expected
        protPreseg = self.protPreseg
        outputStar = 'presegVesiclesCentered_pre.star'
        outputVesiclesPattern = 'emd_10439_tid_%i.mrc'
        outputVesiclesSegPattern = 'emd_10439_tid_%i_seg.mrc'
        self.assertTrue(exists(protPreseg._getExtraPath(outputStar)))
        for i in range(self.nVesicles):
            self.assertTrue(exists(protPreseg._getExtraPath('segs', outputVesiclesPattern % i)))
            self.assertTrue(exists(protPreseg._getExtraPath('segs', outputVesiclesSegPattern % i)))

        # Check the generated outputs
        binnedSamplingRate = 2 * self.samplingRate
        vesicle0Size = (186, 178, 135)
        vesicle1Size = (282, 180, 143)
        vesicle2Size = (133, 133, 139)
        vesicleSizeList = [vesicle0Size, vesicle1Size, vesicle2Size]
        setOfVesiclesTomomasks = getattr(protPreseg, presegOutputs.segmentations.name, None)
        setOfVesiclesSubtomograms = getattr(protPreseg, presegOutputs.vesicles.name, None)
        self.assertSetSize(setOfVesiclesTomomasks, self.nVesicles)
        self.assertSetSize(setOfVesiclesSubtomograms, self.nVesicles)
        self.assertEqual(setOfVesiclesTomomasks.getSamplingRate(), binnedSamplingRate)
        self.assertEqual(setOfVesiclesSubtomograms.getSamplingRate(), binnedSamplingRate)
        for vesicle, vesicleMask in zip(setOfVesiclesSubtomograms, setOfVesiclesTomomasks):
            self.assertEqual(vesicle.getSamplingRate(), binnedSamplingRate)
            self.assertEqual(vesicleMask.getSamplingRate(), binnedSamplingRate)
            # Preseg tid assignment to the output varies from one execution to another
            self.assertTrue(vesicle.getDimensions() in vesicleSizeList)
            self.assertTrue(vesicleMask.getDimensions() in vesicleSizeList)
        return protPreseg

    @classmethod
    def _runGraphs(cls):
        print(magentaStr("\n==> Running graphs:"))
        protGraphs = cls.newProtocol(
            ProtPySegGraphs,
            presegFrom=FROM_SCIPION,
            inSegProt=cls.protPreseg,
            vesiclePkgSize=cls.vesiclesPackagesSize,
            maxLen=330,
            binThreads=cls.nVesicles,
            numberOfThreads=1
        )

        protGraphs.setObjLabel('Graphs')
        protGraphs = cls.launchProtocol(protGraphs)
        return protGraphs

    def testGraphs(self):
        # Check that resulting files are created as expected
        # Vesicles were packed into packages of 3 (n of threads), so 1 star files should have
        # been generated
        self.assertTrue(exists(self.ProtGraphs._getExtraPath(OUT_STARS_DIR, 'graphs_%03d.star' % 1)))
        # By default, the Disperse program intermediate results directories aren't kept
        self.assertTrue(not glob(self.ProtGraphs._getExtraPath('disperse_*')))

    @classmethod
    def _runFils(cls):
        print(magentaStr("\n==> Running fils:"))
        protFils = cls.newProtocol(
            ProtPySegFils,
            graphsFrom=FROM_SCIPION,
            inGraphsProt=cls.ProtGraphs,
            segLabelS=MEMBRANE,
            segLabelT=MEMBRANE_OUTER_SURROUNDINGS,
            maxEucDistT=30,
            maxGeoLenT=60,
            maxSinuT=2,
            gRgEud='1 30',
            gRgLen='1 60',
            gRgSin='0 2',
            numberOfThreads=1 + cls.nVesicles
        )

        protFils.setObjLabel('Fils')
        protFils = cls.launchProtocol(protFils)
        return protFils

    def testFils(self):
        # Check that resulting files are created as expected
        xmlFiles = ['mb_sources.xml', 'no_mb_targets.xml']
        filssFilePattern = 'emd_10439_tid_%i_fil_mb_sources_to_no_mb_targets_net'
        graphsFilesPerVesicle = [filssFilePattern + '.pkl',
                                 filssFilePattern + '.vtp',
                                 filssFilePattern + '_edges.vtp',
                                 filssFilePattern + '_edges2.vtp',
                                 filssFilePattern + '_skel.vtp']

        [self.assertTrue(exists(self.ProtFils._getExtraPath(file))) for file in xmlFiles]
        for i in range(self.nVesicles):
            self.assertTrue(exists(self.ProtFils._getExtraPath(OUT_STARS_DIR, 'fils_%03d.star' % (i + 1))))
            for file in graphsFilesPerVesicle:
                self.assertTrue(exists(self.ProtFils._getExtraPath(FILS_FILES, 'outDir_%03d' % i, file % i)))

    @classmethod
    def _runPicking(cls):
        print(magentaStr("\n==> Running picking:"))
        protPicking = cls.newProtocol(
            ProtPySegPicking,
            inFilsProt=cls.ProtFils,
            side=MEMBRANE_OUTER_SURROUNDINGS,
            cont=PROJECTIONS,
            numberOfThreads=1 + cls.nVesicles
        )

        protPicking.setObjLabel('Picking')
        protPicking = cls.launchProtocol(protPicking)
        return protPicking

    def testPicking(self):
        outputCoordinates = getattr(self.ProtPicking, pickingOutputs.coordinates.name)

        # Check that resulting files are created as expected
        xmlFile = 'mb_ext.xml'
        pickingFilePattern = 'emd_10439_tid_%i_fil_mb_sources_to_no_mb_targets_net_mb_ext'
        graphsFilesPerVesicle = [pickingFilePattern + '.vtp',
                                 pickingFilePattern + '_peak.vtp',
                                 pickingFilePattern + '_surf.vtp']

        self.assertTrue(exists(self.ProtPicking._getExtraPath(xmlFile)))
        for i in range(self.nVesicles):
            for file in graphsFilesPerVesicle:
                self.assertTrue(exists(self.ProtPicking._getExtraPath(file % i)))
                self.assertTrue(exists(self.ProtPicking._getExtraPath(OUT_STARS_DIR, 'picking_%03d_parts.star' % (i + 1))))

        # Output coordinates must have an attribute named _groupId, which stores the corresponding
        # vesicles index. In this case there are 3 vesicles, so it should be contained in range(2)
        testVesicleInd = list(range(self.nVesicles))
        self.assertEqual(outputCoordinates.getBoxSize(), 20)
        for coord3d in outputCoordinates:
            self.assertTrue(int(coord3d.getGroupId()) in testVesicleInd)


