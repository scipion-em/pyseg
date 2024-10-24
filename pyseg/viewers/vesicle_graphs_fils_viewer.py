# **************************************************************************
# *
# * Authors:     Scipion Team
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
import glob
from os.path import join

from pyseg.constants import FILS_FILES
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.utils.path import removeBaseExt
from tomoviz.viewers.viewer_vtk import VtkPlot
from tomoviz.viewers.viewer_triangulations import guiThread

FROM_GRAPHS = 0
FROM_FILS = 1
FROM_PICKING = 2


class VesicleViewerDialog(ToolbarListDialog):
    """
    This class allows to  call a MembraneAnnotator subprocess from a list of Tomograms.
    """

    def __init__(self, parent, vtiPath, areTomoMasks=False, **kwargs):
        self.vtiPath = vtiPath
        self.provider = kwargs.get("provider", None)
        self.prot = kwargs.get('prot', None)
        self.source = kwargs.get('source', None)
        self.areTomoMasks = areTomoMasks
        ToolbarListDialog.__init__(self, parent,
                                   "Vesicle Visualization Object Manager",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.launchVesicleViewer,
                                   allowSelect=False,
                                   **kwargs)

    def launchVesicleViewer(self, inVesicle):
        print("\n==> Running Vesicle Viewer:")
        fileName = inVesicle.getVolName() if self.areTomoMasks else inVesicle.getFileName()
        vesicleBaseName = removeBaseExt(fileName)
        vtiName = join(self.vtiPath, vesicleBaseName + '.vti')
        args = {'vti_file': vtiName}
        if self.source == FROM_GRAPHS:
            args['graph_file'] = glob.glob(self.prot._getExtraPath(vesicleBaseName + '*_edges_2.vtp'))[0]
        elif self.source == FROM_FILS:
            args['net_file'] = glob.glob(self.prot._getExtraPath(FILS_FILES, '*', vesicleBaseName + '*_net.vtp'))[0]
        elif self.source == FROM_PICKING:
            args['peaks_file'] = glob.glob(self.prot._getExtraPath(vesicleBaseName + '*_peak.vtp'))[0]
            args['surf_file'] = glob.glob(self.prot._getExtraPath(vesicleBaseName + '*_surf.vtp'))[0]

        guiThread(VtkPlot, 'initializePlot', **args)
