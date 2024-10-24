# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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

import pyworkflow.viewer as pwviewer
import pwem.viewers.views as vi

from .vesicle_graphs_fils_viewer import VesicleViewerDialog, FROM_GRAPHS, FROM_FILS, FROM_PICKING
from .vesicle_visualization_tree import VesicleViewerProvider
from ..protocols import ProtPySegGraphs, ProtPySegFils, ProtPySegPicking, ProtPySegPreSegParticles


class TomoViz4PysegDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    using pyvista
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        ProtPySegGraphs,
        ProtPySegFils,
        ProtPySegPicking
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        vesicleSubTomos = []
        cls = type(obj)
        areTomoMasks = False
        if hasattr(obj, 'inSegProt'):  # Accessing from graphs
            vtiPath = obj._getExtraPath()
            inObj = obj.inSegProt.get()
            if type(inObj) is ProtPySegPreSegParticles:
                vesicleSubTomos = inObj.vesicles
            else:
                vesicleSubTomos = inObj
                areTomoMasks = True
        elif hasattr(obj, 'inGraphsProt'):  # Accessing from fils
            vtiPath = obj.inGraphsProt.get()._getExtraPath()
            inObj = obj.inGraphsProt.get().inSegProt.get()
            if type(inObj) is ProtPySegPreSegParticles:
                vesicleSubTomos = inObj.vesicles
            else:
                vesicleSubTomos = inObj
                areTomoMasks = True
        elif hasattr(obj, 'inFilsProt'):  # Accessing from picking
            vtiPath = obj.inFilsProt.get().inGraphsProt.get()._getExtraPath()
            vesicleSubTomos = obj.inFilsProt.get().inGraphsProt.get().inSegProt.get().vesicles

        vesicleFileList = [vesicle.clone() for vesicle in vesicleSubTomos.iterItems()]
        vesicleProvider = VesicleViewerProvider(vesicleFileList, None, None, areTomoMasks=areTomoMasks)

        if issubclass(cls, ProtPySegGraphs):
            source = FROM_GRAPHS
        elif issubclass(cls, ProtPySegFils):
            source = FROM_FILS
        elif issubclass(cls, ProtPySegPicking):
            source = FROM_PICKING
        VesicleViewerDialog(self._tkRoot, vtiPath=vtiPath, provider=vesicleProvider,
                            source=source, prot=obj, lockGui=False, areTomoMasks=areTomoMasks)
        return views
