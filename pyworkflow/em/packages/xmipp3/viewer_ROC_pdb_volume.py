# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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


from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.protocol.params import LabelParam
from protocol_ROC_pdb_volume import XmippProtROCPdbVolume
import matplotlib.pyplot as plt
import numpy as np


        
class XmippProtROCPdbVolumeViewer(ProtocolViewer):
    """ Visualization of ROC pdb volume protocol results """

    _label = 'viewer ROC pdb vol'
    _targets = [XmippProtROCPdbVolume]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')

        form.addParam('doShowROC', LabelParam, label="Show ROC")

    def _getVisualizeDict(self):
        self.protocol._defineFileNames()  # Load filename templates
        return {'doShowROC': self.showROC }


    def showROC(self, param=None):

        fprFn = self.protocol._getFileName('fpr')
        fpr = np.loadtxt(fprFn, delimiter=',')

        tprFn = self.protocol._getFileName('tpr')
        tpr = np.loadtxt(tprFn, delimiter=',')

        plt.plot(fpr, tpr)
        plt.title("ROC curve")
        plt.show()