# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:    Amaya Jimenez Moreno (ajimenez@cnb.csic.es)
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
from pyworkflow import VERSION_1_2
from pyworkflow.protocol.params import PointerParam
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from shutil import copy
from xmipp3 import Image, DT_DOUBLE
import numpy as np
from sklearn import metrics

class XmippProtROCPdbVolume(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'ROC pdb volume'
    _lastUpdateVersion = VERSION_1_2
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)



    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'fpr': self._getExtraPath('fpr.txt'),
                  'tpr': self._getExtraPath('tpr.txt'),
                  'thresholds': self._getExtraPath('thresholds.txt')
                  }
        self._updateFilenamesDict(myDict)

    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True)

        form.addParam('inputPdb', PointerParam, pointerClass='PdbFile',
                      label="Input Pdb", important=True)


    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._defineFileNames()
        self._insertFunctionStep('calculateROC')

    def calculateROC(self):

        volume = self.inputPdb.get().getVolume()
        fnPdbVol = self._getExtraPath('pdbVolume')
        if self.inputPdb.get().getVolume() is None:
            params = ' -i %s --sampling %f -o %s --centerPDB --size %d ' % \
                     (self.inputPdb.get().getFileName(),
                      self.inputVolume.get().getSamplingRate(),
                      fnPdbVol, self.inputVolume.get().getDim()[0])
            self.runJob('xmipp_volume_from_pdb', params)
        else:
            copy(volume.getFileName(), fnPdbVol)

        fnPdbVol = self._getExtraPath('pdbVolume.vol')
        fnOrigVol = self.inputVolume.get().getFileName()
        params = "--i1 %s --i2 %s --apply %s --least_squares --local" % \
               (fnPdbVol, fnOrigVol, fnOrigVol)
        self.runJob("xmipp_volume_align", params)

        fnPdbVolThres = self._getExtraPath('pdbVolumeThres.vol')
        params = ' -i %s -o %s --select above %f --substitute value 1' % \
                 (fnPdbVol, fnPdbVolThres, 0.02)
        self.runJob('xmipp_transform_threshold', params)

        volPdb = Image(fnPdbVolThres)
        volPdb.convert2DataType(DT_DOUBLE)
        volOrig = Image(fnOrigVol)
        volOrig.convert2DataType(DT_DOUBLE)

        Xdim, Ydim, Zdim, Ndim = volPdb.getDimensions()
        labels = np.zeros((Xdim, Ydim, Zdim), dtype=np.int)
        labels[:,:,:] = volPdb.getData()
        labels = np.reshape(labels, Xdim*Ydim*Zdim)
        scores = np.zeros((Xdim,Ydim,Zdim),dtype=np.float64)
        scores[:,:,:] = volOrig.getData()
        scores = np.reshape(scores, Xdim * Ydim * Zdim)

        fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
        np.savetxt(self._getFileName('fpr'), fpr, delimiter=',')
        np.savetxt(self._getFileName('tpr'), tpr, delimiter=',')
        np.savetxt(self._getFileName('thresholds'), thresholds, delimiter=',')



    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        summary.append("Calculating ROC curve for PDB %s and volume %s."
                       % (self.inputPdb.get().getFileName(),
                          self.inputVolume.get().getFileName()))
        return summary

    def _validate(self):
        pass


