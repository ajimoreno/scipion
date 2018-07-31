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
from pyworkflow.protocol.params import (PointerParam, IntParam,
                                        BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.em.data import Volume


class XmippProtVolumeGain(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'volume gain'
    _lastUpdateVersion = VERSION_1_2
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      help='Select a volume for equalizing level amplitude values.')

        form.addParam('mask', PointerParam, pointerClass='VolumeMask',
                      label="Binary Mask", important=True,
                      help='The mask determines which points are specimen'
                      ' and which are not.')

        form.addParam('useMonores', BooleanParam, default=True, label='Use MonoRes?',
                      help='Select if to use monoRes volume to perform histogram matching')

        form.addParam('monores', PointerParam, pointerClass='Volume',
                      label="MonoRes volume", condition="useMonores==True",
                      help='The monoRes volume to perform histogram matching by '
                           'frequency bands taking into account only pixel '
                           'values relevant to each band.')
        #
        # form.addParam('sampling', FloatParam, label='Sampling rate (A)',
        #               important = True, help='Provide the sampling rate in A of the volume')

        form.addParam('boxSize', IntParam, label='Box size (px)',
                      expertLevel=LEVEL_ADVANCED, default=5,
                      help='Size of the box to calculate every histogram')

        form.addParam('bandPass', BooleanParam, default=True, label='Band pass?',
                      expertLevel=LEVEL_ADVANCED,
                      help='Histogram matching by frequency bands or with the '
                           'complete frequencies.')

        form.addParam('nBands', IntParam, default=5, label='Number of bands',
                      expertLevel=LEVEL_ADVANCED,
                      condition='bandPass==True or useMonores==True',
                      help='Choose the number of bands in the band pass filter.')

        form.addParam('iter', IntParam, default=1, label='Iterations',
                      expertLevel=LEVEL_ADVANCED,
                      help='Number of iterations.')

        form.addParam('sigma', FloatParam, default=20.0, label='Sigma',
                      expertLevel=LEVEL_ADVANCED,
                      help='Sigma value for gaussian weighting for combining the resolution bands.')


    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
            # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep', )
        self._insertFunctionStep('volumeGainStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        """ Read the input volume.
        """
        self.volFn = self.inputVolume.get().getFileName()
        self.maskFn = self.mask.get().getFileName()
        if self.useMonores == True:
            self.monores = self.monores.get().getFileName()
        self.monoInputVol = self._getExtraPath('monogenicVol.vol')
        self.outputVol = self._getExtraPath('outputVol.vol')

    def volumeGainStep(self):

        # params = ' -i %s' % self.volFn
        # params += ' -o %s' % self.monoInputVol
        # params += ' --monogenic'
        # self.runJob('xmipp_transform_filter', params)

        params = ' -i %s ' % self.volFn #self.monoInputVol
        params += ' --mask %s ' % self.maskFn
        if self.useMonores==True:
            params += ' --mono %s ' % self.monores
        if self.bandPass == True or self.useMonores==True:
            params += ' --bandpass %i'% self.nBands.get()
        params += ' --sampling %f ' % self.inputVolume.get().getSamplingRate()
        params += ' --sigma %f ' % self.sigma.get()
        params += ' --boxSize %f ' % self.boxSize.get()
        params += ' -o %s ' %self.outputVol
        self.runJob('xmipp_volume_gain', params)

    def createOutputStep(self):
        # import time
        # time.sleep(20)
        volume=Volume()
        volume.setFileName(self._getExtraPath('outputVol.vol'))
        volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputVolume, volume)
            

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        summary.append("Equalizing level amplitudes in volume %s."
                       % self.inputVolume.get().getFileName())
        return summary

    def _validate(self):
        pass


