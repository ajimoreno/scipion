# **************************************************************************
# *
# * Authors:  Amaya Jimenez (ajimenez@cnb.csic.es)
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

from pyworkflow.em.metadata.utils import getSize
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtAnalysis2D
from pyworkflow import VERSION_1_2
from convert import writeSetOfParticles
import pyworkflow.em.metadata as md
from os import remove
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles
from scipy.interpolate import UnivariateSpline


class XmippProtParticlePolishing(ProtAnalysis2D):
    """Protocol to make particle polishing"""
    _label = 'particle polishing'
    _lastUpdateVersion = VERSION_1_2

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputMovieParticles', params.PointerParam,
                      pointerClass='SetOfMovieParticles',
                      label="Input movie particles",
                      important=True,
                      help='Select the input movie particles')
        form.addParam('numberOfFrames', params.IntParam,
                      label="Number of movie frames",
                      important=True,
                      help='Number of frames in the input movie particle set')
        form.addParam('invertContrast', params.BooleanParam,
                      label='Invert contrast?',
                      default = False,
                      help='Decide if to invert contrast or not in the output particles')

        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        self.inputParticles = self._getExtraPath('inputParticles.xmd')

        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('processingStep')
        if (self.invertContrast):
            self._insertFunctionStep('invertContrastStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------

    def convertStep(self):
        writeSetOfParticles(self.inputMovieParticles.get(), self.inputParticles)

    def processingStep(self):

        numberOfParticles = long(getSize(self.inputParticles))/long(self.numberOfFrames)
        mdStackOut = md.MetaData()

        for j in range(numberOfParticles):

            args = '-i %s ' % (self.inputParticles)
            args += ' --query select \"itemId%' + '%d==%d\" ' % (numberOfParticles,j)
            args += ' -o %s' %(self._getExtraPath('movie%d.xmd'%j))
            self.runJob('xmipp_metadata_utilities', args, numberOfMpi=1)

            args = '-i %s ' %(self._getExtraPath('movie%d.xmd'%j))
            args += ' -o %s' % (self._getExtraPath('alignParticle%d.xmd' % j))
            self.runJob('xmipp_movie_alignment_correlation', args, numberOfMpi=1)

            mdAlignFn = self._getExtraPath('alignParticle%d.xmd' % j)
            mdAlign = md.MetaData(mdAlignFn)
            fnStack = self._getExtraPath('alignParticle%d.stk' % j)
            list_shifts_x=[]
            list_shifts_y=[]
            count=0
            for i, row in enumerate(md.iterRows(mdAlign)):
                count+=1
                list_shifts_x.append(row.getValue(md.MDL_SHIFT_X))
                list_shifts_y.append(row.getValue(md.MDL_SHIFT_Y))
            spl_shifts_x = UnivariateSpline(range(count), list_shifts_x)
            spl_shifts_y = UnivariateSpline(range(count), list_shifts_y)
            spl_shifts_x.set_smoothing_factor(len(list_shifts_x)*10)
            spl_shifts_y.set_smoothing_factor(len(list_shifts_y)*10)
            smooth_shifts_x = spl_shifts_x(range(count))
            smooth_shifts_y = spl_shifts_y(range(count))
            print("list_shifts_x", list_shifts_x)
            print("smooth_shifts_x", smooth_shifts_x)
            print("list_shifts_y", list_shifts_y)
            print("smooth_shifts_y", smooth_shifts_y)

            for i, row in enumerate(md.iterRows(mdAlign)):
                args = '-i %s' % row.getValue(md.MDL_IMAGE)
                args += ' --shift %d %d ' % (smooth_shifts_x[i], smooth_shifts_y[i])
                args += ' -o %06d@%s' % (i + 1, fnStack)
                self.runJob('xmipp_transform_geometry', args, numberOfMpi=1)

            total=0
            fnStackOut = self._getExtraPath('average.stk')
            for i, row in enumerate(md.iterRows(mdAlign)):
                total+=1
                if i==0:
                    args = '-i %06d@%s' % (1, fnStack)
                    args += ' --plus 0'
                    args += ' -o %06d@%s'%(j+1, fnStackOut)
                else:
                    args = '-i %06d@%s'%(i+1, fnStack)
                    args += ' --plus %06d@%s' % (j+1, fnStackOut)
                    args += ' -o %06d@%s'%(j+1, fnStackOut)
                self.runJob('xmipp_image_operate', args, numberOfMpi=1)

            args = ' -i %06d@%s'%(j+1, fnStackOut)
            args += ' --divide %d' % total
            args += ' -o %06d@%s'%(j+1, fnStackOut)
            self.runJob('xmipp_image_operate', args, numberOfMpi=1)

            row = md.Row()
            row.setValue(md.MDL_IMAGE, '%06d@%s'%(j+1, fnStackOut))
            row.addToMd(mdStackOut)
            remove(mdAlignFn)
            remove(fnStack)

        mdStackOut.write(self._getExtraPath('average.xmd'), md.MD_APPEND)


    def invertContrastStep(self):

        args = ' -i %s' % self._getExtraPath('average.xmd')
        args += ' -o %s' % self._getExtraPath('output_images.stk')
        args += ' --save_metadata_stack %s' % self._getExtraPath('output_images.xmd')
        args += ' --keep_input_columns --mult -1'
        self.runJob('xmipp_image_operate', args, numberOfMpi=1)


    def createOutputStep(self):
        if not self.invertContrast:
            fnOut = self._getExtraPath('average.xmd')
        else:
            fnOut = self._getExtraPath('output_images.xmd')
        outMovieSet = self._createSetOfParticles()
        outMovieSet.setSamplingRate(self.inputMovieParticles.get().getSamplingRate())
        readSetOfParticles(fnOut, outMovieSet)

        result = {'outputParticles': outMovieSet}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputMovieParticles, outMovieSet)