# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *              Javier Vargas (jvargas@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains the XmippProtExtractParticles protocol
"""

from glob import glob
from os.path import exists

import xmipp
from pyworkflow.object import String, Float
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER, ORIGINAL
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, STATUS_FINISHED
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, IntParam, 
                                        BooleanParam, RelationParam, Positive)
from pyworkflow.em.protocol import  ProtExtractParticles
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.constants import RELATION_CTF
from pyworkflow.utils.path import removeBaseExt, replaceBaseExt

from convert import writeSetOfCoordinates, readSetOfParticles, micrographToCTFParam
from xmipp3 import XmippProtocol

# Rejection method constants
REJECT_NONE = 0
REJECT_MAXZSCORE = 1
REJECT_PERCENTAGE = 2

                
class XmippProtExtractParticles(ProtExtractParticles, XmippProtocol):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'extract particles'
    
    def __init__(self, **args):
        ProtExtractParticles.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputCoordinates', PointerParam, label="Coordinates", important=True,
                      pointerClass='SetOfCoordinates',
                      help='Select the SetOfCoordinates ')
        
        form.addParam('downsampleType', EnumParam, choices=['same as picking', 'other', 'original'], 
                      default=0, important=True, label='Downsampling type', display=EnumParam.DISPLAY_COMBO, 
                      help='Select the downsampling type.')
        form.addParam('downFactor', FloatParam, default=2, condition='downsampleType==1',
                      label='Downsampling factor',
                      help='This factor is always referred to the original sampling rate. '
                      'You may use independent downsampling factors for extracting the '
                      'particles, picking them and estimating the CTF. All downsampling '
                      'factors are always referred to the original sampling rate, and '
                      'the differences are correctly handled by Xmipp.')        

        form.addParam('inputMicrographs', PointerParam, label="Micrographs", 
                      condition='downsampleType != 0',
                      pointerClass='SetOfMicrographs',
                      help='Select the original SetOfMicrographs')

        form.addParam('ctfRelations', RelationParam, allowsNull=True,
                      relationName=RELATION_CTF, attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input micrographs. \n'
                           'CTF estimation is need if you want to do phase flipping or \n'
                           'you want to associate CTF information to the particles.')

        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed particles, '
                      'actual particles may be smaller than this.')

        form.addParam('doSort', BooleanParam, default=False,
                      label='Perform sort by statistics',
                      help='Perform sort by statistics to add zscore info to particles.')

        form.addParam('rejectionMethod', EnumParam, choices=['None','MaxZscore', 'Percentage'],
                      default=REJECT_NONE, display=EnumParam.DISPLAY_COMBO, condition='doSort==True',
                      label='Automatic particle rejection',
                      help='How to automatically reject particles. It can be none (no rejection),'
                      ' maxZscore (reject a particle if its Zscore is larger than this value), '
                      'Percentage (reject a given percentage in each one of the screening criteria).',
                      expertLevel=LEVEL_ADVANCED)

        form.addParam('maxZscore', IntParam, default=3, expertLevel=LEVEL_ADVANCED,
                      condition='doSort==True and rejectionMethod==%d' % REJECT_MAXZSCORE,
                      label='Maximum Zscore',
                      help='Maximum Zscore above which particles are rejected.')
        
        form.addParam('percentage', IntParam, default=5, expertLevel=LEVEL_ADVANCED, 
                      condition='rejectionMethod==%d' % REJECT_PERCENTAGE,
                      label='Percentage (%)',
                      help='Percentage of particles to reject')
        
        form.addSection(label='Preprocess')
        form.addParam('doRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5, 
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.')
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        
        form.addParam('doFlip', BooleanParam, default=True,
                      label='Phase flipping (Recommended)', 
                      help='Use the information from the CTF to compensate for phase reversals.')
        
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], default=2, 
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ')
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is 0, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
        
        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        """for each micrograph insert the steps to preprocess it
        """       
        # Set sampling rate and inputMics according to downsample type
        self.inputCoords = self.inputCoordinates.get() 
        
        self.samplingInput = self.inputCoords.getMicrographs().getSamplingRate()
        
        if self.downsampleType.get() == SAME_AS_PICKING:
            # If 'same as picking' get samplingRate from input micrographs  
            self.inputMics = self.inputCoords.getMicrographs()
            self.samplingFinal = self.samplingInput
        else:
            self.inputMics = self.inputMicrographs.get()
            self.samplingOriginal = self.inputMics.getSamplingRate()
            if self.downsampleType.get() == ORIGINAL:
                # If 'original' get sampling rate from original micrographs
                self.samplingFinal = self.samplingOriginal
            else:
                # IF 'other' multiply the original sampling rate by the factor provided
                self.samplingFinal = self.samplingOriginal*self.downFactor.get()
                
        # Write pos files for each micrograph
        firstStepId = self._insertFunctionStep('writePosFilesStep')
                
        #if self.doFlip.get():
        ctfSet = self.ctfRelations.get()
           
        # For each micrograph insert the steps
        #run in parallel
        
        deps = []
        for mic in self.inputMics:
            localDeps = [firstStepId]
            micrographToExtract = mic.getFileName()
            micName = removeBaseExt(mic.getFileName())
            micId = mic.getObjId()
                                            
            #if self.doFlip.get():
            if self.ctfRelations.hasValue():
                mic.setCTF(ctfSet[micId])       

            # If downsample type is 'other' perform a downsample
            if self.downsampleType == OTHER:
                fnDownsampled = self._getTmpPath(micName+"_downsampled.xmp")
                downFactor = self.downFactor.get()
                args = "-i %(micrographToExtract)s -o %(fnDownsampled)s --step %(downFactor)f --method fourier"
                localDeps=[self._insertRunJobStep("xmipp_transform_downsample", args % locals(),prerequisites=localDeps)]
                micrographToExtract = fnDownsampled
            # If remove dust 
            if self.doRemoveDust:
                fnNoDust = self._getTmpPath(micName+"_noDust.xmp")
                
                thresholdDust = self.thresholdDust.get() #TODO: remove this extra variable
                args=" -i %(micrographToExtract)s -o %(fnNoDust)s --bad_pixels outliers %(thresholdDust)f"
                localDeps=[self._insertRunJobStep("xmipp_transform_filter", args % locals(),prerequisites=localDeps)]
                micrographToExtract = fnNoDust
                
                        
            #self._insertFunctionStep('getCTF', micId, micName, micrographToExtract)
            micName = removeBaseExt(mic.getFileName())
            #FIXME: Check only if mic has CTF when implemented ok
            #if self.doFlip or mic.hasCTF():
            fnCTF=None
            if self.ctfRelations.hasValue():
                # If the micrograph doesn't come from Xmipp, we need to write
                # a Xmipp ctfparam file to perform the phase flip on the micrograph                     
                fnCTF = micrographToCTFParam(mic, self._getTmpPath("%s.ctfParam" % micName))
                # Insert step to flip micrograph
                if self.doFlip:
                    localDeps = [self._insertFunctionStep('flipMicrographStep', 
                                                      micName, fnCTF, micrographToExtract,
                                                      prerequisites=localDeps)]
                    micrographToExtract = self._getTmpPath(micName +"_flipped.xmp")
            else:
                fnCTF = None        
            # Actually extract
            deps.append(self._insertFunctionStep('extractParticlesStep', micId, micName, 
                                              fnCTF, micrographToExtract, prerequisites=localDeps))
        # TODO: Delete temporary files
                        
        # Insert step to create output objects      
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
                
    #--------------------------- STEPS functions --------------------------------------------
    def writePosFilesStep(self):
        """ Write the pos file for each micrograph on metadata format. """
        #self.posFiles = writeSetOfCoordinates(self._getExtraPath(), self.inputCoords)
        writeSetOfCoordinates(self._getExtraPath(), self.inputCoords)
               
    def flipMicrographStep(self, micName, fnCTF, micrographToExtract):
        """ Flip micrograph. """           
        fnFlipped = self._getTmpPath(micName +"_flipped.xmp")

        args = " -i %(micrographToExtract)s --ctf %(fnCTF)s -o %(fnFlipped)s --downsampling %(downFactor)f"
        # xmipp_ctf_phase_flip expects the sampling rate of the micrographs, that has been used
        # to calculate the CTF, in the ctfparam file. If its no given, the program asumes that is equal to 1;
        # therefore, downsampling factor must be equal to sampling rate of the final mics.
        downFactor = self.samplingFinal
        self.runJob("xmipp_ctf_phase_flip", args % locals())
        
    def extractParticlesStep(self, micId, micName, fnCTF, micrographToExtract):
        """ Extract particles from one micrograph """
        #If flip selected and exists CTF model use the flip output
#        if self.doFlip and self.fnCTF:
#            micrographToExtract = self._getTmpPath(micName +"_flipped.xmp")
                
        outputRoot = str(self._getExtraPath(micName))
        #fnPosFile = self.getConvertedInput('inputCoords').getMicrographCoordFile(micId)
        fnPosFile =  self._getExtraPath(micName + ".pos")

        # If it has coordinates extract the particles      
        particlesMd = 'particles@%s' % fnPosFile
        
        boxSize = self.boxSize.get()
        
        #if fnPosFile is not None and xmipp.existsBlockInMetaDataFile(particlesMd):
        if exists(fnPosFile):
            args = "-i %(micrographToExtract)s --pos %(particlesMd)s -o %(outputRoot)s --Xdim %(boxSize)d" % locals()
            if self.downsampleType.get() != SAME_AS_PICKING:
                args += " --downsampling %f" % (self.samplingFinal/self.samplingInput)
            if self.doInvert:
                args += " --invert"
            if fnCTF:
                args += " --ctfparam " + fnCTF
            self.runJob("xmipp_micrograph_scissor", args)
            # Normalize 
            if self.doNormalize:
                self.runNormalize(outputRoot + '.stk', self.normType.get(), self.backRadius.get())          
                               
            if self.downsampleType.get() == OTHER:
                selfile = outputRoot + ".xmd"
                md = xmipp.MetaData(selfile)
                downsamplingFactor = self.samplingFinal/self.samplingInput
                md.operate("Xcoor=Xcoor*%f" % downsamplingFactor)
                md.operate("Ycoor=Ycoor*%f" % downsamplingFactor)
                md.write(selfile)
                
    def runNormalize(self, stack, normType, bgRadius):
        program = "xmipp_transform_normalize"
        args = "-i %(stack)s "
        
        if bgRadius <= 0:
            particleSize = xmipp.MetaDataInfo(stack)[0]
            bgRadius = int(particleSize/2)
        
        if normType=="OldXmipp":
            args += "--method OldXmipp"
        elif normType=="NewXmipp":
            args += "--method NewXmipp --background circle %(bgRadius)d"
        else:
            args += "--method Ramp --background circle %(bgRadius)d"
        self.runJob(program, args % locals())
        
    def createOutputStep(self):
        # Create the SetOfImages object on the database
        #imgSet = XmippSetOfParticles(self._getPath('images.xmd'))
                  
        #Create images.xmd metadata
        fnImages = self._getPath('images.xmd')
        imgsXmd = xmipp.MetaData() 
        posFiles = glob(self._getExtraPath('*.pos')) 
        for posFn in posFiles:
            xmdFn = self._getExtraPath(replaceBaseExt(posFn, "xmd"))
            md = xmipp.MetaData(xmdFn)
            mdPos = xmipp.MetaData('particles@%s' % posFn)
            mdPos.merge(md) 
            #imgSet.appendFromMd(mdPos)
            imgsXmd.unionAll(mdPos)
   
        #TODO: CHECK WITH JAVI
        #imgsXmd.sort(xmipp.MDL_IMAGE)
        imgsXmd.write(fnImages)

        # IF selected run xmipp_image_sort_by_statistics to add zscore info to images.xmd
        if self.doSort:
            args="-i %(fnImages)s --addToInput"
            if self.rejectionMethod == REJECT_MAXZSCORE:
                maxZscore = self.maxZscore.get()
                args += " --zcut " + str(maxZscore)
            elif self.rejectionMethod == REJECT_PERCENTAGE:
                percentage = self.percentage.get()
                args += " --percent " + str(percentage)

            self.runJob("xmipp_image_sort_by_statistics", args % locals())
        # Create output SetOfParticles
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputMics)
        if self.doFlip:
            imgSet.setIsPhaseFlipped(True)
        
        #imgSet.setHasCTF(self.fnCTF is not None)
        if self.downsampleType == OTHER:
            imgSet.setSamplingRate(self.inputMics.getSamplingRate()*self.downFactor.get())
        imgSet.setCoordinates(self.inputCoords)
        
        # Create a temporary set to read from the metadata file
        # and later create the good one with the coordinates 
        # properly set. We need this because the .update is not
        # working in the mapper when new attributes are added.
        imgSet.setHasCTF(self.ctfRelations.hasValue())
        auxSet = SetOfParticles(filename=':memory:')
        auxSet.copyInfo(imgSet)
        readSetOfParticles(fnImages, auxSet)
        # For each particle retrieve micId from SetOFCoordinates and set it on the CTFModel
        for img in auxSet:
            #FIXME: This can be slow to make a query to grab the coord, maybe use zip(imgSet, coordSet)???
            coord = self.inputCoords[img.getObjId()]
            ctfModel = img.getCTF()
            if ctfModel is not None:
                ctfModel.setObjId(coord.getMicId())
                ##img.setCTF(ctfModel)####JM
            img.setMicId(coord.getMicId())
            img.setCoordinate(coord)
            imgSet.append(img)
            
        self._storeMethodsInfo(fnImages)
        self._defineOutputs(outputParticles=imgSet)
        self._defineSourceRelation(self.inputCoords, imgSet)
        #TODO: pass CTF relation from input micrographs to imgSet
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # doFlip can only be True if CTF information is available on picked micrographs
        if self.doFlip and not self.ctfRelations.hasValue():
            validateMsgs.append('Phase flipping cannot be performed unless CTF information is provided.')
        return validateMsgs
    
    def _citations(self):
        return ['Vargas2013b']
        
    def _summary(self):
        downsampleTypeText = {
                              ORIGINAL:'Original micrographs',
                              SAME_AS_PICKING:'Same as picking',
                              OTHER: 'Other downsampling factor'}
        summary = []
        summary.append("Downsample type: %s" % downsampleTypeText.get(self.downsampleType.get()))
        if self.downsampleType == OTHER:
            summary.append("Downsampling factor: %.2f" % self.downFactor)
        summary.append("Particle box size: %d" % self.boxSize)
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            summary.append("Particles extracted: %d" % (self.outputParticles.getSize()))
            
        return summary
    
    def _methods(self):
        methodsMsgs = []

        if self.getStatus() == STATUS_FINISHED:
            msg = "A total of %d particles of size %d were extracted" % (self.getOutput().getSize(), self.boxSize)

            if self.downsampleType == ORIGINAL:
                msg += " from original micrographs."

            if self.downsampleType == OTHER:
                msg += " from original micrographs with downsampling factor of %.2f." % self.downFactor

            if self.downsampleType == SAME_AS_PICKING:
                msg += "."

            msg += self.methodsVar.get('')

            methodsMsgs.append(msg)

            if self.doInvert:
                methodsMsgs.append("Inverted contrast on images.")

            if self.doNormalize.get():
                methodsMsgs.append("Normalization performed of type %s." % (self.getEnumText('normType')))

            if self.doRemoveDust.get():
                methodsMsgs.append("Removed dust over a threshold of %s." % (self.thresholdDust))

        return methodsMsgs

    #--------------------------- UTILS functions --------------------------------------------
    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates or to the 
        Selected micrographs if Same as Picking not chosen. """
        if self.downsampleType == SAME_AS_PICKING:
            return self.inputCoordinates.get().getMicrographs()
        else:
            return self.inputMicrographs.get()
    
    def getImgIdFromCoord(self, coordId):
        """ Get the image id from the related coordinate id. """
        '%s:%06d'
        parts = coordId.split(':')
        imgFn = self._getExtraPath(replaceBaseExt(parts[0], "stk")) 
        
        return '%06d@%s' %(int(parts[1]), imgFn)
    
    def _storeMethodsInfo(self, fnImages):
        """ Store some information when the protocol finishes. """
        md = xmipp.MetaData(fnImages)
        total = md.size() 
        md.removeDisabled()
        zScoreMax = md.getValue(xmipp.MDL_ZSCORE, md.lastObject())
        numEnabled = md.size()
        numRejected = total - numEnabled

        msg = ""

        if self.doSort:
            if self.rejectionMethod != REJECT_NONE:
                msg = " %d of them were rejected with Zscore greater than %.2f." % (numRejected, zScoreMax)

        if self.doFlip:
            msg += "\nPhase flipping was performed."

        self.methodsVar.set(msg)

    def getCoords(self):
        if self.inputCoordinates.hasValue():
            return self.inputCoordinates.get()
        else:
            return None

    def getOutput(self):
        if (self.hasAttribute('outputParticles') and
            self.outputParticles.hasValue()):
            return self.outputParticles
        else:
            return None