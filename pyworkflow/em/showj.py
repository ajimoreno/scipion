# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
Configuration definition for launching the ShowJ visualization utility.
The name of ShowJ is due historical reasons. This data visualization tools
was first called xmipp_show. And later was re-written using Java and 
became xmipp_showj.
"""
import os
from os.path import join
from collections import OrderedDict

import pyworkflow as pw
from pyworkflow.utils import runJob


#------------------------ Showj constants ---------------------------
PATH = 'path'
DATASET = 'dataset'

TABLE_NAME = 'blockComboBox'
COLS_CONFIG = 'tableLayoutConfiguration'
COLS_CONFIG_DEFAULT = 'defaultColumnsLayoutProperties'
LABEL_SELECTED = 'labelsToRenderComboBox'

MODE = 'mode'
MODE_GALLERY = 'gallery'
MODE_TABLE = 'table'
MODE_VOL_ASTEX = 'volume_astex'
MODE_VOL_CHIMERA = 'volume_chimera'
RENDER = 'render'
ORDER = 'order'
ZOOM = 'zoom'

GOTO = 'goto'
ROWS = 'rows'
COLS = 'cols'
ALLOW_RENDER = 'allowRender'
MANUAL_ADJUST = 'colRowMode'

VOL_SELECTED = 'volumesToRenderComboBox'
VOL_TYPE = 'typeVolume'
VOL_VIEW = 'resliceComboBox'

IMG_DIMS = 'imageDimensions'
IMG_ZOOM = 'zoom'
IMG_ZOOM_DEFAULT = 'defaultZoom'
IMG_MIRRORY = 'mirrorY'
IMG_APPLY_TRANSFORM = 'applyTransformMatrix'
IMG_ONLY_SHIFTS = 'onlyShifts'
IMG_WRAP = 'wrap'
IMG_MAX_WIDTH = 'imageMaxWidth'
IMG_MIN_WIDTH = 'imageMinWidth'
IMG_MAX_HEIGHT = 'imageMaxHeight'
IMG_MIN_HEIGHT  = 'imageMinHeight'

PROJECT_NAME = 'projectName'
PROJECT_PATH = 'projectPath'
OBJECT_ID = 'objectId'


class ColumnsConfig():
    """ Store the configuration of the columsn for a given table in a dataset.
    The order of the columns will be stored and configuration for each columns.
    For each column, we store properties:
    - visible
    - allowSetVisible
    - renderable
    - allowSetRenderable
    - editable
    - allowSetEditable
    - renderFunc
    - renderFuncExtra
    """
    def __init__(self, ds, table, allowRender=True, defaultColumnsLayoutProperties=None):
        
        self._columnsDict = OrderedDict() 
         
        for col in table.iterColumns():
            
            self._columnsDict[col.getName()] = ColumnProperties(col, ds, allowRender, defaultColumnsLayoutProperties[col.getName()] if defaultColumnsLayoutProperties != None else {})
        
        
    def getRenderableColumns(self):
        """ Return a list with the name of renderable columns. """
        columns = [col.getLabel() for col in self._columnsDict.values() if col.isRenderable()]
        return columns
    
    def hasEnableColumn(self):
        for columnLayout in self._columnsDict.values():
            if "enable" == columnLayout.label:
                return True
        return False;
    
    def getColumnProperty(self, colName, propName):
        """ Get some property value of a given column. """
        col = self._columnsDict[colName]
        return getattr(col, propName)
    
    def configColumn(self, colName, **kwargs):
        """ Configure properties of a given column. """
        col = self._columnsDict[colName]
        for k, v in kwargs.iteritems():
            setattr(col, k, v)
            
    def printColumns(self):
        for col in self._columnsDict.values():
            print "column: ", col.getLabel()
            print "  values: ", col.getValues()
        
            
class ColumnProperties():
    """ Store some properties to customize how each column
    will be display in the table. 
    """
    def __init__(self, col, ds, allowRender, defaultColumnLayoutProperties):
        self._column = col        
        self.columnType = ds.getTypeOfColumn(col.getName())
        
        self.visible = not (self.columnType == 'id')
        self.allowSetVisible = True 
        
        self.editable = (self.columnType == 'text')
        self.allowSetEditable = self.editable
        
        self.renderable = False
        
        for k in defaultColumnLayoutProperties:
            if k == 'renderable':
                self.renderable = True
            else:
                self.renderable = False

#        self.renderable = False
            
        self.allowSetRenderable = (self.columnType == 'image' and allowRender)

        self.renderFunc = "get_image"
        self.extraRenderFunc = ""
        
    def getLabel(self):
        return self._column.getName()
    
    def getColumnType(self):
        return self.columnType
    
    def isRenderable(self):
        return self.renderable or self.allowSetRenderable
        
    def setValues(self, defaultColumnLayoutProperties):
        for key in defaultColumnLayoutProperties:
            setattr(self, key, defaultColumnLayoutProperties[key])
    
    def getValues(self):
        return {"visible":self.visible,
                "allowSetVisible":self.allowSetVisible,
                "editable":self.editable,
                "allowSetEditable":self.allowSetEditable,
                "renderable":self.renderable,
                "allowSetRenderable":self.allowSetRenderable,
                "renderFunc":self.renderFunc,
                "extraRenderFunc":self.extraRenderFunc,
                'columnType': self.columnType
                }
        
def getArchitecture():
    import platform
    arch = platform.architecture()[0]
    for a in ['32', '64']:
        if a in arch:
            return a
    return 'NO_ARCH' 
    
def getJavaIJappArguments(memory, appName, appArgs):
    """ Build the command line arguments to launch 
    a Java application based on ImageJ. 
    memory: the amount of memory passed to the JVM.
    appName: the qualified name of the java application.
    appArgs: the arguments specific to the application.
    """ 
    if len(memory) == 0:
        memory = "1g"
        print "No memory size provided. Using default: " + memory
    
    imagej_home = join(os.environ['XMIPP_HOME'], "external", "imagej")
    lib = join(os.environ['XMIPP_HOME'], "lib")
    javaLib = join(os.environ['XMIPP_HOME'], 'java', 'lib')
    plugins_dir = os.path.join(imagej_home, "plugins")
    arch = getArchitecture()
    args = "-Xmx%(memory)s -d%(arch)s -Djava.library.path=%(lib)s -Dplugins.dir=%(plugins_dir)s -cp %(imagej_home)s/*:%(javaLib)s/* %(appName)s %(appArgs)s" % locals()

    return args
    
def runJavaIJapp(memory, appName, args, batchMode=True, env=None):
    args = getJavaIJappArguments(memory, appName, args)
    runJob(None, "java", args, runInBackground=batchMode, env=env)
    
