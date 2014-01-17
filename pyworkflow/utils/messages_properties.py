# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
This module define the text used in the application.
"""

class Message():
    # Example Usage: 
        #       MyMessage = Message()
        #       print MyMessage.label

    # Header List
    VIEW_PROJECTS = 'Projects'
    VIEW_PROTOCOLS = 'Protocols'
    VIEW_DATA = 'Data'
    VIEW_HOSTS = 'Hosts'
    
    # Projects Template
    LABEL_PROJECTS='Projects'
    LABEL_CREATE_PROJECT='Create Project'
    TITLE_CREATE_PROJECT='Enter the project name'
    TITLE_CREATE_PROJECT_NAME='Project Name: '
    MESSAGE_CREATE_PROJECT='Are you sure to *DELETE* the project and all its *DATA*?'
    LABEL_DELETE_PROJECT='Delete Project'
    TITLE_DELETE_PROJECT='Confirm project deletion'
    LABEL_MODIFIED='Modified: '
    
    # Project Content Template
    LABEL_PROJECT='Project '
    #-- Protocol Treeview --
    LABEL_WORKFLOW='Workflow View: '
    LABEL_PROTTREE_NONE='No athletes.'
    #-- Toolbar --
    LABEL_NEW='New'
    LABEL_NEW_ACTION='New     '
    LABEL_EDIT='Edit'
    LABEL_EDIT_ACTION='Edit     '
    LABEL_COPY='Copy'
    LABEL_COPY_ACTION='Copy   '
    LABEL_DELETE='Delete'
    LABEL_DELETE_ACTION='Delete    '
    LABEL_BROWSE='Browse'
    LABEL_BROWSE_ACTION='Browse '
    LABEL_STOP='Stop'
    LABEL_STOP_ACTION='Stop execution'
    LABEL_ANALYZE='Analyze Results'
    LABEL_TREE='Tree'
    LABEL_LIST='List'
    LABEL_REFRESH='Refresh'
    LABEL_DEFAULT = 'Default'
    LABEL_CONTINUE = 'Continue'
    LABEL_CONTINUE_ACTION = 'Approve continue'
    #-- Tabs --
    LABEL_DATA='Data'
    LABEL_SUMMARY='Summary'
    LABEL_INPUT='Input'
    LABEL_OUTPUT='Output'
    NO_INFO_SUMMARY='No summary information.'
    NO_SAVE_SETTINGS='Error try to save settings. '
    
    # Protocol Form Template
    TITLE_NAME_RUN=' Protocol Run: '
    TITLE_RUN='Run'
    TITLE_RUN_NAME='Run name'
    LABEL_COMMENT='Describe your run here...'
    TITLE_RUN_MODE='Run mode'
    LABEL_RUN_MODE_RESUME='resume'
    LABEL_RUN_MODE_RESTART='restart'
    TITLE_EXPERT='Expert Level'
    LABEL_EXPERT_NORMAL='Normal'
    LABEL_EXPERT_ADVANCE='Advanced'
    LABEL_EXPERT_EXPERT='Expert'
    TITLE_EXEC='Execution'
    TITLE_EXEC_HOST='Execution host'
    TITLE_THREADS='Threads'
    TITLE_MPI='MPI'
    TITLE_QUEUE='Launch to queue?'
    TITLE_BROWSE_DATA='Protocol data'
    LABEL_QUEUE_YES='Yes'
    LABEL_QUEUE_NO='No'
    LABEL_PARAM_YES='Yes'
    LABEL_PARAM_NO='No'
    LABEL_BUTTON_CLOSE='Close'
    LABEL_BUTTON_SAVE='Save'
    LABEL_BUTTON_EXEC='Execute'
    LABEL_BUTTON_VIS='Visualize'
    LABEL_BUTTON_RETURN='Save'
    
    TITLE_LAUNCHED='Success'
    LABEL_LAUNCHED='The protocol was launched successfuly.'
    LABEL_FOUND_ERROR='Errors found'
    TITLE_SAVED_FORM='Success'
    LABEL_SAVED_FORM='The protocol was saved successfuly.'
    TITLE_DELETE_FORM='Confirm DELETE'
    LABEL_DELETE_FORM='*ALL DATA* related to this _protocol run_ will be *DELETED*.\nDo you really want to continue?'    
    TITLE_STOP_FORM='Confirm STOP'
    LABEL_STOP_FORM='Do you really want to *STOP* this run?'
    
    NO_VIEWER_FOUND='There is not viewer for protocol:'
    
    # Hosts Template
    LABEL_HOSTS='Hosts'
    LABEL_HOSTS_ACTION=' Hosts '
    TITLE_HOST_DELETE_FORM='Confirm DELETE'
    LABEL_HOST_DELETE_FORM='*ALL CONFIGURATION* related to this _host_ will be *DELETED*.\nDo you really want to continue?'    
    BUTTON_TEST='Test configuration'
    TITLE_HOST_SAVE_FORM='Action SAVE'
    LABEL_HOST_SAVE_FORM_SUCESS='Host sucessfully saved. \n'
    LABEL_HOST_SAVE_FORM_FAIL='Error saving host. \n'
    LABEL_HOST_SAVE_FORM='*ALL CONFIGURATION* related to this _host_ will be *DELETED*.\nDo you really want to continue?'    
    LABEL_HOST_CONECTION_SUCESS='Connection SUCCEEDED to host'
    LABEL_HOST_CONECTION_FAIL='Connection FAILED to host'
    LABEL_HOST_CONECTION_TEST='Connection TEST'

    
    
class Icon():
    # Project Content Template
    RUNS_TREE = 'fa-sitemap.png'
    RUNS_LIST = 'fa-bars.png'
    ACTION_NEW= 'new_object.gif',
    ACTION_EDIT= 'fa-pencil.png'
    ACTION_COPY=  'fa-files-o.png'
    ACTION_DELETE=  'fa-trash-o.png'
    ACTION_REFRESH=  'fa-refresh.png'
    ACTION_STEPS=  'fa-folder-open.png'
    ACTION_TREE=  None
    ACTION_LIST=  'fa-bars.png'
    ACTION_STOP= 'fa-stop.png'
    ACTION_CONTINUE= 'fa-play-circle-o.png'
    ACTION_RESULTS= 'fa-eye.png'
    
    
    #Host template
    BUTTON_CLOSE='dialog_close.png'
    BUTTON_SAVE='filesave.png'
    
