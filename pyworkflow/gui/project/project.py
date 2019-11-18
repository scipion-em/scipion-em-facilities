#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software: you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation, either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program.  If not, see <https://www.gnu.org/licenses/>.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
Main Project window implementation.
It is composed by three panels:
1. Left: protocol tree.
2. Right upper: VIEWS (Data/Protocols)
3. Summary/Details
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import threading
import shlex
import subprocess
import socketserver
import tempfile

import pyworkflow as pw
import pyworkflow.utils as pwutils
from pyworkflow.gui.project.utils import OS
from pyworkflow.project import MenuConfig, ProjectSettings
from pyworkflow.gui import Message, Icon
from pyworkflow.gui.browser import FileBrowserWindow
# Usage commented.
# from pyworkflow.em.viewers import EmPlotter
# Moved to Scipion-app
#from pyworkflow.gui.plugin_manager import PluginManager
from pyworkflow.gui.plotter import Plotter
from pyworkflow.gui.text import _open_cmd, openTextFileEditor
from pyworkflow.webservices import ProjectWorkflowNotifier, WorkflowRepository

from .labels import LabelsDialog
# Import possible Object commands to be handled
from .base import ProjectBaseWindow, VIEW_PROTOCOLS, VIEW_PROJECTS


class ProjectWindow(ProjectBaseWindow):
    """ Main window for working in a Project. """
    _OBJECT_COMMANDS = {}

    def __init__(self, path, master=None):
        # Load global configuration
        self.projName = os.path.basename(path)
        try:
            projTitle = '%s (%s on %s)' % (self.projName,
                                           pwutils.getLocalUserName(),
                                           pwutils.getLocalHostName())
        except Exception:
            projTitle = self.projName 

        self.projPath = path
        self.project = self.loadProject()

        # TODO: put the menu part more nicely. From here:
        menu = MenuConfig()

        projMenu = menu.addSubMenu('Project')
        projMenu.addSubMenu('Browse files', 'browse',
                            icon='fa-folder-open.gif')
        projMenu.addSubMenu('Remove temporary files', 'delete',
                            icon='fa-trash-o.gif')
        projMenu.addSubMenu('Manage project labels', 'labels',
                            icon=Icon.TAGS)
        projMenu.addSubMenu('Toggle color mode', 'color_mode',
                            shortCut="Ctrl+t", icon=Icon.ACTION_VISUALIZE)
        projMenu.addSubMenu('Select all protocols', 'select all',
                            shortCut="Ctrl+a")
        projMenu.addSubMenu('Find protocol to add', 'find protocol',
                            shortCut="Ctrl+f")
        projMenu.addSubMenu('', '')  # add separator
        projMenu.addSubMenu('Import workflow', 'load_workflow',
                            icon='fa-download.gif')
        projMenu.addSubMenu('Search workflow', 'search_workflow',
                            icon = 'fa-search.gif')
        projMenu.addSubMenu('Export tree graph', 'export_tree')
        projMenu.addSubMenu('', '')  # add separator
        projMenu.addSubMenu('Notes', 'notes', icon='fa-pencil.gif')
        projMenu.addSubMenu('', '')  # add separator
        projMenu.addSubMenu('Exit', 'exit', icon='fa-sign-out.gif')

        helpMenu = menu.addSubMenu('Help')
        helpMenu.addSubMenu('Online help', 'online_help',
                            icon='fa-external-link.gif')
        helpMenu.addSubMenu('About', 'about',
                            icon='fa-question-circle.gif')
        helpMenu.addSubMenu('Contact support', 'contact_us',
                            icon='fa-question-circle.gif')

        self.menuCfg = menu
        # TODO: up to here

        if self.project.openedAsReadOnly():
            self.projName += "<READ ONLY>"

        # Notify about the workflow in this project
        self.selectedProtocol = None
        self.showGraph = False
        Plotter.setBackend('TkAgg')
        ProjectBaseWindow.__init__(self, projTitle, master,
                                   minsize=(90,50), icon=Icon.SCIPION_ICON_PROJ)

        OS.handler().maximizeWindow(self.root)

        self.switchView(VIEW_PROTOCOLS)

        self.initProjectTCPServer()  # Socket thread to communicate with clients

        ProjectWorkflowNotifier(self.project).notifyWorkflow()

    def createHeaderFrame(self, parent):
        """Create the header and add the view selection frame at the right."""
        header = ProjectBaseWindow.createHeaderFrame(self, parent)
        self.addViewList(header)
        return header

    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        self.settings.write()
        
    def _onClosing(self):
        try:
            if not self.project.openedAsReadOnly():
                self.saveSettings()
        except Exception as ex:
            print("%s %s" % (Message.NO_SAVE_SETTINGS, str(ex)))
        ProjectBaseWindow._onClosing(self)
     
    def loadProject(self):
        proj = pw.project.Project(pw.Config.getDomain(), self.projPath)
        proj.load()

        # Check if we have settings.sqlite, generate if not
        settingsPath = os.path.join(proj.path, proj.settingsPath)
        if os.path.exists(settingsPath):
            self.settings = proj.getSettings()
        else:
            print('Warning: settings.sqlite not found! '
                  'Creating default settings..')
            self.settings = proj.createSettings()

        self.generalCfg = self.settings.getConfig()
        self.protCfg = proj.getCurrentProtocolView()

        return proj

    # The next functions are callbacks from the menu options.
    # See how it is done in pyworkflow/gui/gui.py:Window._addMenuChilds()
    #
    def onBrowseFiles(self):
        # Project -> Browse files
        FileBrowserWindow("Browse Project files",
                          self, self.project.getPath(''), 
                          selectButton=None  # we will select nothing
                          ).show()

    def onNotes(self):
        if not all(var in os.environ for var in ['SCIPION_NOTES_PROGRAM',
                                                 'SCIPION_NOTES_FILE',
                                                 'SCIPION_NOTES_ARGS']):
            return self.showError("Missing variables SCIPION_NOTES_* under\n"
                                  "[VARIABLES] section in the configuration file\n"
                                  "~/.config/scipion/scipion.conf")
        args = []
        # Program name
        program = os.environ.get('SCIPION_NOTES_PROGRAM', None)
        notesFile = self.project.getPath('Logs', os.environ['SCIPION_NOTES_FILE'])

        if program:
            args.append(program)
            # Custom arguments
            if os.environ.get('SCIPION_NOTES_ARGS', None):
                args.append(os.environ['SCIPION_NOTES_ARGS'])
            args.append(notesFile)
            subprocess.Popen(args)  #nonblocking
        else:
            # if no program has been selected
            # xdg-open will try to guess but
            # if the file does not exist it
            # will return an error so If the file does
            # not exist I will create an empty one
            # 'a' will avoid accidental truncation
            if not os.path.exists(notesFile):
                open(notesFile, 'a').close()
            openTextFileEditor(notesFile)

    def onRemoveTemporaryFiles(self):
        # Project -> Remove temporary files
        tmpPath = os.path.join(self.project.path, self.project.tmpPath)
        n = 0
        try:
            for fname in os.listdir(tmpPath):
                fpath = "%s/%s" % (tmpPath, fname)
                if os.path.isfile(fpath):
                    os.remove(fpath)
                    n += 1
                # TODO: think what to do with directories. Delete? Report?
            self.showInfo("Deleted content of %s -- %d file(s)." % (tmpPath, n))
        except Exception as e:
            self.showError(str(e))
        
    def _loadWorkflow(self, obj):
        try:
            self.project.loadProtocols(obj.getPath())
            self.getViewWidget().updateRunsGraph(True, reorganize=True)
        except Exception as ex:
            self.showError(str(ex))
            
    def onImportWorkflow(self):
        FileBrowserWindow("Select workflow .json file",
                          self, self.project.getPath(''),
                          onSelect=self._loadWorkflow,
                          selectButton='Import'
                          ).show()

    def onSearchWorkflow(self):
        WorkflowRepository().search()

    def onExportTreeGraph(self):
        runsGraph = self.project.getRunsGraph(refresh=True)
        useId = not pwutils.envVarOn('SCIPION_TREE_NAME')
        dotStr = runsGraph.printDot(useId=useId)
        with tempfile.NamedTemporaryFile(suffix='.gv') as dotFile:
            dotFile.write(dotStr)
            dotFile.flush()
            openTextFileEditor(dotFile.name)

        if useId:
            print("\nexport SCIPION_TREE_NAME=1 # to use names instead of ids")
        else:
            print("\nexport SCIPION_TREE_NAME=0 # to use ids instead of names")

    def onManageProjectLabels(self):
        self.manageLabels()

    def onToggleColorMode(self):
        self.getViewWidget()._toggleColorScheme(None)

    def onSelectAllProtocols(self):
        self.getViewWidget()._selectAllProtocols(None)

    def onFindProtocolToAdd(self):
        self.getViewWidget()._findProtocol(None)

    def manageLabels(self):
        return LabelsDialog(self.root,
                            self.project.settings.getLabels(),
                            allowSelect=True)

    def initProjectTCPServer(self):
        server = ProjectTCPServer((self.project.address, self.project.port),
                                  ProjectTCPRequestHandler)
        server.project = self.project
        server.window = self
        server_thread = threading.Thread(name="projectTCPserver", target=server.serve_forever)
        # Exit the server thread when the main thread terminates
        server_thread.daemon = True
        server_thread.start()

    # Seems it is not used and should be in scipion-em
    # def schedulePlot(self, path, *args):
    #     self.enqueue(lambda: EmPlotter.createFromFile(path, *args).show())

    @classmethod
    def registerObjectCommand(cls, cmd, func):
        """ Register an object command to be handled when receiving the
        action from showj. """
        cls._OBJECT_COMMANDS[cmd] = func

    def runObjectCommand(self, cmd, inputStrId, objStrId):
        try:
            objId = int(objStrId)
            project = self.project

            if os.path.isfile(inputStrId) and os.path.exists(inputStrId):
                from pwem.utils import loadSetFromDb
                inputObj = loadSetFromDb(inputStrId)
            else:
                inputId = int(inputStrId)
                inputObj = project.mapper.selectById(inputId)

            func = self._OBJECT_COMMANDS.get(cmd, None)

            if func is None:
                print("Error, command '%s' not found. " % cmd)
            else:
                def myfunc():
                    func(inputObj, objId)
                    inputObj.close()
                self.enqueue(myfunc)

        except Exception as ex:
            print("There was an error executing object command !!!:")
            print(ex)
    
    def recalculateCTF(self, inputObjId, sqliteFile):
        """ Load the project and launch the protocol to
        create the subset.
        """
        # Retrieve project, input protocol and object from db
        project = self.project
        inputObj = project.mapper.selectById(int(inputObjId))
        parentProtId = inputObj.getObjParentId()
        parentProt = project.mapper.selectById(parentProtId)
        protDep = project._getProtocolsDependencies([parentProt])
        if protDep:
            prot = project.copyProtocol(parentProt)
            prot.continueRun.set(parentProt)
        else:
            prot = parentProt
            prot.isFirstTime.set(True)
        
        # Define the input params of the new protocol
        prot.recalculate.set(True)
        prot.sqliteFile.set(sqliteFile)
        # Launch the protocol
        self.getViewWidget().executeProtocol(prot)


class ProjectManagerWindow(ProjectBaseWindow):
    """ Windows to manage all projects. """
    def __init__(self, **kwargs):
        # Load global configuration
        settings = ProjectSettings()

        # TODO: put the menu part more nicely. From here:
        menu = MenuConfig()

        fileMenu = menu.addSubMenu('File')
        fileMenu.addSubMenu('Browse files', 'browse', icon='fa-folder-open.gif')
        fileMenu.addSubMenu('Exit', 'exit', icon='fa-sign-out.gif')

        confMenu = menu.addSubMenu('Configuration')
        confMenu.addSubMenu('General', 'general')
        confMenu.addSubMenu('Hosts', 'hosts')
        confMenu.addSubMenu('Protocols', 'protocols')
        # Moved to scipion app: confMenu.addSubMenu('Plugins', 'plugins')
        confMenu.addSubMenu('User', 'user')

        helpMenu = menu.addSubMenu('Help')
        helpMenu.addSubMenu('Online help', 'online_help', icon='fa-external-link.gif')
        helpMenu.addSubMenu('About', 'about', icon='fa-question-circle.gif')

        self.menuCfg = menu
        self.generalCfg = settings.getConfig()

        try:
            title = '%s (%s on %s)' % (Message.LABEL_PROJECTS, 
                                       pwutils.getLocalUserName(),
                                       pwutils.getLocalHostName())
        except Exception:
            title = Message.LABEL_PROJECTS
        
        ProjectBaseWindow.__init__(self, title, minsize=(750, 500),
                                   icon=Icon.SCIPION_ICON_PROJS, **kwargs)
        self.manager = pw.project.Manager()
        self.switchView(VIEW_PROJECTS)

    #
    # The next functions are callbacks from the menu options.
    # See how it is done in pyworkflow/gui/gui.py:Window._addMenuChilds()
    #
    def onBrowseFiles(self):
        # File -> Browse files
        FileBrowserWindow("Browse files", self,
                          pw.Config.SCIPION_USER_DATA,
                          selectButton=None).show()

    def onGeneral(self):
        # Config -> General
        self._openConfigFile(pw.Config.SCIPION_CONFIG)

    @staticmethod
    def _openConfigFile(configFile):
        """ Open an Scipion configuration file, if the user have one defined,
        also open that one with the defined text editor.
        """
        _open_cmd(pw.getConfigPath(configFile))


    @staticmethod
    def onHosts():
        # Config -> Hosts
        ProjectManagerWindow._openConfigFile(pw.Config.SCIPION_HOSTS)

    @staticmethod
    def onProtocols():
        ProjectManagerWindow._openConfigFile(pw.Config.SCIPION_PROTOCOLS)

    @staticmethod
    def onUser():
        ProjectManagerWindow._openConfigFile(pw.Config.SCIPION_LOCAL_CONFIG)

    # Moved to scipion-app
    # def onPlugins(self):
    #     # Config -> Plugins
    #     PluginManager("Plugin Manager", self, pw.Config.SCIPION_USER_DATA,
    #                   selectButton=None).show()


class ProjectTCPServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    pass


class ProjectTCPRequestHandler(socketserver.BaseRequestHandler):

    def handle(self):
        try:
            project = self.server.project
            window = self.server.window
            msg = self.request.recv(1024)
            tokens = shlex.split(msg)
            if msg.startswith('run protocol'):
                protocolName = tokens[2]
                protocolClass = pw.Config.getDomain().getProtocols()[protocolName]
                # Create the new protocol instance and set the input values
                protocol = project.newProtocol(protocolClass)

                for token in tokens[3:]:
                    #print token
                    param, value = token.split('=')
                    attr = getattr(protocol, param, None)
                    if param == 'label':
                        protocol.setObjLabel(value)
                    elif attr.isPointer():
                        obj = project.getObject(int(value))
                        attr.set(obj)
                    elif value:
                        attr.set(value)
                #project.launchProtocol(protocol)
                # We need to enqueue the action of execute a new protocol
                # to be run in the same GUI thread and avoid concurrent
                # access to the project sqlite database
                window.getViewWidget().executeProtocol(protocol)
            if msg.startswith('run function'):
                functionName = tokens[2]
                functionPointer = getattr(window, functionName)
                functionPointer(*tokens[3:])
            else:
                answer = 'no answer available'
                self.request.sendall(answer + '\n')
        except Exception as e:
            print(e)
            import traceback
            traceback.print_stack()

