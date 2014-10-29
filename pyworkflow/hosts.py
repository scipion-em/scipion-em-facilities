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
This modules contains classes to store information about
execution hosts.
"""

import os
from ConfigParser import ConfigParser

from pyworkflow.object import *
from pyworkflow.mapper import SqliteMapper, XmlMapper


class HostMapper(SqliteMapper):
    def __init__(self, filename, dictClasses=None):
        if dictClasses is None:
            dictClasses = globals()
        SqliteMapper.__init__(self, filename, dictClasses)
        
    def selectByLabel(self, objLabel):
        hostsList = self.selectAll()
        for host in hostsList:
            if host.label == objLabel:
                return host
        return None  
        
class HostConfig(OrderedObject):
    """ Main store the configuration for execution hosts. """
    
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.label = String()
        self.hostName = String()
        self.userName = String()
        self.password = String()
        self.hostPath = String()
        self.mpiCommand = String()
        self.queueSystem = QueueSystemConfig()
    
    def getLabel(self):
        return self.label.get()
    
    def getHostName(self):
        return self.hostName.get()
    
    def getUserName(self):
        return self.userName.get()
    
    def getPassword(self):
        return self.password.get()
    
    def getHostPath(self):
        return self.hostPath.get()
    
    def getSubmitCommand(self):
        return self.queueSystem.submitCommand.get()
    
    def getCancelCommand(self):
        return self.queueSystem.cancelCommand.get()
    
    def isQueueMandatory(self):
        return self.queueSystem.mandatory.get()
    
    def getSubmitTemplate(self):
        return self.queueSystem.getSubmitTemplate()
    
    def getMpiCommand(self):
        return self.mpiCommand.get()
    
    def getQueueSystem(self):
        return self.queueSystem
    
    def setLabel(self, label):
        self.label.set(label)
    
    def setHostName(self, hostName):
        self.hostName.set(hostName)
    
    def setUserName(self, userName):
        self.userName.set(userName)
    
    def setPassword(self, password):
        self.password.set(password)
    
    def setHostPath(self, hostPath):
        self.hostPath.set(hostPath)
    
    def setMpiCommand(self, mpiCommand):
        self.mpiCommand.set(mpiCommand)  
        
    def setQueueSystem(self, queueSystem):
        self.queueSystem = queueSystem    
  
    def addQueueSystem(self, queueConf=None):
        # Read from users' config file.
        cp = ConfigParser()
        cp.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
        SCIPION_CONFIG = queueConf or os.environ['SCIPION_CONFIG']
        # Also mentioned in /scipion . Maybe we could do better.
    
        try:
            assert cp.read(SCIPION_CONFIG) != [], 'Missing file %s' % SCIPION_CONFIG
    
            # Helper functions (to write less)
            def get(var): return cp.get('HOSTS', var).replace('%_(', '%(')
            def isOn(var): return var.lower() in ['true', 'yes', '1']
    
            self.mpiCommand.set(get('PARALLEL_COMMAND'))
            self.queueSystem = QueueSystemConfig()
            queueSys = self.queueSystem
            #queueSys = QueueSystemConfig()
            queueSys.name.set(get('NAME'))
            queueSys.mandatory.set(isOn(get('MANDATORY')))
            queueSys.submitCommand.set(get('SUBMIT_COMMAND'))
            queueSys.submitTemplate.set(get('SUBMIT_TEMPLATE'))
            queueSys.cancelCommand.set(get('CANCEL_COMMAND'))
            queueSys.checkCommand.set(get('CHECK_COMMAND'))
    
            queue = QueueConfig()
            queue.maxCores.set(get('MAX_CORES'))
            queue.allowMPI.set(isOn(get('ALLOW_MPI')))
            queue.allowThreads.set(isOn(get('ALLOW_THREADS')))
        except Exception as e:
            sys.exit('Failed to read settings. The reported error was:\n  %s\n'
                     'To solve it, delete %s and run again.' % (e, SCIPION_CONFIG))
    
        queueSys.queues = List()
        queueSys.queues.append(queue)      
        
class QueueSystemConfig(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args) 
        self.name = String()
        self.mandatory = Boolean()
        self.queues = None # List for queue configurations
        self.submitCommand = String()
        self.checkCommand = String()
        self.cancelCommand = String()
        self.submitTemplate = String()
        
    def hasValue(self):
        return self.name.hasValue() and len(self.queues)
    
    def getName(self):
        return self.name.get()
    
    def getMandatory(self):
        return self.mandatory.get()
    
    def getSubmitTemplate(self):
        return self.submitTemplate.get()
    
    def getSubmitCommand(self):
        return self.submitCommand.get()
    
    def getCheckCommand(self):
        return self.checkCommand.get()
    
    def getCancelCommand(self):
        return self.cancelCommand.get()
    
    def getQueues(self):
        return self.queues
    
    def setName(self, name):
        self.name.set(name)
    
    def setMandatory(self, mandatory):
        self.mandatory.set(mandatory)
    
    def setSubmitTemplate(self, submitTemplate):
        self.submitTemplate.set(submitTemplate)
    
    def setSubmitCommand(self, submitCommand):
        self.submitCommand.set(submitCommand)
    
    def setCheckCommand(self, checkCommand):
        self.checkCommand.set(checkCommand)
    
    def setCancelCommand(self, cancelCommand):
        self.cancelCommand.set(cancelCommand)
    
    def setQueues(self, queues):
        self.queues = queues
        
    def getQueueConfig(self, objId):
        if objId is not None and self.queues is not None:
            for queueConfig in self.queues:
                if objId == queueConfig.getObjId():
                    return queueConfig
        return None
        
        
class QueueConfig(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args) 
        self.name = String('default')
        self.maxCores = Integer()
        self.allowMPI = Boolean()
        self.allowThreads = Boolean()
        self.maxHours = Integer()
        
    def getName(self):
        return self.name.get()
    
    def getMaxCores(self):
        return self.maxCores.get()
    
    def getAllowMPI(self):
        return self.allowMPI.get()
    
    def getAllowThreads(self):
        return self.allowThreads.get()
    
    def getMaxHours(self):
        return self.maxHours.get()
    
    def setName(self, name):
        self.name.set(name)
    
    def setMaxCores(self, maxCores):
        self.maxCores.set(maxCores)
        
    def setAllowMPI(self, allowMPI):
        self.allowMPI.set(allowMPI)
    
    def setAllowThreads(self, allowThreads):
        self.allowThreads.set(allowThreads)
    
    def setMaxHours(self, maxHours):
        self.maxHours.set(maxHours)