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
This module contains some sqlite basic tools to handle Databases.
"""

from sqlite3 import dbapi2 as sqlite

from pyworkflow.utils import envVarOn


class SqliteDb():
    """Class to handle a Sqlite database.
    It will create connection, execute queries and commands.
    """
    OPEN_CONNECTIONS = {} # Store all conections made
    
    def __init__(self):
        self._reuseConnections = False
        
    def _createConnection(self, dbName, timeout):
        """Establish db connection"""
        self._dbName = dbName
        if self._reuseConnections and dbName in self.OPEN_CONNECTIONS:
            self.connection = self.OPEN_CONNECTIONS[dbName]
        else:
            self.connection = sqlite.Connection(dbName, timeout, check_same_thread=False)
            self.connection.row_factory = sqlite.Row
            self.OPEN_CONNECTIONS[dbName] = self.connection
            
        self.cursor = self.connection.cursor()
        # Define some shortcuts functions
        if envVarOn('SCIPION_DEBUG_SQLITE'):
            self.executeCommand = self._debugExecute
        else:
            self.executeCommand = self.cursor.execute
        self.commit = self.connection.commit
        
    def getDbName(self):
        return self._dbName
    
    def close(self):
        self.connection.close()
        if self._dbName in self.OPEN_CONNECTIONS:
            del self.OPEN_CONNECTIONS[self._dbName]
        
    def _debugExecute(self, *args):
        try:
            return self.cursor.execute(*args)
        except Exception, ex:
            print ">>>> FAILED cursor.execute on db: '%s'" % self._dbName
            print "COMMAND: ", args[0]
            print "ARGUMENTS: ", args[1:]
            raise ex
            
        
        #return self.cursor.fetchone()
    
    def _iterResults(self):
        row = self.cursor.fetchone()
        while row is not None:
            yield row
            row = self.cursor.fetchone()
        
    def _results(self, iterate=False):
        """ Return the results to which cursor, point to. 
        If iterates=True, iterate yielding each result independenly"""
        if not iterate:
            return self.cursor.fetchall()
        else:
            return self._iterResults()
        
    def getTables(self, tablePattern=None):
        """ Return the table names existing in the Database.
        If  tablePattern is not None, only tables matching 
        the pattern will be returned.
        """
        self.executeCommand("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';")
        return [str(row['name']) for row in self._iterResults()]
        