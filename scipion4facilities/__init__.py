# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
# *
# *    Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
"""
This modules contains classes useful for cryo-EM facilities
"""

import os

import pwem

from .constants import *

_logo = "scipion_icon.gif"
_references = ["delaRosaTrevin201693"]



class Plugin(pwem.Plugin):
    _url = URL

    @classmethod
    def _defineVariables(cls):
        pass

    @classmethod
    def _defineEmVar(cls, varName, defaultValue):
        """ Shortcut method to define variables prepending EM_ROOT if variable is not absolute"""
        pass

    @classmethod
    def getHTMLTemplatePath(self):
        """ Returns the path of the customized template at
        config/execution.summary.html or the standard scipion HTML template"""
        # Try if there is a customized template
        template = os.path.join(pwem.Config.SCIPION_CONFIG, 'execution.summary.html')

        if not os.path.exists(template):
            pluginDir = os.path.dirname(os.path.realpath(__file__))
            template = os.path.join(pluginDir, 'execution.summary.template.html')
        else:
            print("Customized HTML template found at %s." % template)
        return template

    @classmethod
    def defineBinaries(cls, env):
        pass
