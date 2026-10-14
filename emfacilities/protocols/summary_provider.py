# -*- coding: utf-8 -*-
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.object as pwobj
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.protocol import getUpdatedProtocol

from pwem.protocols import ProtImportImages


class SummaryProvider(TreeProvider):
    """
    Creates and organizes a hierarchical summary of protocol execution results, allowing users to inspect
    generated outputs and associated metadata in a structured tree view.

    AI Generated:

    Summary Provider (SummaryProvider) — User Manual
        Overview

        The SummaryProvider class is responsible for presenting the outputs generated during a protocol
        execution in an organized and navigable structure. Its primary objective is to help users inspect
        the results of previously executed protocols, understand the relationships between produced data
        objects, and review acquisition-related metadata when available.

        In practical workflows, this component acts as an intermediate representation layer between
        protocol execution and graphical visualization. It allows execution outputs to be displayed as a
        tree where parent protocols and their generated datasets can be explored interactively. This is
        especially useful in large processing pipelines where multiple protocols generate interconnected
        outputs that must be reviewed together.

        General Workflow

        The provider collects information from all protocols connected as inputs to the current workflow.
        Each protocol is analyzed to identify the datasets and outputs it produced. These outputs are then
        organized into a tree-like structure where protocols appear as parent entries and their generated
        objects are represented as child elements.

        The resulting hierarchy allows users to inspect both the identity of the producing protocol and
        the characteristics of the generated datasets. Typical displayed information includes output names,
        object sizes, and associations between related processing steps.

        Output Organization

        The generated summary is designed to avoid duplication and maintain a consistent representation of
        workflow outputs. When several protocols reference related objects, the provider ensures that the
        same item is not redundantly inserted multiple times into the visualization tree.

        Each object in the hierarchy contains descriptive information suitable for graphical interfaces or
        reporting tools. This organization helps users quickly understand the structure of a processing
        pipeline and identify the outputs produced at each stage.

        Acquisition Metadata

        When image import protocols are detected, the provider also extracts acquisition parameters
        associated with the imported data. These parameters commonly include microscope voltage,
        magnification, spherical aberration, pixel size, and dose per frame.

        This information is particularly valuable in cryo-EM workflows because acquisition conditions
        strongly influence downstream processing quality and interpretation. By exposing these parameters
        alongside the processing outputs, the provider allows users to verify that datasets were imported
        with the expected experimental settings.

        Biological and Practical Relevance

        In biological imaging workflows, maintaining clear visibility of generated datasets is essential
        for reproducibility and quality control. Processing pipelines often involve many intermediate
        objects, and the ability to inspect relationships between protocols and outputs helps users detect
        inconsistencies, missing data, or unintended workflow branches.

        The inclusion of acquisition metadata further improves traceability by preserving the connection
        between experimental conditions and computational results. This becomes especially important when
        comparing datasets acquired under different microscope settings or evaluating the origin of
        reconstruction variability.

        Final Perspective

        The SummaryProvider class serves as an organizational component that improves transparency and
        usability in complex processing workflows. By presenting protocol outputs and acquisition metadata
        in a structured hierarchy, it enables users to review, inspect, and validate workflow results in
        an efficient and biologically meaningful manner.
    """
    def __init__(self, protocol):
        TreeProvider.__init__(self)
        self.protocol = protocol
        self.getColumns = lambda: [('Name', 300), ('Output', 150),
                                   ('Number', 100)]
        self._parentDict = {}
        self.acquisition = []
        self.refreshObjects()

    def getObjects(self):
        return self._objects

    def refreshObjects(self):
        objects = []
        objIds = []  # need to store ids too to avoid duplication in runs table

        def addObj(objId, name, output='', size='', parent=None):
            if objId not in objIds:
                obj = pwobj.Object(objId=objId)
                obj.name = name
                obj.output = output
                obj.outSize = size
                obj._objParent = parent
                objIds.append(objId)
                objects.append(obj)
                return obj
            else:
                return None

        prots = [getUpdatedProtocol(p) for p in self.protocol.getInputProtocols()]

        for prot in prots:
            pobj = addObj(prot.getObjId(),
                          '%s (id=%s)' % (prot.getRunName(), prot.strId()))
            for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                outSet.load()
                outSet.loadAllProperties()
                # outSetId needs to be compound id to avoid duplicate ids
                outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
                addObj(outSetId, '', outName, outSet.getSize(), pobj)
                outSet.close()
                # Store acquisition parameters in case of the import protocol
                # NOTE by Yaiza: we force the string containing the Å to be unicode
                # because this is the encoding used when generating report in report_html.py
                if isinstance(prot, ProtImportImages):
                    self.acquisition = [("Microscope Voltage (kV): ",
                                         prot.voltage.get()),
                                        ("Spherical aberration (mm): ",
                                         prot.sphericalAberration.get()),
                                        ("Magnification: ",
                                         prot.magnification.get()),
                                        (u"Pixel Size (Å/px): ",
                                         round(outSet.getSamplingRate(), 2))
                                        ]
                    if prot.dosePerFrame.get() is not None:
                        self.acquisition.append((u"Dose per frame (e/Å²):",
                                                 prot.dosePerFrame.get()))

        self._objects = objects

    def getObjectInfo(self, obj):
        info = {'key': obj.strId(),
                'parent': obj._objParent,
                'text': obj.name,
                'values': (obj.output, obj.outSize),
                'open': True
                }

        return info
