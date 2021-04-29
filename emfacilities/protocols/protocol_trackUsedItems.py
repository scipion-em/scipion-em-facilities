# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
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

from pwem.protocols.protocol import EMProtocol
from pwem.objects import Micrograph
from pyworkflow.protocol import params
import copy, os
import networkx as nx

class UsedItemsTracker(EMProtocol):
  """
  This protocol will track the items (micrographs, classes2D,...) that has been used in a scipion protocol to
  generate a final volume
  """
  _label = 'Track used items'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    form.addParam('inputVolumes', params.PointerParam,
                  pointerClass='Volume', allowsNull=False,
                  label="Input volumes",
                  help='Select the sets of volumes to track their used items')

    form.addSection(label='Tracking')
    form.addParam('trackParticles', params.BooleanParam, default=True,
                  label='Track used particles?',
                  help='The particles used in the final volume reconstruction will be tracked')
    form.addParam('inputParticles', params.PointerParam,
                  pointerClass='SetOfParticles', expertLevel=params.LEVEL_ADVANCED,
                  label="Used Particles", allowsNull=True,
                  help='Set of Particles containing the used particles. If None, the set of particles'
                       'closest to the volume will be used. It will also be used to find the used mics, ctfs'
                       'and coordinates.')
    
    form.addParam('trackMics', params.BooleanParam, default=True,
                  label='Track used micrographs?',
                  help='The micrographs with picked particles used in the final volume reconstruction '
                       'will be tracked and scored according to the number of used particles in them')
    form.addParam('inputMicrographs', params.PointerParam,
                  pointerClass='SetOfMicrographs', expertLevel=params.LEVEL_ADVANCED,
                  label="Input micrographs", condition='trackMics', allowsNull=True,
                  help='Set of micrographs where to track the used items. If None, the first protocol'
                       'with output micrographs will be used. Each of these micrographs will be scored'
                       'depending on the number of used particles in them')

    form.addParam('trackCTFs', params.BooleanParam, default=True,
                  label='Track used CTFs?', condition='trackMics',
                  help='The CTFs from micrographs with picked particles used in the final volume reconstruction '
                       'will be tracked and scored according to the number of used particles in them')
    form.addParam('inputCTFs', params.PointerParam,
                  pointerClass='SetOfCTFs', expertLevel=params.LEVEL_ADVANCED,
                  label="Input CTFs", condition='trackCTFs and trackMics', allowsNull=True,
                  help='Set of CTFs where to track the used items. If None, the first protocol'
                       'with output CTFs will be used')

    form.addParam('trackCoordinates', params.BooleanParam, default=True,
                  label='Track used coordinates?', condition='trackMics',
                  help='The coordinates used in the final volume reconstruction will be tracked')
    form.addParam('inputCoordinates', params.PointerParam,
                  pointerClass='SetOfCoordinates', expertLevel=params.LEVEL_ADVANCED,
                  label="Input Coordinates", condition='trackCoordinates and trackMics', allowsNull=True,
                  help='Set of Coordinates where to track the used items. If None, the first protocol'
                       'with output Coordinates will be used')

    '''form.addParam('trackClasses2D', params.BooleanParam, default=True,
                  label='Track used Classes2D?',
                  help='The Classes2D used in the final volume reconstruction will be tracked')
    form.addParam('inputClasses2D', params.PointerParam,
                  pointerClass='SetOfClasses2D', expertLevel=params.LEVEL_ADVANCED,
                  label="Input Classes2D", condition='trackClasses2D', allowsNull=True,
                  help='Set of Classes2D where to track the used items. If None, the first protocol'
                       'with output Classes2D will be used')'''


  # --------------------------- INSERT steps functions ----------------------
  def _insertAllSteps(self):
      self._insertFunctionStep('getTracksStep')
      self._insertFunctionStep('getUsedItemsStep')

  def getTracksStep(self):
    self.project = self.getProject()

    self.inpDic, self.outDic = self.buildInputsDic(self.project), self.buildOutputsDic(self.project)
    self.protDic = self.project.getProtocolsDict()
    self.outGraph = self.generateOutputsGraphRec(self.getObjId())

  def getUsedItemsStep(self):
    particleSets = self.getInputParticles()
    particlesSet = self.joinSetsOfParticles(particleSets)
    partSamplingRate = particlesSet.getSamplingRate()
    firstParticlesCodes = self.getFurthestCodesFromNode(self.outGraph, self.particlesCodes[0],
                                                        nodeConditions=[('objType', 'SetOfParticles'),
                                                                        ('samplingRate', partSamplingRate)])
    firstParticlesSets = self.getCodesSets(firstParticlesCodes)
    print('Original particles code: ', firstParticlesCodes)
    notParticlesSet = self.getUnused(particleSets, firstParticlesSets)

    if self.trackParticles.get():
      print('\nTracking used particles')
      self._defineOutputs(usedParticles=particlesSet)
      self._defineOutputs(unusedParticles=notParticlesSet)

    if self.trackMics.get():
      print('\nTracking used micrographs')
      #todo: join the sets of particles or join the created outs
      particlesSet = self.joinSetsOfParticles(particleSets)
      micSets = self.getInputMicrographs()

      micDic = self.buildMicDic(micSets)
      micCountDic = self.countParticlesPerMic(particlesSet, micDic)
      outFile = self._getExtraPath('micScores.tsv')
      self.createMicScoreOutput(micCountDic, outFile)
      print('Number of used particles per micrographs written in ', outFile)

      outputMics = self.buildSetOfMics(micDic, micSets)
      self._defineOutputs(usedMicrographs=outputMics)

      if self.trackCTFs.get():
        print('\nTracking used CTFs')
        ctfSets = self.getInputCTFs()

        ctfSet = self.filterFilenamesSets(ctfSets)
        outputCTF = self.buildSetOfCTF(ctfSet, ctfSets)
        self._defineOutputs(usedCTFs=outputCTF)

      if self.trackCoordinates.get():
        print('\nTracking used coordinates')
        coordSets = self.getInputCoordinates()
        downSamplingFactor = self.getDownsampling(outputMics, particlesSet)

        coords = self.getParticleCoordinates(particlesSet)
        coords = self.scaleCoords(coords, downSamplingFactor)
        outputCoordinates = self.buildSetOfCoordinates(outputMics, coords, coordSets, '_used')
        self._defineOutputs(usedCoordinates=outputCoordinates)

        notCoords = self.getParticleCoordinates(notParticlesSet)
        notCoords = self.scaleCoords(notCoords, downSamplingFactor)
        notOutputCoordinates = self.buildSetOfCoordinates(outputMics, notCoords, coordSets, '_unused')
        self._defineOutputs(unusedCoordinates=notOutputCoordinates)

  # --------------------------- STEPS functions -----------------------------
  def getInputParticles(self):
    if self.inputParticles.get() is None:
      self.particlesCodes = self.getFurthestCodesFromNode(self.outGraph, 'root',
                                                          nodeConditions=[('objType', 'SetOfParticles')])
      particleSets = self.getCodesSets(self.particlesCodes)
    else:
      self.particlesCodes = \
        ['{}.{}'.format(self.inputParticles.getObjValue().strId(), self.inputParticles.getExtended())]
      particleSets = [self.inputParticles.get()]
    return particleSets

  def getInputMicrographs(self):
    if self.inputMicrographs.get() is None:
      micCodes = self.getClosestCodesFromNode(self.outGraph, 'root',
                                              nodeConditions=[('objType', 'SetOfMicrographs')])
      micSets = self.getCodesSets(micCodes)
    else:
      micSets = [self.inputMicrographs.get()]
    return micSets

  def getInputCTFs(self):
    if self.inputCTFs.get() is None:
      ctfCodes = self.getFurthestCodesFromNode(self.outGraph, 'root',
                                              nodeConditions=[('objType', 'SetOfCTF')])
      ctfSets = self.getCodesSets(ctfCodes)
    else:
      ctfSets = [self.inputCTFs.get()]
    return ctfSets

  def getInputCoordinates(self):
    if self.inputCoordinates.get() is None:
      coordCodes = self.getClosestCodesFromNode(self.outGraph, 'root',
                                                nodeConditions=[('objType', 'SetOfCoordinates')])
      coordSets = self.getCodesSets(coordCodes)
    else:
      coordSets = [self.inputCoordinates.get()]
    return coordSets

  def getLongestPath(self, outGraph, node1, node2):
    if nx.has_path(outGraph, node1, node2):
      allPaths = nx.all_simple_paths(outGraph, node1, node2)
    elif nx.has_path(outGraph, node2, node1):
      allPaths = nx.all_simple_paths(outGraph, node2, node1)
    else:
      return None
    return max(allPaths, key=len)

  def validateConditions(self, node, conditions):
    for condition in conditions:
      if not node[condition[0]] == condition[1]:
        return False
    return True

  def getFurthestCodesFromNode(self, outGraph, oriNode, nodeConditions):
    '''In a directed graph representing the project outputs, retrieves the output codes (protId.atrName)
    of the outputs of type setType which are furthest from the project root, following the longest possible path'''
    maxDist, furthestCodes = 0, []
    for nodeKey in outGraph.nodes:
      if nodeKey not in ['root', oriNode] and self.validateConditions(outGraph.nodes[nodeKey], nodeConditions):
        longPath = self.getLongestPath(outGraph, oriNode, nodeKey)
        if longPath is not None:
          dist = len(longPath)
          if dist > maxDist:
            maxDist = dist
            furthestCodes = [nodeKey]
          elif dist == maxDist:
            furthestCodes.append(nodeKey)
    return furthestCodes

  def getClosestCodesFromNode(self, outGraph, oriNode, nodeConditions):
    '''In a directed graph representing the project outputs, retrieves the output codes (protId.atrName)
    of the outputs of type setType which are closest from the project root, following the longest possible path'''
    minDist, closestCodes = 10000, []
    for nodeKey in outGraph.nodes:
      if not nodeKey in ['root', oriNode] and self.validateConditions(outGraph.nodes[nodeKey], nodeConditions):
        dist = len(self.getLongestPath(outGraph, oriNode, nodeKey))
        if dist < minDist:
          minDist = dist
          closestCodes = [nodeKey]
        elif dist == minDist:
          closestCodes.append(nodeKey)
    return closestCodes

  def createScipionSet(self, particles, suffix=''):
    scipionSet = self._createSetOfParticles(suffix)
    for part in particles:
      scipionSet.append(part)
    return scipionSet

  def joinParticleSetsIds(self, particleSets):
    ids = set([])
    for particles in particleSets:
      for part in particles:
        ids.add(part.getObjId())
    return ids

  def getUnused(self, usedParticlesSets, firstParticlesSets):
    usedParticlesIds = self.joinParticleSetsIds(usedParticlesSets)
    notParticles, notUsedIds = [], []

    for particleSet in firstParticlesSets:
      for part in particleSet:
        if not part.getObjId() in usedParticlesIds and not part.getObjId() in notUsedIds:
          cPart = part.clone()
          notParticles.append(cPart)
          notUsedIds.append(cPart.getObjId())

    notParticlesSet = self.createScipionSet(notParticles, '_unused')
    notParticlesSet.copyInfo(firstParticlesSets[0])
    return notParticlesSet

  def joinSetsOfParticles(self, particlesSets):
    if len(particlesSets) == 1:
      return particlesSets[0]
    elif len(particlesSets) > 1:
      #todo: efectively join sets of particles
      print('Warning: several sets of used particles detected, outputing just ', particlesSets[0])
      return particlesSets[0]

  def getDownsampling(self, iniSet, finalSet):
    return finalSet.getSamplingRate() / iniSet.getSamplingRate()

  def scaleCoords(self, coords, factor):
    for co in coords:
      co.scale(factor)
    return coords

  def buildSetOfCoordinates(self, outputMics, coords, coordSets, suffix=''):
    outputCoordinates = self._createSetOfCoordinates(outputMics, suffix)
    outputCoordinates.copyInfo(coordSets[0])
    for coord in coords:
      outputCoordinates.append(coord.clone())

    return outputCoordinates

  def buildSetOfCTF(self, ctfSet, ctfSets):
    outputCTF = self._createSetOfCTF()
    outputCTF.copyInfo(ctfSets[0])
    for ctf in ctfSet:
      mic = ctf.getMicrograph()
      ctf.copyObjId(mic)
      try:
        outputCTF.append(ctf)
      except:
        pass
    return outputCTF

  def buildSetOfMics(self, micDic, micSets):
    outputMics = self._createSetOfMicrographs()
    outputMics.copyInfo(micSets[0])
    for micId in micDic:
      mic = Micrograph(micDic[micId])
      mic.setObjId(micId)
      outputMics.append(mic)
    return outputMics

  def createMicScoreOutput(self, micCountDic, outFile):
    with open(outFile, 'w') as f:
      for micFile in micCountDic:
        f.write('{}\t{}\n'.format(micFile, micCountDic[micFile]))

  def countParticlesPerMic(self, particlesSet, micDic):
    micCountDic = {}
    for part in particlesSet:
      micId = part.getMicId()
      if micId in micDic:
        micFile = micDic[micId]
        if micFile in micCountDic:
          micCountDic[micFile] += 1
        else:
          micCountDic[micFile] = 1

    micCountDic = dict(sorted(micCountDic.items()))
    return micCountDic

  def filterFilenamesSets(self, ctfSets):
    newCtfSet, filebases = [], []
    for ctfSet in ctfSets:
      for ctf in ctfSet:
        ctfFile = ctf.getPsdFile()
        filebase = os.path.basename(ctfFile)
        if not filebase in filebases:
          newCtfSet.append(ctf.clone())
          filebases.append(filebase)
    return newCtfSet

  def getNodesIds(self, nodes):
    labs = []
    for nod in nodes:
      if nod.run != None:
        labs.append(nod.run.getObjId())
    return labs

  def getParentsDic(self, graph):
    '''Dictionary of the form {protId: [parentProtsIds]}'''
    pDic = {}
    for node in graph.getNodes():
      parentNodes = node.getParents()
      if parentNodes:
        pDic[node.run.getObjId()] = self.getNodesIds(parentNodes)
    return pDic

  def getChildrenDic(self, graph):
    '''Dictionary of the form {protId: [childrenProtsIds]}'''
    pDic = {}
    for node in graph.getNodes():
      childrenNodes = node.getChilds()
      if childrenNodes:
        pDic[node.run.getObjId()] = self.getNodesIds(childrenNodes)
    return pDic

  def buildInputsDic(self, proj):
    '''Dictionary of protocol inputs of the form:
       {protId: {key: pointer}} '''
    inpDic = {}
    for prot in proj.getRuns():
      protId = prot.getObjId()
      inpDic[protId] = {}
      for inp in prot.iterInputAttributes():
        key, ipointer = inp
        inpDic[protId][key] = ipointer
    return inpDic

  def buildOutputsDic(self, proj):
    '''Dictionary of protocol outputs of the form:
       {protId: {key: pointer}} '''
    outDic = {}
    for prot in proj.getRuns():
      protId = prot.getObjId()
      outDic[protId] = {}
      for out in prot.iterOutputAttributes():
        key, oValue = out
        outDic[protId][key] = oValue
    return outDic

  def generateOutputsGraphRec(self, curProtId, outGraph=nx.DiGraph(), prevCode=None):
    '''Generates a directed graph from a project with the protocol outputs as nodes, starting from the
    current protocol Id (curProtId) and moving upwards. Therefore, the resulting outputs graphs contains
    only those outputs necessary for the generation of the input of this protocol (default volume)'''
    inpKeys = self.inpDic[curProtId].keys()
    if len(inpKeys) == 0:
      if not 'root' in outGraph:
        outGraph.add_node('root')
      outGraph.add_edge('root', prevCode)
      return outGraph

    for k in inpKeys:
      outCodes = self.protDic[curProtId][k]
      if not type(outCodes) == list:
        outCodes = [outCodes]

      for outCode in outCodes:
        protId, outKey, outCode = self.manageCodes(outCode)
        outGraph, nextProId = self.addNode(outCode, outGraph)
        if prevCode is not None:
          outGraph.add_edge(outCode, prevCode)

        outGraph = self.generateOutputsGraphRec(nextProId, outGraph, outCode)
    return outGraph

  def manageCodes(self, outCode):
    if outCode.split('.')[1] == '':
      # Special case of subset with entire protocol as input
      parentProtId = int(outCode.split('.')[0])
      outKey = list(self.outDic[parentProtId])[0]
      outCode = str(parentProtId) + '.' + outKey
    else:
      parentProtId, outKey = int(outCode.split('.')[0]), outCode.split('.')[1]
    return parentProtId, outKey, outCode

  def addNode(self, outCode, outGraph):
    protId, outKey, outCode = self.manageCodes(outCode)

    outGraph.add_node(outCode, code=outCode, protId=protId, atrName=outKey,
                      obj=self.outDic[protId][outKey],
                      objType=self.outDic[protId][outKey].__str__().split()[0],
                      samplingRate=self.tryGetSamplingRate(self.outDic[protId][outKey]))

    return outGraph, protId

  def tryGetSamplingRate(self, scipionSet):
    try:
      return scipionSet.getSamplingRate()
    except:
      return scipionSet.getMicrographs().getSamplingRate()

  def getCodeType(self, code):
    if len(code.split('.')) == 2:
      protId, key = code.split('.')
    elif len(code.split('.')) == 3:
      protId, key, idx = code.split('.')
    else:
      print(code)

    outObj = self.outDic[int(protId)][key]
    objType = outObj.__str__().split()[0]
    #print('Code: {}\nType: {}'.format(code, objType))
    return objType

  def getUpperSetCodes(self, allOutCodes, setType = None):
    '''Select the set code which contain a setType object at Upper level (counting from earlier protocols)'''
    upper = -10000
    upperCodes = set([])
    for codeTrack in allOutCodes:
      for i in range(-1, -len(codeTrack), -1):
        if i < upper: break

        code = codeTrack[i]
        codeType = self.getCodeType(code)
        if codeType == setType:
          if i == upper:
            upperCodes.add(code)
          elif i > upper:
            upperCodes, upper = {code}, i

    return upperCodes

  def getLowerSetCodes(self, allOutCodes, setType = None):
    '''Select the set code which contain a setType object at lower level (counting from earlier protocols)'''
    lower = 10000
    lowerCodes = set([])
    for codeTrack in allOutCodes:
      for i in range(-1, -len(codeTrack), -1):
        if i > lower: break

        code = codeTrack[i]
        codeType = self.getCodeType(code)
        if codeType == setType:
          if i == lower:
            lowerCodes.add(code)
          elif i < lower:
            lowerCodes, lower = {code}, i

    return lowerCodes

  def filterTracks(self, tracks, codes):
    '''Filter the tracks which do not contain the specified codes'''
    ntracks = []
    for t in tracks:
      for code in codes:
        if code in t:
          ntracks.append(t)
          break
    return ntracks

  def getCodesSets(self, codes):
    sets = []
    for code in codes:
      if len(code.split('.')) == 2:
        protId, key = code.split('.')
      elif len(code.split('.')) == 3:
        protId, key, idx = code.split('.')

      sets.append(self.outDic[int(protId)][key])
    return sets

  def buildMicDic(self, micSets):
    micDic = {}
    for micSet in micSets:
      for mic in micSet:
        micDic[mic.getObjId()] = mic.getFileName()
    return micDic

  def getParticleCoordinates(self, particlesSet):
    coords = []
    for part in particlesSet:
      nPart = part.clone()
      newCoord = nPart.getCoordinate()
      newCoord.setObjId(nPart.getObjId())
      coords.append(newCoord)

    return coords

