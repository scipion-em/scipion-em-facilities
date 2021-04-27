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
                  label="Input Particles", condition='trackParticles', allowsNull=True,
                  help='Set of Particles where to track the used items. If None, the set of particles'
                       'closer to the volume will be used. It will also be used to find the used mics, ctfs'
                       'and coordinates.')
    
    form.addParam('trackMics', params.BooleanParam, default=True,
                  label='Track used micrographs?',
                  help='The micrographs with picked particles used in the final volume reconstruction '
                       'will be tracked and scored according to the number of used particles in them')
    form.addParam('inputMicrographs', params.PointerParam,
                  pointerClass='SetOfMicrographs', expertLevel=params.LEVEL_ADVANCED,
                  label="Input micrographs", condition='trackMics', allowsNull=True,
                  help='Set of micrographs where to track the used items. If None, the first protocol'
                       'with output micrographs will be used')

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
    self.allTracks = self.generateOutputsTrackRec(self.getObjId(), [], [])

  def getUsedItemsStep(self):
    if self.inputParticles.get() is None:
      self.particlesTracks, self.particlesCodes = self.getLowerSetTracks(self.allTracks, setType='SetOfParticles')
      particleSets = self.getCodesSets(self.particlesCodes)
    else:
      partCode = '{}.{}'.format(self.inputParticles.getObjValue().strId(), self.inputParticles.getExtended())
      self.particlesTracks, self.particlesCodes = self.getLowerSetTracks(self.allTracks, exactCode=partCode)
      particleSets = [self.inputParticles.get()]

    if self.trackParticles.get():
      print('\nTracking used particles')
      particlesSet = self.joinSetsOfParticles(particleSets)
      self._defineOutputs(usedParticles=particlesSet)

    if self.trackMics.get():
      print('\nTracking used micrographs')
      #todo: join the sets of particles or join the created outs
      particlesSet = particleSets[0]
      if self.inputMicrographs.get() is None:
        micTracks, micCodes = self.getUpperSetTracks(self.particlesTracks, setType='SetOfMicrographs')
        micSets = self.getCodesSets(micCodes)
      else:
        micSets = [self.inputMicrographs.get()]

      micDic = self.buildMicDic(micSets)
      micCountDic = self.countParticlesPerMic(particlesSet, micDic)
      outFile = self._getExtraPath('micScores.tsv')
      self.createMicScoreOutput(micCountDic, outFile)
      print('Number of used particles per micrographs written in ', outFile)

      outputMics = self.buildSetOfMics(micDic, micSets)
      self._defineOutputs(usedMicrographs=outputMics)

      if self.trackCTFs.get():
        print('\nTracking used CTFs')
        if self.inputCTFs.get() is None:
          ctfTracks, ctfCodes = self.getUpperSetTracks(self.particlesTracks, setType='SetOfCTF')
          ctfSets = self.getCodesSets(ctfCodes)
        else:
          ctfSets = [self.inputCTFs.get()]

        ctfSet = self.filterFilenamesSets(ctfSets)
        outputCTF = self.buildSetOfCTF(ctfSet, ctfSets)
        self._defineOutputs(usedCTFs=outputCTF)

      if self.trackCoordinates.get():
        print('\nTracking used coordinates')
        if self.inputCoordinates.get() is None:
          coordTracks, coordCodes = self.getUpperSetTracks(self.particlesTracks, setType='SetOfCoordinates')
          coordSets = self.getCodesSets(coordCodes)
        else:
          coordSets = [self.inputCoordinates.get()]

        coords = self.getParticleCoordinates(particlesSet)
        downSamplingFactor = self.getDownsampling(outputMics, particlesSet)
        coords = self.scaleCoords(coords, downSamplingFactor)

        outputCoordinates = self.buildSetOfCoordinates(outputMics, coords, coordSets)
        self._defineOutputs(usedCoordinates=outputCoordinates)

  # --------------------------- STEPS functions -----------------------------
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

  def buildSetOfCoordinates(self, outputMics, coords, coordSets):
    outputCoordinates = self._createSetOfCoordinates(outputMics)
    outputCoordinates.copyInfo(coordSets[0])
    for coord in coords:
        outputCoordinates.append(coord)

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
        print('Este ctf malo: ', ctf)
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

  def generateOutputsTrackRec(self, curProtId, curOutCodes, allOutCodes):
    inpKeys = self.inpDic[curProtId].keys()
    if len(inpKeys) == 0:
      allOutCodes.append(curOutCodes)
      return allOutCodes
    for k in inpKeys:
      outCodes = self.protDic[curProtId][k]
      if not type(outCodes) == list:
        outCodes = [outCodes]
      for outCode in outCodes:
        if outCode.split('.')[1] == '':
          # Special case of subset with entire protocol as input
          parentProtId = int(outCode.split('.')[0])
          outKey = list(self.outDic[parentProtId])[0]
          outCode = str(parentProtId) + '.' + outKey
        else:
          parentProtId, outKey = int(outCode.split('.')[0]), outCode.split('.')[1]

        nextOutCodes = copy.deepcopy(curOutCodes)
        nextOutCodes.append(copy.deepcopy(outCode))
        allOutCodes = self.generateOutputsTrackRec(parentProtId, nextOutCodes, allOutCodes)

    return allOutCodes

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

  def getLowerSetTracks(self, allOutCodes, setType = None, exactCode = None):
    '''Select the tracks which contain a setType object at lower level (counting from earlier protocols)'''
    lower=10000
    lowerTracks, lowerCodes = [], set([])
    for codeTrack in allOutCodes:
      for i in range(-len(codeTrack), 0):
        if i > lower:
          continue

        code = codeTrack[i]
        codeType = self.getCodeType(code)
        if codeType == setType or code == exactCode:
          if i == lower:
            lowerTracks.append(codeTrack)
            lowerCodes.add(code)
          elif i < lower:
            lowerTracks = [codeTrack]
            lowerCodes = {code}
            lower=i

    return lowerTracks, lowerCodes

  def getUpperSetTracks(self, allOutCodes, setType = None, exactCode = None):
    '''Select the tracks which contain a setType object at Upper level (counting from earlier protocols)'''
    upper = -10000
    upperTracks, upperCodes = [], set([])
    for codeTrack in allOutCodes:
      for i in range(-len(codeTrack), 0):
        if i < upper:
          continue

        code = codeTrack[i]
        codeType = self.getCodeType(code)
        if codeType == setType or code == exactCode:
          if i == upper:
            upperTracks.append(codeTrack)
            upperCodes.add(code)
          elif i > upper:
            upperTracks = [codeTrack]
            upperCodes = {code}
            upper = i

    return upperTracks, upperCodes

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

