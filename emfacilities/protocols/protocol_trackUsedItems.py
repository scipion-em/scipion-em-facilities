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

import pyworkflow.utils as pwutils
from pyworkflow.protocol import params, constants
from pwem.protocols.protocol import EMProtocol
from pwem.objects import Micrograph, Particle, Class3D, Class2D
from pwem.emlib.image import ImageHandler as ih
from pwem.emlib import MDL_XCOOR, MDL_YCOOR, MDL_MICROGRAPH_ID

import os
from joblib import delayed, Parallel
import networkx as nx
import numpy as np
from scipy.spatial.distance import cdist
from xmipp3.convert import (writeSetOfCoordinates, writeCoordsListToPosFname,
                            readPosCoordinates, readSetOfCoordsFromPosFnames)


class UsedItemsTracker(EMProtocol):
  """
  This protocol will track the items (micrographs, classes2D,...) that has been used in a scipion protocol to
  generate a final volume. If the ids have been maintained, it will also track the not used items.
  """
  _label = 'Track used items'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = constants.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    form.addParam('inputVolumes', params.PointerParam,
                  pointerClass='Volume', allowsNull=False,
                  label="Input volumes",
                  help='Select the sets of volumes to track their used items')
    form.addParam('saveJPG', params.BooleanParam, default=False,
                  label="Save images as JPG",
                  help='Save micrographs and particles as jpg')

    form.addSection(label='Tracking')
    group = form.addGroup('Particles')
    group.addParam('trackParticles', params.BooleanParam, default=True,
                  label='Track used particles?',
                  help='The particles used in the final volume reconstruction will be tracked')
    group.addParam('originalParticles', params.PointerParam,
                  pointerClass='SetOfParticles', expertLevel=params.LEVEL_ADVANCED,
                  label="Original Particles", allowsNull=True,
                  help='Set of Particles containing the complete set of original particles, including the used and '
                       'filtered particles. If None, the set of particles with same sampling rate and furthest to the '
                       'used particles will be used.')

    group = form.addGroup('Micrographs')
    group.addParam('trackMics', params.BooleanParam, default=True,
                  label='Track used micrographs?',
                  help='The micrographs with picked particles used in the final volume reconstruction '
                       'will be tracked and scored according to the number of used particles in them')
    group.addParam('originalMicrographs', params.PointerParam,
                  pointerClass='SetOfMicrographs', expertLevel=params.LEVEL_ADVANCED,
                  label="Original micrographs", condition='trackMics', allowsNull=True,
                  help='Set of micrographs containing the complete set of original mics, including the used and '
                       'filtered mics. If None, the first protocol with output micrographs will be used. '
                       'Each of these micrographs will be scored depending on the number of used particles in them')

    group = form.addGroup('CTFs')
    group.addParam('trackCTFs', params.BooleanParam, default=True,
                  label='Track used CTFs?', condition='trackMics',
                  help='The CTFs from micrographs with picked particles used in the final volume reconstruction '
                       'will be tracked and scored according to the number of used particles in them')

    group = form.addGroup('Coordinates')
    group.addParam('trackCoordinates', params.BooleanParam, default=True,
                  label='Track used coordinates?', condition='trackMics',
                  help='The coordinates used in the final volume reconstruction will be tracked')
    group.addParam('originalCoordinates', params.PointerParam,
                  pointerClass='SetOfCoordinates', expertLevel=params.LEVEL_ADVANCED,
                  label="Original Coordinates", condition='trackCoordinates and trackMics', allowsNull=True,
                  help='Set of Coordinates where to get the coordinates parameters. If None, the first protocol '
                       'with output Coordinates will be used')
    group.addParam('pickNoise', params.BooleanParam, default=False, expertLevel=params.LEVEL_ADVANCED,
                  label='Pick negative coordinates', condition='trackCoordinates and trackMics',
                  help='Picks noise from the micrographs as negative particle examples')
    group.addParam('extractNoiseNumber', params.IntParam,
                  default=-1, expertLevel=params.LEVEL_ADVANCED,
                  label='Number of noise particles', condition="pickNoise",
                  help='Number of noise particles to extract from each micrograph. '
                       'Set to -1 for extracting the same amount of noise '
                       'particles as the number true particles for that micrograph')

    group = form.addGroup('Classes')
    group.addParam('trackClasses2D', params.BooleanParam, default=True,
                  label='Track used Classes2D?',
                  help='The Classes2D used in the final volume reconstruction will be tracked')
    group.addParam('inputClasses2D', params.PointerParam,
                  pointerClass='SetOfClasses2D', expertLevel=params.LEVEL_ADVANCED,
                  label="Used Classes2D", condition='trackClasses2D', allowsNull=True,
                  help='Set of Classes2D where to track the used items. If None, the newer protocol'
                       'with output Classes2D will be used')

    group.addParam('trackClasses3D', params.BooleanParam, default=True,
                  label='Track used Classes3D?',
                  help='The Classes3D used in the final volume reconstruction will be tracked')
    group.addParam('inputClasses3D', params.PointerParam,
                  pointerClass='SetOfClasses3D', expertLevel=params.LEVEL_ADVANCED,
                  label="Used Classes3D", condition='trackClasses3D', allowsNull=True,
                  help='Set of Classes3D where to track the used items. If None, the newer protocol'
                       'with output Classes3D will be used')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions ----------------------
  def _insertAllSteps(self):
    allSteps = []
    graphStep = self._insertFunctionStep('getOutputGraphStep')
    usedStep = self._insertFunctionStep('getUsedItemsStep', prerequisites=[graphStep])
    if self.trackParticles.get():
      allSteps.append(self._insertFunctionStep('trackParticlesStep', prerequisites=[usedStep]))
    if self.trackMics.get():
      allSteps.append(self._insertFunctionStep('trackMicrographsStep', prerequisites=[usedStep]))

    if self.trackClasses2D.get():
      allSteps.append(self._insertFunctionStep('trackClasses2DStep', prerequisites=[usedStep]))

    if self.trackClasses2D.get():
      allSteps.append(self._insertFunctionStep('trackClasses3DStep', prerequisites=[usedStep]))

    if self.saveJPG.get():
      self._insertFunctionStep('saveJPGsStep', prerequisites=allSteps)

  def getOutputGraphStep(self):
    self.project = self.getProject()

    self.inpDic, self.outDic = self.buildInputsDic(self.project), self.buildOutputsDic(self.project)
    self.protDic = self.project.getProtocolsDict()
    self.outGraph = self.generateOutputsGraphRec(self.getObjId())

  def getUsedItemsStep(self):
    #Used particles: further from root
    particleSets, self.particlesCodes = self.getUsedParticles()
    self.particlesSet = self.joinSets(particleSets)

    #original particles: further from used with same sampling rate
    originalParticlesSets, originalParticlesCodes = self.getOriginalParticles()
    print('Used particles code: ', self.particlesCodes)
    print('Original particles code: ', originalParticlesCodes)
    #NotUsed: filtered from original to used
    self.notParticlesSet = self.getUnused(particleSets, originalParticlesSets)

  def trackParticlesStep(self):
      print('\nTracking used particles')
      self._defineOutputs(usedParticles=self.particlesSet)
      self._defineOutputs(notUsedParticles=self.notParticlesSet)

  def trackMicrographsStep(self):
      print('\nTracking used micrographs')
      #Original mics: closest from root
      micSets = self.getOriginalMicrographs()
      #Counting particles per mic
      self.micDic = self.buildMicDic(micSets)
      self.micCountDic, _ = self.countParticlesPerMic(self.particlesSet)
      self.createMicScoreOutput(self.micCountDic)

      self.outputMics = self.buildSetOfMics(micSets[0])
      self._defineOutputs(usedMicrographs=self.outputMics)

      if self.trackCTFs.get():
        # Used CTFs: further from root
        print('\nTracking used CTFs')
        ctfSets = self.getUsedCTFs()
        self.ctfSet = self.joinSets(ctfSets)
        #Generating psd files with image handler
        self.psdDic = self.generatePSDs()
        self.createMicScoreOutput(self.micCountDic, ctfSet=self.ctfSet)

        self.outputCTF = self.buildSetOfCTF(self.ctfSet, ctfSets)
        self._defineOutputs(usedCTFs=self.outputCTF)

      if self.trackCoordinates.get():
        # Original coords: closest from root
        print('\nTracking used coordinates')
        coordSets = self.getOriginalCoordinates()
        downSamplingFactor = self.getDownsampling(self.outputMics, self.particlesSet)
        self.boxSize = coordSets[0].getBoxSize()
        
        #Getting coordinates from particles
        coords = self.getParticleCoordinates(self.particlesSet)
        coords = self.scaleCoords(coords, downSamplingFactor)
        outputCoordinates = self.buildSetOfCoordinates(self.outputMics, coords, coordSets, '_used')
        self._defineOutputs(usedCoordinates=outputCoordinates)

        notCoords = self.getParticleCoordinates(self.notParticlesSet)
        notCoords = self.scaleCoords(notCoords, downSamplingFactor)
        notOutputCoordinates = self.buildSetOfCoordinates(self.outputMics, notCoords, coordSets, '_notUsed')
        self._defineOutputs(notUsedCoordinates=notOutputCoordinates)

        if self.pickNoise.get():
          #Picking noise particles as negative examples using xmipp software
          allCoordsSet = self.getAllCoordSet(coords, notCoords, downSamplingFactor, coordSets)
          allCoordsDir, noiseCoordsDir = self._getTmpPath('all_coords'), self._getExtraPath('noiseCoords')
          os.mkdir(allCoordsDir), os.mkdir(noiseCoordsDir)
          writeSetOfCoordinates(allCoordsDir, allCoordsSet)

          self.getNoiseCoordinates(allCoordsDir, noiseCoordsDir, self.extractNoiseNumber.get())
          noiseCoordinates = readSetOfCoordsFromPosFnames(noiseCoordsDir,
                                                     sqliteOutName=self._getExtraPath('noiseCoords.sqlite'),
                                                     setOfInputCoords=outputCoordinates)
          self._defineOutputs(noiseCoordinates=noiseCoordinates)

  def trackClasses2DStep(self):
      # Used classes2D: further from root
      print('\nTracking used classes2D')
      classes2DSets, classes2DCodes = self.getUsedClasses(dd='2D')
      sampRate2D = classes2DSets[0].getFirstItem().getFirstItem().getSamplingRate()
      #Getting used particles in class2D generation as the newer with same sampling rate
      usedParticlesSets2D = self.getClassesUsedParticles(sampRate2D)

      classes2DDic = self.buildClassesDic(classes2DSets)
      self.classes2DCountDic = self.countParticlesPerClass(usedParticlesSets2D[0], classes2DDic)
      self.usedClasses2D, notUsedClasses2D = self.getUsedClasses2D(classes2DSets, usedParticlesSets2D[0])
      self.writeUsedClasses()

      self._defineOutputs(usedClasses2D=self.usedClasses2D, notUsedClasses2D=notUsedClasses2D)

  def trackClasses3DStep(self):
    if self.trackClasses3D.get():
      # Used Classes3D: further from root
      print('\nTracking used classes3D')
      classes3DSets, classes3DCodes = self.getUsedClasses('3D')
      sampRate3D = classes3DSets[0].getFirstItem().getFirstItem().getSamplingRate()
      # Getting used particles in class3D generation as the newer with same sampling rate
      usedParticlesSets3D = self.getClassesUsedParticles(sampRate3D)

      classes3DDic = self.buildClassesDic(classes3DSets)
      self.classes3DCountDic = self.countParticlesPerClass(usedParticlesSets3D[0], classes3DDic)
      usedClasses3D, notUsedClasses3D = self.getUsedClasses3D(classes3DSets, usedParticlesSets3D[0])

      self._defineOutputs(usedClasses3D=usedClasses3D, notUsedClasses3D=notUsedClasses3D)

  def saveJPGsStep(self):
    if self.trackParticles.get():
      self.saveSetAsJPG(self.particlesSet, self._getExtraPath('usedParticles'))

    if self.trackMics.get():
      self.saveSetAsJPG(self.outputMics, self._getExtraPath('usedMics'))
      self.createMicScoreOutput(self.micCountDic, jpg=True)

      if self.trackCTFs.get():
        self.savePSDsAsJPG(self._getExtraPath('usedPSDs'))
        self.createMicScoreOutput(self.micCountDic, ctfSet=self.ctfSet, jpg=True)

    if self.trackClasses2D.get():
      self.saveSetAsJPG(self.usedClasses2D, self._getExtraPath('usedClasses2D'))
      self.writeUsedClasses(jpg=True)

  # --------------------------- UTILS functions -----------------------------
  def getUsedParticles(self):
    particlesCodes = self.getFurthestCodesFromNode(self.outGraph, 'root',
                                                   nodeConditions=[('objType', 'SetOfParticles')])
    particleSets = self.getCodesSets(particlesCodes)
    return particleSets, particlesCodes

  def getOriginalParticles(self):
    if self.originalParticles.get() is None:
      partSamplingRate = self.particlesSet.getSamplingRate()
      originalParticlesCodes = self.getFurthestCodesFromNode(self.outGraph, self.particlesCodes[0],
                                                             nodeConditions=[('objType', 'SetOfParticles'),
                                                                             ('samplingRate', partSamplingRate)])
      particleSets = self.getCodesSets(originalParticlesCodes)
    else:
      originalParticlesCodes = \
        ['{}.{}'.format(self.originalParticles.getObjValue().strId(), self.originalParticles.getExtended())]
      particleSets = [self.originalParticles.get()]
    return particleSets, originalParticlesCodes

  def getUnused(self, usedParticlesSets, originalParticlesSets):
    usedParticlesIds = self.joinParticleSetsIds(usedParticlesSets)
    notParticles, notUsedIds = [], []

    for particleSet in originalParticlesSets:
      for part in particleSet:
        if not part.getObjId() in usedParticlesIds and not part.getObjId() in notUsedIds:
          cPart = part.clone()
          notParticles.append(cPart)
          notUsedIds.append(cPart.getObjId())

    notParticlesSet = self.createScipionSet(notParticles, '_unused')
    notParticlesSet.copyInfo(originalParticlesSets[0])
    return notParticlesSet

  def getOriginalMicrographs(self):
    if self.originalMicrographs.get() is None:
      micCodes = self.getClosestCodesFromNode(self.outGraph, 'root',
                                              nodeConditions=[('objType', 'SetOfMicrographs')])
      print('original mics codes: ', micCodes)
      micSets = self.getCodesSets(micCodes)
    else:
      micSets = [self.originalMicrographs.get()]
    return micSets

  def getUsedCTFs(self):
    ctfCodes = self.getFurthestCodesFromNode(self.outGraph, 'root',
                                             nodeConditions=[('objType', 'SetOfCTF')])
    print('Used CTF codes: ', ctfCodes)
    ctfSets = self.getCodesSets(ctfCodes)
    return ctfSets

  def getOriginalCTFs(self):
    if self.originalCTFs.get() is None:
      ctfCodes = self.getClosestCodesFromNode(self.outGraph, 'root',
                                              nodeConditions=[('objType', 'SetOfCTF')])
      print('Original CTF codes: ', ctfCodes)
      ctfSets = self.getCodesSets(ctfCodes)
    else:
      ctfSets = [self.originalCTFs.get()]
    return ctfSets

  def getOriginalCoordinates(self):
    if self.originalCoordinates.get() is None:
      coordCodes = self.getClosestCodesFromNode(self.outGraph, 'root',
                                                nodeConditions=[('objType', 'SetOfCoordinates')])
      print('Original Coordinate codes: ', coordCodes)
      coordSets = self.getCodesSets(coordCodes)
    else:
      coordSets = [self.originalCoordinates.get()]
    return coordSets

  def getUsedClasses(self, dd='3D'):
    if (dd=='3D' and self.inputClasses3D.get() is None) or \
            (dd=='2D' and self.inputClasses2D.get() is None):
      clCodes = self.getFurthestCodesFromNode(self.outGraph, 'root',
                                              nodeConditions=[('objType', 'SetOfClasses{}'.format(dd))])
      print('Used classes{} codes: '.format(dd), clCodes)
      clSets = self.getCodesSets(clCodes)
    else:
      if dd=='3D':
        clSets = [self.inputClasses3D.get()]
        clCodes = [self.inputClasses3D.getUniqueId()]
      elif dd=='2D':
        clSets = [self.inputClasses2D.get()]
        clCodes = [self.inputClasses2D.getUniqueId()]
    return clSets, clCodes

  def getUsedClasses2D(self, classes2DSets, oriParticlesSet):
    '''Returns the used and not used 2Dclasses sets based on the count of used particles per
    class (classes2DCountDic)'''
    usedClasses2D = self._createSetOfClasses2D(oriParticlesSet, suffix='_used')
    notUsedClasses2D = self._createSetOfClasses2D(oriParticlesSet, suffix='_notUsed')
    for cSet in classes2DSets:
      for cl2D in cSet:
        newClass = Class2D()
        newClass.copyInfo(cl2D)
        newClass.setObjId(cl2D.getObjId())

        if cl2D.getObjId() in self.classes2DCountDic:
          usedClasses2D = self.updateClassSet(usedClasses2D, newClass, cl2D)
        else:
          notUsedClasses2D = self.updateClassSet(notUsedClasses2D, newClass, cl2D)
    return usedClasses2D, notUsedClasses2D

  def getUsedClasses3D(self, classes3DSets, oriParticlesSet, writeFile='usedClasses3D.tsv'):
    '''Returns the used and not used 3Dclasses sets based on the count of used particles per
    class (classes3DCountDic)'''
    usedClasses3D = self._createSetOfClasses3D(oriParticlesSet, suffix='_used')
    notUsedClasses3D = self._createSetOfClasses3D(oriParticlesSet, suffix='_notUsed')
    f = open(self._getExtraPath(writeFile),'w')
    f.write('Class3D Id\tNumber of particles\n')
    for cSet in classes3DSets:
      for cl3D in cSet:
        newClass = Class3D()
        newClass.copyInfo(cl3D)
        newClass.setObjId(cl3D.getObjId())

        if cl3D.getObjId() in self.classes3DCountDic:
          f.write('{}\t{}\n'.format(cl3D.getObjId(), self.classes3DCountDic[cl3D.getObjId()]))
          usedClasses3D = self.updateClassSet(usedClasses3D, newClass, cl3D)
        else:
          notUsedClasses3D = self.updateClassSet(notUsedClasses3D, newClass, cl3D)
    f.close()
    return usedClasses3D, notUsedClasses3D

  def getClassesUsedParticles(self, sampRate):
    usedParticlesCodes = self.getFurthestCodesFromNode(self.outGraph, 'root',
                                                           nodeConditions=[('objType', 'SetOfParticles'),
                                                                           ('samplingRate', sampRate)])
    print('Classes particles codes: ', usedParticlesCodes)
    particleSets = self.getCodesSets(usedParticlesCodes)
    return particleSets

  def generatePSDs(self):
    outDir = self._getExtraPath('computedPSDs')
    os.mkdir(outDir)
    psdDic, argsList = {}, []
    print('Computing micrographs PSD')
    for micId, micFile in self.micDic.items():
      micBase = os.path.basename(micFile)
      micExt = pwutils.getExt(micFile)
      psdFile = micBase.replace(micExt, '_psd.psd')
      psdPath = os.path.join(outDir, psdFile)
      psdDic[micFile] = psdPath

      argsList.append((micFile, psdPath))

    Parallel(n_jobs=self.numberOfThreads.get(), backend="multiprocessing", verbose=1)(
      delayed(computeMicPSD)(*arg) for arg in argsList)

    return psdDic

  def getNoiseCoordinates(self, coordsDir, outDir, noiseNumber=-1):
    micsBaseToFullName = {}
    for micId, micPath in self.micDic.items():
      if micPath.endswith(".mrc") or micPath.endswith(".tif"):
        baseName, _ = os.path.splitext(os.path.basename(micPath))
        micsBaseToFullName[baseName] = micPath

    I = ih().read(micPath)
    xDim, yDim, _, _ = I.getDimensions()

    argsList = []
    for posName in os.listdir(coordsDir):
      if posName.endswith(".pos"):
        baseName, _ = os.path.splitext(os.path.basename(posName))
        micName = micsBaseToFullName[baseName]
        posName = os.path.join(coordsDir, posName)
        argsList.append((baseName, micName, posName, (xDim, yDim), outDir,
                         noiseNumber, self.boxSize))
    Parallel(n_jobs=self.numberOfThreads.get(), backend="multiprocessing", verbose=1)(
      delayed(pickNoiseOneMic)(*arg) for arg in argsList)

  def getMicsDir(self):
    for mic in self.outputMics:
      dir = os.path.dirname(mic.getFileName())
      return dir

  def getAllCoordSet(self, coords, notCoords, downSamplingFactor, coordSets):
    allCoords = coords + notCoords
    allCoords = self.scaleCoords(allCoords, downSamplingFactor)
    allCoordsSet = self.buildSetOfCoordinates(self.outputMics, allCoords, coordSets, '_all')
    return allCoordsSet

  def writeUsedClasses(self, writeFile='usedClasses2D.tsv', jpg=False):
    with open(self._getExtraPath(writeFile), 'w') as f:
      if not jpg:
        f.write('Class2D Id\tNumber of particles\n')
      else:
        f.write('Class2D representative\tNumber of particles\n')
      for cl2dId in self.classes2DCountDic:
        if not jpg:
          f.write('{}\t{}\n'.format(cl2dId, self.classes2DCountDic[cl2dId]))
        else:
          f.write('{}\t{}\n'.format('class2D_{}.jpg'.format(cl2dId), self.classes2DCountDic[cl2dId]))

  def savePSDsAsJPG(self, outDir):
    os.mkdir(outDir)
    for micFile, psdFile in self.psdDic.items():
      jpgPath = outDir + '/' + os.path.basename(pwutils.replaceExt(psdFile, 'jpg'))
      self.exportAsJpg(psdFile, jpgPath)

  def saveSetAsJPG(self, scipionSet, outDir):
    os.mkdir(outDir)
    for item in scipionSet:
      cItem = item.clone()
      if isinstance(cItem, Particle):
        itemId = cItem.getObjId()
        itemPath = cItem.getLocation()
        itemName = cItem.getBaseName()
        jpgPath = outDir + '/' + pwutils.replaceExt(itemName, 'jpg').replace('.jpg', '_{}.jpg'.format(itemId))

      elif isinstance(cItem, Micrograph):
        itemPath = cItem.getLocation()
        itemName = cItem.getBaseName()
        jpgPath = outDir + '/' + pwutils.replaceExt(itemName, 'jpg')

      elif isinstance(cItem, Class2D):
        itemPath = cItem.getRepresentative().getLocation()
        jpgPath = outDir + '/' + 'class2D_{}.jpg'.format(cItem.getObjId())

      self.exportAsJpg(itemPath, jpgPath)

  def convertPSDFile(self, psdFile):
    psd = ih().read(psdFile)
    psd.convertPSD()

    psdName = psdFile.split('/')[-1]
    outFile = self._getTmpPath(psdName)
    psd.write(outFile)
    return outFile

  def exportAsJpg(self, itemPath, jpgPath):
    ih().convert(itemPath, jpgPath)

  def updateClassSet(self, outClasses, newClass, cl3D):
    '''Updated the set of classes "outClassess" with a "newClass", conformed by the particles
    of a previous class "cl3D"'''
    outClasses.append(newClass)
    enabledClass = outClasses[newClass.getObjId()]
    enabledClass.enableAppend()
    for part in cl3D:
      enabledClass.append(part)
    enabledClass.setSamplingRate(part.getSamplingRate())
    outClasses.update(enabledClass)
    return outClasses

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

  def joinSets(self, sets):
    if len(sets) == 1:
      return sets[0]
    elif len(sets) > 1:
      # todo: efectively join sets of particles
      print('Warning: several sets of used particles detected, outputing just ', sets[0])
      return sets[0]

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

  def buildSetOfMics(self, micSet):
    outputMics = self._createSetOfMicrographs()
    outputMics.copyInfo(micSet)
    for micId in self.micDic:
      mic = Micrograph(self.micDic[micId])
      mic.setObjId(micId)
      outputMics.append(mic)
    return outputMics

  def createMicScoreOutput(self, micCountDic, ctfSet=None, jpg=False, outName='micScores.tsv'):
    outFile = self._getExtraPath(outName)
    with open(outFile, 'w') as f:
      if ctfSet is None:
        f.write('Micrograph filename\tNumber of particles\n')
        for micFile in micCountDic:
          if jpg:
            micName = pwutils.replaceExt(micFile, 'jpg').split('/')[-1]
          else:
            micName = micFile.split('/')[-1]
          f.write('{}\t{}\n'.format(micName, micCountDic[micFile]))

      else:
        f.write('Micrograph filename\tNumber of particles\tPSD file\t'
                'DefocusU\tDefocusV\tDefocusAngle\n')
        ctfDic = self.getMic2CTFDic(ctfSet)
        for micFile in micCountDic:
          psdPath = self.psdDic[micFile]
          if jpg:
            micName = pwutils.replaceExt(micFile, 'jpg').split('/')[-1]
            psdName = pwutils.replaceExt(psdPath, 'jpg').split('/')[-1]
          else:
            micName, psdName = micFile.split('/')[-1], psdPath.split('/')[-1]

          if micFile in ctfDic:
            dfU, dfV = ctfDic[micFile].getDefocusU(), ctfDic[micFile].getDefocusV()
            dfA = ctfDic[micFile].getDefocusAngle()
          else:
            dfU, dfV, dfA = '-', '-', '-'
          f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.
                  format(micName, micCountDic[micFile], psdName, dfU, dfV, dfA))

  def getMic2CTFDic(self, ctfSet):
    dic = {}
    for ctf in ctfSet:
      dic[ctf.getMicrograph().getFileName()] = ctf.clone()
    return dic

  def countParticlesPerMic(self, particlesSet):
    micCountDic, nMicDic = {}, {}

    for part in particlesSet:
      micId = part.getMicId()
      if micId in self.micDic:
        micFile = self.micDic[micId]
        if micFile in micCountDic:
          micCountDic[micFile] += 1
        else:
          micCountDic[micFile] = 1
          nMicDic[micId] = micFile

    for micId, micFile in self.micDic.items():
      if not micFile in micCountDic:
        micCountDic[micFile] = 0

    micCountDic = dict(sorted(micCountDic.items()))
    return micCountDic, nMicDic

  def countParticlesPerClass(self, particlesSet, classDic):
    classCountDic = {}
    for part in particlesSet:
      partId = part.getObjId()
      if partId in classDic:
        classId = classDic[partId]
        if classId in classCountDic:
          classCountDic[classId] += 1
        else:
          classCountDic[classId] = 1
    return classCountDic

  def getLongestPath(self, outGraph, node1, node2, dir='both'):
    '''Return the longest path between two nodes in a directed graph'''
    if nx.has_path(outGraph, node1, node2) and dir in ['both', 'down']:
      allPaths = nx.all_simple_paths(outGraph, node1, node2)
      return max(allPaths, key=len)
    if nx.has_path(outGraph, node2, node1) and dir in ['both', 'up']:
      allPaths = nx.all_simple_paths(outGraph, node2, node1)
      return max(allPaths, key=len)
    return None

  def validateConditions(self, node, conditions):
    '''Validates that a node has the specified conditions'''
    for condition in conditions:
      if not node[condition[0]] == condition[1]:
        return False
    return True

  def checkPathExists(self, outGraph, oriNode, targetNode, dir='both'):
    if nx.has_path(outGraph, oriNode, targetNode) and dir in ['both', 'down']:
      return True
    elif nx.has_path(outGraph, targetNode, oriNode) and dir in ['both', 'up']:
      return True
    else:
      return False

  def getGraphNeighbors(self, oriNode):
    for nodeKey in self.outGraph.successors(oriNode):
      yield nodeKey
    for nodeKey in self.outGraph.predecessors(oriNode):
      yield nodeKey

  def getNeighborsCodes(self, oriNode, nodeConditions, dir='both'):
    '''In a directed graph representing the project outputs, retrieves the output codes (protId.atrName)
        of the outputs of type setType which are neighbors from oriNode, following the longest possible path'''
    neighCodes = []
    for nodeKey in self.getGraphNeighbors(oriNode):
      if nodeKey != 'root' and self.validateConditions(self.outGraph.nodes[nodeKey], nodeConditions):
        if self.checkPathExists(self.outGraph, oriNode, nodeKey, dir):
          neighCodes.append(nodeKey)
    return neighCodes

  def getFurthestCodesFromNode(self, outGraph, oriNode, nodeConditions, dir='both'):
    '''In a directed graph representing the project outputs, retrieves the output codes (protId.atrName)
    of the outputs of type setType which are furthest from oriNode, following the longest possible path'''
    maxDist, furthestCodes = 0, []
    for nodeKey in outGraph.nodes:
      if nodeKey != 'root' and nodeKey.split('.')[0] != oriNode.split('.')[0] and \
              self.validateConditions(outGraph.nodes[nodeKey], nodeConditions):
        longPath = self.getLongestPath(outGraph, oriNode, nodeKey, dir)
        if longPath is not None:
          dist = len(longPath)
          if dist > maxDist:
            maxDist = dist
            furthestCodes = [nodeKey]
          elif dist == maxDist:
            furthestCodes.append(nodeKey)
    return furthestCodes

  def getClosestCodesFromNode(self, outGraph, oriNode, nodeConditions, dir='both'):
    '''In a directed graph representing the project outputs, retrieves the output codes (protId.atrName)
    of the outputs of type setType which are closest from oriNode, following the longest possible path'''
    minDist, closestCodes = 10000, []
    for nodeKey in outGraph.nodes:
      if nodeKey != 'root' and nodeKey.split('.')[0] != oriNode.split('.')[0] and \
              self.validateConditions(outGraph.nodes[nodeKey], nodeConditions):
        dist = len(self.getLongestPath(outGraph, oriNode, nodeKey, dir))
        if dist < minDist:
          minDist = dist
          closestCodes = [nodeKey]
        elif dist == minDist:
          closestCodes.append(nodeKey)
    return closestCodes

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
    # print('Code: {}\nType: {}'.format(code, objType))
    return objType

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

  def buildClassesDic(self, classSets):
    '''Returns a dictionary: {partId: classID}'''
    classDic = {}
    for classSet in classSets:
      for classObj in classSet:
        for part in classObj:
          classDic[part.getObjId()] = classObj.getObjId()
    return classDic

  def getParticleCoordinates(self, particlesSet):
    coords = []
    for part in particlesSet:
      nPart = part.clone()
      newCoord = nPart.getCoordinate()
      newCoord.setObjId(nPart.getObjId())
      coords.append(newCoord)

    return coords

  def getCompleteSet(self, nodeKey):
    '''A set is considered to be complete if there was no discarded item
    (the maximum id equals the number of items)'''
    protId, outKey, outCode = self.manageCodes(nodeKey)
    scipionSet = self.outDic[protId][outKey]
    try:
      for item in scipionSet:
        if item.getObjId() > len(scipionSet):
          return False
      return True
    except:
      return False

  def searchCompleteSet(self, curCodes):
    for curCode in curCodes:
      if self.getCompleteSet(curCode):
        return curCode
    return None

  def searchOriginalClassesSet(self, clCodes, dd='3D'):
    '''Iterably looks to upward neighbors from used set to search for a complete set (max(ids) == len(ids))'''
    conditions = [('objType', 'SetOfClasses{}'.format(dd))] #

    completeCode = self.searchCompleteSet(clCodes)
    curCodes = clCodes
    while completeCode == None and len(curCodes)>0:
      neighCodes = []
      for curCode in curCodes:
        neighCodes += self.getNeighborsCodes(curCode, conditions, dir='up')

      curCodes = neighCodes.copy()
      completeCode = self.searchCompleteSet(curCodes)

    oriClassesSets = self.getCodesSets([completeCode])
    return oriClassesSets[0], completeCode


def computeMicPSD(micFile, psdPath):
  imageMic = ih().read(micFile)
  dimX, dimY, _, _ = imageMic.getDimensions()
  psd1 = imageMic.computePSD(0.4, dimX, dimY, 1)
  psd1.convertPSD()
  psd1.write(psdPath)

def pickNoiseOneMic(baseName, mic_fname, posName, mic_shape, outputRoot,
                    extractNoiseNumber, boxSize):
  """ Pick noise from one micrograph
  """

  #print("pick noise one mic %s %s %s %s %s %s %s" % (baseName, mic_fname, posName, mic_shape,
                                                     #outputRoot, extractNoiseNumber, boxSize))

  coords_in_mic_list = readPosCoordsFromFName(posName)

  extractNoiseNumber = extractNoiseNumber if extractNoiseNumber != -1 else len(coords_in_mic_list)

  min_required_distance = boxSize // 2
  currentCoords = np.array(coords_in_mic_list)
  if currentCoords.shape[0] > 0:
    good_new_coords = []
    n_good_coords = 0
    for iterNum in range(9999):
      randomCoordsX = np.random.randint(boxSize // 2, mic_shape[0] - boxSize // 2, size=extractNoiseNumber * 2)
      randomCoordsY = np.random.randint(boxSize // 2, mic_shape[1] - boxSize // 2, size=extractNoiseNumber * 2)
      randomCoords = np.stack([randomCoordsX, randomCoordsY], axis=1)
      del randomCoordsX, randomCoordsY
      dist_mat = cdist(randomCoords, currentCoords)  # dist_mat: i:random_particles   j: currentCoords
      dist_mat = np.min(dist_mat, axis=1)
      g_n_c = randomCoords[dist_mat >= min_required_distance]
      n_good_coords += g_n_c.shape[0]
      good_new_coords.append(g_n_c)
      if n_good_coords >= extractNoiseNumber:
        break
    good_new_coords = np.concatenate(good_new_coords)[:extractNoiseNumber]
    writeCoordsListToPosFname(mic_fname, good_new_coords, outputRoot)

def readPosCoordsFromFName(fname, returnAlsoMicId=False):
  mData= readPosCoordinates(fname)
  coords=[]
  mdId=None
  micId=None
  for mdId in mData:
    x = int( mData.getValue( MDL_XCOOR, mdId) )
    y = int( mData.getValue( MDL_YCOOR, mdId) )
    coords.append((x,y) )
  print("N coords: %d"%(len(coords) ))
  if returnAlsoMicId:
    if mdId:
      micId= mData.getValue( MDL_MICROGRAPH_ID, mdId)
    return coords, micId
  else:
    return coords
