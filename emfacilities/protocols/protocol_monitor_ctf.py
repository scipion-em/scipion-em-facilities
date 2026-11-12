# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import os
import sys
from math import isinf
import sqlite3 as lite
import datetime
import math
import pytz
from configparser import ConfigParser
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.constants import STATUS_RUNNING
from pyworkflow.protocol import getUpdatedProtocol

from .protocol_monitor import ProtMonitor, Monitor

PHASE_SHIFT = 'phaseShift'
TIME_STAMP = 'timeStamp'
DEFOCUS_U = 'defocusU'
RESOLUTION = 'resolution'
CTF_LOG_SQLITE = 'ctf_log.sqlite'


class ProtMonitorCTF(ProtMonitor):
    """
    Monitors the quality and stability of Contrast Transfer Function
    estimation during cryo-EM image processing workflows. The protocol
    supervises defocus behavior, astigmatism levels, phase shift
    estimation, fit quality, and related acquisition parameters in order
    to help users identify problematic micrographs or unstable imaging
    conditions early in data collection and processing.

    AI Generated:

    CTF Monitor (ProtMonitorCTF) — User Manual
        Overview

        The CTF Monitor protocol is designed to supervise the evolution
        of CTF estimation results while cryo-EM processing is running.
        Its primary objective is to provide continuous quality control
        during data acquisition or preprocessing so that imaging problems,
        microscope instabilities, or processing failures can be detected
        before they compromise an entire dataset.

        In practical cryo-EM workflows, CTF estimation represents one of
        the most important indicators of data quality. Parameters such as
        defocus, astigmatism, resolution estimation, and phase shift are
        directly related to microscope alignment, specimen stability, and
        imaging conditions. By continuously monitoring these values, users
        can rapidly identify deviations that may indicate ice contamination,
        charging effects, poor microscope tuning, stage instability, or
        processing inconsistencies.

        Biological Importance of CTF Monitoring

        For biological interpretation, reliable CTF estimation is essential
        because downstream reconstruction quality strongly depends on the
        accuracy of these measurements. Incorrect defocus estimation can
        reduce achievable resolution, while excessive astigmatism may
        indicate optical misalignment or beam instability. Monitoring these
        values in real time allows facility operators and researchers to
        intervene early and preserve data quality.

        Continuous supervision becomes especially valuable during automated
        acquisition sessions that may last many hours or even days. In such
        situations, gradual microscope drift or environmental instability
        may remain unnoticed without automated monitoring tools.

        Defocus Supervision

        One of the central functions of the protocol is the supervision of
        defocus values. Users can define acceptable minimum and maximum
        defocus ranges according to the experimental design. Defocus values
        outside these limits may indicate acquisition problems or unexpected
        microscope behavior.

        In biological workflows, the appropriate defocus range depends on
        the desired balance between contrast and high-resolution information.
        Higher defocus values improve particle visibility but reduce the
        preservation of fine structural details. Lower defocus values improve
        high-resolution transfer but may produce images with weaker contrast.
        Monitoring ensures that the acquisition remains within the intended
        operating range throughout the session.

        Astigmatism Control

        The protocol also evaluates astigmatism levels derived from the
        difference between principal defocus directions. Elevated astigmatism
        often reflects imperfect microscope alignment, contamination, or
        optical instability.

        From a practical perspective, moderate astigmatism may still allow
        acceptable reconstructions, but persistent or increasing astigmatism
        frequently predicts reduced map quality. Early detection enables
        corrective microscope alignment before additional data are collected.

        Resolution and Fit Quality

        In addition to defocus supervision, the monitor tracks estimated CTF
        resolution and fit quality metrics. These measurements provide a
        broader assessment of how reliably the experimental power spectrum is
        being modeled.

        Poor resolution estimates or unstable fit quality values may suggest
        motion contamination, crystalline ice, low signal-to-noise ratio, or
        acquisition defects. In high-throughput facilities, these indicators
        are especially useful for detecting systematic problems affecting
        large groups of micrographs.

        Phase Shift Monitoring

        When phase plates are used, the protocol can also supervise phase
        shift evolution. Stable phase shift values are important for obtaining
        reproducible contrast conditions and consistent image interpretation.
        Unexpected variations may indicate phase plate contamination, charging,
        or instability during acquisition.

        Data Logging and Visualization

        The protocol maintains a continuous record of monitored parameters,
        enabling users to inspect trends over time and correlate imaging
        quality with acquisition conditions. This historical information is
        valuable both for troubleshooting and for documenting microscope
        performance across multiple sessions.

        The monitoring system can support both local visualization and remote
        reporting environments, making it suitable for automated facilities,
        distributed collaborations, and long unattended acquisition sessions.

        Alarm and Notification System

        A major feature of the protocol is the ability to generate warnings
        whenever monitored values exceed user-defined thresholds. These alerts
        help users respond quickly to problematic acquisition conditions before
        substantial processing time or microscope time is lost.

        Typical warning scenarios include excessive astigmatism, unexpectedly
        high defocus values, unusually low defocus values, or unstable quality
        metrics. In facility environments, automated notifications can provide
        continuous operational supervision even when users are not directly
        observing the processing session.

        Practical Recommendations

        In routine cryo-EM practice, users should define threshold values that
        match the biological target and imaging strategy. High-resolution
        single-particle projects generally require tighter astigmatism and
        defocus constraints than exploratory screening datasets.

        During long automated acquisitions, it is advisable to periodically
        review the monitored trends instead of relying solely on isolated
        warnings. Progressive deterioration often becomes visible in temporal
        trends before absolute thresholds are exceeded.

        For challenging specimens, occasional outliers may be expected and do
        not necessarily invalidate the dataset. However, persistent deviations
        usually indicate underlying microscope or sample preparation problems
        that should be investigated.

        Final Perspective

        The CTF Monitor protocol serves as an automated quality assurance layer
        for cryo-EM preprocessing workflows. By continuously supervising the
        physical and computational parameters associated with CTF estimation,
        it helps researchers maintain acquisition stability, detect technical
        issues early, and maximize the reliability of downstream structural
        analysis. For modern automated cryo-EM facilities, such monitoring is
        increasingly essential for efficient and reproducible data collection.
    """
    _label = 'ctf monitor'
    _lastUpdateVersion = VERSION_1_1

    # -------------------------- DEFINE param functions ---------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')

        form.addParam('inputProtocol', params.PointerParam,
                      label="Input protocols", important=True,
                      pointerClass='ProtCTFMicrographs',
                      help="this protocol will be monitorized")
        form.addParam('samplingInterval', params.IntParam, default=60,
                      label="Sampling Interval (sec)",
                      help="Take one sample each SamplinInteval seconds")

        form.addParam('maxDefocus', params.FloatParam, default=40000,
                      label="Raise Alarm if maximum defocus (A) >",
                      help="Raise alarm if defocus is greater than given value")
        form.addParam('minDefocus', params.FloatParam,default=1000,
                      label="Raise Alarm if minimum defocus (A) <",
                      help="Raise alarm if defocus is smaller than given value")
        form.addParam('astigmatism', params.FloatParam,default=2000,
                      label="Raise Alarm if astigmatism (A) >",
                      help="Raise alarm if astigmatism is greater than given value")

        form.addParam('monitorTime', params.FloatParam, default=300,
                      label="Total Logging time (min)",
                      help="Log during this interval")

        ProtMonitor._sendMailParams(self, form)

    # -------------------------- STEPS functions ------------------------------
    def monitorStep(self):

        self.createMonitor().loop()

    def createMonitor(self):

        ctfProt = self.inputProtocol.get()
        ctfProt.setProject(self.getProject())

        ctfMonitor = MonitorCTF(ctfProt,
                                workingDir=self.workingDir.get(),
                                samplingInterval=self.samplingInterval.get(),
                                monitorTime=self.monitorTime.get(),
                                email=self.createEmailNotifier(),
                                stdout=True,
                                minDefocus=self.minDefocus.get(),
                                maxDefocus=self.maxDefocus.get(),
                                astigmatism=self.astigmatism.get())
        return ctfMonitor

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        # TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        return ['Monitor CTF defocus']


class MonitorCTF(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocol, influx=False, **kwargs):
        Monitor.__init__(self, **kwargs)

        # The CTF protocol to monitor
        self.protocol = protocol

        self.maxDefocus = kwargs['maxDefocus']
        self.minDefocus = kwargs['minDefocus']
        self.astigmatism = kwargs['astigmatism']
        self._dataBase = kwargs.get('dbName', CTF_LOG_SQLITE)
        self._tableName = kwargs.get('tableName', 'log')
        self.readCTFs = set()

        self.conn = lite.connect(os.path.join(self.workingDir, self._dataBase),
                                 isolation_level=None)
        self.influx = influx
        if self.influx:
            # get results as a list of dictionaries
            self.conn.row_factory = \
                lambda c, r: dict([(col[0], r[idx])
                                   for idx, col in enumerate(c.description)])
            # read timezone and offset
            from emfacilities.constants import (SECRETSFILE, 
                                                EMFACILITIES_HOME_VARNAME)
            _path = os.getenv(EMFACILITIES_HOME_VARNAME)
            secretsfile = os.path.join(_path, SECRETSFILE)
            confParser = ConfigParser()
            confParser.read(secretsfile)
            self.timeDelta = int(confParser.get('influx', 'timeDelta'))
            self.timeZone = confParser.get('influx', 'timeZone')
        self.cur = self.conn.cursor()

    def warning(self, msg):
        self.notify("Scipion CTF Monitor WARNING", msg)

    def initLoop(self):
        self._createTable()

    def step(self):
        prot = getUpdatedProtocol(self.protocol)
        # Create set of processed CTF from CTF protocol
        if hasattr(prot, 'outputCTF'):
            CTFset = prot.outputCTF.getIdSet()
        else:
            return False
        # find difference
        sys.stdout.flush()
        diffSet = CTFset - self.readCTFs
        setOfCTFs = prot.outputCTF
        astigmatism = self.astigmatism

        for ctfID in diffSet:
            ctf = setOfCTFs[ctfID]
            defocusU = ctf.getDefocusU()
            defocusV = ctf.getDefocusV()

            # Defocus angle
            defocusAngle = ctf.getDefocusAngle()
            if defocusAngle > 360 or defocusAngle < -360:
                defocusAngle = 0

            # Astigmatism
            astig = abs(defocusU - defocusV)

            # Resolution
            resolution = ctf.getResolution()
            if isinf(resolution):
                resolution = 0.
            
            # Fit quality
            fitQuality = ctf.getFitQuality()
            if fitQuality is None or isinf(fitQuality): 
                fitQuality = 0.

            # PhaseShift
            phaseShift = ctf.getPhaseShift() if ctf.hasPhaseShift() else 0.

            psdPath = os.path.abspath(ctf.getPsdFile())
            micPath = os.path.abspath(ctf.getMicrograph().getFileName())
            shiftPlot = (getattr(ctf.getMicrograph(), 'plotCart', None)
                         or getattr(ctf.getMicrograph(), 'plotGlobal', None))
            if shiftPlot is not None:
                shiftPlotPath = os.path.abspath(shiftPlot.getFileName())
            else:
                shiftPlotPath = ""

            if defocusU < defocusV:
                aux = defocusV
                defocusV = defocusU
                defocusU = aux
                # TODO: check if this is always true
                defocusAngle = 180. - defocusAngle
                print("ERROR: defocusU should be greater than defocusV")

            ctfCreationTime = ctf.getObjCreation()

            # get CTFs with this ids a fill table
            # do not forget to compute astigmatism
            defocus = math.sqrt(defocusV*defocusV + defocusU * defocusU)
            sql = """INSERT INTO %s(timestamp
                                    , ctfID
                                    , defocusU
                                    , defocusV
                                    , defocus
                                    , astigmatism
                                    , ratio
                                    , resolution
                                    , fitQuality
                                    , phaseShift
                                    , micPath
                                    , psdPath
                                    , shiftPlotPath)
                     VALUES("%s",%d,%f,%f,%f,%f,%f,%f,%f,%f,"%s","%s","%s");""" %\
                  (self._tableName, ctfCreationTime, ctfID, defocusU,
                   defocusV, defocus, astig, defocusU / defocusV, resolution,
                   fitQuality, phaseShift,  micPath, psdPath, shiftPlotPath)
            try:
                self.cur.execute(sql)
            except Exception as e:
                print("ERROR: saving one data point (CTF monitor). I continue")
                print(e)
                print(sql)

            if abs(defocusU - defocusV) > astigmatism:
                self.warning("Astigmatism (defocusU - defocusV)  = %f."
                             % abs(defocusU - defocusV))

            if defocusU > self.maxDefocus:
                self.warning("DefocusU (%f) is larger than defocus "
                             "maximum (%f)" % (defocusU, self.maxDefocus))
                self.maxDefocus = defocusU

            if defocusV < self.minDefocus:
                self.warning("DefocusV (%f) is smaller than defocus "
                             "minumum (%f)" % (defocusV, self.maxDefocus))
                self.minDefocus = defocusV

        self.readCTFs.update(diffSet)
        # Finish when protocol is not longer running
        return prot.getStatus() != STATUS_RUNNING

    def _createTable(self):
        self.cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now', 'localtime')),
                                ctfID INTEGER,
                                defocusU FLOAT,
                                defocusV FLOAT,
                                defocus FLOAT,
                                astigmatism FLOAT,
                                ratio FLOAT,
                                resolution FLOAT,
                                fitQuality FLOAT,
                                phaseShift FLOAT,
                                micPath STRING,
                                psdPath STRING,
                                shiftPlotPath STRING)
                                """ % self._tableName)

    def getData(self, lastId=-1):
        if self.influx:
            return self.getDataInflux(lastId)
        else:
            return self.getDataHtml()


    def getDataInflux(self, lastId=-1):
        "return query as a list of dictionaries"
        try:
            command = "select * from %s where id > %d " \
                      "order by id" % (self._tableName, lastId)
            self.cur.execute(command)
        except Exception as e:
            print("MonitorCTF, ERROR reading data from db: %s" %
                  os.path.join(self.workingDir, self._dataBase))
        # As we are using a row factory, fetchall returns a list of
        # dictionaries, each item in list(each dictionary)
        # represents a row of the table
        listOfDictionaries = self.cur.fetchall()
        for item in listOfDictionaries:
            local = pytz.timezone(self.timeZone)
            # convert dates from scipion to datetime.datetime
            for d in listOfDictionaries:
                datum = d['timestamp']
                              # string -> date time
                              # oposite -> strftime
                if isinstance(datum, str):
                    naive = datetime.datetime.strptime(datum, "%Y-%m-%d %H:%M:%S")
                    if self.timeDelta != 0:
                        naive = naive + datetime.timedelta(hours=self.timeDelta)
                elif isinstance(datum, datetime.datetime):
                    continue
                else:
                    raise Exception('Error: (CTF:getData()) Can not convert timestamp')
                local_dt = local.localize(naive, is_dst=None)
                d['timestamp'] = local_dt.astimezone(pytz.utc)


        return listOfDictionaries

    def getDataHtml(self):
        """Fill a dictionary for each label in self.labeldisk.
        The key is the label name. The value a list with
        data read from the database"""
        def get(name):
            try:
                self.cur.execute("select %s from %s" % (name, self._tableName))
            except Exception as e:
                print("MonitorCTF, ERROR reading data from db: %s" %
                      os.path.join(self.workingDir, self._dataBase))
            return [r[0] for r in self.cur.fetchall()]

        # TODO: this multiple calls to get may raise a
        # race condition. That is data may be added to the DB
        # between two get call and therefore the lengths of
        # the list can be different. IT would be better to do a single
        # call to the database and retrieve all at the same time
        # see getDataInflux for details
        data = {
            DEFOCUS_U: get('defocusU'),
            'defocusV': get('defocusV'),
            'astigmatism': get('astigmatism'),
            'ratio': get('ratio'),
            'idValues': get('ctfID'),
            RESOLUTION: get('resolution'),
            'fitQuality': get('fitQuality'),
            PHASE_SHIFT: get('phaseShift'),
            'imgMicPath': get('micPath'),
            'imgPsdPath': get('psdPath'),
            'imgShiftPath': get('shiftPlotPath'),
            TIME_STAMP: get("strftime('%s', timestamp) * 1000")
         }
        # conn.close()
        return data


