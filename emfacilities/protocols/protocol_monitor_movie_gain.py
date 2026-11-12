# **************************************************************************
# *
# * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STATUS_RUNNING
from pyworkflow import VERSION_1_1

from .protocol_monitor import ProtMonitor, Monitor


class ProtMonitorMovieGain(ProtMonitor):
    """
    Monitors movie gain estimation quality during cryo-EM preprocessing workflows by supervising statistical
    indicators associated with detector gain normalization. The protocol is designed to help users identify abnormal
    gain behavior, unstable detector responses, or calibration problems that could compromise downstream image
    processing and reconstruction quality.

    AI Generated:

    Movie Gain Monitor (ProtMonitorMovieGain) — User Manual
        Overview

        The Movie Gain Monitor protocol provides continuous supervision of gain estimation results generated during
        cryo-EM preprocessing workflows. Its primary purpose is to detect abnormal detector behavior and identify
        problematic gain references before these issues propagate into alignment, CTF estimation, particle extraction,
        or high-resolution reconstruction stages.

        In cryo-EM data collection, gain normalization is a critical preprocessing operation because it compensates
        for detector response variations across the camera surface. Errors or instabilities in gain estimation may
        introduce structured artifacts, increase image noise, or bias quantitative analyses. Monitoring gain quality
        therefore becomes especially important in automated acquisition environments and high-throughput facilities.

        General Workflow

        The protocol supervises an existing movie gain estimation procedure and periodically evaluates the numerical
        statistics associated with the generated gain references. During execution, the monitoring system continuously
        checks whether these statistical indicators remain within biologically and technically acceptable ranges.

        Users typically connect the monitor to an active gain estimation workflow so that supervision occurs in parallel
        with data processing. This approach allows rapid identification of problems while acquisition or preprocessing
        is still ongoing.

        Statistical Quality Indicators

        The monitoring process evaluates several complementary statistical measurements that reflect the stability and
        uniformity of the gain reference. These indicators are designed to reveal excessive variability, abnormal
        intensity distributions, or extreme outlier behavior.

        One monitored parameter is the residual gain standard deviation. Elevated values may indicate unstable detector
        calibration, acquisition inconsistencies, or excessive noise in the estimated gain reference. In practical
        workflows, unusually large deviations can correlate with detector instability or poor normalization quality.

        Additional indicators compare different regions of the intensity distribution to evaluate whether the gain
        values remain balanced and biologically plausible. Large percentile ratios or strong deviations between extreme
        values and the main intensity population may reveal detector defects, corrupted normalization patterns, or
        unstable acquisition conditions.

        Alarm and Warning System

        The protocol supports configurable thresholds for all monitored statistics. When one or more measurements exceed
        the selected limits, warning notifications are generated automatically. These alerts allow users and facility
        operators to intervene rapidly before problematic gain references affect large portions of the dataset.

        In automated cryo-EM facilities, this type of supervision is especially valuable during unattended overnight
        acquisitions or long collection sessions where detector issues might otherwise remain unnoticed for many hours.

        Biological and Practical Relevance

        Although gain monitoring is primarily a technical quality-control procedure, it has direct consequences for
        biological interpretation. Poor gain normalization can degrade image contrast, reduce alignment accuracy, and
        limit the achievable resolution of reconstructions. Even subtle detector artifacts may propagate through
        multiple processing stages and become difficult to correct later.

        From a practical perspective, early detection of detector instability can save substantial microscope time and
        computational resources. Detecting abnormal gain behavior during acquisition is far more efficient than
        discovering normalization problems after extensive downstream processing has already been completed.

        Data Visualization and Reporting

        The protocol provides structured monitoring data suitable for both lightweight visualization environments and
        more advanced dashboard infrastructures. Statistical measurements may be displayed as evolving time series,
        allowing users to identify gradual deterioration trends, sudden instabilities, or persistent detector anomalies.

        In facility environments, these monitoring outputs can contribute to broader operational supervision systems
        that combine detector quality metrics with microscope status and computational resource monitoring.

        Practical Recommendations

        In routine cryo-EM workflows, threshold values should initially remain moderately permissive until the normal
        behavior of the detector is well characterized. Once stable operational ranges are understood, tighter
        thresholds can improve sensitivity to subtle instabilities.

        Users should pay particular attention to persistent increases in variability or recurring warnings affecting
        multiple movies. Isolated deviations may occasionally occur during acquisition, but sustained abnormalities
        often indicate detector calibration issues or environmental instability requiring intervention.

        Gain monitoring is especially recommended for high-throughput acquisition pipelines, remote data collection,
        and automated facility operation where manual supervision is limited.

        Final Perspective

        For cryo-EM users and facility operators, gain normalization quality is a foundational requirement for reliable
        image processing and reproducible structural interpretation. Continuous supervision of gain estimation behavior
        helps ensure detector stability, reduces the risk of hidden preprocessing artifacts, and supports more robust
        downstream reconstruction workflows. By providing automated quality assessment and early warning capabilities,
        the protocol contributes to more reliable and efficient cryo-EM data collection operations.
    """
    _label = 'movie gain monitor'
    _version = VERSION_1_1

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        # ProtMonitor._defineParams(self, form)
        form.addSection(label='Input')

        form.addParam('inputProtocol', params.PointerParam,
                      label="Input protocol", important=True,
                      pointerClass='XmippProtMovieGain',
                      help="This protocol will be monitorized")

        form.addParam('stddevValue', params.FloatParam, default=0.04,
                      label="Raise Alarm if residual gain standard "
                            "deviation >",
                      help="Raise alarm if residual gain standard deviation "
                           "is greater than given value")
        form.addParam('ratio1Value', params.FloatParam, default=1.15,
                      label="Raise Alarm if the ratio between the 97.5 and "
                            "2.5 percentiles >",
                      help="Raise alarm if the ratio between the 97.5 and "
                           "2.5 percentiles is greater than given value")
        form.addParam('ratio2Value', params.FloatParam, default=4.5,
                      label="Raise Alarm if the ratio between the maximum "
                            "gain value and the 97.5 percentile >",
                      help="Raise alarm if the ratio between the maximum "
                           "gain value and the 97.5 percentile is greater "
                           "than given value")

        form.addParam('monitorTime', params.FloatParam, default=300,
                      label="Total Logging time (min)",
                      help="Log during this interval")

        form.addParam('samplingInterval', params.IntParam, default=60,
                      label="Sampling Interval (sec)",
                      pointerClass='EMProtocol',
                      help="Take one sample each SamplinInteval seconds")
        ProtMonitor._sendMailParams(self, form)

    # -------------------------- STEPS functions ------------------------------
    def monitorStep(self):
        self.createMonitor().loop()

    def createMonitor(self):

        movieGainProt = self.inputProtocol.get()
        movieGainProt.setProject(self.getProject())

        movieGainMonitor = MonitorMovieGain(movieGainProt,
                                            workingDir=self.workingDir.get(),
                                            samplingInterval=self.samplingInterval.get(),
                                            monitorTime=self.monitorTime.get(),
                                            email=self.createEmailNotifier(),
                                            stdout=True,
                                            stddevValue=self.stddevValue.get(),
                                            ratio1Value=self.ratio1Value.get(),
                                            ratio2Value=self.ratio2Value.get())
        return movieGainMonitor

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        # TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        prot = self.inputProtocol.get()
        fnWarning = prot._getPath("warningsMonitor.txt")
        if not os.path.exists(fnWarning):
            summary = ["Monitor movie gain, no warnings yet."]
        else:
            fhWarning = open(fnWarning, "r")
            summary = []
            for line in fhWarning.readlines():
                summary.append(line.rstrip())
            fhWarning.close()
        return summary


class MonitorMovieGain(Monitor):
    """ This will be monitoring a movie gain estimation protocol.
    It will internally handle a database to store produced
    movie gain values.
    """
    def __init__(self, protocol, influx=False, **kwargs):
        Monitor.__init__(self, **kwargs)

        # The movieGain protocol to monitor
        self.protocol = protocol

        self.stddevValue = kwargs['stddevValue']
        self.ratio1Value = kwargs['ratio1Value']
        self.ratio2Value = kwargs['ratio2Value']
        self.influx = influx

    def warning(self, msg):
        self.notify("Scipion Movie Gain Monitor WARNING", msg)

    def initLoop(self):
        pass

    def step(self):
        prot = self.protocol
        fnSummary = prot._getPath("summaryForMonitor.txt")
        fnWarning = prot._getPath("warningsMonitor.txt")
        if not os.path.exists(fnSummary) or os.path.getsize(fnSummary) < 1:
            return False
        fhSummary = open(fnSummary, "r")
        if not os.path.exists(fnWarning):
            fhWarning = open(fnWarning, "w")
        else:
            fhWarning = open(fnWarning, "a")

        line = fhSummary.readlines()[-1]
        stddev, perc25, perc975, maxVal = map(float, line.split()[1:])
        movie_name = map(str, line.split()[0])
        values = line.split()

        if float(values[1]) > self.stddevValue:
            self.warning("Residual gain standard deviation is %f."
                         % stddev)
            fhWarning.write("%s: Residual gain standard deviation is %f.\n"
                            % (movie_name, stddev))

        if (perc975 / perc25) > self.ratio1Value:
            self.warning("The ratio between the 97.5 and 2.5 "
                         "percentiles is %f."
                         % (perc975 / perc25))
            fhWarning.write("%s: The ratio between the 97.5 and 2.5 "
                            "percentiles is %f.\n"
                            % (movie_name, (perc975 / perc25)))

        if (maxVal / perc975) > self.ratio2Value:
            self.warning("The ratio between the maximum gain value "
                         "and the 97.5 percentile is %f."
                         % (maxVal / perc975))
            fhWarning.write("%s: The ratio between the maximum gain value "
                            "and the 97.5 percentile is %f.\n"
                            % (movie_name, (maxVal / perc975)))
        fhSummary.close()
        fhWarning.close()
        return prot.getStatus() != STATUS_RUNNING

    def getData(self, lastId=-1):
        if self.influx:
            return self.getDataInflux(lastId)
        else:
            return self.getDataHtml()


    def getDataInflux(self, lastId=None):
        """retuen data as a list of dictionaries"""
        prot = self.protocol
        fnSummary = prot._getPath("summaryForMonitor.txt")
        data = []
        if not os.path.exists(fnSummary) < 1:
            fhSummary = open(fnSummary, "r")
            for idx, line in enumerate(fhSummary.readlines(), 1):  # star idx in 1
                                             # in other protocols index start in 1
                idx = int(idx)
                if idx > lastId:  # where id > last_id
                    row={}
                    stddev, perc25, perc975, maxVal = map(float, line.split()[1:])
                    row['idx'] = idx
                    row['stddev'] = stddev
                    row['ratio1'] = perc975 / perc25
                    row['ratio2'] = maxVal / perc975
                    data.append(row)
            fhSummary.close()
        # movie_000001: 0.016680 0.976350 1.028174 17.423561
        return data


    def getDataHtml(self):
        idValues = []
        stddevValues = []
        ratio1Values = []
        ratio2Values = []

        prot = self.protocol
        fnSummary = prot._getPath("summaryForMonitor.txt")
        if not os.path.exists(fnSummary) < 1:
            fhSummary = open(fnSummary, "r")
            for idx, line in enumerate(fhSummary.readlines()):
                stddev, perc25, perc975, maxVal = map(float, line.split()[1:])
                idValues.append(idx)
                stddevValues.append(stddev)
                ratio1Values.append(perc975 / perc25)
                ratio2Values.append(maxVal / perc975)
            fhSummary.close()

        data = {
            'idValues': idValues,
            'standard_deviation': stddevValues,
            'ratio1': ratio1Values,
            'ratio2': ratio2Values
        }
        return data

