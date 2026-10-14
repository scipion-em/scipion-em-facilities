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

import sys
import time

import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol


class ProtMonitor(EMProtocol):
    """
    Provides the foundational framework for monitoring streaming cryo-EM
    protocols within Scipion environments. The protocol is intended to
    supervise ongoing processing tasks, collect operational information over
    time, and notify users when relevant conditions or events occur during
    execution.

    AI Generated:

    Monitor Framework (ProtMonitor) — User Manual
        Overview

        The Monitor framework provides a general infrastructure for supervising
        long-running or streaming cryo-EM workflows inside Scipion. Its primary
        purpose is to continuously observe the state of one or more active
        protocols while data acquisition or processing is still ongoing. This
        allows users and facility operators to detect problems early, follow
        processing evolution in real time, and react rapidly when unexpected
        conditions appear.

        In practical cryo-EM environments, monitoring becomes especially
        important during automated acquisition pipelines, high-throughput
        screening sessions, facility services, and unattended overnight
        processing. Rather than waiting until a workflow is completed, the
        monitoring system enables continuous inspection of processing quality
        and computational activity as new data are generated.

        General Monitoring Philosophy

        The framework is designed around the concept of periodic observation.
        At regular intervals, monitoring routines inspect the status of running
        protocols and evaluate whether predefined conditions have been reached.
        These observations may include processing quality indicators, resource
        consumption, reconstruction statistics, or any other workflow-specific
        metric considered important for the experiment.

        From a biological and operational perspective, continuous monitoring
        improves reliability and reproducibility by allowing users to identify
        acquisition problems, unstable processing behavior, or data-quality
        degradation before large amounts of microscope time are wasted.

        Streaming Workflow Integration

        The monitoring system is specifically intended for streaming workflows.
        In streaming cryo-EM processing, data are continuously generated while
        downstream analysis proceeds simultaneously. Under these conditions,
        users often need immediate feedback rather than delayed post-processing
        reports.

        The framework therefore operates naturally alongside acquisition and
        online processing protocols. It repeatedly checks protocol outputs
        during execution and continues until the monitored workflow finishes or
        until the configured monitoring duration expires.

        This behavior is particularly valuable in automated facility pipelines,
        where operators may supervise multiple microscope sessions at the same
        time and require centralized monitoring of all ongoing computations.

        Input Protocol Selection

        The framework allows one or several protocols to be monitored
        simultaneously. These protocols may correspond to motion correction,
        CTF estimation, particle extraction, classification, reconstruction,
        or any other streaming-compatible cryo-EM task.

        In biological practice, selecting the appropriate targets for
        monitoring depends on the critical points of the workflow. Early-stage
        monitoring is often focused on microscope stability and image quality,
        whereas later monitoring may focus on reconstruction resolution,
        classification consistency, or computational throughput.

        Monitoring Frequency and Duration

        The monitoring interval determines how frequently the system evaluates
        the monitored protocols. Short intervals provide near real-time
        feedback and are useful during microscope setup, troubleshooting, or
        critical acquisition periods. Longer intervals reduce computational
        overhead and are generally sufficient for stable production runs.

        The total monitoring duration defines how long supervision remains
        active. In extended cryo-EM sessions lasting multiple days, prolonged
        monitoring allows the framework to accompany the entire acquisition and
        processing workflow without user intervention.

        Notification System

        One of the most important aspects of the framework is its notification
        capability. When relevant events occur, the system can inform users
        through configurable notification channels. Notifications may be sent
        directly to the terminal output or through email delivery systems.

        Email notifications are especially valuable in unattended workflows,
        overnight processing, or remote facility operations. They allow users
        to receive warnings or updates without continuously supervising the
        processing interface.

        In practical terms, notifications can indicate abnormal processing
        behavior, quality degradation, threshold violations, or the completion
        of important processing milestones. This early-warning capability helps
        prevent unnecessary data collection and supports rapid decision-making
        during experiments.

        Extensible Monitoring Architecture

        The framework is intentionally generic and extensible. Specialized
        monitors can build on this infrastructure to supervise highly specific
        cryo-EM metrics such as motion correction quality, CTF estimation
        behavior, gain stability, reconstruction convergence, or hardware
        utilization.

        This extensibility allows facilities and advanced users to adapt the
        monitoring system to their own acquisition pipelines and scientific
        requirements. Different monitoring strategies can therefore coexist
        within the same Scipion ecosystem while sharing a common execution and
        notification model.

        Operational Reliability

        Continuous monitoring contributes significantly to operational
        robustness in large cryo-EM installations. Automated supervision helps
        identify failures in data transfer, unstable processing conditions,
        abnormal acquisition behavior, or quality-control deviations before
        they propagate through the workflow.

        In facility-scale environments, this capability reduces downtime,
        improves resource utilization, and increases confidence in unattended
        processing pipelines. For individual users, monitoring provides an
        additional layer of security and transparency during demanding
        experiments.

        Final Perspective

        In modern cryo-EM practice, monitoring is not simply a technical
        convenience but an essential component of reliable streaming workflows.
        Real-time supervision allows users to maintain awareness of acquisition
        quality, computational stability, and processing evolution throughout
        the entire experiment.

        By combining periodic inspection, configurable notifications, and
        extensible monitoring logic, the framework provides a flexible
        foundation for building robust online supervision systems adapted to
        the needs of both individual laboratories and large cryo-EM facilities.
    """
    _label = None

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputProtocols', params.MultiPointerParam,
                      label="Input protocols", important=True,
                      pointerClass='EMProtocol',
                      help="this protocol/s will be monitorized")

        form.addParam('samplingInterval', params.IntParam, default=60,
                      label="Sampling Interval (sec)",
                      help="Take one sample each *samplingInterval* seconds")

        form.addParam('monitorTime', params.FloatParam, default=34560,
                      label="Total Logging time (min)",
                      help="Log during this interval. 21 days by default")

    def _sendMailParams(self, form):
        g = form.addGroup('Email settings')

        g.addParam('doMail', params.BooleanParam,
                   label="Enable Email notification?", default=False,
                   help="Allow monitors to notify via email.")

        g.addParam('emailFrom', params.StringParam, condition='doMail',
                   default="from@from.fakeadress.com",
                   label='From',
                   help='Provide the sender address for notifications.')

        g.addParam('emailTo', params.StringParam, condition='doMail',
                   default="to@to.fakeadress.com",
                   label='To',
                   help='Provide the destination address for notifications.')

        g.addParam('smtp', params.StringParam, condition='doMail',
                   default="smtp.fakeadress.com",
                   label='SMTP Mail server',
                   help='Provide the address of SMTP mail server.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # -------------------------- STEPS functions ------------------------------
    def monitorStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []

    def getInputProtocols(self):
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)
        return protocols

    def sendEMail(self, emailSubject, emailMessage):
        # Import smtplib for the actual sending function
        import smtplib
        # Import the email modules we'll need
        from email.mime.text import MIMEText

        msg = MIMEText(emailMessage)

        msg['Subject'] = emailSubject
        msg['From'] = self.emailFrom.get()
        msg['To'] = self.emailTo.get()

        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP(self.smtp.get())
        s.sendmail(self.emailFrom.get(), self.emailTo.get(), msg.as_string())
        s.quit()

    def createEmailNotifier(self):
        if getattr(self, 'doMail', False):
            email = EmailNotifier(self.smtp.get(),
                                  self.emailFrom.get(),
                                  self.emailTo.get())
        else:
            email = None

        return email

    @classmethod
    def worksInStreaming(cls):
        # A monitor protocol always work in streaming
        return True


class Monitor:
    def __init__(self, **kwargs):
        # Where to store any data from this monitor
        self.workingDir = kwargs['workingDir']
        self.samplingInterval = kwargs.get('samplingInterval', None)
        self.monitorTime = kwargs.get('monitorTime', None)

        self._notifiers = []

        if kwargs.get('email', None) is not None:
            self._notifiers.append(kwargs['email'])

        if 'stdout' in kwargs:
            self._notifiers.append(PrintNotifier())

    def notify(self, title, message):
        for n in self._notifiers:
            if n: 
                n.notify(title, message)

    def info(self, message):
        self.notify("INFO", message)

    def initLoop(self):
        """ To be defined in subclasses. """
        pass

    def loop(self):
        self.initLoop()
        timeout = time.time() + 60. * self.monitorTime   # interval minutes from now

        while True:
            finished = self.step()
            if (time.time() > timeout) or finished:
                break
            time.sleep(self.samplingInterval)

    def step(self):
        """ To be defined in subclasses. """
        pass

    def addNotifier(self, notifier):
        self._notifiers.append(notifier)


class EmailNotifier:
    def __init__(self, smtpServer, emailFrom, emailTo):
        self._smtpServer = smtpServer
        self._emailFrom = emailFrom
        self._emailTo = emailTo

    def notify(self, title, message):
        # Import smtplib for the actual sending function
        import smtplib
        # Import the email modules we'll need
        from email.mime.text import MIMEText

        msg = MIMEText(message)

        msg['Subject'] = title
        msg['From'] = self._emailFrom
        msg['To'] = self._emailTo

        try:
            # Send the message via our own SMTP server, but don't include the
            # envelope header.
            s = smtplib.SMTP(self._smtpServer)
            s.sendmail(self._emailFrom, self._emailTo, msg.as_string())
            s.quit()
        except Exception as ex:
            from traceback import print_exc
            print("Some error happened while trying to send email warning.")
            print(" > Error:")
            print_exc()
            print(" > Message:")
            print(msg.as_string())


class PrintNotifier:
    def notify(self, title, message):
        print(title, message)
        sys.stdout.flush()
