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

import os.path

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1

from pwem.protocols import ProtCTFMicrographs, ProtAlignMovies
from pwem import Domain
import subprocess

from .report_influx import ReportInflux
from .report_html import ReportHtml
from .protocol_monitor import ProtMonitor, Monitor
from .protocol_monitor_ctf import MonitorCTF
from .protocol_monitor_movie_gain import MonitorMovieGain
from .protocol_monitor_system import MonitorSystem
from pyworkflow import BETA, UPDATED, NEW, PROD


class ProtMonitorSummary(ProtMonitor):
    """
    Provides an integrated monitoring and reporting environment for key stages of cryo-EM data processing workflows,
    including movie alignment, CTF estimation, movie gain evaluation, and computational resource supervision. The
    protocol is intended to give facility operators, platform administrators, and cryo-EM users a unified overview
    of both data quality and system health during ongoing processing campaigns.

    AI Generated:

    Monitor Summary (ProtMonitorSummary) — User Manual
        Overview

        The Monitor Summary protocol combines several monitoring utilities into a single reporting framework designed
        for continuous supervision of cryo-EM processing pipelines. Its primary goal is to help users detect data
        quality problems, computational bottlenecks, and hardware instability while acquisitions or processing jobs
        are still running. By consolidating multiple monitoring sources into a unified report, the protocol supports
        rapid decision making and early intervention before significant processing time or microscope resources are lost.

        In practical cryo-EM environments, monitoring is especially important during automated or unattended workflows.
        Long acquisition sessions, large movie collections, and GPU-intensive processing can produce failures that may
        remain unnoticed for hours without active supervision. This protocol addresses that need by generating live
        summaries of processing behavior and system activity.

        General Workflow

        The protocol operates as a centralized supervisor that coordinates several independent monitoring tasks. These
        tasks may include movie gain analysis, CTF quality inspection, CPU and memory tracking, GPU usage monitoring,
        network activity observation, and disk throughput analysis. Each monitoring component contributes information
        to a common reporting interface.

        Users typically connect the protocol to existing processing steps already present in a Scipion workflow. The
        monitoring process then runs periodically during execution and updates reports as new information becomes
        available. This design allows the monitoring system to evolve continuously alongside the processing pipeline.

        Movie Gain Monitoring

        One important function of the protocol is the supervision of movie gain estimation quality. Gain references are
        essential for detector normalization, and problems in gain estimation may introduce structured artifacts or
        reduce overall reconstruction quality. The monitoring system evaluates statistical indicators associated with
        gain behavior and can raise alarms when abnormal variability or intensity distributions are detected.

        From a biological and practical perspective, sudden changes in gain statistics may indicate detector instability,
        acquisition inconsistencies, or calibration problems. Early detection helps prevent the accumulation of unusable
        datasets and supports more stable downstream image processing.

        CTF Quality Supervision

        The protocol also supervises contrast transfer function estimation results. Monitoring defocus ranges and
        astigmatism values provides a rapid overview of acquisition quality and microscope stability. Excessively large
        defocus values, unusually small defocus ranges, or strong astigmatism can indicate focusing problems, optical
        instability, contamination, or acquisition setup issues.

        In routine cryo-EM practice, continuous monitoring of CTF parameters is particularly valuable during overnight
        acquisition sessions or automated facility operation. Detecting problematic trends early may allow microscope
        operators to intervene before large portions of a dataset become compromised.

        System Resource Monitoring

        Beyond image quality assessment, the protocol supervises computational infrastructure during processing.
        Monitoring CPU utilization, memory allocation, swap activity, GPU load, network traffic, and disk input/output
        helps identify resource saturation and infrastructure limitations that may slow or destabilize workflows.

        GPU supervision is especially important in modern cryo-EM processing because many alignment and reconstruction
        tasks rely heavily on accelerator hardware. Monitoring GPU usage and memory consumption allows users to verify
        that computational resources are being used efficiently and helps detect overloaded or malfunctioning devices.

        Network and disk activity monitoring are particularly relevant in shared facilities or distributed processing
        environments where large cryo-EM datasets are continuously transferred between storage systems and compute
        nodes. Sustained bottlenecks in these areas may strongly affect processing throughput.

        Alarm and Notification System

        The protocol supports configurable warning thresholds for several monitored quantities. Users can define limits
        for resource utilization and quality indicators so that abnormal conditions trigger notifications. This approach
        allows facilities to implement proactive monitoring policies adapted to their own hardware and acquisition
        standards.

        In production environments, automated alarms are especially useful for minimizing downtime and avoiding wasted
        microscope time. Notifications can help users react rapidly to hardware overload, unstable GPU behavior,
        problematic gain estimation, or deteriorating image quality.

        HTML and Dashboard Reporting

        A central feature of the protocol is the generation of continuously updated reports summarizing all monitored
        information. Reports may be produced as standard HTML summaries or integrated into external visualization
        infrastructures based on Grafana and InfluxDB.

        HTML reports are suitable for lightweight deployment and local monitoring environments. They provide an
        accessible overview of acquisition and processing activity that can easily be shared within a laboratory or
        facility. Dashboard-oriented infrastructures offer a more scalable solution for larger installations requiring
        long-term monitoring, centralized visualization, or multi-user access.

        Publication and Remote Access

        The protocol supports automatic publication of generated reports through external commands. This capability is
        useful for facilities that maintain remote dashboards, institutional web portals, or centralized monitoring
        servers. Reports may therefore remain accessible even when processing occurs on remote clusters or isolated
        compute systems.

        From an operational perspective, centralized report publication simplifies supervision of multiple ongoing
        projects and allows facility staff to monitor workflows without direct access to the processing nodes.

        Practical Recommendations

        In most cryo-EM facilities, it is advisable to enable both image-quality monitoring and system-resource
        supervision simultaneously. Problems in computational infrastructure frequently correlate with processing
        instability, while image-quality metrics provide direct biological feedback about acquisition performance.

        GPU monitoring should generally be activated whenever acceleration hardware is used extensively. Similarly,
        disk and network supervision become increasingly important in high-throughput acquisition facilities where
        storage bandwidth may become a limiting factor.

        Thresholds for alarms should initially remain conservative and later be adapted to the normal behavior of the
        local microscope and computational infrastructure. Excessively strict thresholds may produce unnecessary alarms,
        while overly permissive values may delay the detection of important problems.

        Final Perspective

        For cryo-EM users and facility administrators, workflow monitoring is not only a technical convenience but an
        essential operational safeguard. Continuous supervision of image quality, detector behavior, and computational
        infrastructure improves reliability, reduces wasted acquisition time, and supports more reproducible biological
        results. By combining processing supervision with infrastructure monitoring into a unified reporting framework,
        the protocol provides a practical foundation for stable and efficient cryo-EM operations.
    """
    _label = 'monitor summary'
    _lastUpdateVersion = VERSION_1_1
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self.reportDir = ''
        self.reportPath = ''

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        form.addSection('MovieGain Monitor')
        form.addParam('stddevValue', params.FloatParam, default=0.04,
                      label="Raise Alarm if residual gain standard "
                            "deviation >",
                      help="Raise alarm if residual gain standard deviation "
                           "is greater than given value")
        form.addParam('ratio1Value', params.FloatParam, default=1.15,
                      label="Raise Alarm if the ratio between the 97.5 "
                            "and 2.5 percentiles >",
                      help="Raise alarm if the ratio between the 97.5 "
                           "and 2.5 percentiles is greater than given value")
        form.addParam('ratio2Value', params.FloatParam, default=4.5,
                      label="Raise Alarm if the ratio between the maximum "
                            "gain value and the 97.5 percentile >",
                      help="Raise alarm if the ratio between the maximum "
                           "gain value and the 97.5 percentile is greater "
                           "than given value")

        form.addSection('CTF Monitor')
        form.addParam('maxDefocus', params.FloatParam, default=40000,
                      label="Raise Alarm if maximum defocus (A) >",
                      help="Raise alarm if defocus is greater than given "
                           "value")
        form.addParam('minDefocus', params.FloatParam, default=1000,
                      label="Raise Alarm if minimum defocus (A) <",
                      help="Raise alarm if defocus is smaller than given "
                           "value")
        form.addParam('astigmatism', params.FloatParam, default=1000,
                      label="Raise Alarm if astigmatism >",
                      help="Raise alarm if astigmatism (defocusU-defocusV)is greater than given "
                           "value")



        form.addSection('System Monitor')
        form.addParam('cpuAlert', params.FloatParam, default=101,
                      label="Raise Alarm if CPU > XX%",
                      help="Raise alarm if memory allocated is greater "
                           "than given percentage")

        form.addParam('memAlert', params.FloatParam, default=101,
                      label="Raise Alarm if Memory > XX%",
                      help="Raise alarm if cpu allocated is greater "
                           "than given percentage")
        form.addParam('swapAlert', params.FloatParam, default=101,
                      label="Raise Alarm if Swap > XX%",
                      help="Raise alarm if swap allocated is greater "
                           "than given percentage")

        group = form.addGroup('GPU')
        group.addParam('doGpu', params.BooleanParam, default=False,
                       label="Check GPU",
                       help="Set to true if you want to monitor the GPU")
        group.addParam('gpusToUse', params.StringParam, default='0',
                       label='Which GPUs to use:', condition='doGpu',
                       help='Provide a list of GPUs '
                            '(e.g. "0 1 2 3"). Default is to monitor GPU 0 only')
        group = form.addGroup('NETWORK')
        group.addParam('doNetwork', params.BooleanParam, default=False,
                       label="Check Network",
                       help="Set to true if you want to monitor the Network")
        group.addParam('netInterfaces', params.EnumParam,
                       choices=MonitorSystem.getNifsNameList(),
                       default=1,  # usually 0 is the loopback
                       label="Interface", condition='doNetwork',
                       help="Name of the network interface to be checked")

        group = form.addGroup('Disk')
        group.addParam('doDiskIO', params.BooleanParam, default=False,
                       label="Check Disk IO",
                       help="Set to true if you want to monitor the Disk "
                            "Acces")

        form.addSection('Mail settings')
        ProtMonitor._sendMailParams(self, form)

        form.addSection('HTML Report')

        form.addParam("doInflux", params.BooleanParam,
                      label="use grafana/influx",
                      default=False,
                      help="Use grafana+influx vs apache for reports")
        form.addParam('publishCmd', params.StringParam, default='',
                      label="Publish command",
                      help="Specify a command to publish the template. "
                           "You can use the special token %(REPORT_FOLDER)s "
                           "that will be replaced with the report folder. "
                           "For example: \n"
                           "rsync -avL %(REPORT_FOLDER)s "
                           "scipion@webserver:public_html/")

    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # --------------------------- STEPS functions ----------------------------
    def monitorStep(self):
        # create monitors
        movieGainMonitor = self.createMovieGainMonitor()
        ctfMonitor = self.createCtfMonitor()
        sysMonitor = self.createSystemMonitor()
        reportHtml = self.createHtmlReport(ctfMonitor, sysMonitor,
                                           movieGainMonitor)

        monitor = Monitor(workingDir=self.workingDir.get(),
                          samplingInterval=self.samplingInterval.get(),
                          monitorTime=self.monitorTime.get())

        def initAll():
            if ctfMonitor is not None:
                ctfMonitor.initLoop()
            if movieGainMonitor is not None:
                movieGainMonitor.initLoop()
            sysMonitor.initLoop()

        def stepAll():
            finished = False
            try:
                if ctfMonitor is not None:
                    # Call ctf monitor step
                    ctfMonitor.step()

                if movieGainMonitor is not None:
                    # Call movie gain step
                    movieGainMonitor.step()

                # sysmonitor watches all input protocols so
                # when sysmonitor done all protocols done
                sysMonitorFinished = sysMonitor.step()
                htmlFinished = reportHtml.generate(finished)
                if sysMonitorFinished and htmlFinished:
                    finished = True
                    reportHtml.generate(finished)

            except Exception as ex:
                print("An error happened:")
                import traceback
                traceback.print_exc()

            return finished

        monitor.initLoop = initAll
        monitor.step = stepAll

        monitor.loop()

    def createReportDir(self):
        self.reportDir = os.path.abspath(self._getExtraPath(self.getProject().getShortName()))
        self.reportPath = os.path.join(self.reportDir, 'index.html')
        # create report dir
        pwutils.makePath(self.reportDir)

        pathRepoSummary = self._getPath("pathRepo.txt")
        pathRepoSummary = open(pathRepoSummary, "w")
        pathRepoSummary.write("HTML path to summary: " + self.reportPath)
        pathRepoSummary.close()


    def _getAlignProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtAlignMovies):
                return prot
        return None

    def _getCtfProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtCTFMicrographs):
                return prot
        return None

    def _getMovieGainProtocol(self):
        XmippProtMovieGain = Domain.importFromPlugin('xmipp3.protocols',
                                                     'XmippProtMovieGain')

        if XmippProtMovieGain is None:
            return None

        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if prot.getClassName() == XmippProtMovieGain.__name__:
                return prot
        return None

    def createMovieGainMonitor(self):
        movieGainProt = self._getMovieGainProtocol()

        if movieGainProt is None:
            return None

        movieGainProt.setProject(self.getProject())

        movieGainMonitor = MonitorMovieGain(
                movieGainProt,
                influx=self.doInflux,
                workingDir=self.workingDir.get(),
                samplingInterval=self.samplingInterval.get(),
                monitorTime=self.monitorTime.get(),
                email=self.createEmailNotifier(),
                stdout=True,
                stddevValue=self.stddevValue.get(),
                ratio1Value=self.ratio1Value.get(),
                ratio2Value=self.ratio2Value.get())
        return movieGainMonitor

    def createCtfMonitor(self):
        ctfProt = self._getCtfProtocol()

        if ctfProt is None:
            return None

        ctfProt.setProject(self.getProject())

        ctfMonitor = MonitorCTF(ctfProt,
                                influx=self.doInflux,
                                workingDir=self.workingDir.get(),
                                samplingInterval=self.samplingInterval.get(),
                                monitorTime=self.monitorTime.get(),
                                email=self.createEmailNotifier(),
                                stdout=True,
                                minDefocus=self.minDefocus.get(),
                                maxDefocus=self.maxDefocus.get(),
                                astigmatism=self.astigmatism.get())
        return ctfMonitor

    def createSystemMonitor(self):
        protocols = self.getInputProtocols()

        sysMon = MonitorSystem(protocols,
                               influx=self.doInflux,
                               workingDir=self.workingDir.get(),
                               samplingInterval=self.samplingInterval.get(),
                               monitorTime=self.monitorTime.get(),
                               email=self.createEmailNotifier(),
                               stdout=True,
                               cpuAlert=self.cpuAlert.get(),
                               memAlert=self.memAlert.get(),
                               swapAlert=self.swapAlert.get(),
                               doGpu=self.doGpu.get(),
                               gpusToUse=self.gpusToUse.get(),
                               doNetwork=self.doNetwork.get(),
                               doDiskIO=self.doDiskIO.get(),
                               nif=MonitorSystem.getNifsNameList()[
                                   self.netInterfaces.get()])

        return sysMon

    def getReportPath(self):
        return self.reportPath

    def createHtmlReport(self, ctfMonitor=None, sysMonitor=None,
                         movieGainMonitor=None):
        ctfMonitor = ctfMonitor or self.createCtfMonitor()
        sysMonitor = sysMonitor or self.createSystemMonitor()
        movieGainMonitor = movieGainMonitor or self.createMovieGainMonitor()
        self.createReportDir()
        if self.doInflux:
            htmlReport = ReportInflux(self, ctfMonitor, sysMonitor, movieGainMonitor,
                                    self.publishCmd.get(),
                                    refreshSecs=self.samplingInterval.get())
        else:
            htmlReport = ReportHtml(self, ctfMonitor, sysMonitor, movieGainMonitor,
                                self.publishCmd.get(),
                                refreshSecs=self.samplingInterval.get())
            htmlReport.setUp()

        return htmlReport
    def _summary(self):
        summary = []
        pathRepoSummary = self._getPath("pathRepo.txt")
        if not os.path.exists(pathRepoSummary):
            summary.append("No summary file yet.")
        else:
            pathRepoSummary = open(pathRepoSummary, "r")
            for line in pathRepoSummary.readlines():
                summary.append(line.rstrip())
            pathRepoSummary.close()
        return summary


    def validate(self):
        errors = []
        if self.publishCmd.get() != '':
            self.reportDir = os.path.abspath(
                self._getExtraPath(self.getProject().getShortName()))
            pwutils.makePath(self.reportDir)

            cmd = str(self.publishCmd) % {'REPORT_FOLDER': self.reportDir}
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            output, err = p.communicate()
            if err.decode("utf-8") != '':
                errors.append('The publish command {} is wrong, please check it{}'.format(cmd, err.decode("utf-8")))


        return errors
