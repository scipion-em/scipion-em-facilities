# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     R. Marabini (roberto@cnb.csic.es)
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
from os.path import basename
from datetime import datetime, timedelta
from configparser import ConfigParser, NoOptionError, NoSectionError
import urllib3
import base64
import time
from .transport import Connect

import pyworkflow.utils as pwutils
from pwem.emlib.image import ImageHandler

from emfacilities.protocols import SummaryProvider
from influxdb import InfluxDBClient


# --------------------- CONSTANTS -----------------------------------
# These constants are the keys used in the ctfMonitor function
# getData() to store paths to micrographs, psd files and shift plots
from .. import Plugin

MIC_PATH = 'imgMicPath'
PSD_PATH = 'imgPsdPath'
SHIFT_PATH = 'imgShiftPath'
# These constants are the name of the folders where thumbnails
# for the html report will be stored. They are also the keys to
# used in the execution.summary.template.html to read data (where
# they will need to be changed if they're changed here)
CONFILE = 'monitor.conf'

class ReportInflux:
    """ Create a report using influxdb with a summary of the processing.
    The report will be updated with a given frequency.
    """
    def __init__(self, protocol, ctfMonitor,
                 sysMonitor, movieGainMonitor,
                 **kwargs):
        """
        :param protocol:
        :param ctfMonitor:
        :param sysMonitor:
        :param movieGainMonitor:
        :param kwargs:     refreshSecs -> update each refreshSecs seconds
        """

        # The CTF protocol to monitor
        self.protocol = protocol
        self.micThumbSymlinks = False
        self.reportPath = protocol.reportPath
        self.reportDir = protocol.reportDir
        self.provider = SummaryProvider(protocol)
        self.ctfMonitor = ctfMonitor
        self.sysMonitor = sysMonitor
        self.movieGainMonitor = movieGainMonitor
        self.lastThumbIndex = 0
        self.thumbsReady = 0
        self.refreshSecs = kwargs.get('refreshSecs', 60)
        if self.refreshSecs < 10:
            self.refreshSecs = 10
        self.ih = ImageHandler()

        # Create a config file to store data we are going to
        # need in other iteration. So far is only used by influx
        # but it certainly  be used by hte standard html
        # for example to decide which files need to be tranfered
        self.confFileName = self.protocol._getTmpPath(CONFILE)
        if not os.path.isfile(self.confFileName):
            self.confParser = self.createStatusConfigFile(self.confFileName,
                                                          self.protocol.getProject().getShortName())
        else:
            self.confParser = ConfigParser()
            self.confParser.read(self.confFileName)

        # Clean project name from special characters
        self.projectName = slugify(self.protocol.getProject().getShortName())

        # Parse secrets file and get authentication info
        hostinflux, passwordInflux, port, ssl, usernameInflux, verify_ssl = self.parseSecrets()

        # Connect to influx
        self.connectToInflux(hostinflux, passwordInflux, port, ssl, usernameInflux, verify_ssl)
    
    @staticmethod
    def createStatusConfigFile(configFileName, projectName):
        
        # Create the configuration file as it doesn't exist yet
        confParser = ConfigParser()
        confParser.add_section("project")
        confParser.set("project", "projectName",
                            projectName)
        confParser.set("project", "properties", "1")
        confParser.set("project", "acquisition", "1")
        confParser.set("project", "summary", "1")
        confParser.add_section('ctf')
        confParser.set("ctf", "lastId", "0")
        confParser.add_section("gain")
        confParser.set("gain", "lastId", "0")
        confParser.add_section("system")
        confParser.set("system", "lastId", "0")

        with open(configFileName, 'w') as confFile:
            confParser.write(confFile)

        return confParser

    def connectToInflux(self, hostinflux, passwordInflux, port, ssl, usernameInflux, verify_ssl):
        try:
            # since I am using a self generated certificate
            # the certificate can not be verified against a
            # certificate authority. So I disable the warnings
            urllib3.disable_warnings()

            #  open conection to influx
            self.client = InfluxDBClient(host=hostinflux,
                                         port=port,
                                         username=deCrypt(usernameInflux),
                                         password=deCrypt(passwordInflux),
                                         ssl=ssl,
                                         verify_ssl=verify_ssl)
            # we can only create the database here if scipion_writer
            # has admin privileges
            # If you attempt to create a database that already exists,
            # InfluxDB does nothing and does not return an error.
            # self.client.create_database("scipion")
            self.client.switch_database(self.dataBaseName)

            # delete meassurement
            # project names may contain forbiden character
            # IF this is a problem we will need to slugify the projName
            # self.client.drop_measurement(self.projectName)
            self.client.delete_series(measurement=self.projectName)

            # replication -> number of copies of the DB stored in the cluster
            # 12w -> delete data after 12 weeks
            self.client.create_retention_policy("ret_pol", "12w",
                                                replication=1,
                                                default=True)

        except Exception as e:
            print("Error:", e)

    def parseSecrets(self):
        # Read secrets file to connect to influx and send file using paramiko
        # open connection to database
        # I put this import here so users with no database
        # can run traditional html report
        try:
        
            secretsParser = ConfigParser()
            from emfacilities.constants import (SECRETSFILE,
                                                EMFACILITIES_HOME_VARNAME)
            _path = Plugin.getHome(EMFACILITIES_HOME_VARNAME)
            secretsfile = os.path.join(_path, SECRETSFILE)
            secretsParser.read(secretsfile)
            # influx data
            dataBase = secretsParser.get('influx', 'dataBase')
            passwordInflux = secretsParser.get('influx', 'passwordInflux')
            usernameInflux = secretsParser.get('influx', 'usernameInflux')
            hostinflux = secretsParser.get('influx', 'hostinflux')
            port = int(secretsParser.get('influx', 'port'))
            ssl = secretsParser.get('influx', 'ssl')
            ssl = ssl.lower() in ['true', 1, 't', 'y']
            verify_ssl = secretsParser.get('influx', 'verify_ssl')
            verify_ssl = verify_ssl.lower() in ['true', 1, 't', 'y']
            self.timeDelta = int(secretsParser.get('influx', 'TimeDelta'))
            self.dataBaseName = dataBase
            # paramiko data
            self.usernameParamiko = secretsParser.get('paramiko', 'usernameParamiko')
            self.passwordParamiko = secretsParser.get('paramiko', 'passwordParamiko')
            self.keyfilepath = secretsParser.get('paramiko', 'keyfilepath')
            self.keyfiletype = secretsParser.get('paramiko', 'keyfiletype')
            self.remote_path = secretsParser.get('paramiko', 'remote_path')
            self.hostparamiko = secretsParser.get('paramiko', 'hostparamiko')

        except NoSectionError as e:
            print("Can not process secrets file %s: %s" % (secretsfile, str(e)))
            raise e

        return hostinflux, passwordInflux, port, ssl, usernameInflux, verify_ssl

    #def __del__(self):
    #    if self.influxsb:
    #        self.client.close()

    def generate(self, finished):

        project = self.protocol.getProject()

        # Project Properties Section
        # Do not delete this variables. We are using them
        # in an eval command
        self.projectName = project.getShortName()
        startTime = pwutils.dateStr(project.getCreationTime(), secs=True),
        tnow = datetime.now()
        _now = project.getCreationTime()
        dateStr = pwutils.prettyTime(dt=tnow, secs=True),
        projectDuration = pwutils.prettyDelta(tnow - project.getCreationTime()),
        projectStatus = "FINISHED" if finished else "RUNNING",
        scipionVersion = os.environ['SCIPION_VERSION'],

        # create structure with information related with the properties section
        # the primmary key is time, projectionName, section
        # if we insert two points with these three values idential
        # the second updates the first
        pointsDict = {}  # dictionary for data points
        pointsDict['measurement'] = self.projectName
        tags = {}
        tags['section'] = 'properties'
        pointsDict['tags'] = tags
        fields = {}
        firstTime = self.confParser.getint("project", "properties")

        fieldKeys  = {'dateStr': 4, 'projectDuration': 3, 'projectStatus': 2}
        fieldNames = {'dateStr': "<b>Last Update</b>",
                      'projectDuration': '<b>Duration</b>',
                      'projectStatus': '<b>Status</b>'}
        if firstTime:
            fieldKeys['startTime'] = 5
            fieldKeys['scipionVersion'] = 1
            fieldNames['startTime'] = '<b>Start Time</b>'
            fieldNames['scipionVersion'] = '<b>Scipion Version</b>'

        for metric, delta in fieldKeys.items():
            localNow = _now + timedelta(seconds=delta)
            pointsDict['time'] = localNow.strftime('%Y-%m-%dT%H:%M:%SZ')
            fields['metric'] = fieldNames[metric]
            # str is need because all values must have the same type
            fields['valueStr'] = str(eval(metric)[0])
            pointsDict['fields'] = fields
            self.client.write_points([pointsDict])
        self.confParser.set("project", "properties", "0")
        with open(self.confFileName, 'w') as confFile:
            self.confParser.write(confFile)

        # acquisition section
        self.provider.refreshObjects()
        pointsDict = {}  # dictionary for data points
        pointsDict['measurement'] = self.projectName
        tags = {}
        tags['section'] = 'acquisition'
        pointsDict['tags'] = tags

        fields = {}
        firstTime = self.confParser.getint("project", "acquisition")

        if firstTime:
            delta = 0
            for metricName, value in self.provider.acquisition:
                localNow = _now + timedelta(seconds=delta)
                pointsDict['time'] = localNow.strftime('%Y-%m-%dT%H:%M:%SZ')
                fields['metric'] = "<b>" + metricName + "</b>"
                fields['valueNum'] = float(value)
                pointsDict['fields'] = fields
                self.client.write_points([pointsDict])
                delta = delta -1
                # update first time only if some date has been read.
                # do not place this upside the loop
                self.confParser.set("project", "acquisition", "0")
                with open(self.confFileName, 'w') as confFile:
                    self.confParser.write(confFile)

        # send summary section
        pointsDict = {}  # dictionary for data points
        pointsDict['measurement'] = self.projectName
        tags = {}
        tags['section'] = 'summary'
        pointsDict['tags'] = tags
        localNow = _now
        for obj in self.provider.getObjects():
            fields = {}
            localNow = localNow + timedelta(seconds=-1)
            pointsDict['time'] = localNow.strftime('%Y-%m-%dT%H:%M:%SZ')
            # If it's a protocol
            isProtocol = True if obj.name else False

            if isProtocol:
                fields['protocolName'] = '<b>' + obj.name + '</b>'
                fields['output'] = ""
                fields['size'] = ""
            else:
                fields['protocolName'] = ""
                fields['output'] = obj.output
                fields['size'] = str(obj.outSize)
            pointsDict['fields'] = fields
            self.client.write_points([pointsDict])

        # Ctf monitor chart data
        last_id = self.confParser.getint("ctf", "lastId")
        listDictionaryCTF = {} if self.ctfMonitor is None else \
            self.ctfMonitor.getData(lastId=last_id)
        if listDictionaryCTF:
            pointsDict = {}  # dictionary for data points
            pointsDict['measurement'] = self.projectName
            tags = {}
            tags['section'] = 'ctf'
            for ctf in listDictionaryCTF:
                fields = {}
                for key in ctf:
                    if key == 'timestamp':
                        pointsDict['time'] = ctf['timestamp']
                    # id must be a tag since tw CTF may have the same timestamp
                    # and time + tags must be unique
                    elif key == 'id':
                        tags['id'] = ctf['id']
                    # elif key == 'ctfID':
                    #     fields['ctfID'] = "<b>" + str(ctf[key]) + "</b>"
                    elif key == 'shiftPlotPath':
                        temp = os.path.join('/public/img/scipionbox', self.projectName, basename(ctf[key]))
                        popUpStr = """<a href = "%s" target = "_blank"> <img src="%s"  alt="%s" width="128px" height="128px"> </a>""" % \
                                   (temp, temp, basename(ctf[key]))
                        fields[key] = popUpStr
                        fields[key + 'Local'] = ctf[key]
                    elif key == 'psdPath':
                        # convert to pǹg
                        baseName = basename(ctf[key])
                        baseNamePng = pwutils.replaceBaseExt(ctf[key], "png")
                        temp = os.path.join('/public/img/scipionbox', self.projectName, baseName)
                        temp = pwutils.replaceExt(temp, "png")
                        popUpStr = """<a href = "%s" target = "_blank"> <img src="%s"  alt="%s" width="128px" height="128px"> </a>""" % \
                                   (temp, temp, baseNamePng)
                        fields[key] = popUpStr
                        fields[key + 'Local'] = ctf[key]
                        fields[key + 'LocalPng'] = pwutils.replaceExt(ctf[key], "png")
                    elif key == 'micPath':
                        # convert to pǹg
                        baseName = basename(ctf[key])
                        baseNamePng = pwutils.replaceBaseExt(ctf[key], "png")
                        temp = os.path.join('/public/img/scipionbox', self.projectName, baseName)
                        temp = pwutils.replaceExt(temp, "png")
                        popUpStr = """<a href = "%s" target = "_blank"> <img src="%s"  alt="%s" width="128px" height="128px"> </a>""" % \
                                   (temp, temp, baseNamePng)
                        fields[key] = popUpStr
                        fields[key + 'Local'] = ctf[key]
                        fields[key + 'LocalPng'] = pwutils.replaceExt(ctf[key], "png")
                    else:
                        fields[key] = ctf[key]
                # while be use to control image creation
                fields["transferImage"] = False
                pointsDict['fields'] = fields
                pointsDict['tags'] = tags
                self.client.write_points([pointsDict])
                last_id += 1
            self.confParser.set("ctf", "lastId", str(last_id))
            with open(self.confFileName, 'w') as confFile:
                self.confParser.write(confFile)

        # GAIN Section
        last_id = self.confParser.getint("gain", "lastId")
        listDictionaryGain = {} if self.movieGainMonitor is None else \
            self.movieGainMonitor.getData(lastId=last_id)
        if listDictionaryGain:
            pointsDict = {}  # dictionary for data points
            pointsDict['measurement'] = self.projectName
            tags = {}
            tags['section'] = 'gain'
            pointsDict['tags'] = tags
            counter = 1
            tnow = datetime.now()
            for gain in listDictionaryGain:
                fields = {}
                for key in gain:
                    fields[key] = gain[key]
                pointsDict['fields'] = fields
                # gain has no time information so I just
                # put current time
                localNow = tnow + timedelta(seconds=counter)
                pointsDict['time'] = localNow # .strftime('%Y-%m-%dT%H:%M:%SZ')
                self.client.write_points([pointsDict])
                last_id += 1
            self.confParser.set("gain", "lastId", str(last_id))
            with open(self.confFileName, 'w') as confFile:
                self.confParser.write(confFile)

        # SYSTEM data
        last_id = self.confParser.getint("system", "lastId")
        listDictionarySystem = {} if self.sysMonitor is None else \
            self.sysMonitor.getData(lastId=last_id)
        if listDictionarySystem:
            pointsDict = {}  # dictionary for data points
            pointsDict['measurement'] = self.projectName
            tags = {}
            tags['section'] = 'system'
            pointsDict['tags'] = tags
            for system in listDictionarySystem:
                fields = {}
                for key in system:
                    if key=='timestamp':
                        pointsDict['time'] = system['timestamp']
                    elif key=='id':
                        continue
                    else:
                        fields[key] = system[key]
                pointsDict['fields'] = fields
                self.client.write_points([pointsDict])
                last_id += 1
            self.confParser.set("system", "lastId", str(last_id))
            with open(self.confFileName, 'w') as confFile:
                self.confParser.write(confFile)
        self.transferFiles()
        return last_id # reportFinished


    def transferFiles(self):
        # get images that need to  be transfered
        # in grupos of 10 images
        start_time = time.time()
        query = '''select * 
                   from "%s"  
                   where section = 'ctf'
                     and transferImage = false
                   order by time desc
                   limit 10'''% self.projectName

        result = self.client.query(query)
        if len(result) < 1:
            return True
 
        connect = Connect(host=self.hostparamiko,
                          port=22,
                          username= deCrypt(self.usernameParamiko),
                          password=None,
                          keyfilepath=deCrypt(self.keyfilepath),
                          keyfiletype=deCrypt(self.keyfiletype),
                          remote_path=self.remote_path,
                          projectName=self.projectName)

        while True:
            result = self.client.query(query)
            source = []
            target = []
            newPoints = []
            points = result.get_points()
            for point in points:
                tags = {}
                tags['section'] = 'ctf'
                pointsDict = {}
                source.append(point['shiftPlotPathLocal']) # shift plot
                source.append(point['psdPathLocalPng'])    # psd image
                source.append(point['micPathLocalPng'])    # micrograph

                # default value 512
                try:
                    X, Y, Z, N = self.ih.getDimensions(point['psdPathLocal'])
                except:
                    print("file %s does not exist"% point['psdPathLocal'])
                if X > 512:
                    scaleFactor = X // 512
                else:
                    scaleFactor = 1

                self.ih.computeThumbnail(point['psdPathLocal'], point['psdPathLocalPng'],
                                    scaleFactor=scaleFactor, flipOnY=True)

                X, Y, Z, N = self.ih.getDimensions(point['micPathLocal'])
                if X > 512:
                    scaleFactor = X // 512
                else:
                    scaleFactor = 1
                self.ih.computeThumbnail(point['micPathLocal'], point['micPathLocalPng'],
                                    scaleFactor=scaleFactor, flipOnY=True)

                # add files to transfer list
                target.append(os.path.join(self.projectName,
                                           basename(point['shiftPlotPathLocal'])))
                target.append(os.path.join(self.projectName,
                                           basename(point['psdPathLocalPng'])))
                target.append(os.path.join(self.projectName,
                                           basename(point['micPathLocalPng'])))
                #
                point['transferImage'] = True
                pointsDict['time'] = point['time']
                tags['id'] = point['id']
                del point['time']
                del point['section']
                del point['id']
                pointsDict['fields'] = point
                pointsDict['tags'] = tags
                pointsDict['measurement'] = self.projectName
                newPoints.append(pointsDict)
            try:
                connect.put(source, target)
            except Exception as e:
                print("Error tranfering files: ", e)
                return False # do not mark files are read
            # update influx
            try:
                self.client.write_points(points=newPoints)
            except Exception as e:
                print("Error updating database: ", e)
                return False # do not mark files are read

            elapsed_time = time.time() - start_time
            if elapsed_time > self.refreshSecs:
                break
            elif len(result) == 0:
                break
        connect.close()
        return True

def slugify(text):
    """
    Remove character that can not be used in databases
    as keys
    """
    non_url_safe = ['"', '#', '$', '%', '&', '+', #'-',
                    ',', '/', ':', ';', '=', '?',
                    '@', '[', '\\', ']', '^', '`',
                    '{', '|', '}', '~', "'"]

    non_safe = [c for c in text if c in non_url_safe]
    if non_safe:
        for c in non_safe:
            text = text.replace(c, '')
    # Strip leading, trailing and multiple whitespace, convert remaining whitespace to _

    text = u'_'.join(text.split())
    return text


def enCrypt(message):
    """Totally naive encryption routine that will not
    stop a hacker
    """
    message_bytes = message.encode('ascii')
    base64_bytes = base64.b64encode(message_bytes)
    return base64_bytes.decode('ascii')


def deCrypt(base64_message):
    """
    Decodes the string created by enCrypt
    """
    base64_bytes = base64_message.encode('ascii')
    message_bytes = base64.b64decode(base64_bytes)
    message = message_bytes.decode('ascii')
    return message
