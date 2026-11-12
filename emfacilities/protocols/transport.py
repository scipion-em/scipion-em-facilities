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

import paramiko
import os


class Connect:
    def __init__(self, host, port, username, password, keyfilepath, keyfiletype,
                 remote_path, projectName):
        """
        Establishes and manages secure remote file transfer sessions through SFTP, allowing projects and
        associated files to be stored, organized, and synchronized on external systems.

        AI Generated:

        Connect (Connect) — User Manual
            Overview

            The Connect class provides a secure communication layer for transferring files between a local
            environment and a remote storage system using the SFTP protocol. Its main objective is to support
            reliable exchange of project data while maintaining authenticated and encrypted communication with
            remote servers.

            In practical scientific and computational workflows, remote file transfer becomes essential when
            processing data across distributed infrastructures, high performance computing facilities, shared
            servers, or institutional repositories. This component simplifies the management of these remote
            interactions by handling authentication, connection establishment, directory preparation, and file
            upload operations within a unified interface.

            Authentication and Secure Access

            The connection mechanism supports multiple authentication strategies in order to adapt to different
            institutional security requirements. Users may authenticate either with standard credentials or
            with private key authentication using supported cryptographic key formats.

            Key-based authentication is especially important in automated workflows and cluster environments
            where secure unattended access is required. This approach improves security while also enabling
            reproducible and automated data management pipelines without repeated manual intervention.

            Remote Workspace Organization

            Once connected, the class prepares a dedicated project directory within the remote storage
            location. This ensures that transferred data are grouped consistently according to project identity
            and remain organized across repeated executions or collaborative workflows.

            Such organization is particularly valuable in scientific environments where multiple experiments,
            datasets, or processing runs may coexist on the same remote infrastructure. Maintaining isolated
            project directories reduces confusion and improves traceability of generated results.

            File Transfer Operations

            The primary operational purpose of the class is the transfer of local files into the remote
            workspace. This capability supports workflows where processed outputs, reports, metadata, or
            intermediate results must be archived remotely or shared with collaborators and computational
            facilities.

            In biological and imaging workflows, remote transfers are frequently used for moving microscopy
            datasets, reconstruction results, or processing outputs between acquisition systems and analysis
            platforms. Efficient remote synchronization becomes especially important when working with large
            cryo-EM datasets or distributed computational infrastructures.

            Reliability and Resource Management

            The connection management process includes safeguards to ensure that communication resources are
            properly released when operations finish or when unexpected situations occur. This minimizes the
            risk of leaving incomplete sessions active on remote systems and improves the stability of
            long-running automated workflows.

            By isolating transfer operations within a dedicated communication component, the class also helps
            maintain cleaner workflow organization and simplifies integration into larger processing systems.

            Final Perspective

            The Connect class serves as a secure bridge between local workflows and remote computational or
            storage environments. By combining authenticated communication, remote workspace management, and
            reliable file transfer capabilities, it enables efficient movement of scientific data across
            distributed infrastructures while supporting reproducibility, collaboration, and organized project
            management.
        """
        self.sftp = None
        key = None
        self.transport = None
        try:
            if keyfilepath is not None:
                # Get private key used to authenticate user.
                if keyfiletype == 'DSA':
                    # The private key is a DSA type key.
                    key = paramiko.DSSKey.from_private_key_file(keyfilepath)
                else:
                    # The private key is a RSA type key.
                    key = paramiko.RSAKey.from_private_key_file(keyfilepath)

            # Create Transport object using supplied method of authentication

            self.transport = paramiko.Transport((host, port))
            self.transport.connect(None, username, password, key)

            self.sftp = paramiko.SFTPClient.from_transport(self.transport)

        except Exception as e:
            print('An error occurred creating SFTP client: %s: %s' % (e.__class__, e))
            if self.sftp is not None:
                self.sftp.close()
            if self.transport is not None:
                self.transport.close()
        self.remote_path = remote_path
        directory = os.path.join(self.remote_path, projectName)
        try:
            self.sftp.chdir(directory)  # Test if directory exists
        except IOError:
            self.sftp.mkdir(directory)  # Create directory

    def put(self,listLocalPaths, listRemotePaths):
        try:
            for local, remote in zip(listLocalPaths, listRemotePaths):
                remote = os.path.join(self.remote_path, remote)
                print(local, "-->", remote)
                self.sftp.put(local, remote, confirm=True)
        except IOError:
            pass
        except OSError:
            pass
        except Exception as e:
            print(str(e))
        return 0

    def close(self):
        if self.sftp is not None:
            self.sftp.close()
        if self.transport is not None:
            self.transport.close()
