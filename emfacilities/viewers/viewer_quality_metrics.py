from emfacilities.protocols import ProtQualityMetrics
from emfacilities.protocols.protocol_quality_metrics import MONITOR_FN

from pwem.viewers import ObjectView, EmProtocolViewer
from pwem.viewers.showj import *

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
import pyworkflow.protocol.params as params

import matplotlib.pyplot as plt
import os
import subprocess
import sys
import tempfile
import webbrowser
import time
import socket

class QualityMetricsViewer(EmProtocolViewer):
    """ This viewer is intended to visualize the selection made by
        the Miffi - categorize micrographs protocol.
    """
    _label = 'viewer Quality Metrics'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtQualityMetrics]


    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('visualizeQualityMonitor', params.LabelParam,
                      label="Visualize Live Quality Monitor",
                      help="Visualize a Streamlit Quality Monitor.")

    def _getVisualizeDict(self):
        return {
                 'visualizeQualityMonitor': self._visualizeQualityMonitor
                }

    def _visualizeQualityMonitor(self, e=None):
        self._showStreamLitMonitor(os.path.join(self.protocol._getExtraPath(), MONITOR_FN))

    def _showStreamLitMonitor(self, monitorPath):
        print(monitorPath)
        if os.path.exists(monitorPath):
            env = os.environ.copy()
            env["STREAMLIT_MONITOR_PATH"] = monitorPath

            # Aquí llamas directamente al script guardado en tu plugin
            script_path = os.path.join(os.path.dirname(__file__), 'quality_dashboard.py')
            # Launch Streamlit app on a custom port (e.g., 8503 to avoid conflicts)
            port = find_free_port()

            subprocess.Popen(
                [sys.executable, "-m", "streamlit", "run", script_path,
                 "--server.headless", "true",
                 f"--server.port={port}",
                 "--server.address=0.0.0.0"],  # This makes it externally accessible
                env=env
            )
            # Wait a bit for the server to start
            time.sleep(2)
            # Open the browser automatically
            webbrowser.open(f"http://localhost:{port}")

            hostname = socket.gethostname()
            ip_address = socket.gethostbyname(hostname)
            print(f"Streamlit available at: http://{ip_address}:{port}")
        else:
            self.error("Non CSV metadata files was found.")


def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))  # Let OS pick an available port
        return s.getsockname()[1]