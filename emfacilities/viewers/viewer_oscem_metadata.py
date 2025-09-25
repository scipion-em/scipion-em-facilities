import subprocess
import time
from os.path import abspath
import pyworkflow.viewer as pwviewer
from pwem import Plugin
from emfacilities import STRM_ENV_NAME
from emfacilities.protocols import ProtOSCEM

class OscemViewer(pwviewer.Viewer):
    _label = 'OSC-EM Viewer'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [ProtOSCEM]

    def _visualize(self, obj, **kwargs):
        view = OSCEMView(obj)
        view._tkParent = self.getTkRoot()
        return [view]

class OSCEMView(pwviewer.View):
    def __init__(self, protocol):
        self.protocol = protocol
        pyFile = abspath(self.generateReport())
        cmd = Plugin.getCondaActivationCmd()
        cmd += f"conda activate {STRM_ENV_NAME} && "
        cmd += f"streamlit run {pyFile}"

        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            process.terminate()

    def generateReport(self):
        code = f"""
import streamlit as st
import yaml 
from pathlib import Path
from PIL import Image
import base64
from io import BytesIO

Image.MAX_IMAGE_PIXELS = 100_000_000
{genDataStruct()}
# --- AUXILIAR FUNCTIONS ---
def image_to_base64(image_path):
    img = Image.open(image_path)
    buffered = BytesIO()
    img.save(buffered, format="JPEG")
    return base64.b64encode(buffered.getvalue()).decode()

def get_prefix(level, is_last):
    if level == 0:
        return ""
    elif is_last:
        return "    " * (level - 1) + "└── "
    else:
        return "    " * (level - 1) + "├── "

def format_key(key: str) -> str:
    return " ".join([w.capitalize() for w in key.split("_")])

def display_node(key, value, level=0, is_last=True):

    # --- TOOLTIPS ---
    tooltips = {{
        'astigmatism': 'Astigmatism was calculated using the defocus ratio method.',
        'images_classes_3d': 'Images show central section.',
        'defocus_mic_examples': 'Micrographs are shown in increasing order of defocus',
        'particles_mic_examples': 'Micrographs are shown in decreasing order of particle number'
    }}
    tooltip = tooltips.get(key.lower(), "")

    # --- Personalized titles ---
    if key.lower() == "classes2d":
        formatted_key = "Classes 2D"
    elif key.lower() == "classes3d":
        formatted_key = "Classes 3D"
    elif key.lower() == "ctfs":
        formatted_key = "CTFs"
    else:
        formatted_key = format_key(key)

    # --- Level 0 (Processing) ---
    if level == 0:
        st.markdown(
            f"<div style='background:#dce6f0; padding:16px; margin-bottom:20px; "
            f"border-radius:12px; box-shadow:3px 3px 8px rgba(0,0,0,0.15);'>"
            f"<div style='font-size:24px; font-weight:700; color:#2c3e50; margin-bottom:12px'>{{formatted_key}}:</div>",
            unsafe_allow_html=True
        )

        if isinstance(value, dict):
            for i, (subkey, subval) in enumerate(value.items()):
                display_node(subkey, subval, level + 1, i == len(value) - 1)

        st.markdown("</div>", unsafe_allow_html=True)
        return

    # --- Level 1 (movies, CTFs, etc) ---
    if level == 1:
        st.markdown(
            f"<div style='background:#f5f7fa; padding:12px; margin-bottom:16px; "
            f"border-radius:10px; box-shadow:2px 2px 6px rgba(0,0,0,0.1);'>"
            f"<div style='font-size:20px; font-weight:600; color:#2c3e50; margin-bottom:8px'>{{formatted_key}}</div>",
            unsafe_allow_html=True
        )

        if isinstance(value, dict):
            for i, (subkey, subval) in enumerate(value.items()):
                display_node(subkey, subval, level + 1, i == len(value) - 1)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                display_node(f"{{key}} [{{i}}]", item, level + 1, i == len(value) - 1)
        else:
            st.markdown(f"<div style='margin-left:1rem'>{{value}}</div>", unsafe_allow_html=True)

        st.markdown("</div>", unsafe_allow_html=True)
        return
        
    # --- Levels > 1 ---
    prefix = get_prefix(level, is_last)

    if isinstance(value, dict):
        st.markdown(
            f"<div style='font-family:monospace; white-space:pre; line-height:1.8; margin-bottom:6px;'>{{prefix}}<strong>{{formatted_key}}</strong></div>",
            unsafe_allow_html=True
        )
        for i, (subkey, subval) in enumerate(value.items()):
            display_node(subkey, subval, level + 1, i == len(value) - 1)

    elif isinstance(value, list):
        if all(isinstance(item, (int, float, str)) for item in value):
            value_text = ", ".join(map(str, value))
            st.markdown(
                f"<div style='font-family:monospace; white-space:pre; line-height:1.8; margin-bottom:6px;'>{{prefix}}<strong>{{formatted_key}}:</strong> {{value_text}}</div>",
                unsafe_allow_html=True
            )
        else:
            st.markdown(
                f"<div style='font-family:monospace; white-space:pre; line-height:1.8; margin-bottom:6px;'>{{prefix}}<strong>{{formatted_key}}</strong></div>",
                unsafe_allow_html=True
            )
            for i, item in enumerate(value):
                display_node(f"{{key}} [{{i}}]", item, level + 1, i == len(value) - 1)

    elif isinstance(value, str) and value.lower().endswith(".jpg"):
        img_path = Path('{self.protocol._getExtraPath()}') / value
        if img_path.exists():
            img_base64 = image_to_base64(img_path)
            img_html = f"<img src='data:image/jpeg;base64,{{img_base64}}' class='zoom-img' alt='{{value}}'>"
        else:
            img_html = f"<span style='color:red;'>Image not found: {{img_path}}</span>"

        st.markdown(
            f"<div style='font-family:monospace; white-space:pre; line-height:1.8; margin-bottom:6px;'>{{prefix}}<strong>{{formatted_key}}:</strong> {{img_html}}</div>",
            unsafe_allow_html=True
        )

    else:
        line = f"{{prefix}}<strong>{{formatted_key}}:</strong> {{value}}"
        if tooltip:
            line += f" 🔍 <span style='color:gray;font-size:14px'>{{tooltip}}</span>"
        st.markdown(
            f"<div style='font-family:monospace; white-space:pre; line-height:1.8; margin-bottom:6px;'>{{line}}</div>",
            unsafe_allow_html=True
        )

# --- MAIN ---
yaml_path = '{self.protocol._getExtraPath('Processing_metadata.yaml')}'
st.markdown(
    "<div style='font-size:36px; font-weight:700; color:#2c3e50; margin-bottom:24px;'>OSC-EM Metadata Viewer</div>",
    unsafe_allow_html=True
)
if Path(yaml_path).exists():
    with open(yaml_path) as f:
        yaml_data = yaml.safe_load(f)
    for i, (key, value) in enumerate(yaml_data.items()):
        display_node(key, value, level=0, is_last=(i == len(yaml_data) - 1))
# else:
#     st.error("YAML file not found.")
"""
        genReport = self.protocol._getExtraPath('genReport.py')
        with open(genReport, "w") as pyFile:
            pyFile.write(code)
        return genReport


def genDataStruct():
    stMarkdown = f"""
st.markdown(\"\"\"
    <style>
    @import url("https://fonts.googleapis.com/css2?family=Roboto:wght@400;600&display=swap");
    
    .stApp {{
        background-color: #e6ecf1;
        font-family: "Roboto", sans-serif;
        color: #333333;
    }}
    
    .block-container {{
        padding-left: 2rem;
        padding-right: 2rem;
    }}
    
    .tree-label {{
        font-family: monospace;
        white-space: pre;
        margin-bottom: 6px;
        line-height: 1.5;
        font-size: 15px;
        padding: 4px 8px;
        border-radius: 6px;
        box-shadow: 1px 1px 4px rgba(0,0,0,0.1);
        transition: transform 0.2s ease, box-shadow 0.2s ease;
        position: relative;   
        overflow: visible;
    }}
    .tree-label:hover {{
        transform: translateY(-2px);
        box-shadow: 2px 4px 8px rgba(0,0,0,0.15);
    }}
    
    .zoom-img {{
        width: 200px;
        transition: transform 0.3s ease, box-shadow 0.3s ease;
        border: 2px solid #ccc;
        border-radius: 10px;
        margin-top: 5px;
        box-shadow: 2px 2px 8px rgba(0,0,0,0.15);
        transform-origin: center center;
        position: relative;
        z-index: 10;   
    }}
    
    .zoom-img:hover {{
        transform: scale(2.5);
        box-shadow: 4px 6px 12px rgba(0,0,0,0.25);
        z-index: 999;
    }}
    
    h1 {{
        font-family: "Roboto", sans-serif;
        font-weight: 600;
        color: #2c3e50;
    }}
    </style>
    \"\"\", unsafe_allow_html=True)
"""
    return stMarkdown