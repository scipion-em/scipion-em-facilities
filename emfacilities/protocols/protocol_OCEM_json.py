import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol
from pyworkflow.object import String


class ProtocolOSCEM(EMProtocol):
    """ !!!!!!!!!!!!!!!!!!!!!!
    """
    _label = 'OSCEMJSon'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputProtocols', params.PointerParam,
                      label="Input", important=True,
                      pointerClass='EMProtocol',  #'CTFModel'],
                      help="input")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.generateJson)

    # -------------------------- STEPS functions ------------------------------
    def generateJson(self):
        protocol = self.inputProtocols.get()
        protocol.alignFrame0.get()

        protocol.get('alignFrame0', String()).get()


    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []

    def getInputProtocols(self):  # multipointer por el for
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)
        return protocols
