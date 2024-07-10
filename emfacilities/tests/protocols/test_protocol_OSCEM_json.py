import json
from os.path import join, exists, abspath

from pwem.protocols import ProtImportMovies
from pyworkflow.tests import BaseTest, tests, DataSet
from pyworkflow.utils import magentaStr
from xmipp3.protocols import XmippProtMovieGain
from ...protocols import ProtOSCEM
from ...protocols.protocol_OSCEM_json import INPUT_MOVIES


class TestOscemJson(BaseTest):
    """
    """
    sampling_rate = 0.495

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('OSCEM_jsons')
        cls.protimportmovies, cls.importedmovies = cls.runImportMovies()

    @classmethod
    def runImportMovies(cls):
        prot = cls.newProtocol(ProtImportMovies,
                               filesPath=cls.dataset.getFile('movies_dir'),
                               filesPattern='*.tiff',
                               samplingRate=cls.sampling_rate)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMovies', None)
        return prot, output

    def runMovieGain(cls):
        prot = cls.newProtocol(XmippProtMovieGain,
                               inputMovies=cls.importedmovies)

        cls.launchProtocol(prot)
        # output=getattr(prot, 'SetOffMovies', None)
        # cls.assertIsNotNone(output, 'No set of movies was generated')
        return prot

    def test_only_import(self):
        print(magentaStr(f"\n==> Running import movies test:"))
        test_data = {
            "Import_movies": {
                "Microscope_voltage": 300.0,
                "Spherical_aberration": 2.7,
                "Amplitud_contrast": 0.1,
                "Sampling_rate": 0.495,
                "Pixel_size": 2.475,
                "Number_movies": 5
            }
        }


        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies)
        self.launchProtocol(prot)
        # recuperar resultados
        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, current_test_dict in test_data.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in current_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                self.assertEqual(current_test_value, current_value)
