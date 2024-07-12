import json
from os.path import join, exists, abspath

from pwem.protocols import ProtImportMovies, ProtAlignMovies
from pyworkflow.tests import BaseTest, tests, DataSet
from pyworkflow.utils import magentaStr
from xmipp3.protocols import XmippProtMovieGain, XmippProtFlexAlign, XmippProtMovieMaxShift, XmippProtCTFMicrographs
from ...protocols import ProtOSCEM
from ...protocols.protocol_OSCEM_json import INPUT_MOVIES


class TestOscemJson(BaseTest):
    """
    """
    sampling_rate = 0.495
    max_movie_shift = 20.0
    ctf_down_factor = 2.0
    high_res = 0.5
    test_data = {
        "Import_movies": {
            "Microscope_voltage": 300.0,
            "Spherical_aberration": 2.7,
            "Amplitud_contrast": 0.1,
            "Sampling_rate": 0.495,
            "Pixel_size": 2.475,
            "Gain_image": True,
            "Number_movies": 10
        },
        "Movie_alignment": {
            "Method": "XmippProtFlexAlign",
            "Binning_factor": 1.0,
            "Maximum_resolution": 30.0,
            "Frames_aligned": {
                "Frame0": 1,
                "FrameN": 0
            },
            "Ouput_avg_shift": 12.085274343980803,
            "Output_max_shift": 34.31086742604614
        },
        "Movie_maxshift": {
            "Discarded_movies": 5,
            "Max_frame_shift": 5.0,
            "Max_movie_shift": 20.0,
            "Rejection_type": "By frame or movie",
            "Ouput_avg_shift": 10.655840604909368,
            "Output_max_shift": 28.289413400687277
        },
        "CTF_estimation": {
            "Defocus": {
                "Output_max_defocus": 9238.84,
                "Output_min_defocus": 1509.145,
                "Ouput_avg_defocus": 5364.463750000001
            },
            "Resolution": {
                "Output_max_resolution": 3.062813,
                "Output_min_resolution": 2.45025,
                "Ouput_avg_resolution": 2.793877375
            }
        }
    }


    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('OSCEM_jsons')
        cls.protimportmovies, cls.importedmovies = cls.runImportMovies()
        cls.protmoviegain, cls.gainmovies = cls.runMovieGain()
        cls.protalignmovie, cls.alignedmovies = cls.runMovieAlign()
        cls.protmaxshift, cls.maxshiftmicro = cls.runMaxShift()
        # cls.protCTF, cls.CTFmicro = cls.runCTFestimation()

    @classmethod
    def runImportMovies(cls):

        prot = cls.newProtocol(ProtImportMovies,
                               filesPath=cls.dataset.getFile('movies_dir'),
                               filesPattern='*.tiff',
                               samplingRate=cls.sampling_rate,
                               gainFile=cls.dataset.getFile('gain_im'))

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMovies', None)
        return prot, output

    @classmethod
    def runMovieGain(cls):
        prot = cls.newProtocol(XmippProtMovieGain,
                               inputMovies=cls.importedmovies)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMovies', None)
        return prot, output

    @classmethod
    def runMovieAlign(cls):
        prot = cls.newProtocol(XmippProtFlexAlign,
                               inputMovies=cls.gainmovies)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMovies', None)
        return prot, output

    @classmethod
    def runMaxShift(cls):
        prot = cls.newProtocol(XmippProtMovieMaxShift,
                               inputMovies=cls.alignedmovies,
                               maxMovieShift=cls.max_movie_shift)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMicrographs', None)
        return prot, output

    @classmethod
    def runCTFestimation(cls):
        prot = cls.newProtocol(XmippProtCTFMicrographs,
                               inputMicrographs=cls.maxshiftmicro,
                               AutoDownsampling=False,
                               ctfDownFactor=cls.ctf_down_factor,
                               highRes=cls.high_res)

        cls.launchProtocol(prot)
        import os
        fname = "/home/lsanchez/debugs/test.txt"
        if os.path.exists(fname):
            os.remove(fname)
        fjj = open(fname, "a+")
        fjj.write('-------->onDebugMode PID {}'.format(os.getpid()))
        fjj.close()
        print('-------->onDebugMode PID {}'.format(os.getpid()))
        import time
        time.sleep(10)

        output = getattr(prot, 'outputMovies', None)
        return prot, output

    def test_only_import(self):
        print(magentaStr(f"\n==> Running import movies test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies)
        self.launchProtocol(prot)

        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, current_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in current_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                self.assertEqual(current_test_value, current_value)

    def test_import_and_movie_align(self):
        print(magentaStr(f"\n==> Running import movies and movie alignment test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Movie_alignment": self.test_data["Movie_alignment"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie)
        self.launchProtocol(prot)

        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, current_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in current_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Ouput_avg_shift" or key_in == "Output_max_shift":
                    # these values change decimals each time alignment is run
                    self.assertEqual(round(current_test_value), round(current_value))
                else:
                    self.assertEqual(current_test_value, current_value)

    def test_import_and_max_shift(self):
        print(magentaStr(f"\n==> Running import movies and movie max shift test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Movie_maxshift": self.test_data["Movie_maxshift"]
                            }

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                maxShift=self.protmaxshift)
        self.launchProtocol(prot)
        #
        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, current_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in current_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Ouput_avg_shift" or key_in == "Output_max_shift":
                    # these values change decimals each time alignment is run
                    self.assertEqual(round(current_test_value), round(current_value))
                else:
                    self.assertEqual(current_test_value, current_value)

    # def test_import_and_CTF(self):
    #     print(magentaStr(f"\n==> Running import movies and CTF estimation test:"))
    #     test_data_import = {"Import_movies": self.test_data["Import_movies"],
    #                         "CTF_estimation": self.test_data["CTF_estimation"]
    #                         }
    #
    #     prot = self.newProtocol(ProtOSCEM,
    #                             inputType=INPUT_MOVIES,
    #                             importMovies=self.protimportmovies,
    #                             CTF=self.protCTF)
    #     self.launchProtocol(prot)
    #     # recuperar resultados -> meter en funcion y meto el protocolo de entrada
    #     file_path = prot.getOutFile()
    #     self.assertTrue(exists(file_path))
    #
    #     with open(abspath(file_path), 'r') as json_file:
    #         load_json = json.load(json_file)
    #
    #     for key, current_test_dict in test_data_import.items():
    #         current_dict = load_json.get(key, None)
    #         self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
    #         for key_in, current_test_value in current_test_dict.items():
    #             current_value = current_dict.get(key_in, None)
    #             self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
    #             if key_in == "Ouput_avg_shift" or key_in == "Output_max_shift":
    #                 # these values change decimals each time alignment is run
    #                 self.assertEqual(round(current_test_value), round(current_value))
    #             else:
    #                 self.assertEqual(current_test_value, current_value)

    def test_5(self):
        print(magentaStr(f"\n==> Running import movies, movie alignment and max shift test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Movie_alignment": self.test_data["Movie_alignment"],
                            "Movie_maxshift": self.test_data["Movie_maxshift"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift)
        self.launchProtocol(prot)

        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, current_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in current_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Ouput_avg_shift" or key_in == "Output_max_shift":
                    # these values change decimals each time alignment is run
                    self.assertEqual(round(current_test_value), round(current_value))
                else:
                    self.assertEqual(current_test_value, current_value)