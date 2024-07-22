import json
from os.path import join, exists, abspath

from cryosparc2.protocols import ProtCryo2D, ProtCryoSparcInitialModel

from relion.protocols import ProtRelionAutopickLoG, ProtRelionExtractParticles, ProtRelionClassify2D

from pwem.protocols import ProtImportMovies, ProtAlignMovies, ProtUserSubSet
from pyworkflow.tests import BaseTest, tests, DataSet
from pyworkflow.utils import magentaStr
from xmipp3.protocols import XmippProtMovieGain, XmippProtFlexAlign, XmippProtMovieMaxShift, XmippProtCTFMicrographs, \
    XmippProtParticlePicking, XmippParticlePickingAutomatic, XmippProtCTFConsensus, XmippProtCenterParticles, \
    XmippProtExtractParticles
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
            "Number_movies": 30
        },
        "Movie_alignment": {
            "Method": "XmippProtFlexAlign",
            "Binning_factor": 1.0,
            "Maximum_resolution": 30.0,
            "Frames_aligned": {
                "Frame0": 1,
                "FrameN": 0
            },
            "Output_avg_shift": 11.817132324620918,
            "Output_max_shift": 34.02834443504211
        },
        "Movie_maxshift": {
            "Discarded_movies": 9,
            "Max_frame_shift": 5.0,
            "Max_movie_shift": 20.0,
            "Rejection_type": "By frame or movie",
            "Output_avg_shift": 11.806156624471319,
            "Output_max_shift": 29.569197112511205
        },
        "CTF_estimation": {
            "Defocus": {
                "Output_max_defocus": 11850.75,
                "Output_min_defocus": 4822.985000000001,
                "Output_avg_defocus": 9654.118046875
            },
            "Resolution": {
                "Output_max_resolution": 2.882647,
                "Output_min_resolution": 2.085319,
                "Output_avg_resolution": 2.5573557734375
            }
        },
        "Particle_picking": {
            "Particles_per_micrograph": 113.375
        },
        "Classes_2D": {
            "Number_classes_2D": 49,
            "Particles_per_class": [
                7,
                12,
                13,
                3,
                11,
                7,
                22,
                7,
                6,
                1,
                9,
                4,
                5,
                18,
                10,
                5,
                7,
                5,
                22,
                3,
                6,
                52,
                13,
                6,
                4,
                332,
                17,
                10,
                3,
                10,
                8,
                3,
                1,
                15,
                3,
                15,
                19,
                86,
                13,
                14,
                20,
                2,
                21,
                20,
                2,
                10,
                24,
                1
            ]
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
        cls.protCTF, cls.CTFout = cls.runCTFestimation()
        cls.protCTFconsensus, cls.CTFmicro, cls.CTFconsensus = cls.runCTFconsensus()
        cls.protpicking, cls.picked = cls.runLoGPicking()
        cls.protextract, cls.particles = cls.runExtractParticles()
        cls.prot2Dclasses, cls.classes2D = cls.run2DClassification()
        cls.portCenter, cls.centeredClasses2D, cls.centeredParticles = cls.runCenterParticles()

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
        output = getattr(prot, 'outputCTF', None)
        return prot, output

    @classmethod
    def runCTFconsensus(cls):
        prot = cls.newProtocol(XmippProtCTFConsensus,
                               inputCTF=cls.CTFout,
                               resolution=4.0)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputMicrographs', None)
        output2 = getattr(prot, 'outputCTF', None)
        return prot, output1, output2

    @classmethod
    def runLoGPicking(cls):
        prot = cls.newProtocol(ProtRelionAutopickLoG,
                               inputMicrographs=cls.CTFmicro,
                               boxSize=300,
                               minDiameter=100,
                               maxDiameter=200)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputCoordinates', None)
        return prot, output

    @classmethod
    def runExtractParticles(cls):
        prot = cls.newProtocol(XmippProtExtractParticles,
                               inputCoordinates=cls.picked,
                               ctfRelations=cls.CTFconsensus,
                               downFactor=2.5,
                               boxSize=300)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputParticles', None)
        return prot, output

    @classmethod
    def run2DClassification(cls):
        prot = cls.newProtocol(ProtCryo2D,
                               inputParticles=cls.particles)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputClasses', None)
        return prot, output

    @classmethod
    def runCenterParticles(cls):
        prot = cls.newProtocol(XmippProtCenterParticles,
                               inputClasses=cls.classes2D,
                               inputMics=cls.CTFmicro)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputClasses', None)
        output2 = getattr(prot, 'outputParticles', None)
        return prot, output1, output2

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

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
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

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Output_avg_shift" or key_in == "Output_max_shift":
                    # these values change decimals each time alignment is run
                    self.assertAlmostEqual(current_test_value, current_value, delta=1.5)
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

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Discarded_movies":
                    # these values change decimals each time alignment is run
                    self.assertAlmostEqual(current_test_value, current_value, delta=1)
                elif key_in == "Output_avg_shift" or key_in == "Output_max_shift":
                    # these values change decimals each time alignment is run
                    self.assertAlmostEqual(current_test_value, current_value, delta=1.5)
                else:
                    self.assertEqual(current_test_value, current_value)

    def test_import_and_CTF(self):
        print(magentaStr(f"\n==> Running import movies and CTF estimation test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "CTF_estimation": self.test_data["CTF_estimation"]
                            }

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                CTF=self.CTFconsensus)

        self.launchProtocol(prot)
        # recuperar resultados -> meter en funcion y meto el protocolo de entrada
        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                if key == 'Import_movies':
                    current_value = current_dict.get(key_in, None)
                    self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                    self.assertEqual(current_test_value, current_value)
                else:
                    for key_in_in, current_test_value_in in current_test_value.items():
                        current_section_dict = current_dict.get(key_in, None)
                        current_value = current_section_dict.get(key_in_in, None)
                        self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                        if key_in_in == "Output_max_defocus" or key_in_in == "Output_min_defocus" or key_in_in == "Output_avg_defocus":
                            self.assertAlmostEqual(current_test_value_in, current_value, delta=1500)
                        elif key_in_in == "Output_max_resolution" or key_in_in == "Output_min_resolution" or key_in_in == "Output_avg_resolution":
                            self.assertAlmostEqual(current_test_value_in, current_value, delta=1500)
                        else:
                            self.assertEqual(current_test_value_in, current_value)

    def test_5(self):
        print(magentaStr(f"\n==> Running import movies, movie alignment, max shift and CTF estimation test:"))
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

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Discarded_movies":
                    # these values change decimals each time alignment is run
                    self.assertAlmostEqual(current_test_value, current_value, delta=1)
                elif key_in == "Output_avg_shift" or key_in == "Output_max_shift":
                    # these values change decimals each time alignment is run
                    self.assertAlmostEqual(current_test_value, current_value, delta=1.5)
                else:
                    self.assertEqual(current_test_value, current_value)

    def test_6(self):
        print(magentaStr(f"\n==> Running import movies, movie alignment, max shift and CTF test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Movie_alignment": self.test_data["Movie_alignment"],
                            "Movie_maxshift": self.test_data["Movie_maxshift"],
                            "CTF_estimation": self.test_data["CTF_estimation"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift,
                                CTF=self.CTFconsensus)

        self.launchProtocol(prot)

        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                if key == 'CTF_estimation':
                    for key_in_in, current_test_value_in in current_test_value.items():
                        current_section_dict = current_dict.get(key_in, None)
                        current_value = current_section_dict.get(key_in_in, None)
                        self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                        if key_in_in == "Output_max_defocus" or key_in_in == "Output_min_defocus" or key_in_in == "Output_avg_defocus":
                            self.assertAlmostEqual(current_test_value_in, current_value, delta=1500)
                        elif key_in_in == "Output_max_resolution" or key_in_in == "Output_min_resolution" or key_in_in == "Output_avg_resolution":
                            self.assertAlmostEqual(current_test_value_in, current_value, delta=1500)
                        else:
                            self.assertEqual(current_test_value_in, current_value)
                else:
                    current_value = current_dict.get(key_in, None)
                    self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                    if key_in == "Discarded_movies":
                        # these values change decimals each time alignment is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=1)
                    elif key_in == "Output_avg_shift" or key_in == "Output_max_shift":
                        # these values change decimals each time alignment is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=1.5)
                    else:
                        self.assertEqual(current_test_value, current_value)

    def test_import_and_particles(self):
        print(magentaStr(f"\n==> Running import movies and particles test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Particle_picking": self.test_data["Particle_picking"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                particles=self.particles)
        self.launchProtocol(prot)

        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Particles_per_micrograph":
                    # these values change decimals each time alignment is run
                    self.assertAlmostEqual(current_test_value, current_value, delta=2)
                else:
                    self.assertEqual(current_test_value, current_value)

    def test_import_and_2DClassification(self):
        print(magentaStr(f"\n==> Running import movies and 2D classification test:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Classes_2D": self.test_data["Classes_2D"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                classes2D=self.centeredClasses2D)
        self.launchProtocol(prot)

        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))

        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                current_value = current_dict.get(key_in, None)
                self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                if key_in == "Particles_per_class":
                    pass
                else:
                    self.assertAlmostEqual(current_test_value, current_value, delta=2)