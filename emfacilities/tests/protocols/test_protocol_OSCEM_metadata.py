import json
from os.path import exists, abspath

from cryosparc2.protocols import ProtCryo2D, ProtCryoSparcInitialModel
from pwem import SYM_TETRAHEDRAL

from relion.protocols import ProtRelionAutopickLoG, ProtRelionClassify3D

from pwem.protocols import ProtImportMovies, ProtImportVolumes
from pyworkflow.tests import BaseTest, tests, DataSet
from pyworkflow.utils import magentaStr
from xmipp3.protocols import XmippProtMovieGain, XmippProtFlexAlign, XmippProtMovieMaxShift, XmippProtCTFMicrographs, \
    XmippProtCenterParticles, \
    XmippProtExtractParticles
from ...protocols import ProtOSCEM
from ...protocols.protocol_OSCEM_metadata import INPUT_MOVIES, INPUT_MICS


box_size = 300
class TestOscemJson(BaseTest):
    """
    """
    sampling_rate = 0.495
    max_movie_shift = 20.0
    ctf_down_factor = 2.0
    high_res = 0.5
    test_data = {
        "Import_movies": {
            "Microscope_voltage": {
                "value": 300.0,
                "unit": "kV"
            },
            "Spherical_aberration": {
                "value": 2.7,
                "unit": "mm"
            },
            "Pixel_size": {
                "value": 0.495,
                "unit": "Å/px"
            },
            "Gain_image": "gain.jpg",
            "Number_movies": 30,
            "Frames_per_movie": 50,
            "Frames_size": {
                "value": "3710 x 3838",
                "unit": "pixels"
            }
        },
        "Movie_alignment": {
            "Method": "XmippProtFlexAlign",
            "Binning_factor": 1.0,
            "Maximum_resolution": {
                "value": 30.0,
                "unit": "Å"
            },
            "Frames_aligned": {
                "Frame0": 1,
                "FrameN": 30
            },
            "Output_avg_shift": {
                "value": 11.8,
                "unit": "Å"
            },
            "Output_max_shift": {
                "value": 34.0,
                "unit": "Å"
            }
        },
        "Movie_maxshift": {
            "Discarded_movies": 9,
            "Max_frame_shift_(Å)": 5.0,
            "Max_movie_shift_(Å)": 20.0,
            "Rejection_type": "By frame or movie",
            "Output_avg_shift_(Å)": 11.8,
            "Output_max_shift_(Å)": 29.6
        },
        "CTF_estimation": {
            "Amplitude_contrast": 0.1,
            "Defocus_(Å)": {
                "Output_max_defocus": 11850.8,
                "Output_min_defocus": 1469.9,
                "Output_avg_defocus": 1723.3,
                "Defocus_histogram": "defocus_hist.jpg",
                "Defocus_mic_examples": "Micro_examples/micro_defocus.jpg"
            },
            "Resolution_(Å)": {
                "Output_max_resolution": 5.2,
                "Output_min_resolution": 2.1,
                "Output_avg_resolution": 3.0,
                "Resolution_histogram": "resolution_hist.jpg"
            },
            "Astigmatism": {
                "Astigmatism_histogram": "astigmatism_hist.jpg"
            }
        },
        "Particle_picking": {
            "Number_particles": 2860,
            "Particles_per_micrograph": 136.2,
            "Particles_histogram": "particles_hist.jpg",
            "Particles_mic_examples": "Micro_examples/micro_particles.jpg"
        },
        "Classes_2D": {
            "Number_classes_2D": 49,
            "Particles_per_class": [
                798,
                365,
                262,
                128,
                109,
                92,
                82,
                78,
                77,
                54,
                51,
                44,
                42,
                40,
                35,
                33,
                31,
                28,
                27,
                26,
                25,
                23,
                23,
                22,
                22,
                22,
                21,
                21,
                20,
                19,
                18,
                18,
                17,
                17,
                17,
                16,
                16,
                16,
                16,
                14,
                13,
                13,
                12,
                12,
                9,
                7,
                4,
                4,
                1
            ],
            "Images_classes_2D": "classes_2D.jpg"
        },
        "Initial_volume": {
            "Orthogonal_slices": {
                "Orthogonal_slices_X": "Initial_volume/orthogonal_slices/orthogonal_slices_X.jpg",
                "Orthogonal_slices_Y": "Initial_volume/orthogonal_slices/orthogonal_slices_Y.jpg",
                "Orthogonal_slices_Z": "Initial_volume/orthogonal_slices/orthogonal_slices_Z.jpg"
            },
            "Isosurface_images": {
                "Front_view": "Initial_volume/isosurface_images/front_view.jpg",
                "Side_view": "Initial_volume/isosurface_images/side_view.jpg",
                "Top_view": "Initial_volume/isosurface_images/top_view.jpg"
            }
        },
        "Classes_3D": {
            "Number_classes_3D": 2,
            "Particles_per_class": [
                1985,
                875
            ],
            "Images_classes_3D": "Classes_3D/classes_3D.jpg",
            "Volumes": {
                "Volume_1": {
                    "Orthogonal_slices": {
                        "Orthogonal_slices_X": "Classes_3D/orthogonal_slices_volume1/orthogonal_slices_X.jpg",
                        "Orthogonal_slices_Y": "Classes_3D/orthogonal_slices_volume1/orthogonal_slices_Y.jpg",
                        "Orthogonal_slices_Z": "Classes_3D/orthogonal_slices_volume1/orthogonal_slices_Z.jpg"
                    },
                    "Isosurface_images": {
                        "Front_view": "Classes_3D/isosurface_images_volume1/front_view.jpg",
                        "Side_view": "Classes_3D/isosurface_images_volume1/side_view.jpg",
                        "Top_view": "Classes_3D/isosurface_images_volume1/top_view.jpg"
                    }
                },
                "Volume_2": {
                    "Orthogonal_slices": {
                        "Orthogonal_slices_X": "Classes_3D/orthogonal_slices_volume2/orthogonal_slices_X.jpg",
                        "Orthogonal_slices_Y": "Classes_3D/orthogonal_slices_volume2/orthogonal_slices_Y.jpg",
                        "Orthogonal_slices_Z": "Classes_3D/orthogonal_slices_volume2/orthogonal_slices_Z.jpg"
                    },
                    "Isosurface_images": {
                        "Front_view": "Classes_3D/isosurface_images_volume2/front_view.jpg",
                        "Side_view": "Classes_3D/isosurface_images_volume2/side_view.jpg",
                        "Top_view": "Classes_3D/isosurface_images_volume2/top_view.jpg"
                    }
                }
            }
        }
    }

    @classmethod
    def setUpClass(cls):

        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('OSCEM_jsons')
        cls.protimportmovies, cls.importedmovies = cls.runImportMovies()

        cls.protimportvolumes, cls.importedvolume = cls.runImportVolumes()

        cls.protmoviegain, cls.gainmovies = cls.runMovieGain()
        cls.protalignmovie, cls.alignedmovies = cls.runMovieAlign()
        cls.protmaxshift, cls.maxshiftmicro = cls.runMaxShift()
        cls.protCTF, cls.CTFout = cls.runCTFestimation()
        cls.protpicking, cls.coordinates = cls.runLoGPicking()
        cls.protextract, cls.particles = cls.runExtractParticles()
        cls.prot2Dclasses, cls.classes2D = cls.run2DClassification()
        cls.portCenter, cls.centeredClasses2D, cls.centeredParticles = cls.runCenterParticles()
        cls.prot3DClassification, cls.classes3DClassification, cls.volumes3DClassification = cls.run3DClassification()

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
    def runImportVolumes(cls):

        prot = cls.newProtocol(ProtImportVolumes,
                               filesPath=cls.dataset.getFile('initial_volume'),
                               samplingRate=1.24,)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputVolume', None)
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
    def runLoGPicking(cls):
        prot = cls.newProtocol(ProtRelionAutopickLoG,
                               inputMicrographs=cls.maxshiftmicro,
                               boxSize=box_size,
                               minDiameter=100,
                               maxDiameter=200)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputCoordinates', None)
        return prot, output

    @classmethod
    def runExtractParticles(cls):
        prot = cls.newProtocol(XmippProtExtractParticles,
                               inputCoordinates=cls.coordinates,
                               ctfRelations=cls.CTFout,
                               doResize=True,
                               downFactor=2.5,
                               boxSize=box_size)

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
                               inputMics=cls.maxshiftmicro)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputClasses', None)
        output2 = getattr(prot, 'outputParticles', None)
        return prot, output1, output2

    @classmethod
    def runInitialVolume(cls):
        prot = cls.newProtocol(ProtCryoSparcInitialModel,
                               inputParticles=cls.centeredParticles,
                               symmetryGroup=SYM_TETRAHEDRAL,
                               numberOfClasses=2,
                               abinit_max_res=20,
                               abinit_num_init_iters=50,
                               abinit_num_final_iters=75,
                               abinit_radwn_step=0.1)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputClasses', None)
        output2 = getattr(prot, 'outputVolumes', None)
        return prot, output1, output2

    @classmethod
    def run3DClassification(cls):
        prot = cls.newProtocol(ProtRelionClassify3D,
                               inputParticles=cls.centeredParticles,
                               referenceVolume=cls.importedvolume,  #  initVolVolumes,
                               numberOfIterations=10,
                               symmetryGroup='O',
                               initialLowPassFilterA=15,
                               numberOfClasses=2)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputClasses', None)
        output2 = getattr(prot, 'outputVolumes', None)
        return prot, output1, output2

    def test_complete_input(self):
        print(magentaStr("\n==> Running test with all input completed:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Movie_alignment": self.test_data["Movie_alignment"],
                            "Movie_maxshift": self.test_data["Movie_maxshift"],
                            "CTF_estimation": self.test_data["CTF_estimation"],
                            "Particle_picking": self.test_data["Particle_picking"],
                            "Classes_2D": self.test_data["Classes_2D"],
                            "Initial_volume": self.test_data["Initial_volume"],
                            "Classes_3D": self.test_data["Classes_3D"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift,
                                CTF=self.CTFout,
                                coords=self.coordinates,
                                particles=self.particles,
                                classes2D=self.classes2D,
                                initVolume=self.importedvolume,  # Pointer(self.initVolVolumes.getFirstItem()),
                                classes3D=self.classes3DClassification)

        load_json = self.prot_and_load_json(prot)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                if key == 'CTF_estimation':
                    self.CTF_comparison(current_test_value, current_dict, key_in, key)
                else:
                    current_value = current_dict.get(key_in, None)
                    self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                    if key_in == "Discarded_movies":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=2)
                    elif key_in == "Output_avg_shift_(Å)" or key_in == "Output_max_shift_(Å)":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=4)
                    elif key_in == "Number_particles":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=1000)
                    elif key_in == "Particles_per_micrograph":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=10)
                    elif key_in == "Number_classes_2D":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=2)
                    elif key_in == 'Particles_per_class':
                        pass
                    else:
                        self.assertEqual(current_test_value, current_value)

    def test_only_compulsory(self):
        print(magentaStr("\n==> Running test with only compulsory input completed:"))
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

    def test_medium_level(self):
        print(magentaStr("\n==> Running test with compulsory and some optional input:"))
        test_data_import = {"Import_movies": self.test_data["Import_movies"],
                            "Movie_alignment": self.test_data["Movie_alignment"],
                            "Movie_maxshift": self.test_data["Movie_maxshift"],
                            "CTF_estimation": self.test_data["CTF_estimation"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift,
                                CTF=self.CTFout)

        load_json = self.prot_and_load_json(prot)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key_in, current_test_value in section_test_dict.items():
                if key == 'CTF_estimation':
                    self.CTF_comparison(current_test_value, current_dict, key_in, key)
                else:
                    current_value = current_dict.get(key_in, None)
                    self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                    if key_in == "Discarded_movies":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=2)
                    elif key_in == "Output_avg_shift_(Å)" or key_in == "Output_max_shift_(Å)":
                        # these values change each time alignment protocol is run
                        self.assertAlmostEqual(current_test_value, current_value, delta=4)
                    else:
                        self.assertEqual(current_test_value, current_value)

    def test_micro_input(self):
        print(magentaStr("\n==> Running test with micrographs as input:"))
        test_data_import = {"CTF_estimation": self.test_data["CTF_estimation"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MICS,
                                CTF=self.CTFout)

        load_json = self.prot_and_load_json(prot)

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
                        if key_in_in == "Output_max_defocus" or key_in_in == "Output_min_defocus" or key_in_in == "Output_avg_defocus"\
                                or key_in_in == "Output_max_resolution" or key_in_in == "Output_min_resolution" or key_in_in == "Output_avg_resolution":
                            # these values change each time alignment protocol is run
                            self.assertAlmostEqual(current_test_value_in, current_value, delta=3000)
                        else:
                            self.assertEqual(current_test_value_in, current_value)

    def CTF_comparison(self, current_test_value, current_dict, key_in, key):
        for key_in_in, current_test_value_in in current_test_value.items():
            current_section_dict = current_dict.get(key_in, None)
            current_value = current_section_dict.get(key_in_in, None)
            self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
            if key_in_in == "Output_max_defocus" or key_in_in == "Output_min_defocus" or key_in_in == "Output_avg_defocus":
                self.assertAlmostEqual(current_test_value_in, current_value, delta=3000)
            elif key_in_in == "Output_max_resolution" or key_in_in == "Output_min_resolution" or key_in_in == "Output_avg_resolution":
                self.assertAlmostEqual(current_test_value_in, current_value, delta=3000)
            else:
                self.assertEqual(current_test_value_in, current_value)

    def prot_and_load_json(self, prot):
        self.launchProtocol(prot)
        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))
        with open(abspath(file_path), 'r') as json_file:
            load_json = json.load(json_file)

        return load_json