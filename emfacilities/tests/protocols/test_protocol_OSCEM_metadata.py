import json
from os.path import exists, abspath

import yaml

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
    "movies": {
        "gain_image": "gain.jpg",
        "descriptors": [
            {
                "descriptor_name": "XmippProtFlexAlign",
                "binning_factor": 1.0,
                "maximum_resolution": {
                    "value": 30.0,
                    "unit": "Å"
                },
                "frames_aligned": {
                    "frame0": 1,
                    "frameN": 30
                },
                "output_avg_shift": {
                    "value": 11.8,
                    "unit": "Å"
                },
                "output_max_shift": {
                    "value": 34.0,
                    "unit": "Å"
                }
            },
            {
                "descriptor_name": "XmippProtMovieMaxShift",
                "discarded_movies": 9,
                "max_frame_shift": {
                    "value": 5.0,
                    "unit": "Å"
                },
                "max_movie_shift": {
                    "value": 20.0,
                    "unit": "Å"
                },
                "rejection_type": "by frame or movie",
                "output_avg_shift": {
                    "value": 11.8,
                    "unit": "Å"
                },
                "output_max_shift": {
                    "value": 29.6,
                    "unit": "Å"
                },
                "shift_histogram": "shift_hist.jpg"
            }
        ]
    },
    "micrographs": {
        "number_micrographs": 21
    },
    "CTFs": {
        "amplitude_contrast": 0.1,
        "defocus": {
            "output_min_defocus": {
                "value": 1469.9,
                "unit": "Å"
            },
            "output_max_defocus": {
                "value": 11850.8,
                "unit": "Å"
            },
            "output_avg_defocus": {
                "value": 1723.3,
                "unit": "Å"
            },
            "defocus_histogram": "defocus_hist.jpg",
            "micrograph_examples": "Micro_examples/micro_defocus.jpg"
        },
        "resolution": {
            "output_min_resolution": {
                "value": 2.1,
                "unit": "Å"
            },
            "output_max_resolution": {
                "value": 5.2,
                "unit": "Å"
            },
            "output_avg_resolution": {
                "value": 3.0,
                "unit": "Å"
            },
            "resolution_histogram": "resolution_hist.jpg"
        },
        "astigmatism": {
            "astigmatism_histogram": "astigmatism_hist.jpg"
        }
    },
    "coordinates": {
        "number_particles": 2860,
        "particles_per_micrograph": 136.2,
        "particles_histogram": "particles_hist.jpg",
        "micrograph_examples": "Micro_examples/micro_particles.jpg"
    },
    "classes2D": {
        "number_classes_2D": 49,
        "particles_per_class": [
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
        "images_classes_2D": "classes_2D.jpg"
    },
    "classes3D": {
        "number_classes_3D": 2,
        "particles_per_class": [
            1985,
            875
        ],
        "images_classes_3D": "Classes_3D/classes_3D.jpg",
        "volumes": [
            {
                "orthogonal_slices": {
                    "orthogonal_slices_X": "Classes_3D/orthogonal_slices_volume1/orthogonal_slices_X.jpg",
                    "orthogonal_slices_Y": "Classes_3D/orthogonal_slices_volume1/orthogonal_slices_Y.jpg",
                    "orthogonal_slices_Z": "Classes_3D/orthogonal_slices_volume1/orthogonal_slices_Z.jpg"
                },
                "isosurface_images": {
                    "front_view": "Classes_3D/isosurface_images_volume1/front_view.jpg",
                    "dide_view": "Classes_3D/isosurface_images_volume1/side_view.jpg",
                    "top_view": "Classes_3D/isosurface_images_volume1/top_view.jpg"
                }
            },
            {
                "orthogonal_slices": {
                    "orthogonal_slices_X": "Classes_3D/orthogonal_slices_volume2/orthogonal_slices_X.jpg",
                    "orthogonal_slices_Y": "Classes_3D/orthogonal_slices_volume2/orthogonal_slices_Y.jpg",
                    "orthogonal_slices_Z": "Classes_3D/orthogonal_slices_volume2/orthogonal_slices_Z.jpg"
                },
                "isosurface_images": {
                    "front_view": "Classes_3D/isosurface_images_volume2/front_view.jpg",
                    "dide_view": "Classes_3D/isosurface_images_volume2/side_view.jpg",
                    "top_view": "Classes_3D/isosurface_images_volume2/top_view.jpg"
                }
            }
        ]
    },
    "volumes": [
        {
            "volume_type": "initial volume",
            "orthogonal_slices": {
                "orthogonal_slices_X": "Initial_volume/orthogonal_slices/orthogonal_slices_X.jpg",
                "orthogonal_slices_Y": "Initial_volume/orthogonal_slices/orthogonal_slices_Y.jpg",
                "orthogonal_slices_Z": "Initial_volume/orthogonal_slices/orthogonal_slices_Z.jpg"
            },
            "isosurface_images": {
                "front_view": "Initial_volume/isosurface_images/front_view.jpg",
                "side_view": "Initial_volume/isosurface_images/side_view.jpg",
                "top_view": "Initial_volume/isosurface_images/top_view.jpg"
            }
        }
    ]
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
        test_data_import = {"movies": self.test_data["movies"],
                            "micrographs": self.test_data["micrographs"],
                            "CTFs": self.test_data["CTFs"],
                            "coordinates": self.test_data["coordinates"],
                            "classes2D": self.test_data["classes2D"],
                            "volumes": self.test_data["volumes"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift,
                                micrographs=self.maxshiftmicro,
                                CTF=self.CTFout,
                                particles=self.particles,
                                classes2D=self.classes2D,
                                initVolume=self.importedvolume,  # Pointer(self.initVolVolumes.getFirstItem()),
                                classes3D=self.classes3DClassification)

        load_json = self.prot_and_load_json(prot)

        # for key, section_test_dict in test_data_import.items():
        #     current_dict = load_json.get(key, None)
        #     self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
        #     for key_in, current_test_value in section_test_dict.items():
        #         if key == 'CTF_estimation':
        #             self.CTF_comparison(current_test_value, current_dict, key_in, key)
        #         else:
        #             current_value = current_dict.get(key_in, None)
        #             self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
        #             # Handle cases where the test data now uses value-unit pairs
        #             if isinstance(current_test_value,
        #                           dict) and 'value' in current_test_value and 'unit' in current_test_value:
        #                 # Retrieve value and unit separately
        #                 test_value = current_test_value.get('value', None)
        #                 test_unit = current_test_value.get('unit', None)
        #
        #                 # Extract value and unit from the current dictionary entry
        #                 current_value_data = current_value.get('value', None)
        #                 current_unit_data = current_value.get('unit', None)
        #
        #                 # Check values
        #                 if key_in == "Output_avg_shift" or key_in == "Output_max_shift":
        #                     self.assertAlmostEqual(test_value, current_value_data, delta=4)
        #                 else:
        #                     self.assertEqual(test_value, current_value_data, msg=f'Value mismatch for {key_in}')
        #
        #                 # Check units
        #                 self.assertEqual(test_unit, current_unit_data, msg=f'Unit mismatch for {key_in}')
        #
        #             else:
        #                 # Original comparison logic if it's not a value-unit pair
        #                 if key_in == "Discarded_movies":
        #                     self.assertAlmostEqual(current_test_value, current_value, delta=2)
        #                 elif key_in == "Number_particles":
        #                     self.assertAlmostEqual(current_test_value, current_value, delta=1000)
        #                 elif key_in == "Particles_per_micrograph":
        #                     self.assertAlmostEqual(current_test_value, current_value, delta=10)
        #                 elif key_in == "Number_classes_2D":
        #                     self.assertAlmostEqual(current_test_value, current_value, delta=2)
        #                 elif key_in == 'Particles_per_class':
        #                     pass
        #                 else:
        #                     self.assertEqual(current_test_value, current_value)
        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')

            for key_in, current_test_value in section_test_dict.items():
                if key == 'CTFs':
                    self.CTF_comparison(current_test_value, current_dict, key_in, key)
                else:
                    current_value = current_dict.get(key_in, None)
                    self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')

                    # Handle cases where the test data uses value-unit pairs or dictionaries
                    if isinstance(current_test_value,
                                  dict) and 'value' in current_test_value and 'unit' in current_test_value:
                        # Handle value-unit dictionary comparison
                        test_value = current_test_value.get('value', None)
                        test_unit = current_test_value.get('unit', None)

                        current_value_data = current_value.get('value', None)
                        current_unit_data = current_value.get('unit', None)

                        # Compare the values with tolerance
                        if key_in == "output_avg_shift" or key_in == "output_max_shift":
                            self.assertAlmostEqual(test_value, current_value_data, delta=0.2,
                                                   msg=f"Value mismatch for {key_in}")
                        else:
                            self.assertEqual(test_value, current_value_data, msg=f"Value mismatch for {key_in}")

                        # Compare the units
                        self.assertEqual(test_unit, current_unit_data, msg=f"Unit mismatch for {key_in}")

                    elif isinstance(current_test_value, list):
                        # Handle cases where the values are lists
                        self.assertIsInstance(current_value, list,
                                              msg=f"Expected a list for {key_in}, but got {type(current_value)}")
                        self.assertEqual(len(current_test_value), len(current_value),
                                         msg=f"List lengths do not match for {key_in}")

                        # Iterate through the lists and compare each item
                        for test_item, current_item in zip(current_test_value, current_value):
                            if isinstance(test_item, dict) and 'value' in test_item and 'unit' in test_item:
                                # Compare dictionaries within the lists
                                self.assertAlmostEqual(test_item['value'], current_item['value'], delta=0.2,
                                                       msg=f"Value mismatch in list for {key_in}")
                                self.assertEqual(test_item['unit'], current_item['unit'],
                                                 msg=f"Unit mismatch in list for {key_in}")
                            else:
                                # Compare items directly if they are not dictionaries
                                self.assertEqual(test_item, current_item,
                                                 msg=f"Item mismatch in list for {key_in}")
                    elif isinstance(current_test_value, dict):
                        # General dictionary comparison (not specifically value-unit)
                        for sub_key, sub_value in current_test_value.items():
                            current_sub_value = current_value.get(sub_key, None)
                            self.assertEqual(sub_value, current_sub_value,
                                             msg=f"Mismatch for {key_in} -> {sub_key}")
                    else:
                        # Handle direct value comparisons
                        if key_in == "discarded_movies":
                            self.assertAlmostEqual(current_test_value, current_value, delta=2)
                        elif key_in == "number_particles":
                            self.assertAlmostEqual(current_test_value, current_value, delta=1000)
                        elif key_in == "particles_per_micrograph":
                            self.assertAlmostEqual(current_test_value, current_value, delta=10)
                        elif key_in == "number_classes_2D":
                            self.assertAlmostEqual(current_test_value, current_value, delta=2)
                        elif key_in == 'particles_per_class':
                            pass  # Custom logic for specific cases if needed
                        else:
                            self.assertEqual(current_test_value, current_value)
    def test_only_compulsory(self):
        print(magentaStr("\n==> Running test with only compulsory input completed:"))
        test_data_import = {"movies": self.test_data["movies"]}

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
                # Handle cases where the test data now uses value-unit pairs
                if isinstance(current_test_value,
                              dict) and 'value' in current_test_value and 'unit' in current_test_value:
                    # Retrieve value and unit separately
                    test_value = current_test_value.get('value', None)
                    test_unit = current_test_value.get('unit', None)

                    # Extract value and unit from the current dictionary entry
                    current_value_data = current_value.get('value', None)
                    current_unit_data = current_value.get('unit', None)
                    # Check value
                    self.assertEqual(test_value, current_value_data, msg=f'Value mismatch for {key_in}')
                    # Check units
                    self.assertEqual(test_unit, current_unit_data, msg=f'Unit mismatch for {key_in}')

                else:
                    self.assertEqual(current_test_value, current_value)

    def test_medium_level(self):
        print(magentaStr("\n==> Running test with compulsory and some optional input:"))
        test_data_import = {"movies": self.test_data["movies"],
                            "CTFs": self.test_data["CTFs"]}

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
                if key == 'CTFs':
                    self.CTF_comparison(current_test_value, current_dict, key_in, key)
                else:
                    current_value = current_dict.get(key_in, None)
                    self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                    if isinstance(current_test_value,
                                  dict) and 'value' in current_test_value and 'unit' in current_test_value:
                        # Retrieve value and unit separately
                        test_value = current_test_value.get('value', None)
                        test_unit = current_test_value.get('unit', None)

                        # Extract value and unit from the current dictionary entry
                        current_value_data = current_value.get('value', None)
                        current_unit_data = current_value.get('unit', None)

                        # Check values
                        if key_in == "output_avg_shift" or key_in == "output_max_shift":
                            self.assertAlmostEqual(test_value, current_value_data, delta=4)
                        else:
                            self.assertEqual(test_value, current_value_data, msg=f'Value mismatch for {key_in}')

                        # Check units
                        self.assertEqual(test_unit, current_unit_data, msg=f'Unit mismatch for {key_in}')

                    else:
                        # Original comparison logic if it's not a value-unit pair
                        if key_in == "discarded_movies":
                            self.assertAlmostEqual(current_test_value, current_value, delta=2)
                        else:
                            self.assertEqual(current_test_value, current_value)


    def test_micro_input(self):
        print(magentaStr("\n==> Running test with micrographs as input:"))
        test_data_import = {"CTFs": self.test_data["CTFs"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MICS,
                                CTF=self.CTFout)

        load_json = self.prot_and_load_json(prot)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            for key, section_test_dict in test_data_import.items():
                current_dict = load_json.get(key, None)
                self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
                for key_in, current_test_value in section_test_dict.items():
                    if key == 'movies':
                        current_value = current_dict.get(key_in, None)
                        self.assertIsNotNone(current_value, msg=f'In dictionary {key}, {key_in} is not found')
                        if isinstance(current_test_value,
                                      dict) and 'value' in current_test_value and 'unit' in current_test_value:
                            # Retrieve value and unit separately
                            test_value = current_test_value.get('value', None)
                            test_unit = current_test_value.get('unit', None)
                            # Extract value and unit from the current dictionary entry
                            current_value_data = current_value.get('value', None)
                            current_unit_data = current_value.get('unit', None)
                            # Check values
                            self.assertEqual(test_value, current_value_data, msg=f'Value mismatch for {key_in}')
                            # Check units
                            self.assertEqual(test_unit, current_unit_data, msg=f'Unit mismatch for {key_in}')
                        else:
                            self.assertEqual(current_test_value, current_value)

                    else:
                        current_section_dict = current_dict.get(key_in, None)
                        self.assertIsNotNone(current_section_dict, msg=f'In dictionary {key}, {key_in} is not found')

                        # If current_test_value is a dictionary, loop through its items
                        if isinstance(current_test_value, dict):
                            for key_in_in, current_test_value_in in current_test_value.items():
                                current_value = current_section_dict.get(key_in_in, None)
                                self.assertIsNotNone(current_value,
                                                     msg=f'In dictionary {key}, {key_in_in} is not found')
                                # Handle cases where the test data now uses value-unit pairs
                                if isinstance(current_test_value_in,
                                              dict) and 'value' in current_test_value_in and 'unit' in current_test_value_in:
                                    # Retrieve value and unit separately
                                    test_value = current_test_value_in.get('value', None)
                                    test_unit = current_test_value_in.get('unit', None)

                                    # Extract value and unit from the current dictionary entry
                                    current_value_data = current_value.get('value', None)
                                    current_unit_data = current_value.get('unit', None)

                                    # Check values
                                    if key_in_in in ["output_max_defocus", "output_min_defocus", "ouput_avg_defocus",
                                                     "output_max_resolution", "output_min_resolution",
                                                     "ouput_avg_resolution"]:
                                        # these values change each time alignment protocol is run
                                        self.assertAlmostEqual(test_value, current_value_data, delta=3000)
                                    else:
                                        self.assertEqual(test_value, current_value_data,
                                                         msg=f'Value mismatch for {key_in_in}')

                                    # Check units
                                    self.assertEqual(test_unit, current_unit_data, msg=f'Unit mismatch for {key_in_in}')

                                else:
                                    # Compare simple values directly
                                    self.assertEqual(current_test_value_in, current_value)
                        else:
                            # If current_test_value is not a dictionary, compare it directly
                            self.assertEqual(current_test_value, current_section_dict,
                                             msg=f"Value mismatch for {key_in}")


    def CTF_comparison(self, current_test_value, current_dict, key_in, key):
        current_section_dict = current_dict.get(key_in, None)
        self.assertIsNotNone(current_section_dict, msg=f'In dictionary {key}, {key_in} is not found')

        # Check if current_test_value is a dictionary
        if isinstance(current_test_value, dict):
            for sub_key, sub_value in current_test_value.items():
                # Check if the sub_value has both 'value' and 'unit'
                if isinstance(sub_value, dict) and 'value' in sub_value and 'unit' in sub_value:
                    test_value = sub_value['value']
                    test_unit = sub_value['unit']

                    # Get the corresponding value and unit from the current section
                    current_value = current_section_dict.get(sub_key, {}).get('value', None)
                    current_unit = current_section_dict.get(sub_key, {}).get('unit', None)

                    # Compare the units
                    self.assertEqual(test_unit, current_unit, msg=f"Unit mismatch for {sub_key} in {key}")

                    # Compare the values with appropriate delta
                    if sub_key in ["output_min_defocus", "output_max_defocus", "ouput_avg_defocus",
                                   "output_min_resolution", "output_max_resolution", "ouput_avg_resolution"]:
                        self.assertAlmostEqual(test_value, current_value, delta=3000)
                    else:
                        self.assertEqual(test_value, current_value)

                else:
                    # Handle simple values  that don't have value-unit structure
                    current_value = current_section_dict[sub_key]
                    self.assertEqual(sub_value, current_value, msg=f"Value mismatch for {sub_key} in {key}")

        else:
            # Handle simple float or int values
            self.assertEqual(current_test_value, current_section_dict, msg=f"Value mismatch for {key_in} in {key}")


    def prot_and_load_json(self, prot):
        self.launchProtocol(prot)
        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))
        # with open(abspath(file_path), 'r') as json_file:
        #     load_json = json.load(json_file)
        # Open the YAML file
        with open(abspath(file_path), 'r', encoding='utf-8') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)
            json_string = json.dumps(yaml_data, ensure_ascii=False)
            load_json = json.loads(json_string)

        return load_json