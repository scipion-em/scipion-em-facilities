import json
from collections.abc import Mapping, Sequence
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

    def test_only_compulsory(self):
        print(magentaStr("\n==> Running test with only compulsory input:"))
        test_data_import = {"movies": self.test_data["movies"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift)

        load_json = self.prot_and_load_json(prot)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)


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

        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)


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
        print(load_json)

        # Recursive comparison
        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)


    def test_micro_input(self):
        print(magentaStr("\n==> Running test with micrographs as input:"))
        test_data_import = {"CTFs": self.test_data["CTFs"]}

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MICS,
                                CTF=self.CTFout)

        load_json = self.prot_and_load_json(prot)

        # Recursive comparison
        for key, section_test_dict in test_data_import.items():
            current_dict = load_json.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)

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

    def recursive_compare(self, test_data, current_data, parent_key=""):
        """
        Compares nested dictionaries and lists
        """
        if isinstance(test_data, dict):
            # Compare dictionaries
            for key, test_value in test_data.items():
                self.assertIn(key, current_data, msg=f"Key '{key}' not found in {parent_key}")
                # Recursively call to compare nested values
                self.recursive_compare(
                    test_value, current_data[key], f"{parent_key}.{key}" if parent_key else key
                )

        elif isinstance(test_data, list):
            # Compare lists
            if parent_key.endswith("particles_per_class"):
                length_difference = abs(len(test_data) - len(current_data))
                self.assertLessEqual(
                    length_difference, 2,
                    msg=f"List length mismatch at {parent_key}: {len(test_data)} != {len(current_data)}"
                )

                min_length = min(len(test_data), len(current_data))
                for i in range(min_length):
                    self.assertAlmostEqual(
                        test_data[i], current_data[i], delta=500,
                        msg=f"Value mismatch at {parent_key}[{i}]: {test_data[i]} != {current_data[i]}"
                    )
            else:
                self.assertEqual(
                    len(test_data), len(current_data),
                    msg=f"List length mismatch at {parent_key}: {len(test_data)} != {len(current_data)}"
                )
                for i, (test_item, current_item) in enumerate(zip(test_data, current_data)):
                    self.recursive_compare(
                        test_item, current_item, f"{parent_key}[{i}]"
                    )
        else:
            # Check if we are dealing with a value-unit pair
            if parent_key.endswith(".value") or parent_key.endswith(".unit"):
                # Handle value-unit pairs
                value_unit_key = ".".join(parent_key.split(".")[:-1])
                test_value = test_data if parent_key.endswith(".value") else None
                test_unit = test_data if parent_key.endswith(".unit") else None
                current_value = current_data if parent_key.endswith(".value") else None
                current_unit = current_data if parent_key.endswith(".unit") else None

                # Compare value and unit if both are present
                if "value" in value_unit_key or "unit" in value_unit_key:
                    print(f"Found value and unit at {value_unit_key}")
                    print(f"Comparing value: {test_value} with {current_value}, unit: {test_unit} with {current_unit}")

                    # Handle flexible comparisons for specific keys
                    last_key = value_unit_key.split('.')[-1]
                    if last_key in [
                        "output_max_defocus", "output_min_defocus", "output_avg_defocus",
                        "output_max_resolution", "output_min_resolution", "output_avg_resolution",
                    ]:
                        self.assertAlmostEqual(
                            test_value, current_value, delta=3000,
                            msg=f"Value mismatch at {value_unit_key}: {test_value} != {current_value}"
                        )
                    elif last_key in ["output_avg_shift", "output_max_shift"]:
                        self.assertAlmostEqual(
                            test_value, current_value, delta=1,
                            msg=f"Value mismatch at {value_unit_key}: {test_value} != {current_value}"
                        )
                    else:
                        self.assertEqual(
                            test_value, current_value,
                            msg=f"Value mismatch at {value_unit_key}: {test_value} != {current_value}"
                        )

                    # Compare units
                    if test_unit and current_unit:
                        self.assertEqual(
                            test_unit, current_unit,
                            msg=f"Unit mismatch at {value_unit_key}: {test_unit} != {current_unit}"
                        )
            else:
                key_in = parent_key.split('.')[-1]
                if key_in == "discarded_movies":
                    self.assertAlmostEqual(
                        test_data, current_data, delta=2,
                        msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
                    )
                elif key_in == "number_particles":
                    self.assertAlmostEqual(
                        test_data, current_data, delta=1000,
                        msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
                    )
                elif key_in == "particles_per_micrograph":
                    self.assertAlmostEqual(
                        test_data, current_data, delta=10,
                        msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
                    )
                elif key_in == "number_classes_2D":
                    self.assertAlmostEqual(
                        test_data, current_data, delta=2,
                        msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
                    )
                elif key_in == "number_classes_3D":
                    self.assertAlmostEqual(
                        test_data, current_data, delta=2,
                        msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
                    )
                else:
                    self.assertEqual(
                        test_data, current_data,
                        msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
                    )
