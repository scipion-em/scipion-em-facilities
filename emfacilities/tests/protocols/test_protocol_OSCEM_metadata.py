import json
from collections.abc import Mapping, Sequence
from os.path import exists, abspath

import yaml
from sphire.protocols import SphireProtCRYOLOPicking

from cistem.protocols import CistemProtCTFFind

from cryosparc2.protocols import ProtCryo2D, ProtCryoSparcInitialModel, ProtCryoSparcNewNonUniformRefine3D
from pwem import SYM_TETRAHEDRAL

from relion.protocols import ProtRelionAutopickLoG, ProtRelionClassify3D, ProtRelionPostprocess

from pwem.protocols import ProtImportMovies, ProtImportVolumes
from pyworkflow.tests import BaseTest, tests, DataSet
from pyworkflow.utils import magentaStr
from xmipp3.protocols import XmippProtMovieGain, XmippProtFlexAlign, XmippProtMovieMaxShift, XmippProtCTFMicrographs, \
    XmippProtCenterParticles, \
    XmippProtExtractParticles, XmippProtCreateMask3D
from ...protocols import ProtOSCEM
from ...protocols.protocol_OSCEM_metadata import INPUT_MOVIES, INPUT_MICS


box_size = 300
class TestOscemMetadata(BaseTest):
    """
    """
    sampling_rate = 0.495
    max_movie_shift = 20.0
    ctf_down_factor = 2.0
    high_res = 0.5
    test_data = {
        "processing":
            {"movies": {
            "dose_per_image": {
                "value": 0.64,
                "unit": "e/Å²"
            },
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
                        "value": 34.1,
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
                        "value": 29.7,
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
                        "value": 4199.7,
                        "unit": "Å"
                    },
                    "output_max_defocus": {
                        "value": 11828.1,
                        "unit": "Å"
                    },
                    "output_avg_defocus": {
                        "value": 9424.5,
                        "unit": "Å"
                    },
                    "defocus_histogram": "defocus_hist.jpg",
                    "defocus_micrograph_examples": "Micro_examples/micro_defocus.jpg"
                },
                "resolution": {
                    "output_min_resolution": {
                        "value": 3.9,
                        "unit": "Å"
                    },
                    "output_max_resolution": {
                        "value": 2.8,
                        "unit": "Å"
                    },
                    "output_avg_resolution": {
                        "value": 3.3,
                        "unit": "Å"
                    },
                    "resolution_histogram": "resolution_hist.jpg"
                },
                "astigmatism": {
                    "astigmatism_histogram": "astigmatism_hist.jpg"
                }
            },
            "coordinates": {
                "number_particles": 2937,
                "particles_per_micrograph": 139.9,
                "particles_histogram": "particles_hist.jpg",
                "micrograph_examples": "Micro_examples/micro_particles.jpg"
            },
            "classes2D": {
                "number_classes_2D": 20,
                "particles_per_class": [
                    333,
                    296,
                    235,
                    231,
                    228,
                    212,
                    202,
                    192,
                    182,
                    168,
                    96,
                    91,
                    82,
                    80,
                    76,
                    65,
                    62,
                    59,
                    27,
                    20
                ],
                "images_classes_2D": "classes_2D.jpg"
            },
            "classes3D": {
                "number_classes_3D": 1,
                "particles_per_class": [
                    2937
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
                            "side_view": "Classes_3D/isosurface_images_volume1/side_view.jpg",
                            "top_view": "Classes_3D/isosurface_images_volume1/top_view.jpg"
                        }
                    }
                ]
            },
            "volumes": [
                {
                    "volume_type": "initial volume",
                    "size": "(250, 250, 250)",
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
                },
                {
                    "volume_type": "final volume",
                    "number_particles": 2937,
                    "vol_resolution": {
                        "value": 3.39,
                        "unit": "Å"
                    },
                    "size": "(250, 250, 250)",
                    "orthogonal_slices": {
                        "orthogonal_slices_X": "Final_volume/orthogonal_slices/orthogonal_slices_X.jpg",
                        "orthogonal_slices_Y": "Final_volume/orthogonal_slices/orthogonal_slices_Y.jpg",
                        "orthogonal_slices_Z": "Final_volume/orthogonal_slices/orthogonal_slices_Z.jpg"
                    },
                    "isosurface_images": {
                        "front_view": "Final_volume/isosurface_images/front_view.jpg",
                        "side_view": "Final_volume/isosurface_images/side_view.jpg",
                        "top_view": "Final_volume/isosurface_images/top_view.jpg"
                    }
                },
                {
                    "volume_type": "sharpened volume",
                    "size": "(250, 250, 250)",
                    "orthogonal_slices": {
                        "orthogonal_slices_X": "Sharpened_volume/orthogonal_slices/orthogonal_slices_X.jpg",
                        "orthogonal_slices_Y": "Sharpened_volume/orthogonal_slices/orthogonal_slices_Y.jpg",
                        "orthogonal_slices_Z": "Sharpened_volume/orthogonal_slices/orthogonal_slices_Z.jpg"
                    },
                    "isosurface_images": {
                        "front_view": "Sharpened_volume/isosurface_images/front_view.jpg",
                        "side_view": "Sharpened_volume/isosurface_images/side_view.jpg",
                        "top_view": "Sharpened_volume/isosurface_images/top_view.jpg"
                    }
                }
            ]}
}

    @classmethod
    def setUpClass(cls):

        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('OSCEM_metadata')
        cls.protimportmovies, cls.importedmovies = cls.runImportMovies()
        cls.imported_volume_classes3D, cls.imported_initial_volume = cls.runImportVolumes()
        cls.protmoviegain, cls.gainmovies = cls.runMovieGain()
        cls.protalignmovie, cls.alignedmovies = cls.runMovieAlign()
        cls.protmaxshift, cls.maxshiftmicro = cls.runMaxShift()
        cls.protCTF, cls.CTFout = cls.runCTFestimation()
        cls.protpicking, cls.coordinates = cls.runPicking()
        cls.protextract, cls.particles = cls.runExtractParticles()
        cls.prot2Dclasses, cls.classes2D = cls.run2DClassification()
        cls.protinitvol, cls.classesinitvol, cls.volinitvol = cls.runInitialVolume()
        cls.prot3DClassification, cls.classes3DClassification, cls.volumes3DClassification = cls.run3DClassification()
        cls.protrefine, cls.volrefine, cls.particlesrefine = cls.runRefinement()
        cls.prot3Dmask, cls.mask = cls.run3DMask()
        cls.protpostprocess, cls.volpostprocess = cls.runPostProcess()

    @classmethod
    def runImportMovies(cls):

        prot = cls.newProtocol(ProtImportMovies,
                               filesPath=cls.dataset.getFile('movies_dir'),
                               filesPattern='*.tiff',
                               samplingRate=cls.sampling_rate,
                               dosePerFrame=0.64,
                               gainFile=cls.dataset.getFile('gain_im'))

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMovies', None)
        return prot, output

    @classmethod
    def runImportVolumes(cls):

        prot1 = cls.newProtocol(ProtImportVolumes,
                               filesPath=cls.dataset.getFile('volume_classification3D'),
                               samplingRate=0.99)

        cls.launchProtocol(prot1)
        output1 = getattr(prot1, 'outputVolume', None)

        prot2 = cls.newProtocol(ProtImportVolumes,
                                filesPath=cls.dataset.getFile('volume_init'),
                                samplingRate=0.99)

        cls.launchProtocol(prot2)
        output2 = getattr(prot2, 'outputVolume', None)
        return output1, output2

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
        prot = cls.newProtocol(CistemProtCTFFind,
                               inputMicrographs=cls.maxshiftmicro,
                               windowSize=1024,
                               lowRes=16.5,
                               highRes=1.98)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputCTF', None)
        return prot, output

    @classmethod
    def runPicking(cls):
        prot = cls.newProtocol(SphireProtCRYOLOPicking,
                               inputMicrographs=cls.maxshiftmicro,
                               lowPassFilter=True)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputCoordinates', None)
        return prot, output

    @classmethod
    def runExtractParticles(cls):
        prot = cls.newProtocol(XmippProtExtractParticles,
                               inputCoordinates=cls.coordinates,
                               ctfRelations=cls.CTFout,
                               doResize=True,
                               downFactor=2,
                               boxSize=500)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputParticles', None)
        return prot, output

    @classmethod
    def run2DClassification(cls):
        prot = cls.newProtocol(ProtCryo2D,
                               inputParticles=cls.particles,
                               numberOfClasses=20,
                               class2D_window_inner_A=140,
                               class2D_window_outer_A=150)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputClasses', None)
        return prot, output

    @classmethod
    def runInitialVolume(cls):
        prot = cls.newProtocol(ProtCryoSparcInitialModel,
                               inputParticles=cls.particles,
                               symmetryGroup=SYM_TETRAHEDRAL,
                               numberOfClasses=1,
                               abinit_max_res=20,
                               abinit_init_res=35,
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
                               inputParticles=cls.particles,
                               maskDiameterA=150,
                               referenceVolume=cls.volinitvol,  #  initVolVolumes,
                               numberOfIterations=10,
                               symmetryGroup='O',
                               numberOfClasses=1,
                               useBlush=True,
                               limitResolEStep=7)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputClasses', None)
        output2 = getattr(prot, 'outputVolumes', None)
        return prot, output1, output2

    @classmethod
    def runRefinement(cls):
        prot = cls.newProtocol(ProtCryoSparcNewNonUniformRefine3D,
                               inputParticles=cls.particles,
                               referenceVolume=cls.imported_volume_classes3D,
                               symmetryGroup=SYM_TETRAHEDRAL)

        cls.launchProtocol(prot)
        output1 = getattr(prot, 'outputVolume', None)
        output2 = getattr(prot, 'outputParticles', None)
        return prot, output1, output2

    @classmethod
    def run3DMask(cls):
        prot = cls.newProtocol(XmippProtCreateMask3D,
                               inputVolume=cls.volrefine,
                               threshold=0.1,
                               doMorphological=True,
                               elementSize=5,
                               doSmooth=True)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputMask', None)
        return prot, output

    @classmethod
    def runPostProcess(cls):
        prot = cls.newProtocol(ProtRelionPostprocess,
                               relionInput=False,
                               inputVolume=cls.volrefine,
                               solventMask=cls.mask,
                               skipFscWeighting=True)

        cls.launchProtocol(prot)
        output = getattr(prot, 'outputVolume', None)
        return prot, output

    def test_complete_input(self):
        print(magentaStr("\n==> Running test with all input completed:"))
        test_data_import = {"processing": {
                "movies": self.test_data["processing"]["movies"],
                "micrographs": self.test_data["processing"]["micrographs"],
                "CTFs": self.test_data["processing"]["CTFs"],
                "coordinates": self.test_data["processing"]["coordinates"],
                "classes2D": self.test_data["processing"]["classes2D"],
                "volumes": self.test_data["processing"]["volumes"]}
                }

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift,
                                micrographs=self.maxshiftmicro,
                                CTF=self.CTFout,
                                particles=self.particles,
                                classes2D=self.classes2D,
                                initVolume=self.imported_initial_volume,  # Pointer(self.initVolVolumes.getFirstItem()),
                                classes3D=self.classes3DClassification,
                                finalVolume=self.volrefine,
                                sharpenedVolume=self.volpostprocess)

        load_metadata = self.prot_and_load_metadata(prot)

        for key, section_test_dict in test_data_import.items():
            current_dict = load_metadata.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)


    def test_medium_level(self):
        print(magentaStr("\n==> Running test with compulsory and some optional input:"))
        test_data_import = {"processing": {
                            "movies": self.test_data["processing"]["movies"],
                            "CTFs": self.test_data["processing"]["CTFs"]}
                            }

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MOVIES,
                                importMovies=self.protimportmovies,
                                movieAlignment=self.protalignmovie,
                                maxShift=self.protmaxshift,
                                CTF=self.CTFout)

        load_metadata = self.prot_and_load_metadata(prot)

        # Recursive comparison
        for key, section_test_dict in test_data_import.items():
            current_dict = load_metadata.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)


    def test_micro_input(self):
        print(magentaStr("\n==> Running test with micrographs as input:"))
        test_data_import = {"processing": {
                            "CTFs": self.test_data["processing"]["CTFs"]}
                            }

        prot = self.newProtocol(ProtOSCEM,
                                inputType=INPUT_MICS,
                                CTF=self.CTFout)

        load_metadata = self.prot_and_load_metadata(prot)

        # Recursive comparison
        for key, section_test_dict in test_data_import.items():
            current_dict = load_metadata.get(key, None)
            self.assertIsNotNone(current_dict, msg=f'Dictionary section {key} is not found')
            self.recursive_compare(section_test_dict, current_dict, parent_key=key)

    def prot_and_load_metadata(self, prot):
        self.launchProtocol(prot)
        file_path = prot.getOutFile()
        self.assertTrue(exists(file_path))
        # with open(abspath(file_path), 'r') as json_file:
        #     load_json = json.load(json_file)
        # Open the YAML file
        with open(abspath(file_path), 'r', encoding='utf-8') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)
            metadata_string = json.dumps(yaml_data, ensure_ascii=False)
            load_metadata = json.loads(metadata_string)

        return load_metadata

    def recursive_compare(self, test_data, current_data, parent_key=""):
        """
        Compares nested dictionaries and lists
        """
        if isinstance(test_data, dict):
            self.compare_dicts(test_data, current_data, parent_key)

        elif isinstance(test_data, list):
            self.compare_lists(test_data, current_data, parent_key)

        else:
            # Check if we are dealing with a value-unit pair
            if parent_key.endswith(".value") or parent_key.endswith(".unit"):
                self.compare_value_unit(test_data, current_data, parent_key)

            else:
                self.compare_rest(test_data, current_data, parent_key)

    def compare_dicts(self, test_data, current_data, parent_key):
        """
            Compares two lists dictionaries and raises assertion errors if a key is not found in it,
            then it performs recursive comparison.
            """
        for key, test_value in test_data.items():
            self.assertIn(key, current_data, msg=f"Key '{key}' not found in {parent_key}")
            # Recursively call to compare nested values
            self.recursive_compare(
                test_value, current_data[key], f"{parent_key}.{key}" if parent_key else key
            )

    def compare_lists(self, test_data, current_data, parent_key):
        """
            Compares two lists (test_data and current_data) and raises assertion errors if there are mismatches.
            If the list ends with 'particles_per_class', it checks that the length difference is <= 2
            and values match within a delta of 500. Otherwise, it recursively compares the lists element by element.
            """
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

    def compare_value_unit(self, test_data, current_data, parent_key):
        """
            Compares value-unit pairs (test_data and current_data) and raises assertion errors if there are mismatches.
            """
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

    def compare_rest(self, test_data, current_data, parent_key):
        """
            Compares rest of key-values (test_data and current_data) and raises assertion errors if there are mismatches.
            """
        key_in = parent_key.split('.')[-1]
        if (
                key_in == "number_micrographs"
                or key_in == "discarded_movies"
                or key_in == "number_classes_2D"
                or key_in == "number_classes_3D"
        ):
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
                test_data, current_data, delta=15,
                msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
            )
        elif key_in == "resolution":
            self.assertAlmostEqual(
                test_data, current_data, delta=0.5,
                msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
            )
        else:
            self.assertEqual(
                test_data, current_data,
                msg=f"Value mismatch at {parent_key}: {test_data} != {current_data}"
            )