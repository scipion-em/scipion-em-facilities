import sqlite3
from functools import reduce
from os.path import abspath, join, dirname, splitext

from PIL import Image, ImageDraw, ImageFont, ImageOps, ImageEnhance
from matplotlib import pyplot as plt
import mrcfile

import pyworkflow.protocol.params as params
import xmipp3
from pwem.objects import Class2D, Class3D
from pwem.protocols import EMProtocol
import pyworkflow.utils as pwutils

import numpy as np
import json
import yaml
import os

INPUT_MOVIES = 0
INPUT_MICS = 1
# OUTFILE = 'Processing_metadata.json'
OUTFILE = 'Processing_metadata.yaml'

# input movies attributes:
_dosePerFrame = 'outputMovies._acquisition._dosePerFrame'
_doseInitial = 'outputMovies._acquisition._doseInitial'
_gainFile = 'outputMovies._gainFile'
_darkFile = 'outputMovies._darkFile'
_size = 'outputMovies._size'

# y label of micrographs
hist_ylabel_mic = 'Frequency of Micrographs'
hist_ylabel_frames = 'Frequency of frames'

class ProtOSCEM(EMProtocol):
    """ This is the class for generating the OSCEM metadata json file from Scipion workflow
    """
    _label = 'OSCEM Metadata'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        CONDITION_MOVIES = f'inputType=={INPUT_MOVIES}'

        form.addSection(label='Input')

        form.addParam('inputType', params.EnumParam, label='Input type',
                      choices=['Movies', 'Micrographs'], important=True,
                      help='Select the type of input, either movies or micrographs',
                      default=INPUT_MOVIES)

        form.addParam('importMovies', params.PointerParam,
                      label="Import Movies", important=True,
                      pointerClass='ProtImportMovies',
                      help="Import Movies",
                      condition=CONDITION_MOVIES)

        form.addParam('movieAlignment', params.PointerParam,
                      label="Movie alignment",
                      pointerClass='ProtAlignMovies',
                      help="Movie alignment used",
                      condition=CONDITION_MOVIES,
                      allowsNull=True)

        form.addParam('maxShift', params.PointerParam,
                      label="Max Shift",
                      pointerClass='XmippProtMovieMaxShift',
                      help="Max Shift",
                      condition=CONDITION_MOVIES,
                      allowsNull=True)

        form.addParam('micrographs', params.PointerParam,
                      label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      help="Micrographs",
                      allowsNull=True)

        form.addParam('CTF', params.PointerParam,
                      label="CTF",
                      pointerClass='SetOfCTF',
                      help="CTF micrographs",
                      allowsNull=True)

        form.addParam('particles', params.PointerParam,
                      label="Particles",
                      pointerClass='SetOfParticles',
                      help="Particles obtained when doing particle extraction",
                      allowsNull=True)

        form.addParam('classes2D', params.PointerParam,
                      label="Classes 2D",
                      pointerClass='SetOfClasses2D',
                      help="Set of 2D classes",
                      allowsNull=True)

        form.addParam('initVolume', params.PointerParam,
                      label="Initial volume",
                      pointerClass='Volume',
                      help="Initial volume",
                      allowsNull=True)

        form.addParam('threshold_initVol', params.IntParam, default=-1,
                      label="Initial volume threshold",
                      help='Threshold to obtain isosurface of volume')

        form.addParam('classes3D', params.PointerParam,
                      label="Classes 3D",
                      pointerClass='SetOfClasses3D',
                      help="Set of 3D classes",
                      allowsNull=True)

        form.addParam('threshold_classes3D', params.IntParam, default=-1,
                      label="3D Classes threshold",
                      help='Threshold to obtain isosurface of classes 3D')

        form.addParam('finalVolume', params.PointerParam,
                      label="Final volume",
                      pointerClass='Volume',
                      help="Final volume",
                      allowsNull=True)

        form.addParam('threshold_finalVol', params.IntParam, default=-1,
                      label="Final volume threshold",
                      help='Threshold to obtain isosurface of volume')

        form.addParam('polishedVolume', params.PointerParam,
                      label="Polished volume",
                      pointerClass='Volume',
                      help="Polished volume",
                      allowsNull=True)

        form.addParam('threshold_polishedVol', params.IntParam, default=-1,
                      label="Polished volume threshold",
                      help='Threshold to obtain isosurface of volume')

        form.addParam('enhancedVolume', params.PointerParam,
                      label="Enhanced volume",
                      pointerClass='Volume',
                      help="Enhanced volume",
                      allowsNull=True)

        form.addParam('threshold_enhancedVol', params.IntParam, default=-1,
                      label="Enhanced volume threshold",
                      help='Threshold to obtain isosurface of volume')



    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.generateJson)
        self._insertFunctionStep(self.saveJson)

    # -------------------------- STEPS functions ------------------------------
    def generateJson(self):

        self.processing_json = {}

        if self.inputType.get() == 0:  # movies as input
            ###### IMPORT MOVIES ######
            import_movies = self.import_movies_generation()
            self.processing_json['Movies'] = import_movies

            # if self.movieAlignment.get() is not None:
            #     ###### MOVIE ALGINMENT ######
            #     movie_alignment = self.movie_alignment_generation()
            #     self.processing_json['Movie_alignment'] = movie_alignment
            #
            # if self.maxShift.get() is not None:
            #     ###### MAX SHIFT ######
            #     max_shift = self.max_shift_generation()
            #     self.processing_json['Movie_maxshift'] = max_shift

        if self.micrographs.get() is not None:
            ###### MICROGRAPHS ######
            micrographs = self.micrographs_generation()
            self.processing_json['Micrographs'] = micrographs

        if self.CTF.get() is not None:
            ###### CTF ######
            CTF = self.CTF_generation()
            self.processing_json['CTFs'] = CTF

        if self.particles.get() is not None:
            ###### PARTICLES ######
            particles = self.particles_generation()
            self.processing_json['Coordinates'] = particles

        if self.classes2D.get() is not None:
            ###### CLASSES 2D ######
            classes_2D = self.classes2D_generation()
            self.processing_json['Classes2D'] = classes_2D

        if self.classes3D.get() is not None:
            ###### CLASSES 3D ######
            classes_3D = self.classes3D_generation()
            self.processing_json['Classes_3D'] = classes_3D

        volumes = []
        if self.initVolume.get() is not None:
            ###### INITIAL VOLUME ######
            volume_type = 'initial volume'
            folder_name = 'Initial_volume'
            volume = self.initVolume.get()
            th = int(self.threshold_initVol.get())
            init_volume = self.volume_generation(volume_type, folder_name, volume, th)
            volumes.append(init_volume)
            # self.processing_json['Initial_volume'] = volume
            self.processing_json['Volumes'] = volumes

        if self.finalVolume.get() is not None:
            ###### FINAL VOLUME ######
            volume_type = 'final volume'
            folder_name = 'Final_volume'
            volume = self.finalVolume.get()
            th = int(self.threshold_finalVol.get())
            final_volume = self.volume_generation(volume_type, folder_name, volume, th)
            volumes.append(final_volume)
            self.processing_json['Volumes'] = volumes

        if self.enhancedVolume.get() is not None:
            ###### ENHANCED VOLUME ######
            volume_type = 'enhanced volume'
            folder_name = 'Enhanced_volume'
            volume = self.enhancedVolume.get()
            th = int(self.threshold_enhancedVol.get())
            enhanced_volume = self.volume_generation(volume_type, folder_name, volume, th)
            volumes.append(enhanced_volume)
            self.processing_json['Volumes'] = volumes

        if self.polishedVolume.get() is not None:
            ###### POLISHED VOLUME ######
            volume_type = 'polished volume'
            folder_name = 'Polished_volume'
            volume = self.polishedVolume.get()
            th = int(self.threshold_polishedVol.get())
            polished_volume = self.volume_generation(volume_type, folder_name, volume, th)
            volumes.append(polished_volume)
            self.processing_json['Volumes'] = volumes

        print(json.dumps(self.processing_json, ensure_ascii=False, indent=4))

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []

    # -------------------- UTILS functions -------------------------
    def import_movies_generation(self):
        ImportMoviesProt = self.importMovies.get()
        input_movies = ImportMoviesProt.getObjDict()

        # List of keys to retrieve
        # if dosePerFrame has a value, then dose initial is also retrieved
        # Otherwise, none of them are retrieved.
        if input_movies[_dosePerFrame] is None:
            keys_to_retrieve = [_gainFile, _darkFile]
        else:
            keys_to_retrieve = [_dosePerFrame, _doseInitial, _gainFile, _darkFile]

        # Mapping dictionary for key name changes
        key_mapping = {
            _dosePerFrame: ('Dose_per_image', 'e/Å²'),
            _doseInitial: ('Initial_dose', 'e/Å²'),
            _gainFile: 'Gain_image',
            _darkFile: 'Dark_image'
        }

        # Filter the dictionary and rename the keys
        extra_folder = self._getExtraPath()
        file_keys = [_gainFile, _darkFile]
        for key in file_keys:
            if input_movies[key] is not None:
                with mrcfile.open(input_movies[key], 'r') as mrc:
                    # Read the data from the MRC file
                    data = mrc.data

                # Normalize the data to 8-bit (0-255) range
                min_val = np.min(data)
                max_val = np.max(data)
                normalized_data = 255 * (data - min_val) / (max_val - min_val)
                normalized_data = normalized_data.astype(np.uint8)

                # Apply Histogram Equalization
                # Convert to PIL Image
                image = Image.fromarray(normalized_data)
                image = ImageOps.equalize(image)

                # Save the image as PNG to the specified path
                # Extract the base filename without extension
                base_filename = os.path.splitext(os.path.basename(input_movies[key]))[0]

                png_path = os.path.join(extra_folder, f'{base_filename}.jpg')
                image.save(png_path)

                image_path = os.path.basename(png_path)
                input_movies[key] = os.path.basename(image_path)

        # Creating the dictionary
        import_movies = {}
        for key in keys_to_retrieve:
            if key in input_movies and input_movies[key] is not None and input_movies[key] != 0:
                if key in key_mapping:
                    mapped_key = key_mapping[key]
                    if isinstance(mapped_key, tuple):  # Value-unit pair
                        field, unit = mapped_key
                        import_movies[field] = {
                            "value": input_movies[key],
                            "unit": unit
                        }
                    else:
                        import_movies[mapped_key] = input_movies[key]


        self.number_movies = input_movies[_size]

        # Descriptors
        descriptors = []
        if self.movieAlignment.get() is not None:
            ###### MOVIE ALGINMENT DESCRIPTOR ######

            MovieAlignmentProt = self.movieAlignment.get()
            movie_align = {'Descriptor_name': MovieAlignmentProt.getClassName()}

            ################################ INPUT #############################################
            input_alignment = MovieAlignmentProt.getObjDict()
            # List of keys to retrieve
            keys_to_retrieve = ['binFactor', 'maxResForCorrelation', 'gainRot', 'gainFlip']

            # Mapping dictionary for key name changes
            key_mapping = {
                'binFactor': 'Binning_factor',
                'maxResForCorrelation': 'Maximum_resolution_(Å)',
                'gainRot': 'Rotate_gain_reference',
                'gainFlip': 'Flip_gain_reference'
            }

            # Map values for the schema
            input_movie_align = {}
            for key in keys_to_retrieve:
                if key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0:
                    if key == 'maxResForCorrelation':
                        input_movie_align['Maximum_resolution'] = {
                            'value': input_alignment[key],
                            'unit': 'Å'
                        }
                    else:
                        input_movie_align[key_mapping[key]] = input_alignment[key]

            movie_align.update(input_movie_align)

            # Dictionary for crop offsets
            crop_offset_mapping = {
                'cropOffsetX': 'Crop_offsetX',
                'cropOffsetY': 'Crop_offsetY',
            }

            crop_offsets = {}
            for key, mapped_key in crop_offset_mapping.items():
                if key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0:
                    crop_offsets[mapped_key] = {
                        'value': input_alignment[key],
                        'unit': 'pixels'
                    }
            if crop_offsets:
                movie_align['Crop_offsets'] = crop_offsets

            # Dictionary for crop dimensions
            crop_dim_mapping = {
                'cropDimX': 'Crop_dimsX',
                'cropDimY': 'Crop_dimsY',
            }
            crop_dims = {}
            for key, mapped_key in crop_dim_mapping.items():
                if key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0:
                    crop_dims[mapped_key] = {
                        'value': input_alignment[key],
                        'unit': 'pixels'
                    }
            if crop_dims:
                movie_align['Crop_dims'] = crop_dims

            # Dictionary for frames aligned
            # if input_alignment['alignFrameN'] != 0:
            keys_to_retrieve = ['alignFrame0', 'alignFrameN']
            key_mapping = {
                'alignFrame0': 'Frame0',
                'alignFrameN': 'FrameN',
            }
            frames_aligned = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
                              key in input_alignment and input_alignment[key] is not None}

            if frames_aligned['FrameN'] == 0:
                frames_aligned['FrameN'] = self.number_movies

            movie_align['Frames_aligned'] = frames_aligned

            ############################### OUTPUT #############################################
            # average and max shift
            for a, output in MovieAlignmentProt.iterOutputAttributes():
                if a == 'outputMovies':
                    for index, item in enumerate(output.iterItems()):
                        attributes = item.getAttributes()
                        attributes_dict = dict(attributes)
                        shiftX = attributes_dict.get('_xmipp_ShiftX')
                        shiftY = attributes_dict.get('_xmipp_ShiftY')
                        norm = np.linalg.norm([shiftX, shiftY], axis=0)

                        # Max and Average shift
                        max_norm = np.max(norm)
                        avgXY = np.mean(norm)
                        if index == 0:
                            avg_shift = avgXY
                            max_shift = max_norm
                        else:
                            avg_shift = np.mean([avgXY, avg_shift])
                            max_shift = max(max_shift, max_norm)

                # movie_align.update(output_movie_align)
                movie_align['Output_avg_shift'] = {
                    'value': round(avg_shift, 1),
                    'unit': 'Å'
                }
                movie_align['Output_max_shift'] = {
                    'value': round(max_shift, 1),
                    'unit': 'Å'
                }

            # Append movie_align to descriptors
            descriptors.append(movie_align)

        if self.maxShift.get() is not None:
            ###### MAX SHIFT DESCRIPTOR######
            MaxShiftProt = self.maxShift.get()
            movie_maxshift = {'Descriptor_name': MaxShiftProt.getClassName()}

            ############################### INPUT #############################################
            input_shift = MaxShiftProt.getObjDict()
            # List of keys to retrieve
            keys_to_retrieve = ['outputMoviesDiscarded._size', 'maxFrameShift', 'maxMovieShift', 'rejType']
            # Mapping dictionary for key name changes
            key_mapping = {
                'outputMoviesDiscarded._size': 'Discarded_movies',
                'maxFrameShift': 'Max_frame_shift',
                'maxMovieShift': 'Max_movie_shift',
                'rejType': 'Rejection_type'
            }

            # Filter dictionary and rename keys
            # movie_maxshift = {}
            rejtype_list = ['By frame', 'By whole movie', 'By frame and movie', 'By frame or movie']

            for key in keys_to_retrieve:
                if key == 'rejType':
                    rej_type = rejtype_list[input_shift[key]]
                    movie_maxshift[key_mapping[key]] = rej_type
                elif key in input_shift and input_shift[key] is not None and input_shift[key] != 0:
                    if key in ['maxFrameShift', 'maxMovieShift']:
                        movie_maxshift[key_mapping[key]] = {
                            'value': input_shift[key],
                            'unit': 'Å'
                        }
                    else:
                        movie_maxshift[key_mapping[key]] = input_shift[key]

            ############################### OUTPUT #############################################
            # average and max shift
            shift_list = []
            for a, output in MaxShiftProt.iterOutputAttributes():
                if a == 'outputMovies':
                    for index, item in enumerate(output.iterItems()):
                        attributes = item.getAttributes()
                        attributes_dict = dict(attributes)
                        shiftX = attributes_dict.get('_xmipp_ShiftX')
                        shiftY = attributes_dict.get('_xmipp_ShiftY')
                        norm = np.linalg.norm([shiftX, shiftY], axis=0)

                        # Shift
                        shift_list.append(norm)

                        # Max and average shift
                        max_norm = np.max(norm)
                        avgXY = np.mean(norm)
                        if index == 0:
                            avg_shift = avgXY
                            max_shift = max_norm
                        else:
                            avg_shift = np.mean([avgXY, avg_shift])
                            max_shift = max(max_shift, max_norm)

            flattened_shift_list = np.hstack(shift_list)

            # Histogram generation
            # shift
            plt.close('all')
            plt.clf()
            plt.cla()

            plt.hist(flattened_shift_list, edgecolor='black')
            plt.xlabel('# Shift (Å)')
            plt.ylabel(hist_ylabel_frames)
            plt.title('Shift histogram')
            shift_hist_name = 'shift_hist.jpg'
            shift_hist = self.hist_path(shift_hist_name)
            plt.savefig(shift_hist)

            output_movie_maxshift = {
                'Output_avg_shift': {
                    'value': round(avg_shift, 1),
                    'unit': 'Å'
                },
                'Output_max_shift': {
                    'value': round(max_shift, 1),
                    'unit': 'Å'
                },
                'Shift_histogram': shift_hist_name
            }

            movie_maxshift.update(output_movie_maxshift)

            # Append movie_align to descriptors
            descriptors.append(movie_maxshift)

        import_movies['descriptors'] = descriptors
        return import_movies

    # def movie_alignment_generation(self):
    #     MovieAlignmentProt = self.movieAlignment.get()
    #     movie_align = {'Method': MovieAlignmentProt.getClassName()}
    #
    #     ################################ INPUT #############################################
    #     input_alignment = MovieAlignmentProt.getObjDict()
    #     # List of keys to retrieve
    #     keys_to_retrieve = ['binFactor', 'maxResForCorrelation', 'gainRot', 'gainFlip']
    #
    #     # Mapping dictionary for key name changes
    #     key_mapping = {
    #         'binFactor': 'Binning_factor',
    #         'maxResForCorrelation': 'Maximum_resolution_(Å)',
    #         'gainRot': 'Rotate_gain_reference',
    #         'gainFlip': 'Flip_gain_reference'
    #     }
    #
    #     # Map values for the schema
    #     input_movie_align = {}
    #     for key in keys_to_retrieve:
    #         if key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0:
    #             if key == 'maxResForCorrelation':
    #                 input_movie_align['Maximum_resolution'] = {
    #                     'value': input_alignment[key],
    #                     'unit': 'Å'
    #                 }
    #             else:
    #                 input_movie_align[key_mapping[key]] = input_alignment[key]
    #
    #     movie_align.update(input_movie_align)
    #
    #     # Dictionary for crop offsets
    #     crop_offset_mapping = {
    #         'cropOffsetX': 'Crop_offsetX',
    #         'cropOffsetY': 'Crop_offsetY',
    #     }
    #
    #     crop_offsets = {}
    #     for key, mapped_key in crop_offset_mapping.items():
    #         if key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0:
    #             crop_offsets[mapped_key] = {
    #                 'value': input_alignment[key],
    #                 'unit': 'pixels'
    #             }
    #     if crop_offsets:
    #         movie_align['Crop_offsets'] = crop_offsets
    #
    #
    #     # Dictionary for crop dimensions
    #     crop_dim_mapping = {
    #         'cropDimX': 'Crop_dimsX',
    #         'cropDimY': 'Crop_dimsY',
    #     }
    #     crop_dims = {}
    #     for key, mapped_key in crop_dim_mapping.items():
    #         if key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0:
    #             crop_dims[mapped_key] = {
    #                 'value': input_alignment[key],
    #                 'unit': 'pixels'
    #             }
    #     if crop_dims:
    #         movie_align['Crop_dims'] = crop_dims
    #
    #
    #     # Dictionary for frames aligned
    #     # if input_alignment['alignFrameN'] != 0:
    #     keys_to_retrieve = ['alignFrame0', 'alignFrameN']
    #     key_mapping = {
    #         'alignFrame0': 'Frame0',
    #         'alignFrameN': 'FrameN',
    #     }
    #     frames_aligned = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
    #                       key in input_alignment and input_alignment[key] is not None}
    #
    #     if frames_aligned['FrameN'] == 0:
    #         frames_aligned['FrameN'] = self.number_movies
    #
    #     movie_align['Frames_aligned'] = frames_aligned
    #
    #     ############################### OUTPUT #############################################
    #     # average and max shift
    #     for a, output in MovieAlignmentProt.iterOutputAttributes():
    #         if a == 'outputMovies':
    #             for index, item in enumerate(output.iterItems()):
    #                 attributes = item.getAttributes()
    #                 attributes_dict = dict(attributes)
    #                 shiftX = attributes_dict.get('_xmipp_ShiftX')
    #                 shiftY = attributes_dict.get('_xmipp_ShiftY')
    #                 norm = np.linalg.norm([shiftX, shiftY], axis=0)
    #
    #                 # Max and Average shift
    #                 max_norm = np.max(norm)
    #                 avgXY = np.mean(norm)
    #                 if index == 0:
    #                     avg_shift = avgXY
    #                     max_shift = max_norm
    #                 else:
    #                     avg_shift = np.mean([avgXY, avg_shift])
    #                     max_shift = max(max_shift, max_norm)
    #
    #         # movie_align.update(output_movie_align)
    #         movie_align['Output_avg_shift'] = {
    #             'value': round(avg_shift, 1),
    #             'unit': 'Å'
    #         }
    #         movie_align['Output_max_shift'] = {
    #             'value': round(max_shift, 1),
    #             'unit': 'Å'
    #         }
    #
    #     return movie_align

    # def max_shift_generation(self):
    #     MaxShiftProt = self.maxShift.get()
    #
    #     ############################### INPUT #############################################
    #     input_shift = MaxShiftProt.getObjDict()
    #     # List of keys to retrieve
    #     keys_to_retrieve = ['outputMoviesDiscarded._size', 'maxFrameShift', 'maxMovieShift', 'rejType']
    #     # Mapping dictionary for key name changes
    #     key_mapping = {
    #         'outputMoviesDiscarded._size': 'Discarded_movies',
    #         'maxFrameShift': 'Max_frame_shift',
    #         'maxMovieShift': 'Max_movie_shift',
    #         'rejType': 'Rejection_type'
    #     }
    #
    #     # Filter dictionary and rename keys
    #     movie_maxshift = {}
    #     rejtype_list = ['By frame', 'By whole movie', 'By frame and movie', 'By frame or movie']
    #
    #     for key in keys_to_retrieve:
    #         if key == 'rejType':
    #             rej_type = rejtype_list[input_shift[key]]
    #             movie_maxshift[key_mapping[key]] = rej_type
    #         elif key in input_shift and input_shift[key] is not None and input_shift[key] != 0:
    #             if key in ['maxFrameShift', 'maxMovieShift']:
    #                 movie_maxshift[key_mapping[key]] = {
    #                     'value': input_shift[key],
    #                     'unit': 'Å'
    #                 }
    #             else:
    #                 movie_maxshift[key_mapping[key]] = input_shift[key]
    #
    #
    #     ############################### OUTPUT #############################################
    #     # average and max shift
    #     shift_list = []
    #     for a, output in MaxShiftProt.iterOutputAttributes():
    #         if a == 'outputMovies':
    #             for index, item in enumerate(output.iterItems()):
    #                 attributes = item.getAttributes()
    #                 attributes_dict = dict(attributes)
    #                 shiftX = attributes_dict.get('_xmipp_ShiftX')
    #                 shiftY = attributes_dict.get('_xmipp_ShiftY')
    #                 norm = np.linalg.norm([shiftX, shiftY], axis=0)
    #
    #                 # Shift
    #                 shift_list.append(norm)
    #
    #                 # Max and average shift
    #                 max_norm = np.max(norm)
    #                 avgXY = np.mean(norm)
    #                 if index == 0:
    #                     avg_shift = avgXY
    #                     max_shift = max_norm
    #                 else:
    #                     avg_shift = np.mean([avgXY, avg_shift])
    #                     max_shift = max(max_shift, max_norm)
    #
    #     flattened_shift_list = np.hstack(shift_list)
    #
    #     # Histogram generation
    #     # shift
    #     plt.close('all')
    #     plt.clf()
    #     plt.cla()
    #
    #
    #     plt.hist(flattened_shift_list, edgecolor='black')
    #     plt.xlabel('# Shift (Å)')
    #     plt.ylabel(hist_ylabel_frames)
    #     plt.title('Shift histogram')
    #     shift_hist_name = 'shift_hist.jpg'
    #     shift_hist = self.hist_path(shift_hist_name)
    #     plt.savefig(shift_hist)
    #
    #     output_movie_maxshift = {
    #         'Output_avg_shift': {
    #             'value': round(avg_shift, 1),
    #             'unit': 'Å'
    #         },
    #         'Output_max_shift': {
    #             'value': round(max_shift, 1),
    #             'unit': 'Å'
    #         },
    #         'Shift_histogram': shift_hist_name
    #     }
    #
    #     movie_maxshift.update(output_movie_maxshift)
    #
    #     return movie_maxshift

    def micrographs_generation(self):
        mics = self.micrographs.get()
        size = mics.getSize()

        micrographs = {
            'number_micrographs': size
        }

        return micrographs

    def CTF_generation(self):
        CTFs = self.CTF.get()
        ############################## OUTPUT #############################################
        CTF_estimation = {}
        defocus_list = []
        resolution_list = []
        astigmatism_list = []
        # dictionary to show 3 micrographs (minimum defocus, medium defocus, max defocus)
        dict_defocus = {}

        for index, item in enumerate(CTFs.iterItems()):
            amplitude_contrast = float(item._micObj._acquisition._amplitudeContrast)
            # Min, max  and average defocus and resolution
            defocus = np.mean([float(item._defocusU), float(item._defocusV)])
            resolution = float(item._resolution)
            astigmatism = float(item._defocusRatio)

            defocus_list.append(defocus)
            resolution_list.append(resolution)
            astigmatism_list.append(astigmatism)

            # Dictionary to save defocus and paths to save 3 micrographs (min defocus, max defocus and medium defocus)
            dict_defocus[f'mic_{index}'] = {'defocus': defocus, 'path_to_mic': item._micObj._filename.get()}

            if index == 0:
                max_defocus = defocus
                min_defocus = defocus
                max_resolution = resolution
                min_resolution = resolution
                avg_defocus = defocus
                avg_resolution = resolution
            else:
                max_defocus = max(max_defocus, defocus)
                min_defocus = min(min_defocus, defocus)
                max_resolution = max(max_resolution, resolution)
                min_resolution = min(min_resolution, resolution)
                avg_defocus = np.mean([avg_defocus, defocus])
                avg_resolution = np.mean([avg_resolution, resolution])

        # Sort the dictionary by defocus values
        sorted_micrographs = sorted(dict_defocus.items(), key=lambda x: x[1]['defocus'])
        # Extract paths for lowest, medium, and highest defocus values
        medium_defocus_mic = min(sorted_micrographs, key=lambda x: abs(x[1]['defocus'] - avg_defocus))
        lowest_defocus_path = sorted_micrographs[0][1]['path_to_mic']
        highest_defocus_path = sorted_micrographs[-1][1]['path_to_mic']
        medium_defocus_path = medium_defocus_mic[1]['path_to_mic']
        medium_defocus_value = medium_defocus_mic[1]['defocus']

        # To draw defocus value on image, we make images smaller and then larger
        # since the default font cannot be increased in size
        with mrcfile.open(lowest_defocus_path, permissive=True) as mrc1:
            data1 = mrc1.data
            img1 = Image.fromarray(np.uint8(data1 / np.max(data1) * 255))
            img1 = img1.convert('RGB')
            img1_resized = self.resize_image(img1, scale=0.05)
        with mrcfile.open(medium_defocus_path, permissive=True) as mrc2:
            data2 = mrc2.data
            img2 = Image.fromarray(np.uint8(data2 / np.max(data2) * 255))
            img2 = img2.convert('RGB')
            img2_resized = self.resize_image(img2, scale=0.05)
        with mrcfile.open(highest_defocus_path, permissive=True) as mrc3:
            data3 = mrc3.data
            img3 = Image.fromarray(np.uint8(data3 / np.max(data3) * 255))
            img3 = img3.convert('RGB')
            img3_resized = self.resize_image(img3, scale=0.05)


        font = ImageFont.load_default()
        # Draw text on images
        draw1 = ImageDraw.Draw(img1_resized)
        text1 = f'{min_defocus} A'
        draw1.text((10, 10), text1, fill='#80FF00', font=font)

        draw2 = ImageDraw.Draw(img2_resized)
        text2 = f'{medium_defocus_value} A'
        draw2.text((10, 10), text2, fill='#80FF00', font=font)

        draw3 = ImageDraw.Draw(img3_resized)
        text3 = f'{max_defocus} A'
        draw3.text((10, 10), text3, fill='#80FF00', font=font)

        # Define original sizes
        img1_original_size = img1.size
        img2_original_size = img2.size
        img3_original_size = img3.size
        # Return to original sizes
        img1_final = self.resize_back(img1_resized, img1_original_size)
        img2_final = self.resize_back(img2_resized, img2_original_size)
        img3_final = self.resize_back(img3_resized, img3_original_size)

        # Create collage with the 3 images
        collage_width = img1_final.width + img2_final.width + img3_final.width
        collage_height = max(img1_final.height, img2_final.height, img3_final.height)
        collage = Image.new('RGB', (collage_width, collage_height))
        # Paste images into the collage
        collage.paste(img1_final, (0, 0))
        collage.paste(img2_final, (img1_final.width, 0))
        collage.paste(img3_final, (img1_final.width + img2_final.width, 0))

        # Save the result
        extra_folder = self._getExtraPath()
        micro_folder_name = 'Micro_examples'
        micro_folder_path = join(extra_folder, micro_folder_name)
        os.makedirs(micro_folder_path, exist_ok=True)

        micro_name = 'micro_defocus.jpg'
        micro_path = join(micro_folder_path, micro_name)
        collage.save(micro_path)

        # Histograms generation
        # DEFOCUS
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(defocus_list, bins='auto', edgecolor='black')
        plt.xlabel('# Defocus')
        plt.ylabel(hist_ylabel_mic)
        plt.title('Defocus histogram')
        defocus_hist_name = 'defocus_hist.jpg'
        defocus_hist = self.hist_path(defocus_hist_name)
        plt.savefig(defocus_hist)

        # RESOLUTION
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(resolution_list, bins='auto', edgecolor='black')
        plt.xlabel("Resolution")
        plt.ylabel(hist_ylabel_mic)
        plt.title('Resolution histogram')
        resolution_hist_name = 'resolution_hist.jpg'
        resolution_hist = self.hist_path(resolution_hist_name)
        plt.savefig(resolution_hist)

        # ASTIGMATISM
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(astigmatism_list, bins='auto', edgecolor='black')
        plt.xlabel("Astigmatism")
        plt.ylabel(hist_ylabel_mic)
        plt.title('Astigmatism histogram')
        astigmatism_hist_name = 'astigmatism_hist.jpg'
        astigmatism_hist = self.hist_path(astigmatism_hist_name)
        plt.savefig(astigmatism_hist)

        CTF_estimation = {
            'Amplitude_contrast': amplitude_contrast,
            'Defocus': {
                'Output_min_defocus': {
                    'value': round(min_defocus, 1),
                    'unit': 'Å'
                },
                'Output_max_defocus': {
                    'value': round(max_defocus, 1),
                    'unit': 'Å'
                },
                'Ouput_avg_defocus': {
                    'value': round(avg_defocus, 1),
                    'unit': 'Å'
                },
                'Defocus_histogram': defocus_hist_name,
                'Micrograph_examples': join(micro_folder_name, micro_name)
            },
            'Resolution': {
                'Output_min_resolution': {
                    'value': round(min_resolution, 1),
                    'unit': 'Å'
                },
                'Output_max_resolution': {
                    'value': round(max_resolution, 1),
                    'unit': 'Å'
                },
                'Ouput_avg_resolution': {
                    'value': round(avg_resolution, 1),
                    'unit': 'Å'
                },
                'Resolution_histogram': resolution_hist_name
            },
            'Astigmatism': {
                'Astigmatism_histogram': astigmatism_hist_name
            }
        }

        return CTF_estimation

    def particles_generation(self):
        parts = self.particles.get()
        mic_dict = {}

        for index, item in enumerate(parts.iterItems()):
            micrograph_num = str(item._micId)
            sampling_rate_part = item._samplingRate.get()
            sampling_rate_ctf = item._ctfModel._xmipp_ctfSamplingRate.get()
            scale = sampling_rate_part/sampling_rate_ctf

            # coordinates scaled
            coordinates_scaled = [item._coordinate._x.get() * scale, item._coordinate._y.get() * scale]

            # Key of dictionary is the micrograph ID
            # If the key already exists in the dictionary, increment the count
            if micrograph_num in mic_dict:
                mic_dict[micrograph_num]['particles_num'] += 1
                mic_dict[micrograph_num]['coordinates'].append(coordinates_scaled)

            else:
                # Otherwise, add the key with an initial count of 1
                mic_dict[micrograph_num] = {
                    'particles_num': 1,
                    'coordinates': [coordinates_scaled],
                    'mic_path': item._ctfModel._micObj._filename.get()
                }

        # Calculate the mean particle values
        particle_counts = [data['particles_num'] for data in mic_dict.values()]
        mean_particles_values = np.mean(particle_counts)

        # Plot histogram of particle number
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(particle_counts, bins='auto', edgecolor='black')
        plt.xlabel('# Particles per Micrograph')
        plt.ylabel(hist_ylabel_mic)
        plt.title('Histogram for particle number per micrograph')

        hist_name = 'particles_hist.jpg'
        particles_hist = self.hist_path(hist_name)
        plt.savefig(particles_hist)


        # Obtain 3 micrographs with particles drawn on them
        # Retrieve mic ID of micrograph with higher, medium and lower number of particles:
        mic_closest_to_mean = min(mic_dict.items(), key=lambda x: abs(x[1]['particles_num']- mean_particles_values))[0]
        mic_with_lowest_particles = min(mic_dict, key=lambda k: mic_dict[k]['particles_num'])
        mic_with_highest_particles = max(mic_dict, key=lambda k: mic_dict[k]['particles_num'])

        mic_ids = [int(mic_with_highest_particles), int(mic_closest_to_mean), int(mic_with_lowest_particles)]

        # Dict to store the path and coordinates of the 3 micrographs keeping mic_ids order
        reduced_mics_dict = {mic_id: {"mic_path": mic_dict.get(str(mic_id), {}).get('mic_path', 0), #None,
                                      "coordinates": mic_dict.get(str(mic_id), {}).get('coordinates', 0)
                             #        "particles_num": mic_dict.get(str(mic_id), {}).get('particles_num', 0)
                                       }
                             for mic_id in mic_ids
                                      }

        # Draw particles in images
        images = [] # List to store the images after drawing particles
        for micrograph, values in reduced_mics_dict.items():
            with mrcfile.open(values['mic_path'], permissive=True) as mrc:
                mrc_data = mrc.data
            mrc_normalized = 255 * (mrc_data - np.min(mrc_data)) / (np.max(mrc_data) - np.min(mrc_data))
            mrc_normalized = mrc_normalized.astype(np.uint8)
            image = Image.fromarray(mrc_normalized).convert('RGB')

            W_jpg, H_jpg = image.size
            draw = ImageDraw.Draw(image)
            r = W_jpg / 256

            for coord in values['coordinates']:
                x = coord[0]
                y = coord[1]
                draw.ellipse((x - r, y - r, x + r, y + r), fill=(0, 255, 0))

                # Draw number of particles
                # font = ImageFont.load_default()
                # # Draw text on images
                # draw = ImageDraw.Draw(image)
                # text = str(values['particles_num'])
                # draw.text((10, 10), text, fill='#80FF00', font=font)

            # Append the image to the list
            images.append(image)

        # Create collage
        width, height = images[0].size
        collage_width = 3 * width
        collage_height = height
        collage = Image.new('RGB', (collage_width, collage_height))
        # Paste each image into the collage
        for i, img in enumerate(images):
            collage.paste(img, (i * width, 0))

        extra_folder = self._getExtraPath()
        micro_folder_name = 'Micro_examples'
        micro_folder_path = join(extra_folder, micro_folder_name)

        micro_name = 'micro_particles.jpg'
        micro_path = join(micro_folder_path, micro_name)
        collage.save(micro_path)

        particles = {"Number_particles": sum(particle_counts),
                     "Particles_per_micrograph": round(mean_particles_values, 1),
                     "Particles_histogram": hist_name,
                     'Particles_mic_examples': join(micro_folder_name, micro_name)
                     }

        return particles

    def classes2D_generation(self):
        classes2D = self.classes2D.get()

        # Creation of a list (copy) of classes2D
        list_classes = []
        for cl in classes2D.iterItems():
            new_class = Class2D()
            new_class.copy(cl)

            new_class_repre = new_class.getRepresentative()
            current_class_repre = cl.getRepresentative()
            new_class_repre.setIndex(current_class_repre.getIndex())
            new_class._size.set(cl.getSize())

            list_classes.append(new_class)

        # Sorting list in descending order regarding the number of particles
        sorted_list_classes = sorted(list_classes, key=lambda x: x.getSize(), reverse=True)
        classes = len(sorted_list_classes)  # Number of classes

        img_classes_file = classes2D.getFirstItem().getRepresentative().getFileName()
        # Saving images in .jpg, drawing number of particles on them
        particles_list = []
        img_filenames = []
        with mrcfile.open(img_classes_file, 'r') as mrc:
            data = mrc.data
            for i, class_2D in enumerate(sorted_list_classes):
                particles = class_2D.getSize()
                particles_list.append(particles)
                index = class_2D.getRepresentative().getIndex()

                img = data[index - 1, :, :]

                img_normalized = 255 * (img - np.min(img)) / (np.max(img) - np.min(img))
                img_normalized = img_normalized.astype(np.uint8)
                image = Image.fromarray(img_normalized)
                image = image.convert('RGB')

                # Draw the number of particles on the images
                draw = ImageDraw.Draw(image)
                font = ImageFont.load_default()
                position = (10, 10)
                draw.text(position, str(particles), fill='#80FF00', font=font)

                # Saving images
                new_img_filename = splitext(img_classes_file)[0] + f'_image_{i}.jpg'
                image.save(new_img_filename)
                img_filenames.append(new_img_filename)

            # Creating collage in .jpg with all images ordered in descending order
            images = [Image.open(filename) for filename in img_filenames]

            output_folder = self._getExtraPath()  # Extra path of current protocol
            collage_filename = 'classes_2D.jpg'
            collage_filepath = join(output_folder, collage_filename)
            self.create_collage(images, collage_filepath)

        classes_2D = {"Number_classes_2D": classes, "Particles_per_class": particles_list,
                      "Images_classes_2D": collage_filename}

        return classes_2D

    def volume_generation(self, volume_type, folder_name, volume, th):
        print(f"att: {dir(volume)}")
        print(f" file name: {volume.getFileName()}")

        extra_folder = self._getExtraPath()
        # initial_vol_folder_name = 'Initial_volume'
        folder_path = join(extra_folder, folder_name)
        os.makedirs(folder_path, exist_ok=True)
        # volume = self.initVolume.get()
        volume_file = volume.getFileName()

        # Access to particles.sqlite to obtain particles in volume:
        reference_path = volume.getFileName()
        base_directory = os.path.dirname(os.path.dirname(reference_path))
        sqlite_file = os.path.join(base_directory, "particles.sqlite")

        if not os.path.exists(sqlite_file):
            print(f"SQLite file not found: {sqlite_file}")
        else:
            # Connect to the SQLite database
            conn = sqlite3.connect(sqlite_file)
            cursor = conn.cursor()

            query = "SELECT value FROM properties WHERE key = '_size'"
            cursor.execute(query)

            particle_num = cursor.fetchone()

            if particle_num:
                # Extract the number from the tuple
                size_value = particle_num[0]
                print(f"The value of '_size' is: {size_value}")
            else:
                print("Key '_size' not found in the properties table.")
            conn.close()


        # Getting orthogonal slices in X, Y and Z
        # Folder to store orthogonal slices
        orthogonal_slices_folder = 'orthogonal_slices'
        orthogonal_slices_path = join(folder_path, orthogonal_slices_folder)
        os.makedirs(orthogonal_slices_path, exist_ok=True)

        self.orthogonalSlices(fnRoot=orthogonal_slices_path, map=volume_file)

        # Getting 3 isosurface images
        # Folder to store isosurface images
        isosurface_images_folder = 'isosurface_images'
        isosurface_images_path = join(folder_path, isosurface_images_folder)
        os.makedirs(isosurface_images_path, exist_ok=True)

        # th = int(self.threshold_initVol.get())

        volume_file_abspath = abspath(volume_file)
        front_view_img = 'front_view.jpg'
        side_view_img = 'side_view.jpg'
        top_view_img = 'top_view.jpg'
        self.generate_isosurfaces(isosurface_images_path, volume_file_abspath,
                                  th, front_view_img, side_view_img, top_view_img)

        if 'size_value' in locals():
            volume = {
                'Volume_type': volume_type,
                'number_particles': int(size_value),
                'Orthogonal_slices': {
                    'Orthogonal_slices_X': join(folder_name, orthogonal_slices_folder, "orthogonal_slices_X.jpg"),
                    'Orthogonal_slices_Y': join(folder_name, orthogonal_slices_folder, "orthogonal_slices_Y.jpg"),
                    'Orthogonal_slices_Z': join(folder_name, orthogonal_slices_folder, "orthogonal_slices_Z.jpg")
                },
                'Isosurface_images': {
                    'Front_view': join(folder_name, isosurface_images_folder, front_view_img),
                    'Side_view': join(folder_name, isosurface_images_folder, side_view_img),
                    'Top_view': join(folder_name, isosurface_images_folder, top_view_img)
                }}
        else:
            volume = {
                'Volume_type': volume_type,
                'Orthogonal_slices': {
                'Orthogonal_slices_X': join(folder_name, orthogonal_slices_folder, "orthogonal_slices_X.jpg"),
                'Orthogonal_slices_Y': join(folder_name, orthogonal_slices_folder, "orthogonal_slices_Y.jpg"),
                'Orthogonal_slices_Z': join(folder_name, orthogonal_slices_folder, "orthogonal_slices_Z.jpg")
            },
                'Isosurface_images': {
                    'Front_view': join(folder_name, isosurface_images_folder, front_view_img),
                    'Side_view': join(folder_name, isosurface_images_folder, side_view_img),
                    'Top_view': join(folder_name, isosurface_images_folder, top_view_img)
                }}

        return volume

    def classes3D_generation(self):
        classes3D = self.classes3D.get()

        extra_folder = self._getExtraPath()  # Extra path of current protocol
        classes_3D_folder_name = 'Classes_3D'
        classes3D_folder_path = join(extra_folder, classes_3D_folder_name)
        os.makedirs(classes3D_folder_path, exist_ok=True)

        # Creation of a list (copy) of classes3D
        list_classes = []
        for cl in classes3D.iterItems():
            new_class = Class3D()
            new_class.copy(cl)

            new_class_repre = new_class.getRepresentative()
            current_class_repre = cl.getRepresentative()
            new_class_repre.setIndex(current_class_repre.getIndex())
            new_class._size.set(cl.getSize())

            list_classes.append(new_class)

        # Sorting list in descending order regarding the number of particles
        sorted_list_classes = sorted(list_classes, key=lambda x: x.getSize(), reverse=True)
        classes = len(sorted_list_classes)  # Number of classes

        # Saving images in .jpg, drawing number of particles on them
        particles_list = []
        img_filenames = []

        volumes_list = []  # List to store volume data

        for i, class_3D in enumerate(sorted_list_classes):
            file_name = sorted_list_classes[i].getRepresentative().getFileName()
            # Option 1:
            # img = self.readMap(file_name)
            # data = img.getData()

            # Option 2:
            if ':' in file_name:
                file_name_without_suffix = file_name.split(':')[0]

            with mrcfile.open(file_name_without_suffix, 'r') as mrc:
                data = mrc.data

                ############################
                ########## CLASSES ##########
                ############################

                particles = class_3D.getSize()
                particles_list.append(particles)

                mid_index = data.shape[0] // 2

                img = data[mid_index, :, :]

                img_normalized = 255 * (img - np.min(img)) / (np.max(img) - np.min(img))
                img_normalized = img_normalized.astype(np.uint8)
                image = Image.fromarray(img_normalized)
                image = image.convert('RGB')

                # Draw the number of particles on the images
                draw = ImageDraw.Draw(image)
                font = ImageFont.load_default()
                position = (10, 10)
                draw.text(position, str(particles), fill='#80FF00', font=font)

                # Saving images
                new_img_filename = splitext(file_name)[0] + '.jpg'
                image.save(new_img_filename)
                img_filenames.append(new_img_filename)

                #############################
                ########## VOLUMES ##########
                #############################

                # Getting orthogonal slices in X, Y and Z
                # Folder to store orthogonal slices
                orthogonal_slices_folder = f'orthogonal_slices_volume{i + 1}'
                orthogonal_slices_path = join(classes3D_folder_path, orthogonal_slices_folder)
                os.makedirs(orthogonal_slices_path, exist_ok=True)

                self.orthogonalSlices(fnRoot=orthogonal_slices_path, map=data)

                # Getting 3 isosurface images
                # Folder to store isosurface images
                isosurface_images_folder = f'isosurface_images_volume{i + 1}'
                isosurface_images_path = join(classes3D_folder_path, isosurface_images_folder)
                os.makedirs(isosurface_images_path, exist_ok=True)

                th = int(self.threshold_classes3D.get())

                volume_file_abspath = abspath(file_name_without_suffix)
                front_view_img = 'front_view.jpg'
                side_view_img = 'side_view.jpg'
                top_view_img = 'top_view.jpg'
                self.generate_isosurfaces(isosurface_images_path, volume_file_abspath,
                                          th, front_view_img, side_view_img, top_view_img)

                # Dictionary fill in:
                volume = {
                    "Orthogonal_slices": {
                        "Orthogonal_slices_X": join(classes_3D_folder_name, orthogonal_slices_folder,
                                                    "orthogonal_slices_X.jpg"),
                        "Orthogonal_slices_Y": join(classes_3D_folder_name, orthogonal_slices_folder,
                                                    "orthogonal_slices_Y.jpg"),
                        "Orthogonal_slices_Z": join(classes_3D_folder_name, orthogonal_slices_folder,
                                                    "orthogonal_slices_Z.jpg")},
                    'Isosurface_images': {
                        'Front_view': join(classes_3D_folder_name, isosurface_images_folder, front_view_img),
                        'Side_view': join(classes_3D_folder_name, isosurface_images_folder, side_view_img),
                        'Top_view': join(classes_3D_folder_name, isosurface_images_folder, top_view_img)
                    }}
                # Add this volume to the volumes list
                volumes_list.append(volume)

        # Creating collage in .jpg with all images ordered in descending order
        images = [Image.open(filename) for filename in img_filenames]
        collage_filename = 'classes_3D.jpg'
        collage_filepath = join(classes3D_folder_path, collage_filename)
        self.create_collage(images, collage_filepath)

        classes_3D = {"Number_classes_3D": classes, "Particles_per_class": particles_list,
                      "Images_classes_3D": join(classes_3D_folder_name, collage_filename),
                      "Volumes": volumes_list}
        return classes_3D

    def preprocess_data(self, data):
        """Recursively convert complex objects into YAML-compatible types."""
        if isinstance(data, dict):
            return {key: self.preprocess_data(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [self.preprocess_data(item) for item in data]
        elif isinstance(data, np.ndarray):  # Convert numpy arrays to lists
            return data.tolist()
        elif isinstance(data, (np.float64, np.float32, np.int64, np.int32)):  # Convert numpy numbers
            return data.item()
        elif isinstance(data, bytes):  # Convert binary data to strings if applicable
            return data.decode("utf-8")
        else:
            return data

    def saveJson(self):
        file_path = self.getOutFile()
        # with open(file_path, 'w', encoding='utf-8') as json_file:
        #     json.dump(self.processing_json, json_file, ensure_ascii=False, indent=4)
        # print(f"JSON data successfully saved to {file_path}")
        # Save the data in YAML format
        preprocessed_data = self.preprocess_data(self.processing_json)

        with open(file_path, 'w', encoding='utf-8') as yaml_file:
            yaml.dump(preprocessed_data, yaml_file, allow_unicode=True, sort_keys=False, indent=4)
        print(f"YAML data successfully saved to {file_path}")


        ### JSON TO SEE
        # Now, convert the YAML to JSON and save it as well
        file_path_json = file_path.replace('.yaml', '.json')  # Change extension from .yaml to .json

        # Convert the data back to JSON
        with open(file_path, 'r', encoding='utf-8') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)  # Load YAML into a Python object

        # Save the JSON data to a new file
        with open(file_path_json, 'w', encoding='utf-8') as json_file:
            json.dump(yaml_data, json_file, ensure_ascii=False, indent=4)
        print(f"JSON data successfully saved to {file_path_json}")




    def getOutFile(self):
        return self._getExtraPath(OUTFILE)

    def hist_path(self, file_name):
        folder_path = self._getExtraPath()
        file_path = join(folder_path, file_name)
        return file_path

    def orthogonalSlices(self, fnRoot, map):
        if type(map) is str:
            V = self.readMap(map)
            mV = V.getData()
            mV = np.squeeze(mV)
            Zdim, Ydim, Xdim = mV.shape
        else:
            mV = map
            Zdim, Ydim, Xdim = mV.shape

        min_val = mV.min()
        max_val = mV.max()
        mV = 255 * (mV - min_val) / (max_val - min_val)
        mV = mV.astype(np.uint8)

        slices_X = [mV[:, :, i] for i in range(Xdim)]
        slices_Y = [mV[:, i, :] for i in range(Ydim)]
        slices_Z = [mV[i, :, :] for i in range(Zdim)]

        # Save slices as images
        images_X = self.slices_to_images(slices_X)
        images_Y = self.slices_to_images(slices_Y)
        images_Z = self.slices_to_images(slices_Z)

        collage_X_path = join(fnRoot, 'orthogonal_slices_X.jpg')
        collage_Y_path = join(fnRoot, 'orthogonal_slices_Y.jpg')
        collage_Z_path = join(fnRoot, 'orthogonal_slices_Z.jpg')

        self.create_collage(images_X, collage_X_path)
        self.create_collage(images_Y, collage_Y_path)
        self.create_collage(images_Z, collage_Z_path)

    # Convert slices to images
    def slices_to_images(self, slices):
        images = [Image.fromarray(slice) for slice in slices]
        return images

    def create_collage(self, images, collage_filename):
        img_width, img_height = images[0].size

        # Define the number of rows and columns for the collage
        num_images = len(images)
        num_columns = int(np.ceil(np.sqrt(num_images)))
        num_rows = int(np.ceil(num_images / num_columns))

        # Create a blank canvas for the collage
        collage_width = num_columns * img_width
        collage_height = num_rows * img_height
        collage = Image.new('RGB', (collage_width, collage_height))

        # Paste each image onto the collage canvas in sorted order
        for index, image in enumerate(images):
            row = index // num_columns
            col = index % num_columns
            collage.paste(image, (col * img_width, row * img_height))

        collage.save(collage_filename)

    def readMap(self, fnMap):
        return xmipp3.Image(fnMap)

    def generateChimeraView(self, fnWorkingDir, fnMap, fnView, isMap=True, threshold=0, angX=0, angY=0, angZ=0,
                            bfactor=False, \
                            occupancy=False, otherAttribute=[], rainbow=True, legendMin=None, legendMax=None):
        maxMemToUse = 60000
        maxVoxelsToOpen = 1500

        chimeraScript = \
            """
            windowsize 1300 700
            set bgColor white
            """

        if isMap:
            chimeraScript += \
                """
                volume dataCacheSize %d
                volume voxelLimitForOpen %d
                volume showPlane false
                """ % (maxMemToUse, maxVoxelsToOpen)
            chimeraScript += \
                """
                open %s
                """ % fnMap

            if threshold == -1:
                chimeraScript += \
                    """show #1 models
                    volume #1 color #4e9a06
                    lighting soft
                    """
            else:
                chimeraScript += \
                    """show #1 models
                    volume #1 level %f
                    volume #1 color #4e9a06
                    lighting soft
                    """ % threshold
            #  # volume #1 level %f
            # % threshold
        else:
            chimeraScript += \
                """hide atoms
                show cartoons
                """
            if bfactor:
                chimeraScript += "color bfactor\n"
            if occupancy:
                chimeraScript += "color byattribute occupancy"
                if rainbow:
                    chimeraScript += " palette rainbow\n"
                if legendMin and legendMax:
                    chimeraScript += " palette bluered range %s,%s\n" % (legendMin, legendMax)
                else:
                    chimeraScript += "\n"
            if len(otherAttribute) > 0:
                chimeraScript += "open %s\n" % otherAttribute[0]
                chimeraScript += "color byattribute %s" % otherAttribute[1]
                if legendMin and legendMax:
                    chimeraScript += " palette bluered range %s,%s\n" % (legendMin, legendMax)
                else:
                    chimeraScript += "\n"
            if legendMin and legendMax:
                chimeraScript += "key blue:%s white: red:%s fontSize 15 size 0.025,0.4 pos 0.01,0.3\n" % (
                    legendMin, legendMax)
        chimeraScript += \
            """turn x %f
            turn y %f
            turn z %f
            view all
            save %s
            exit
            """ % (angX, angY, angZ, fnView)
        fnTmp = join(fnWorkingDir, "chimeraScript.cxc")
        fh = open(fnTmp, "w")
        fh.write(chimeraScript)
        fh.close()

        from chimera import Plugin
        args = "--nogui --offscreen chimeraScript.cxc"
        Plugin.runChimeraProgram(Plugin.getProgram(), args, cwd=fnWorkingDir)
        # cleanPath(fnTmp)

    def generate_isosurfaces(self, isosurface_img_path, volume_file_path,
                             th, front_view_img, side_view_img, top_view_img):
        working_path = dirname(volume_file_path)

        # Front_view
        output_path = abspath(join(isosurface_img_path, front_view_img))
        self.generateChimeraView(fnWorkingDir=working_path, fnMap=volume_file_path,
                                 fnView=output_path, threshold=th, angX=0, angY=0, angZ=0)

        # Side view (rotated 90 degrees around Y-axis)
        output_path = abspath(join(isosurface_img_path, side_view_img))
        self.generateChimeraView(fnWorkingDir=working_path, fnMap=volume_file_path,
                                 fnView=output_path, threshold=th, angX=0, angY=90, angZ=0)

        # Top view (rotated 90 degrees around X-axis)
        output_path = abspath(join(isosurface_img_path, top_view_img))
        self.generateChimeraView(fnWorkingDir=working_path, fnMap=volume_file_path,
                                 fnView=output_path, threshold=th, angX=90, angY=0, angZ=0)

    # Resize images to make text appear larger
    def resize_image(self, image, scale):
        width, height = image.size
        new_size = (int(width * scale), int(height * scale))
        return image.resize(new_size, Image.Resampling.LANCZOS)

    # Resize images back to original size
    def resize_back(self, image, original_size):
        width, height = original_size
        enlarge_factor= 1.5
        new_size = (int(width * enlarge_factor), int(height * enlarge_factor))
        return image.resize(new_size, Image.Resampling.LANCZOS)
