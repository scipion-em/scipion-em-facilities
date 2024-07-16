import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import String

import numpy as np
import json
import os

INPUT_MOVIES = 0
INPUT_MICS = 1
OUTFILE = 'Processing_json.json'


class ProtOSCEM(EMProtocol):
    """ This is the class for generating the OSCEM metadata json file from Scipion workflow
    """
    _label = 'OSCEM Json'

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
                      label="Movie alignment", important=True,
                      pointerClass='ProtAlignMovies',
                      help="Movie alignment used",
                      condition=CONDITION_MOVIES,
                      allowsNull=True)

        form.addParam('maxShift', params.PointerParam,
                      label="Max Shift", important=True,
                      pointerClass='XmippProtMovieMaxShift',
                      help="Max Shift",
                      condition=CONDITION_MOVIES,
                      allowsNull=True)

        form.addParam('CTF', params.PointerParam,
                      label="CTF", important=True,
                      pointerClass='SetOfCTF',
                      help="CTF micrographs",
                      allowsNull=True)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.generateJson)
        self._insertFunctionStep(self.saveJson)

    # -------------------------- STEPS functions ------------------------------
    def generateJson(self):

        self.processing_json = {}
        # self.processing_json['Description'] = 'OSCEM json for processing. Version 1.0'  ## !!!!

        if self.inputType.get() == 0:  # movies as input
            ###### IMPORT MOVIES ######
            import_movies = self.import_movies_generation()
            # print('Import movie dict:')
            # print(json.dumps(import_movies, indent=4))
            self.processing_json['Import_movies'] = import_movies

            if self.movieAlignment.get() is not None:
                ###### MOVIE ALGINMENT ######
                movie_alignment = self.movie_alignment_generation()
                # print('Movie alignment dict:')
                # print(json.dumps(movie_alignment, indent=4))
                self.processing_json['Movie_alignment'] = movie_alignment

            if self.maxShift.get() is not None:
                ###### MAX SHIFT ######
                max_shift = self.max_shift_generation()
                # print('Max shift dict:')
                # print(json.dumps(max_shift, indent=4))
                self.processing_json['Movie_maxshift'] = max_shift

        if self.CTF.get() is not None:
            ###### CTF ######
            CTF = self.CTF_generation()
            # print('CTF estimation dict:')
            # print(json.dumps(CTF, indent=4))
            self.processing_json['CTF_estimation'] = CTF

        print(json.dumps(self.processing_json, indent=4))

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
        # if doseperframe has a value, then dose initial is also retrieveD.
        # Otherwise, none of them are retrieved.
        if input_movies['outputMovies._acquisition._dosePerFrame'] is None:
            keys_to_retrieve = ['outputMovies._acquisition._voltage', 'outputMovies._acquisition._sphericalAberration',
                                'outputMovies._acquisition._amplitudeContrast', 'outputMovies._samplingRate',
                                'outputMovies._scannedPixelSize',
                                'outputMovies._gainFile', 'outputMovies._darkFile', 'outputMovies._size']
        else:
            keys_to_retrieve = ['outputMovies._acquisition._voltage', 'outputMovies._acquisition._sphericalAberration',
                                'outputMovies._acquisition._amplitudeContrast', 'outputMovies._samplingRate',
                                'outputMovies._scannedPixelSize',
                                'outputMovies._acquisition._dosePerFrame', 'outputMovies._acquisition._doseInitial',
                                'outputMovies._gainFile',
                                'outputMovies._darkFile', 'outputMovies._size']

        # Mapping dictionary for key name changes
        key_mapping = {
            'outputMovies._acquisition._voltage': 'Microscope_voltage',
            'outputMovies._acquisition._sphericalAberration': 'Spherical_aberration',
            'outputMovies._acquisition._amplitudeContrast': 'Amplitud_contrast',
            'outputMovies._samplingRate': 'Sampling_rate',
            'outputMovies._scannedPixelSize': 'Pixel_size',
            'outputMovies._acquisition._dosePerFrame': 'Dose_per_image',
            'outputMovies._acquisition._doseInitial': 'Initial_dose',
            'outputMovies._gainFile': 'Gain_image',
            'outputMovies._darkFile': 'Dark_image',
            'outputMovies._size': 'Number_movies'
        }

        # Filter the dictionary and rename the keys
        boolean_keys = ['outputMovies._gainFile', 'outputMovies._darkFile']
        import_movies = {}
        for key in keys_to_retrieve:
            if key in input_movies and input_movies[key] is not None and input_movies[key] != 0:
                if key in boolean_keys:
                    import_movies[key_mapping[key]] = bool(input_movies[key])
                else:
                    import_movies[key_mapping[key]] = input_movies[key]

        return import_movies

    def movie_alignment_generation(self):
        MovieAlignmentProt = self.movieAlignment.get()
        movie_align = {'Method': MovieAlignmentProt.getClassName()}

        ################################ INPUT #############################################
        input_alignment = MovieAlignmentProt.getObjDict()
        # List of keys to retrieve
        keys_to_retrieve = ['binFactor', 'maxResForCorrelation', 'gainRot', 'gainFlip']

        ## TOD: apply dose filter, where?
        # Mapping dictionary for key name changes
        key_mapping = {
            'binFactor': 'Binning_factor',
            'maxResForCorrelation': 'Maximum_resolution',
            'gainRot': 'Rotate_gain_reference',
            'gainFlip': 'Flip_gain_reference'
        }

        # Filter the dictionary and rename the keys
        input_movie_align = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
                             key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0}
        movie_align.update(input_movie_align)

        # Dictionary for crop offsets
        keys_to_retrieve = ['cropOffsetX', 'cropOffsetY']
        key_mapping = {
            'cropOffsetX': 'Crop_offsetX',
            'cropOffsetY': 'Crop_offsetX',
        }
        crop_offsets = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
                        key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0}
        if crop_offsets:
            movie_align['Crop_offsets'] = crop_offsets

        # Dictionary for crop dims
        keys_to_retrieve = ['cropDimX', 'cropDimY']
        key_mapping = {
            'cropDimX': 'Crop_dimsX',
            'cropDimY': 'Crop_dimsY',
        }
        crop_dims = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
                     key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0}
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
        if frames_aligned:
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

            output_movie_align = {'Output_avg_shift': avg_shift, 'Output_max_shift': max_shift}
            movie_align.update(output_movie_align)

        return movie_align

    def max_shift_generation(self):
        MaxShiftProt = self.maxShift.get()

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
        movie_maxshift = {}
        for key in keys_to_retrieve:
            if key == 'rejType':
                if input_shift[key] == 0:
                    rej_type = 'By frame'
                if input_shift[key] == 1:
                    rej_type = 'By whole movie'
                if input_shift[key] == 2:
                    rej_type = 'By frame and movie'
                if input_shift[key] == 3:
                    rej_type = 'By frame or movie'
                movie_maxshift[key_mapping[key]] = rej_type
            elif key in input_shift and input_shift[key] is not None and input_shift[key] != 0:
                movie_maxshift[key_mapping[key]] = input_shift[key]

        ############################### OUTPUT #############################################
        # average and max shift
        max_shift = 0
        for a, output in MaxShiftProt.iterOutputAttributes():
            if a == 'outputMovies':
                for index, item in enumerate(output.iterItems()):
                    attributes = item.getAttributes()
                    attributes_dict = dict(attributes)
                    shiftX = attributes_dict.get('_xmipp_ShiftX')
                    shiftY = attributes_dict.get('_xmipp_ShiftY')
                    norm = np.linalg.norm([shiftX, shiftY], axis=0)

                    # Max shift
                    max_norm = np.max(norm)
                    max_shift = max(max_shift, max_norm)

                    # Average shift
                    avgXY = np.mean(norm)
                    if index == 0:
                        avg_shift = avgXY
                    else:
                        avg_shift = np.mean([avgXY, avg_shift])

        output_movie_maxshift = {'Output_avg_shift': avg_shift, 'Output_max_shift': max_shift}
        movie_maxshift.update(output_movie_maxshift)

        return movie_maxshift

    def CTF_generation(self):
        CTFs = self.CTF.get()
        ############################## OUTPUT #############################################
        CTF_estimation = {}
        for index, item in enumerate(CTFs.iterItems()):
            # Min, max  and average defocus and resolution
            defocus = np.mean([float(item._defocusU), float(item._defocusV)])
            resolution = float(item._resolution)

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

        defocus = {'Output_max_defocus': max_defocus, 'Output_min_defocus': min_defocus,
                   'Output_avg_defocus': avg_defocus}
        resolution = {'Output_max_resolution': max_resolution, 'Output_min_resolution': min_resolution,
                      'Output_avg_resolution': avg_resolution}
        CTF_estimation['Defocus'] = defocus
        CTF_estimation['Resolution'] = resolution

        return CTF_estimation

    def saveJson(self):

        file_path = self.getOutFile()
        try:
            with open(file_path, 'w') as json_file:
                json.dump(self.processing_json, json_file, indent=4)
            print(f"JSON data successfully saved to {file_path}")
        except Exception as e:
            print(f"An error occurred while saving JSON data: {e}")

    def getOutFile(self):
        return self._getExtraPath(OUTFILE)