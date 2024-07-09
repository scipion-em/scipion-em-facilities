import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol
from pyworkflow.object import String

import numpy as np
import json
import os


class ProtOSCEM(EMProtocol):
    """ This is the class for generating the OSCEM metadata json file from Scipion workflow
    """
    _label = 'OSCEM Json'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputType',  params.EnumParam, label='Input type', choices=['Movies','Micrographs'], important=True,
                      help='Select the type of input, either movies or micrographs')

        form.addParam('importMovies', params.PointerParam,
                      label="Import Movies", important=True,
                      pointerClass='ProtImportMovies',  # 'CTFModel'],
                      help="Import Movies",
                      condition='inputType==0')  # extender este help

        form.addParam('movieAlignment', params.PointerParam,
                      label="Movie alignment", important=True,
                      pointerClass='ProtAlignMovies',  # 'CTFModel'],
                      help="Movie alignment used",
                      condition='inputType==0',
                      allowsNull=True)  # extender este help

        form.addParam('maxShift', params.PointerParam,
                      label="Max Shift", important=True,
                      pointerClass='XmippProtMovieMaxShift',  # 'CTFModel'],
                      help="Max Shift",
                      condition='inputType==0',
                      allowsNull=True)  # extender este help

        form.addParam('CTF', params.PointerParam,
                      label="CTF", important=True,
                      pointerClass='SetOfCTF',  # 'CTFModel'],
                      help="Max Shift",
                      allowsNull=True)  # extender este help

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.generateJson)
        self._insertFunctionStep(self.saveJson)

    # -------------------------- STEPS functions ------------------------------
    def generateJson(self):

        self.processing_json = {}
        ###### IMPORT MOVIES ######
        import_movies = self.import_movies_generation()
        # print('Import movie dict:')
        # print(json.dumps(import_movies, indent=4))
        self.processing_json['Import_movies'] = import_movies

        ###### MOVIE ALGINMENT ######
        movie_alignment = self.movie_alignment_generation()
        # print('Movie alignment dict:')
        # print(json.dumps(movie_alignment, indent=4))
        self.processing_json['Movie_alignment'] = movie_alignment

        ###### MAX SHIFT ######
        max_shift = self.max_shift_generation()
        # print('Max shift dict:')
        # print(json.dumps(max_shift, indent=4))
        self.processing_json['Movie_maxshift'] = max_shift

        ###### CTF ######
        CTF = self.CTF_generation()
        # print('CTF estimation dict:')
        # print(json.dumps(CTF, indent=4))
        self.processing_json['CTF_estimation'] = CTF

        print(json.dumps(self.processing_json, indent=4))

        # protocol.alignFrame0.get() # assumes alignFrame0 exists as an attribute
        # importmoviesprotocol.get('alignFrame0', String()).get()

    # -------------------------- INFO functions -------------------------------
    def _validate(self):  # si no se cumple no ejecuta protocolo
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []

    # def getInputProtocols(self):  # multipointer por el for
    # protocols = []
    # for protPointer in self.inputProtocols:
    #     prot = protPointer.get()
    #     prot.setProject(self.getProject())
    #      protocols.append(prot)
    # return protocols

    # -------------------- UTILS functions -------------------------
    def import_movies_generation(self):
        ImportMoviesProt = self.importMovies.get()
        # print(ImportMoviesProt)
        input_movies = ImportMoviesProt.getObjDict()

        # print("Input params:")
        # print(inputParams)

        # List of keys to retrieve
        # if doseperframe has a value, then dose initial is also retrieveD. Otherwise, none of them are retrieved.
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
                    # Check if the value is not None or null (empty string in this case)
                    import_movies[key_mapping[key]] = bool(input_movies[key])
                else:
                    import_movies[key_mapping[key]] = input_movies[key]

        return import_movies

    def movie_alignment_generation(self):
        MovieAlignmentProt = self.movieAlignment.get()
        # print(MovieAlignmentProt)
        movie_align = {'Method': MovieAlignmentProt.getClassName()}

        #     ############################ INPUT #############################################

        input_alignment = MovieAlignmentProt.getObjDict()
        # print("Input params:")
        # print(input_alignment)

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
        keys_to_retrieve = ['alignFrame0', 'alignFrameN']
        key_mapping = {
            'alignFrame0': 'Frame0',
            'alignFrameN': 'FrameN',
        }
        frames_aligned = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
                          key in input_alignment and input_alignment[key] is not None}
        if frames_aligned:
            movie_align['Frames_aligned'] = frames_aligned

        #     ############################ OUTPUT #############################################
        # average and max shift
        max_shift = 0
        for a, output in MovieAlignmentProt.iterOutputAttributes():
            # print(f"a: {a}")
            # print(f"output: {output}")
            if a == 'outputMovies':  ## ESTOY COGIENDO MOVIES NO MICROGRAPHS. DUDA!!!
                for index, item in enumerate(output.iterItems()):
                    # print(f"item: {item}")
                    attributes = item.getAttributes()
                    # print(attributes)
                    attributes_dict = dict(attributes)
                    # print(f'attributes dir: {attributes_dict}')
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
                        avg_shift = np.mean([avgXY, avg_shift])  # average with total X and Y

                    # for key, value in attributes:
                    #     print(f"key: {key}")
                    #     print(f"value: {value}")

            output_movie_align = {'Ouput_avg_shift': avg_shift, 'Output_max_shift': max_shift}
            movie_align.update(output_movie_align)

        return movie_align

    def max_shift_generation(self):
        MaxShiftProt = self.maxShift.get()
        # print(MaxShiftProt)

        #     ############################ INPUT #############################################
        input_shift = MaxShiftProt.getObjDict()
        # print("Input params:")
        # print(inputParams)

        # List of keys to retrieve
        keys_to_retrieve = ['outputMoviesDiscarded._size', 'maxFrameShift', 'maxMovieShift', 'rejType']

        # Mapping dictionary for key name changes
        key_mapping = {
            'outputMoviesDiscarded._size': 'Discarded_movies',
            'maxFrameShift': 'Max_frame_shift',
            'maxMovieShift': 'Max_movie_shift',
            'rejType': 'Rejection_type'
        }

        #     # Filter the dictionary and rename the keys
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

        #     ############################ OUTPUT #############################################
        # average and max shift
        max_shift = 0
        for a, output in MaxShiftProt.iterOutputAttributes():
            # print(f"a: {a}")
            # print(f"output: {output}")
            if a == 'outputMovies':  ## ESTOY COGIENDO MOVIES NO MICROGRAPHS. DUDA!!!
                for index, item in enumerate(output.iterItems()):
                    # print(f"item: {item}")
                    attributes = item.getAttributes()
                    attributes_dict = dict(attributes)
                    # print(f'attributes dir: {attributes_dict}')
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
                        avg_shift = np.mean([avgXY, avg_shift])  # average with total X and Y

                    # for key, value in attributes:
                    #     print(f"key: {key}")
                    #     print(f"value: {value}")

        output_movie_maxshift = {'Ouput_avg_shift': avg_shift, 'Output_max_shift': max_shift}
        movie_maxshift.update(output_movie_maxshift)

        return movie_maxshift

    def CTF_generation(self):
        CTFs = self.CTF.get()
        # print('CTFs:')
        # print(CTFs)

        #     ############################ OUTPUT #############################################

        CTF_estimation = {}
        for index, item in enumerate(CTFs.iterItems()):
            # print(f"item: {item}")
            # print(f"type of item: {type(item)}")
            # print(dir(item))  # available attributes
            # print(f"Index: {index}, Item: {item}")
            # print(f"defU: {item._defocusU}")
            # print(f"defV: {item._defocusV}")
            # print(f"Resolution: {item._resolution}")

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
                   'Ouput_avg_defocus': avg_defocus}
        resolution = {'Output_max_resolution': max_resolution, 'Output_min_resolution': min_resolution,
                      'Ouput_avg_resolution': avg_resolution}
        CTF_estimation['Defocus'] = defocus
        CTF_estimation['Resolution'] = resolution

        return CTF_estimation

    def saveJson(self):

        folder_path = self._getExtraPath()
        # print(folder_path)
        file_path = os.path.join(folder_path, 'Processing_json')
        try:
            with open(file_path, 'w') as json_file:
                json.dump(self.processing_json, json_file, indent=4)
            print(f"JSON data successfully saved to {file_path}")
        except Exception as e:
            print(f"An error occurred while saving JSON data: {e}")
