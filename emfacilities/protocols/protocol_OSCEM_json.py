from matplotlib import pyplot as plt

import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pwem.viewers import EmPlotter
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

        form.addParam('particles', params.PointerParam,
                      label="Particles", important=True,
                      pointerClass='SetOfParticles',
                      help="Particles obtained when doing particle extraction",
                      allowsNull=True)

        form.addParam('classes2D', params.PointerParam,
                      label="Classes 2D", important=True,
                      pointerClass='SetOfClasses2D',
                      help="Set of 2D classes",
                      allowsNull=True)

        form.addParam('initVolume', params.PointerParam,
                      label="Initial volume", important=True,
                      pointerClass='Volume',
                      help="Initial volume",
                      allowsNull=True)

        form.addParam('classes3D', params.PointerParam,
                      label="Classes 3D", important=True,
                      pointerClass='SetOfClasses3D',
                      help="Set of 2D classes",
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

        if self.particles.get() is not None:
            ###### PARTICLES ######
            particles = self.particles_generation()
            self.processing_json['Particle_picking'] = particles

        if self.classes2D.get() is not None:
            ###### CLASSES 2D ######
            classes_2D = self.classes2D_generation()
            self.processing_json['Classes_2D'] = classes_2D

        if self.initVolume.get() is not None:
            ###### INITIAL VOLUME ######
            self.initVolume_generation()
            # self.processing_json['Initial_volume'] = volume

        if self.classes3D.get() is not None:
            ###### CLASSES 3D ######
            classes_3D = self.classes3D_generation()
            self.processing_json['Classes_3D'] = classes_3D


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
        defocus_list = []
        resolution_list = []
        astigmatism_list = []

        for index, item in enumerate(CTFs.iterItems()):
            # Min, max  and average defocus and resolution
            defocus = np.mean([float(item._defocusU), float(item._defocusV)])
            resolution = float(item._resolution)
            astigmatism = float(item._defocusRatio)

            defocus_list.append(defocus)
            resolution_list.append(resolution)
            astigmatism_list.append(astigmatism)

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

        # Histograms generation
        numberOfBins = 10
        plotterDefocus = EmPlotter()
        plotterResolution = EmPlotter()
        plotterAstigmatism = EmPlotter()

        # DEFOCUS
        plotterDefocus.createSubPlot("Defocus histogram", "Defocus (A)", "#")
        plotterDefocus.plotHist(defocus_list, nbins=numberOfBins)
        defocus_hist = self.hist_path('defocus_hist')
        plotterDefocus.savefig(defocus_hist)
        plotterDefocus.close()

        # RESOLUTION
        plotterResolution.createSubPlot("Resolution histogram", "Resolution", "#")
        plotterResolution.plotHist(resolution_list, nbins=numberOfBins)
        resolution_hist = self.hist_path('resolution_hist')
        plotterResolution.savefig(resolution_hist)
        plotterResolution.close()

        # ASTIGMATISM
        plotterAstigmatism.createSubPlot("Astigmatism histogram", "Astigmatism", "#")
        plotterAstigmatism.plotHist(astigmatism_list, nbins=numberOfBins)
        astigmatism_hist = self.hist_path('astigmatism_hist')
        plotterAstigmatism.savefig(astigmatism_hist)
        plotterAstigmatism.close()

        return CTF_estimation

    def particles_generation(self):
        parts = self.particles.get()
        # print(parts)
        # mic_numbers = []
        particle_counts = []
        for index, item in enumerate(parts.iterItems()):
            micrograph_num = item._micId
            key_name = f"mic_{micrograph_num}"
            # retrieve the number of particles per micrograph
            if index == 0: # or key_name not in mic_numbers:
                # mic_numbers.append(key_name)
                particle_counts.append(1)
            else:
                # index = mic_numbers.index(key_name)
                particle_counts[index] += 1

        mean_particles_values = np.mean(particle_counts)
        particles = {"Particles_per_micrograph": mean_particles_values}
        # print(particles)

        plt.hist(particle_counts)
        plt.xlabel('# Micrograph')
        plt.ylabel('# Particles')
        plt.title('Histogram for particle number per micrograph')

        particles_hist = self.hist_path('particles_hist')
        plt.savefig(particles_hist)

        return particles

    def classes2D_generation(self):
        classes2D = self.classes2D.get()
        # print(classes2D)
        attrib = classes2D.getAttributes()
        # print(attrib)
        particles_per_class = []
        for index, item in enumerate(classes2D.iterItems()):
            # print(f"item: {item}")
            # print(f" size: {item._size}")
            # print(f"type of item: {type(item)}")
            # print(dir(item)) # available attributes
            # print(f"Index: {index}, Item: {item}")
            number_particles = item._size.get()
            # print(f"particles: {number_particles}")
            particles_per_class.append(number_particles)

        classes = index + 1
        # print(f"classes: {clasases}")
        # print
        classes_2D = {"Number_classes_2D": classes, "Particles_per_class": particles_per_class}

        return classes_2D

    def initVolume_generation(self):
        # volume = self.initVolume.get()
        # print(volume)
        # # attrib = volume.getAttributes()
        # attributes = [attr for attr in dir(volume) if not attr.startswith('__')]
        # print(attributes)
        particles_per_class = []

    def classes3D_generation(self):
        classes3D = self.classes3D.get()
        # print(classes2D)
        attrib = classes3D.getAttributes()
        # print(attrib)
        particles_per_class = []
        for index, item in enumerate(classes3D.iterItems()):
            # print(f"item: {item}")
            # print(f" size: {item._size}")
            # print(f"type of item: {type(item)}")
            # print(dir(item)) # available attributes
            # print(f"Index: {index}, Item: {item}")
            number_particles = item._size.get()
            # print(f"particles: {number_particles}")
            particles_per_class.append(number_particles)
        #
        classes = index + 1
        # print(f"classes: {classes}")
        # print
        classes_3D = {"Number_classes_3D": classes, "Particles_per_class": particles_per_class}

        return classes_3D

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

    def hist_path(self, file_name):
        folder_path = self._getExtraPath()
        file_path = os.path.join(folder_path, file_name)
        return file_path
