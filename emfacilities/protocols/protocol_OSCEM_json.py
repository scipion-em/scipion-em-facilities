import copy

from PIL import Image, ImageDraw, ImageFont
from matplotlib import pyplot as plt
import mrcfile

import pyworkflow.protocol.params as params
from pwem.objects import Class2D
from pwem.protocols import EMProtocol
from pwem.viewers import EmPlotter
from pyworkflow.object import String
import pyworkflow.utils as pwutils

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

        form.addParam('classes3D', params.PointerParam,
                      label="Classes 3D",
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
        # if doseperframe has a value, then dose initial is also retrieved
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

        # Histograms generation
        # DEFOCUS
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(defocus_list, bins='auto', edgecolor='black')
        plt.xlabel('# Defocus (A)')
        plt.ylabel('Frequency of Micrographs')
        plt.title('Defocus histogram')
        defocus_hist_name = 'defocus_hist.png'
        defocus_hist = self.hist_path(defocus_hist_name)
        plt.savefig(defocus_hist)

        # RESOLUTION
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(resolution_list, bins='auto', edgecolor='black')
        plt.xlabel("Resolution")
        plt.ylabel('Frequency of Micrographs')
        plt.title('Resolution histogram')
        resolution_hist_name = 'resolution_hist.png'
        resolution_hist = self.hist_path(resolution_hist_name)
        plt.savefig(resolution_hist)

        # ASTIGMATISM
        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(astigmatism_list, bins='auto', edgecolor='black')
        plt.xlabel("Astigmatism")
        plt.ylabel('Frequency of Micrographs')
        plt.title('Astigmatism histogram')
        astigmatism_hist_name = 'astigmatism_hist.png'
        astigmatism_hist = self.hist_path(astigmatism_hist_name)
        plt.savefig(astigmatism_hist)

        defocus = {'Output_max_defocus': max_defocus, 'Output_min_defocus': min_defocus,
                   'Output_avg_defocus': avg_defocus, 'Defocus_histogram': defocus_hist_name}
        resolution = {'Output_max_resolution': max_resolution, 'Output_min_resolution': min_resolution,
                      'Output_avg_resolution': avg_resolution, 'Resolution_histogram': resolution_hist_name}
        CTF_estimation['Defocus'] = defocus
        CTF_estimation['Resolution'] = resolution
        CTF_estimation['Astigmatism'] = {'Astigmatism_histogram': astigmatism_hist_name}

        return CTF_estimation

    def particles_generation(self):
        parts = self.particles.get()
        mic_numbers = []
        particle_counts = []
        for index, item in enumerate(parts.iterItems()):
            micrograph_num = item._micId
            key_name = f"mic_{micrograph_num}"
            if index == 0 or key_name not in mic_numbers:
                mic_numbers.append(key_name)
                particle_counts.append(1)
            else:
                index = mic_numbers.index(key_name)
                particle_counts[index] += 1

        mean_particles_values = np.mean(particle_counts)

        plt.close('all')
        plt.clf()
        plt.cla()

        plt.hist(particle_counts, bins='auto', edgecolor='black')
        plt.xlabel('# Particles per Micrograph')
        plt.ylabel('Frequency of Micrographs')
        plt.title('Histogram for particle number per micrograph')

        hist_name = 'particles_hist.png'
        particles_hist = self.hist_path(hist_name)
        plt.savefig(particles_hist)

        particles = {"Particles_per_micrograph": mean_particles_values,
                     "Particles_histogram": hist_name}
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
        print(f"Particles: {sorted_list_classes[10].getSize()}")
        print(f"Index: {sorted_list_classes[10].getRepresentative().getIndex()}")
        classes = len(sorted_list_classes)

        img_classes_file = classes2D.getFirstItem().getRepresentative().getFileName()
        print(img_classes_file)

        # Saving images in .png, drawing number of particles on them
        particles_list = []
        img_filenames = []
        with mrcfile.open(img_classes_file, 'r') as mrc:
            data = mrc.data
            # for i, data in enumerate(mrc.data):
            for i, class_2D in enumerate(sorted_list_classes):
                particles = class_2D.getSize()
                particles_list.append(particles)
                index = class_2D.getRepresentative().getIndex()

                img = data[index-1, :, :]

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
                new_img_filename = os.path.splitext(img_classes_file)[0] + f'_image_{i}.png'
                image.save(new_img_filename)
                img_filenames.append(new_img_filename)

            # Creating collage in .png with all images ordered in descending order
            if img_filenames:
                images = images = [Image.open(filename) for filename in img_filenames]
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

                output_folder = self._getExtraPath()  # Extra path of current protocol
                collage_filename = 'classes_2D.png'
                collage_filepath = os.path.join(output_folder, 'classes_2D.png')
                collage.save(collage_filepath)

        classes_2D = {"Number_classes_2D": classes, "Particles_per_class": particles_list,
                      "Images_classes_2D": collage_filename}

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

    def find_file_with_suffix(self, directory, suffix):
        """ Function to get the filename that ends with certain suffix
        """
        for root, _, files in os.walk(directory):
            for file in files:
                if file.endswith(suffix):
                    return os.path.join(root, file)
        return None  # Return None if no matching file is found
