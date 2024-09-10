from os.path import abspath, join, dirname, splitext

from PIL import Image, ImageDraw, ImageFont, ImageOps, ImageEnhance
from matplotlib import pyplot as plt
import mrcfile

import pyworkflow.protocol.params as params
import xmipp3
from pwem.objects import Class2D, Class3D
from pwem.protocols import EMProtocol

import numpy as np
import json
import os

INPUT_MOVIES = 0
INPUT_MICS = 1
OUTFILE = 'Processing_metadata.json'

# input movies attributes:
_voltage = 'outputMovies._acquisition._voltage'
_sphericalAberration = 'outputMovies._acquisition._sphericalAberration'
_amplitudeContrast = 'outputMovies._acquisition._amplitudeContrast'
_samplingRate = 'outputMovies._samplingRate'
_dosePerFrame = 'outputMovies._acquisition._dosePerFrame'
_doseInitial = 'outputMovies._acquisition._doseInitial'
_gainFile = 'outputMovies._gainFile'
_darkFile = 'outputMovies._darkFile'
_size = 'outputMovies._size'
_firstDim = 'outputMovies._firstDim'

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
            self.processing_json['Import_movies'] = import_movies

            if self.movieAlignment.get() is not None:
                ###### MOVIE ALGINMENT ######
                movie_alignment = self.movie_alignment_generation()
                self.processing_json['Movie_alignment'] = movie_alignment

            if self.maxShift.get() is not None:
                ###### MAX SHIFT ######
                max_shift = self.max_shift_generation()
                self.processing_json['Movie_maxshift'] = max_shift

        if self.CTF.get() is not None:
            ###### CTF ######
            CTF = self.CTF_generation()
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
            volume = self.init_volume_generation()
            self.processing_json['Initial_volume'] = volume

        if self.classes3D.get() is not None:
            ###### CLASSES 3D ######
            classes_3D = self.classes3D_generation()
            self.processing_json['Classes_3D'] = classes_3D

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
        # if doseperframe has a value, then dose initial is also retrieved
        # Otherwise, none of them are retrieved.
        if input_movies[_dosePerFrame] is None:
            keys_to_retrieve = [_voltage, _sphericalAberration, _amplitudeContrast, _samplingRate,
                                _gainFile, _darkFile, _size]
        else:
            keys_to_retrieve = [_voltage, _sphericalAberration, _amplitudeContrast, _samplingRate,
                                _dosePerFrame, _doseInitial, _gainFile, _darkFile, _size]

        # Mapping dictionary for key name changes
        key_mapping = {
            _voltage: 'Microscope_voltage_(kV)',
            _sphericalAberration: 'Spherical_aberration_(mm)',
            _amplitudeContrast: 'Amplitud_contrast',
            _samplingRate: 'Pixel_size_(Å/px)',
            _dosePerFrame: 'Dose_per_image_(e/Å²)',
            _doseInitial: 'Initial_dose_(e/Å²)',
            _gainFile: 'Gain_image',
            _darkFile: 'Dark_image',
            _size: 'Number_movies'
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

                png_path = os.path.join(extra_folder, f'{base_filename}.png')
                image.save(png_path)

                image_path = os.path.basename(png_path)
                input_movies[key] = os.path.basename(image_path)

        import_movies = {key_mapping[key]: input_movies[key] for key in keys_to_retrieve if
                         key in input_movies and input_movies[key] is not None and input_movies[key] != 0}

        # Retrieve nº of frames per movie and size of frames:
        dims = input_movies[_firstDim]
        dims_list = dims.split(',')
        dim1 = int(dims_list[0])
        dim2 = int(dims_list[1])
        n_frames = int(dims_list[2])
        frame_dim = f'{dim1} x {dim2}'

        import_movies['Frames_per_movie'] = n_frames
        import_movies['Frames_size_(pixels)'] = frame_dim

        return import_movies

    def movie_alignment_generation(self):
        MovieAlignmentProt = self.movieAlignment.get()
        movie_align = {'Method': MovieAlignmentProt.getClassName()}

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
            movie_align['Crop_offsets_(pixels)'] = crop_offsets

        # Dictionary for crop dims
        keys_to_retrieve = ['cropDimX', 'cropDimY']
        key_mapping = {
            'cropDimX': 'Crop_dimsX',
            'cropDimY': 'Crop_dimsY',
        }
        crop_dims = {key_mapping[key]: input_alignment[key] for key in keys_to_retrieve if
                     key in input_alignment and input_alignment[key] is not None and input_alignment[key] != 0}
        if crop_dims:
            movie_align['Crop_dims_(pixels)'] = crop_dims

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
            frames_aligned['FrameN'] = self.processing_json['Import_movies']["Number_movies"]

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

            output_movie_align = {'Output_avg_shift_(Å)': round(avg_shift, 1),
                                  'Output_max_shift_(Å)': round(max_shift, 1)}
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
            'maxFrameShift': 'Max_frame_shift_(Å)',
            'maxMovieShift': 'Max_movie_shift_(Å)',
            'rejType': 'Rejection_type'
        }

        # Filter dictionary and rename keys
        movie_maxshift = {}
        rejtype_list = ['By frame', 'By whole movie', 'By frame and movie', 'By frame or movie']
        for key in keys_to_retrieve:
            if key == 'rejType':
                rej_type = rejtype_list[input_shift[key]]
                movie_maxshift[key_mapping[key]] = rej_type
            elif key in input_shift and input_shift[key] is not None and input_shift[key] != 0:
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

        print(flattened_shift_list)
        plt.hist(flattened_shift_list, edgecolor='black')
        plt.xlabel('# Shift (Å)')
        plt.ylabel(hist_ylabel_frames)
        plt.title('Shift histogram')
        shift_hist_name = 'shift_hist.png'
        shift_hist = self.hist_path(shift_hist_name)
        plt.savefig(shift_hist)

        output_movie_maxshift = {'Output_avg_shift_(Å)': round(avg_shift, 1),
                                 'Output_max_shift_(Å)': round(max_shift, 1),
                                 'Shift_histogram': shift_hist_name}
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
        plt.xlabel('# Defocus')
        plt.ylabel(hist_ylabel_mic)
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
        plt.ylabel(hist_ylabel_mic)
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
        plt.ylabel(hist_ylabel_mic)
        plt.title('Astigmatism histogram')
        astigmatism_hist_name = 'astigmatism_hist.png'
        astigmatism_hist = self.hist_path(astigmatism_hist_name)
        plt.savefig(astigmatism_hist)

        defocus = {'Output_max_defocus': round(max_defocus, 1), 'Output_min_defocus': round(min_defocus, 1),
                   'Output_avg_defocus': round(avg_defocus, 1), 'Defocus_histogram': defocus_hist_name}
        resolution = {'Output_max_resolution': round(max_resolution, 1),
                      'Output_min_resolution': round(min_resolution, 1),
                      'Output_avg_resolution': round(avg_resolution, 1), 'Resolution_histogram': resolution_hist_name}
        CTF_estimation['Defocus_(Å)'] = defocus
        CTF_estimation['Resolution_(Å)'] = resolution
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
        plt.ylabel(hist_ylabel_mic)
        plt.title('Histogram for particle number per micrograph')

        hist_name = 'particles_hist.png'
        particles_hist = self.hist_path(hist_name)
        plt.savefig(particles_hist)

        particles = {"Number_particles": sum(particle_counts),
                     "Particles_per_micrograph": round(mean_particles_values, 1),
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
        classes = len(sorted_list_classes)  # Number of classes

        img_classes_file = classes2D.getFirstItem().getRepresentative().getFileName()
        # Saving images in .png, drawing number of particles on them
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
                new_img_filename = splitext(img_classes_file)[0] + f'_image_{i}.png'
                image.save(new_img_filename)
                img_filenames.append(new_img_filename)

            # Creating collage in .png with all images ordered in descending order
            images = [Image.open(filename) for filename in img_filenames]

            output_folder = self._getExtraPath()  # Extra path of current protocol
            collage_filename = 'classes_2D.png'
            collage_filepath = join(output_folder, collage_filename)
            self.create_collage(images, collage_filepath)

        classes_2D = {"Number_classes_2D": classes, "Particles_per_class": particles_list,
                      "Images_classes_2D": collage_filename}

        return classes_2D

    def init_volume_generation(self):
        extra_folder = self._getExtraPath()
        initial_vol_folder_name = 'Initial_volume'
        initial_vol_folder_path = join(extra_folder, initial_vol_folder_name)
        os.makedirs(initial_vol_folder_path, exist_ok=True)

        volume = self.initVolume.get()
        volume_file = volume.getFileName()

        # Getting orthogonal slices in X, Y and Z
        # Folder to store orthogonal slices
        orthogonal_slices_folder = 'orthogonal_slices'
        orthogonal_slices_path = join(initial_vol_folder_path, orthogonal_slices_folder)
        os.makedirs(orthogonal_slices_path, exist_ok=True)

        self.orthogonalSlices(fnRoot=orthogonal_slices_path, map=volume_file)

        # Getting 3 isosurface images
        # Folder to store isosurface images
        isosurface_images_folder = 'isosurface_images'
        isosurface_images_path = join(initial_vol_folder_path, isosurface_images_folder)
        os.makedirs(isosurface_images_path, exist_ok=True)

        th = int(self.threshold_initVol.get())

        volume_file_abspath = abspath(volume_file)
        front_view_img = 'front_view.png'
        side_view_img = 'side_view.png'
        top_view_img = 'top_view.png'
        self.generate_isosurfaces(isosurface_images_path, volume_file_abspath,
                                  th, front_view_img, side_view_img, top_view_img)

        init_volume = {'Orthogonal_slices': {
            'Orthogonal_slices_X': join(initial_vol_folder_name, orthogonal_slices_folder, "orthogonal_slices_X.png"),
            'Orthogonal_slices_Y': join(initial_vol_folder_name, orthogonal_slices_folder, "orthogonal_slices_Y.png"),
            'Orthogonal_slices_Z': join(initial_vol_folder_name, orthogonal_slices_folder, "orthogonal_slices_Z.png")
        },
            'Isosurface_images': {
                'Front_view': join(initial_vol_folder_name, isosurface_images_folder, front_view_img),
                'Side_view': join(initial_vol_folder_name, isosurface_images_folder, side_view_img),
                'Top_view': join(initial_vol_folder_name, isosurface_images_folder, top_view_img)
            }}

        return init_volume

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

        # Saving images in .png, drawing number of particles on them
        particles_list = []
        img_filenames = []

        Volumes = {}  # Dictionary to store image path for each volume

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
                new_img_filename = splitext(file_name)[0] + '.png'
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
                front_view_img = 'front_view.png'
                side_view_img = 'side_view.png'
                top_view_img = 'top_view.png'
                self.generate_isosurfaces(isosurface_images_path, volume_file_abspath,
                                          th, front_view_img, side_view_img, top_view_img)

                # Dictionary fill in:
                volume = {
                    "Orthogonal_slices": {
                        "Orthogonal_slices_X": join(classes_3D_folder_name, orthogonal_slices_folder,
                                                    "orthogonal_slices_X.png"),
                        "Orthogonal_slices_Y": join(classes_3D_folder_name, orthogonal_slices_folder,
                                                    "orthogonal_slices_Y.png"),
                        "Orthogonal_slices_Z": join(classes_3D_folder_name, orthogonal_slices_folder,
                                                    "orthogonal_slices_Z.png")},
                    'Isosurface_images': {
                        'Front_view': join(classes_3D_folder_name, isosurface_images_folder, front_view_img),
                        'Side_view': join(classes_3D_folder_name, isosurface_images_folder, side_view_img),
                        'Top_view': join(classes_3D_folder_name, isosurface_images_folder, top_view_img)
                    }}

                Volumes_key = f'Volume_{i + 1}'
                Volumes[Volumes_key] = volume

        # Creating collage in .png with all images ordered in descending order
        images = [Image.open(filename) for filename in img_filenames]
        collage_filename = 'classes_3D.png'
        collage_filepath = join(classes3D_folder_path, collage_filename)
        self.create_collage(images, collage_filepath)

        classes_3D = {"Number_classes_3D": classes, "Particles_per_class": particles_list,
                      "Images_classes_3D": join(classes_3D_folder_name, collage_filename),
                      "Volumes": Volumes}
        return classes_3D

    def saveJson(self):
        file_path = self.getOutFile()
        with open(file_path, 'w', encoding='utf-8') as json_file:
            json.dump(self.processing_json, json_file, ensure_ascii=False, indent=4)
        print(f"JSON data successfully saved to {file_path}")

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

        collage_X_path = join(fnRoot, 'orthogonal_slices_X.png')
        collage_Y_path = join(fnRoot, 'orthogonal_slices_Y.png')
        collage_Z_path = join(fnRoot, 'orthogonal_slices_Z.png')

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