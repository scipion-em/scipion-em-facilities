from pyworkflow.tests import DataSet

DataSet(name='OSCEM_metadata', folder='OSCEM_metadata',
        files={
            'processing_json': 'processing.json',
            'movies_dir': 'movies',
            'gain_im': 'gain.mrc',
            'CTF_sqlite': 'ctfs.sqlite',
            'volume_classification3D': 'volume_class3D.mrc',
            'volume_init': 'initial_volume.mrc'})



