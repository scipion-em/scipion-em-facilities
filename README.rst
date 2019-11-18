
pyworkflow
===========

**pyworkflow** is a simple workflow platform used in scientific applications. It has been initially developed
within the Scipion framework for image processing in Electron Microscopy.
 
**pyworkflow** allows you to define a working Domain by defining the following group of classes:

  * Objects: input/outputs that will be generated by different programs
  * Protocols: special type of Objects that have defined input parameters and will produce some output
  * Viewers: Entities that provide graphical analysis of results.
  * Wizards: Small GUIs that can be develop to help users to select given parameter values.

Development
-------------
We are now fully going toward Python 3!

To install **pyworkflow** for development purposes, one can do:

.. code-block:: bash

    # Create a clean virtual environment
    python -m venv ~/myenv
    source ~/myenv/bin/activate
    git clone git@github.com:scipion-em/scipion-pyworkflow.git
    cd scipion-pyworkflow
    python -m pip install -e .  # Install in the environment as development

Running tests
.............
First make sure that **pyworkflow** is available as a Python module in your
current Python environment. During development, I tend to set the PYTHONPATH:

.. code-block:: bash

    cd scipion-pyworkflow
    # Either you have installed as mentioned above, or modify the PYTHONPATH
    export PYTHONPATH=$PYTHONPATH:$PWD
    # After pyworkflow is accesible as a module, then:
    cd pyworkflow/tests

    python -m unittest discover

    # Simple project GUI can be shown after running tests:
    cd scipion-pyworkflow

    # At the moment you need to specify SCIPION_DOMAIN and SCIPION_VERSION
    export SCIPION_DOMAIN=scipion-pyworkflow/pyworkflow/tests/mock_domain
    export SCIPION_VERSION=3.0.0

    python pyworkflow/apps/pw_project.py TestProtocolOutputs


Installation
............

Using a python virtual environment you might need:
# For virtual env in ubuntu:
# sudo apt-get install python3-dev
# sudo apt-get install python3-tk