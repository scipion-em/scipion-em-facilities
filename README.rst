====
scipion4facilities
====

**scipion4facilities** is a Python module of Scipion framework for cryo-EM facilities


There is a type of protocols called monitors which are used to produce live analysis plots, generate reports or raise alerts when some problems are detected. A monitor example is the CTF-monitor, that checks the computed defocus values for each micrograph as they are generated. CTF-monitor may raise an alert if the defocus values are above or below certain thresholds. A special case of this monitors is the monitor summary which encapsulates the CTF Monitor, the system monitor and the movie gain monitor and continuosly creates a report.

This module contains protocols and utilities related with monitors.

See https://scipion-em.github.io/docs/release-3.0.0/docs/facilities/customize-report.html for more information

-------------
Development
-------------

To install **scipion4facilities** for development purposes, one can do:

::

    # Create a clean virtual environment
    python -m venv ~/myenv
    source ~/myenv/bin/activate
    git clone git@github.com:scipion-em/scipion-em-facilities.git
    cd scipion-em
    python -m pip install -e .  # Install in the environment as development

-------------
Running tests
-------------

First make sure that **scipion4facilities** is available as a Python module in your
current Python environment. During development, I tend to set the PYTHONPATH:

::

    cd scipion-em
    # Either you have installed as mentioned above, or modify the PYTHONPATH
    export PYTHONPATH=$PYTHONPATH:$PWD
    # After pyworkflow is accesible as a module, then:
    cd scipion4facilities/tests

    python -m unittest discover

