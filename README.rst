====
scipion4facilities
====

**scipion4facilities** is a Python module of Scipion framework for cryo-EM facilities


The entire collection is licensed under the terms of the GNU Public License,
version 3 (GPLv3).

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

