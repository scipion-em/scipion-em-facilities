======================
Scipion for facilities
======================

'scipion-em-facilities' plugin allows to use different utils for cryo-EM facilities (like monitors)
within the Scipion framework.

+------------------+------------------+
| stable: |stable| | devel: | |devel| |
+------------------+------------------+

.. |stable| image:: http://scipion-test.cnb.csic.es:9980/badges/facilities_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/facilities_sdevel.svg


Installation
------------

You will need to use `Scipion 3.0 <https://scipion-em.github.io/docs/release-3.0.0/index.html>`_
to be able to run these protocols.

To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion3 installp -p scipion-em-facilities

b) Developer's version

   * 1st. Download repository

   .. code-block::

      git clone https://github.com/scipion-em/scipion-em-facilities.git /path/to/scipion-em-facilities

   * 2nd. Install the plugin

   .. code-block::

      scipion3 installp -p /path/to/scipion-em-facilities --devel

Testing
-------

To check the installation, simply run the following Scipion test:

.. code-block::

  scipion3 test --grep emfacilities --mode modules --run

If `--run` is not passed, it only shows the available tests to check.