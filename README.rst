=========================
Scipion4facilities plugin
=========================

This plugin allows to use different utils for cryo-EM facilities (like monitors)
within the Scipion framework.

+------------------+------------------+
| stable: |stable| | devel: | |devel| |
+------------------+------------------+

.. |stable| image:: http://scipion-test.cnb.csic.es:9980/badges/facilities_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/facilities_sdevel.svg

Please, check the `facilities documentation <https://scipion-em.github.io/docs/docs/facilities/facilities.html>`_ for more details.


Installation
------------

You will need to use `Scipion 3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_
to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-facilities

b) Developer's version

   * download repository

   .. code-block::

      git clone https://github.com/scipion-em/scipion-em-facilities.git

   * install

   .. code-block::

      scipion installp -p /path/to/scipion-em-facilities --devel

Testing
-------

To check the installation, simply run the following Scipion test:

``scipion test --grep scipion4facilities --mode modules --run``

If '--run' is not passed, it only shows the available tests to check.