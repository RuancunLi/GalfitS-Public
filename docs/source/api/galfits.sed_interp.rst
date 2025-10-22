galfits.sed_interp
==================

.. automodule:: galfits.sed_interp
   :members:
   :undoc-members:
   :show-inheritance:

SED Interpolation Functions
----------------------------

Stellar Populations
^^^^^^^^^^^^^^^^^^^

.. autofunction:: get_host_SED
.. autofunction:: get_host_SED_3D
.. autofunction:: get_host_SED_withL

AGN Components
^^^^^^^^^^^^^^

.. autofunction:: get_AGN_SED
.. autofunction:: starburst_disk
.. autofunction:: cat3Dtorus
.. autofunction:: FeII
.. autofunction:: BaC

Nebular Emission
^^^^^^^^^^^^^^^^

.. autofunction:: get_nebular

Dust Emission
^^^^^^^^^^^^^

.. autofunction:: get_cold_dust_sed
.. autofunction:: get_cold_dust_L

Stellar SEDs
^^^^^^^^^^^^

.. autofunction:: get_stellar_sed

Coordinate Transformations
---------------------------

.. autofunction:: sed_to_obse
.. autofunction:: sed_to_rest

Cosmology
---------

.. autofunction:: luminosity_distance
.. autofunction:: kpc_per_arcsec
.. autofunction:: cosmo_age

Spectral Processing
-------------------

.. autofunction:: conv_spec

Attenuation Curves
------------------

.. autofunction:: mycigale
.. autofunction:: cigale
.. autofunction:: calzetti2000
.. autofunction:: ccm89

IGM Absorption
--------------

.. autofunction:: igm_transmission

Star Formation History
-----------------------

.. autofunction:: f_cont_to_logSFR
.. autofunction:: logSFR_to_f_cont

Table Models
------------

.. autoclass:: TableModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autofunction:: read_table_model

Simple Models
-------------

.. autofunction:: spowerlaw
.. autofunction:: brokenpower
.. autofunction:: Bbody

Helper Functions
----------------

.. autofunction:: gaussian_jax
.. autofunction:: gaussianL_jax3D
.. autofunction:: nearest_neighbors_jax
