Installation
============

Prerequisite Packages
---------------------

**JAX**

GalfitS utilizes `JAX <https://github.com/google/jax>`_ as its primary computational backend.

**CPU installation** (recommended for general platforms):

.. code-block:: sh

   pip install -U "jax[cpu]"

**GPU installation (NVIDIA CUDA)**

JAX computations are significantly accelerated when using NVIDIA GPUs.

- If CUDA is already installed locally:

.. code-block:: sh

   pip install -U "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

- Or directly install the CUDA-enabled JAX version from pip:

.. code-block:: sh

   pip install -U "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

In both cases, a compatible NVIDIA GPU driver is required.

**GPU installation (Apple Silicon)**

JAX also supports acceleration using Apple GPUs. Follow Apple's official instructions:
`Apple JAX Installation Guide <https://developer.apple.com/metal/jax/>`_

Other Requirements
------------------

Navigate to the GalfitS directory and install additional dependencies:

.. code-block:: sh

   cd GalfitS
   pip install -r requirements.txt

All SED templates required by GalfitS are stored on Google Drive.

Download GalfitS
---------------

Clone the main branch of the GalfitS repository from GitHub:

.. code-block:: sh

   git clone https://github.com/RuancunLi/GalfitS-Public

Navigate to the cloned repository:

.. code-block:: sh

   cd GalfitS

Install the GalfitS package and set up the alias and environment variable:

.. code-block:: sh

   python setup.py /path/to/your/data

This command will:

- Install the GalfitS package.
- Create an alias named ``galfits`` for easy command-line access.
- Set the ``GS_DATA_PATH`` environment variable in your shell configuration file (e.g., ``.bashrc`` or ``.zshrc``).

To activate these changes, restart your terminal or source your shell configuration:

.. code-block:: sh

   source ~/.bashrc  # or ~/.zshrc