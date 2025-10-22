Release Notes
=============

Version 1.0 (Public Release)
----------------------------

**Release Date:** [Date of release]

This is the initial public release of GalfitS.

**New Features:**

*   **Core Fitting Engine**:
    *   Simultaneous morphological and SED fitting of galaxy components.
    *   Support for multiple imaging bands and optional spectroscopic data.
    *   JAX-based backend for high-performance computing on CPUs and GPUs.

*   **Modeling Capabilities**:
    *   Library of 2D profiles: Sersic, Ferrer, Edge-on Disk, Gaussian Ring, and more.
    *   Fourier modes for non-axisymmetric features (e.g., spiral arms).
    *   Stellar Population Synthesis (SPS) models for galaxy components, supporting various Star Formation Histories (burst, continuous, binned).
    *   AGN models including continuum (power-law, thin disk), emission lines, and torus dust emission.
    *   Foreground star modeling.

*   **Optimization & Sampling**:
    *   Multiple fitting backends:
        *   Gradient-based optimizer (`jaxopt`).
        *   Nested sampling with `dynesty` for Bayesian inference and posterior estimation.
        *   Evolutionary Strategies (`evosax`) for global optimization.
    *   Flexible parameter constraint system via external files.
    *   Support for astrophysical priors (e.g., Mass-Size relation, Mass-Metallicity relation).

*   **Data Handling**:
    *   Automated image cutout, sky subtraction, and mask generation utilities.
    *   Empirical PSF construction tool.
    *   Support for `astropy` FITS files and WCS.

*   **Usability**:
    *   Command-line interface for running fits.
    *   Detailed configuration files (`.lyric`) for defining fitting runs.
    *   Plotting utilities for visualizing results (image residuals, model SEDs, SFHs).

*   **Documentation**:
    *   Initial version of the documentation, including installation, quickstart, and configuration guides.

**Known Issues:**

*   The documentation is still under development. Some sections may be incomplete.
*   Error handling for invalid configuration file parameters could be more verbose.

---

*Future releases will focus on expanding the model library, improving performance, and enhancing the user interface and documentation.*
