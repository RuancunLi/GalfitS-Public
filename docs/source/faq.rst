FAQ
===

Frequently Asked Questions
--------------------------

**Q: What is GalfitS?**

A: GalfitS is a powerful, open-source tool for detailed galaxy decomposition. It allows for simultaneous fitting of a galaxy's morphological components (like bulge, disk, bar) and their spectral energy distributions (SEDs) across multiple wavelength bands. It is built on JAX for high performance.

**Q: What are the main features of GalfitS?**

A: 
- **Multi-band Image Fitting**: Simultaneously model galaxy structure across different wavelengths.
- **SED Modeling**: Fit physical parameters of stellar populations (age, metallicity, dust) and AGNs.
- **Flexible Models**: Use various profiles (Sersic, Ferrer, etc.) and combine them.
- **Performance**: Leverages JAX for JIT compilation and GPU/TPU acceleration.
- **Advanced Optimization**: Includes multiple fitting backends like `jaxopt`, `dynesty` (nested sampling), and evolutionary strategies.
- **Extensible**: Easily add new models or priors.

**Q: How do I install GalfitS?**

A: Please refer to the :doc:`installation` guide for detailed instructions. The key steps are installing JAX with the appropriate backend (CPU or GPU) and then installing the required Python packages.

**Q: I'm getting an error related to JAX. What should I do?**

A: Make sure you have installed the correct version of JAX for your system (CPU, NVIDIA GPU, or Apple Silicon). The installation commands can be found in the :doc:`installation` guide. If you have a GPU, ensure your drivers are up to date.

**Q: How does the configuration file work?**

A: The configuration file (e.g., `quickstart.lyric`) is the heart of a GalfitS run. It's a text file where you define everything:
- The target galaxy (name, coordinates, redshift).
- The input data (images, PSFs, masks, spectra).
- The model components (e.g., a Sersic profile for a bulge, an exponential disk).
- The initial parameters for each component, including their bounds and whether they are fixed or free.

A detailed explanation can be found in the :doc:`config_intro` section.

**Q: Can GalfitS fit spectra as well?**

A: Yes. GalfitS can perform spectrum-only fitting or a combined fitting of images and spectra. This is useful for detailed studies of stellar populations and AGN properties. See the :doc:`config_example` for an example.

**Q: What is the difference between `optimizer`, `dynesty`, and `ES` fitting methods?**

A: They are different algorithms for finding the best-fit parameters:
- **optimizer**: A fast, gradient-based method (like Adam) that finds a single best-fit solution (a local or global minimum of the chi-square). It's good for quick exploration.
- **dynesty**: A nested sampling method. It is more computationally expensive but provides full posterior distributions for all parameters and calculates the Bayesian evidence, which is useful for model comparison.
- **ES (Evolutionary Strategy)**: A population-based optimization method that is robust against getting stuck in local minima. It's a good alternative to gradient-based optimizers.

**Q: How can I contribute to GalfitS?**

A: GalfitS is an open-source project, and contributions are welcome. You can contribute by reporting bugs, suggesting new features, improving the documentation, or submitting pull requests with new code on the `GitHub repository <https://github.com/RuancunLi/GalfitS-Public>`_.
