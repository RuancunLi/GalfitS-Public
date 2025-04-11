Run GalfitS
===========

To execute GalfitS, use the following command in your terminal:

.. code-block:: bash

    galfits config.lyric --work ./result

This command launches the GalfitS fitting process using a configuration file named ``config.lyric`` and directs the output to the ``./result`` directory. Below, you'll find a comprehensive breakdown of all available command-line arguments to customize the execution of GalfitS.

Command-Line Arguments
----------------------

GalfitS provides a wide range of options to control its behavior. These arguments are passed via the command line and allow you to configure the fitting process, output settings, and more. Each argument is detailed below, including its type, default value, and purpose.

### General Options

- **--config** (type: str, default: "")
  Specifies the path to the configuration file for GalfitS. This file defines the settings and parameters used during the fitting process.

  *Example:* ``--config myconfig.lyric``

- **--workplace** (type: str, default: "./")
  Defines the directory where GalfitS will store its output files.

  *Example:* ``--workplace ./output``

- **--notfit** (action: 'store_true')
  If provided, GalfitS will skip the fitting process and only generate the parameter file.

  *Example:* ``--notfit``

- **--savefull_results** (action: 'store_true')
  When specified, GalfitS saves the full results of the fitting process rather than a summary.

  *Example:* ``--savefull_results``

- **--saveimgs** (action: 'store_true')
  If set, GalfitS saves the generated images from the fitting process.

  *Example:* ``--saveimgs``

### Fitting Method Selection

- **--fit_method** (type: str, default: "optimizer")
  Selects the fitting method to use. Available options are:
  - ``chisq``: Chi-square minimization.
  - ``flowmc``: FlowMC (Markov Chain Monte Carlo method).
  - ``dynesty``: Dynesty (nested sampling method).
  - ``optimizer``: General optimizer (default).
  - ``ES``: Evolutionary Strategy.

  *Example:* ``--fit_method dynesty``

### Optimizer Fitting Options

- **--num_steps** (type: int, default: 3000)
  Sets the number of steps for the optimizer fitting method.

  *Example:* ``--num_steps 5000``

- **--learning_rate** (type: float, default: 0.0008)
  Specifies the learning rate for the optimizer fitting method.

  *Example:* ``--learning_rate 0.001``

- **--baysian** (action: 'store_true')
  Enables the use of Bayesian priors for the optimizer fitting method.

  *Example:* ``--baysian``

- **--cal_sigma** (action: 'store_true')
  If specified, calculates the Hessian sigma for the optimizer fitting method.

  *Example:* ``--cal_sigma``

### Chi-Square Fitting Options

- **--chimethod** (type: str, default: "scipy")
  Defines the method for chi-square fitting. Currently, the only option is:
  - ``scipy``: Uses SciPy's optimization routines.

  *Example:* ``--chimethod scipy``

- **--scimethod** (type: str, default: "SLSQP")
  Specifies the SciPy optimization method to use (e.g., "SLSQP", "L-BFGS-B").

  *Example:* ``--scimethod L-BFGS-B``

- **--constrain** (action: 'store_true')
  Applies constraints during chi-square fitting if specified.

  *Example:* ``--constrain``

### FlowMC Fitting Options

- **--nchains** (type: int, default: 30)
  Sets the number of chains for the FlowMC method.

  *Example:* ``--nchains 50``

- **--rstep** (type: float, default: 0.1)
  Defines the parameter step size for FlowMC.

  *Example:* ``--rstep 0.05``

- **--nlocalsteps** (type: int, default: 3000)
  Specifies the number of local steps for FlowMC.

  *Example:* ``--nlocalsteps 4000``

- **--nlooptraining** (type: int, default: 20)
  Sets the number of training loops for FlowMC.

  *Example:* ``--nlooptraining 30``

- **--nloopproduction** (type: int, default: 10)
  Defines the number of production loops for FlowMC.

  *Example:* ``--nloopproduction 15``

- **--sampler** (type: str, default: "GRW")
  Selects the sampler for FlowMC (e.g., "GRW" for Gaussian Random Walk).

  *Example:* ``--sampler GRW``

### Dynesty Fitting Options

- **--nlive** (type: int, default: 80)
  Sets the number of live points for Dynesty.

  *Example:* ``--nlive 100``

- **--maxiters** (type: int, default: 100000)
  Specifies the maximum number of iterations for Dynesty.

  *Example:* ``--maxiters 50000``

- **--dlogz** (type: float, default: 0.02)
  Defines the stopping criterion for Dynesty based on the change in log-evidence.

  *Example:* ``--dlogz 0.01``

- **--sample_method** (type: str, default: "rwalk")
  Selects the sampling method for Dynesty (e.g., "rwalk", "slice").

  *Example:* ``--sample_method slice``

- **--dynamic** (action: 'store_true')
  Enables dynamic nested sampling for Dynesty if specified.

  *Example:* ``--dynamic``

- **--maxbatch** (type: int, default: 10)
  Sets the maximum number of batches for Dynesty's dynamic sampling.

  *Example:* ``--maxbatch 20``

### Evolutionary Strategy (ES) Options

- **--num_generations** (type: int, default: 100Pocket)
  Defines the number of generations for the Evolutionary Strategy method.

  *Example:* ``--num_generations 15000``

- **--popsize** (type: int, default: 20)
  Sets the population size for the ES method.

  *Example:* ``--popsize 30``

### Parameter and Prior Options

- **--parconstrain** (type: str, default: "default")
  Specifies the file path for parameter constraints.

  *Example:* ``--parconstrain constraints.txt``

- **--readpar** (type: str, default: "default")
  Defines the file path to read initial parameters from.

  *Example:* ``--readpar params.txt``

- **--readsummary** (type: str, default: "default")
  Specifies the file path to read a summary file from.

  *Example:* ``--readsummary summary.txt``

- **--priorpath** (type: str, default: None)
  Specifies the path to the astrophysical prior file.

  *Example:* ``--priorpath priors.txt``

- **--weight_spec** (type: float, default: 1.0)
  Sets the weight of the spectrum in the fitting process.

  *Example:* ``--weight_spec 0.8``

### AGN-Related Options

- **--fixagnlinepro** (action: 'store_true')
  If specified, fixes the AGN line profile during fitting.

  *Example:* ``--fixagnlinepro``

- **--ndisagn** (action: 'store_true')
  If specified, prevents displaying the AGN in the output image.

  *Example:* ``--ndisagn``

Usage Example
-------------

Here’s an example of a more customized GalfitS command:

.. code-block:: bash

    galfits config.lyric --work ./output --fit_method dynesty --nlive 120 --maxiters 20000 --saveimgs

This command runs GalfitS with the Dynesty fitting method, using 120 live points, a maximum of 20,000 iterations, and saves the resulting images to the ``./output`` directory.

Conclusion
----------

These command-line arguments give you extensive flexibility to tailor GalfitS to your specific needs, whether you're performing a quick optimization or a detailed Bayesian analysis. Adjust the parameters as necessary to suit your data and research goals.

