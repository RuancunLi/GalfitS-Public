import numpy as np
import dynesty
import emcee

__all__ = ['fitter']

def fitter(object):
    '''
    Fit the galaxy model according to image data.
    '''
    def __init__(self, data=None, model=None):
        '''
        Parameters
        ----------
        data (optional) : `image_atlas`
            Image data.
        model (optional) : `galaxy_model`
            Galaxy model.
        '''
        self.data = data
        self.model = model
        self.parnames = model.get_parnames()  # Get list of parameter names
        self.ndim = len(self.parnames)

    def lnprob(self, params):
        '''
        The probability of the model and data.
        '''
        prior = self.logprior(params)
        if not np.isclose(prior, -1e25):
            return prior + self.loglike(params)
        else:
            return -1e25

    def loglike(self, params):
        '''
        The likelihood of the spectrum and phase data.

        Parameters
        ----------
        params : 1D array
             The parameters of the model.

        Returns
        -------
        lnlike : float
            Likelihood in log scale.
        '''
        chisqList = []
        for loop, d in enumerate(self.data):
            data_info = d.get_data_info()  # This should be a dict; Need to define the content.
            image = d.data  # This is a numpy masked array
            sigma_image = d.get_sigma_image()  # Need to add this function.
            model_image = model[loop].generate_image(params, data_info)
            chisq = -0.5 * np.sum((image - model_image)**2 / sigma_image**2 +
                                  np.log(2 * np.pi * sigma_image**2))
            chisqList.append(chisq)
        return np.sum(chisqList)

    def logprior(self, params):
        '''
        The prior of the model.

        Parameters
        ----------
        params : list
            Parameters from sampler.

        Returns
        -------
        logp : float
            Value of total log prior.
        '''
        logpList = []
        # For physical parameters
        pDict = self.__prior
        for loop, pn in enumerate(self.__parnames):
            if pDict[pn]['type'] == 'uniform':
                logpList.append(logprior_uniform(params[loop], pDict[pn]['bounds']))
            if pDict[pn]['type'] == 'gaussian':
                logpList.append(logprior_gaussian(params[loop],
                                                  center=pDict[pn]['center'],
                                                  stddev=pDict[pn]['stddev']))
        logp = np.sum(logpList)
        return logp

    def ptform(self, u):
        '''
        Transforms samples "u" drawn from unit cube to samples to prior.
        '''
        params = []
        # Model parameters
        for loop, pn in enumerate(self.__parnames):
            prior = self.__prior[pn]
            if prior['type'] == 'uniform':  # Uniform prior
                params.append(ptform_uniform(u[loop], prior))
            else:  # Gaussian prior
                params.append(ptform_gaussian(u[loop], prior))
        params = np.array(params)
        return params


def logprior_uniform(p, bounds):
    '''
    Calculate the uniform prior given the parameter information.

    Parameters
    ----------
    p : float
        The proposed parameter value.
    bounds : tuple
        The lower and upper bound of the uniform distribution.

    Returns
    -------
    The ln prior, either 0 or -inf.
    '''
    if (p >= bounds[0]) & (p <= bounds[1]):
        return 0.
    else:
        return -1e25


def logprior_gaussian(p, center, stddev):
    '''
    Calculate the Gaussian prior probability given the parameter information.

    Parameters
    ----------
    p : float
        The proposed parameter value.
    center : float
        The center of the Gaussian distribution.
    stddev : float
        The standard deviation of the Gaussian distribution.

    Returns
    -------
    The ln prior, Gaussian distribution.
    '''
    return np.log(scp_stats.norm.pdf(p, loc=center, scale=stddev))


def ptform_uniform(u, prior):
    '''
    Transform the Uniform(0, 1) sample to the model parameter value with the
    uniform prior.

    Parameters
    ----------
    u : float
        The random number ~Uniform(0, 1).
    prior : dict
        The prior of the parameter, prior=dict(bound=[a, b]).
    '''
    pmin, pmax = prior['bounds']
    return (pmax - pmin) * (u - 0.5) + (pmax + pmin) * 0.5


def ptform_gaussian(u, prior):
    '''
    Transform the Uniform(0, 1) sample to the model parameter value with the
    Gaussian prior.

    Parameters
    ----------
    u : float
        The random number ~Uniform(0, 1).
    prior : dict
        The prior of the parameter, prior=dict(center=a, stddev=b).
    '''
    return scp_stats.norm.ppf(u, loc=prior['center'], scale=prior['stddev'])
