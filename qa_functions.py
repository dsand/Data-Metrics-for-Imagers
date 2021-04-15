import numpy as np
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

def CumGauss1D(x, mean=0.0, stddev=1.0):
    return 0.5*(1.0+erf((x-mean) / (1.414213562*stddev)))



def measure_bg_from_image(img, sampling=10, value_only=False, gaussfit=True):
    """
    Return background value, and its std deviation, as measured directly
    from pixels in the image
    """
    for ext in img:
        flags = None
        bg_data = ext.data.ravel()
        bg_data = bg_data[::sampling]

        if len(bg_data) > 0:
            if gaussfit:
                # An ogive fit is more robust than a histogram fit
                bg_data = np.sort(bg_data)
                bg = np.median(bg_data)
                bg_std = 0.5*(np.percentile(bg_data, 84.13) -
                              np.percentile(bg_data, 15.87))
                g_init = CumGauss1D(bg, bg_std)
                fit_g = fitting.LevMarLSQFitter()
                g = fit_g(g_init, bg_data, np.linspace(0.,1.,len(bg_data)+1)[1:])
                bg, bg_std = g.mean.value, abs(g.stddev.value)
            else:
                # Sigma-clipping will screw up the stats of course!
                bg_data = sigma_clip(bg_data, sigma=2.0, maxiters=2)
                bg_data = bg_data.data[~bg_data.mask]
                bg = np.median(bg_data)
                bg_std = np.std(bg_data)
        else:
            bg, bg_std = None, None

        if value_only:
            output_list.append(bg)
        else:
            output_list.append((bg, bg_std, len(bg_data)))

    return output_list
