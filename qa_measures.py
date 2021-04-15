import numpy as np





Measurement = namedtuple('Measurement', 'value std samples')



class QA():

    #  Class containing QA Measurements  #

    def measureBg(self, img):

        filter=get_filter(img)
        pixscale=get_pixscale(img)
        exptime=get_exptime(img)
        bg_list = gt.measure_bg_from_image(ad, sampling=100, gaussfit=False)


        if npz is not None:
            if bg_count.value > 0:
                # convert background to counts/arcsec^2/second, but
                # want to preserve values of sci_bg and sci_std
                fak = 1.0 / (exptime * pixscale * pixscale)
                bg_mag = Measurement(npz - 2.5*math.log10(bg_count.value*fak),
                    2.5*math.log10(1 + bg_count.std/bg_count.value),
                                     bg_count.samples)
                # Need to report to FITSstore in electrons
                bg_mag_list.append(bg_mag)
                qastatus = _get_qa_band('bg', ad, bg_mag, bg_band_limits)
            else:
                print("Background is less than or equal to 0 "
                                    "for {}:{}".format(filename,extver))



    def measureCC(self, img):



        for ext in img():
            nom_phot_zpt = ad.nominal_photometric_zeropoint()
            nom_at_ext = ad.nominal_atmospheric_extinction()
            if nom_at_ext is None:
                log.warning("Cannot get atmospheric extinction. Assuming zero.")
                nom_at_ext = 0.0
            exptime = get_exptime()

            all_zp = []
            all_cloud = []
            info_list = []

            try:
                objcat = get_OBJCAT()
            except AttributeError:
                print("No OBJCAT in {}:{}".format(img,extver))
                    all_zp.append(Measurement(None, None, 0))
                    continue

            # Incrementally cull the catalog: remove sources without mags
            good_obj = objcat[~np.logical_or(objcat['MAG_AUTO'] == -999,
                                                 objcat['MAG_AUTO'] > 90)]
            if len(good_obj) == 0:
                print("No magnitudes found")
                all_zp.append(Measurement(None, None, 0))
                continue

            # Remove sources without reference mags
            good_obj = good_obj[~np.logical_or.reduce(
                [good_obj['REF_MAG'] == -999, np.isnan(good_obj['REF_MAG']),
                 np.isnan(good_obj['REF_MAG_ERR'])])]
            if len(good_obj) == 0:
                print("No reference magnitudes")
                all_zp.append(Measurement(None, None, 0))
                continue

            # Sources must be free of SExtractor flags and unsaturated, and
            # <2% of pixels be otherwise flagged (typically bad/non-linear)
            good_obj = good_obj[np.logical_and.reduce([good_obj['FLAGS'] == 0,
                    good_obj['NIMAFLAGS_ISO'] < 0.02*good_obj['ISOAREA_IMAGE'],
                    good_obj['IMAFLAGS_ISO'] & DQ.saturated == 0])]

            zps = good_obj['REF_MAG'] - nom_at_ext - (good_obj['MAG_AUTO'] +
                                                     2.5*math.log10(exptime))
            zperrs = np.sqrt(good_obj['REF_MAG_ERR']**2 +
                             good_obj['MAGERR_AUTO']**2)

            # There shouldn't be any NaN left
            assert sum(np.logical_or(np.isnan(zps), np.isnan(zperrs))) == 0

                # TODO: weight instead?
            # Trim out where zeropoint error > err_threshold
            if len([z for z in zps if z is not None]) <= 5:
                # 5 sources or less. Beggars are not choosers.
                ok = zperrs<0.2
            else:
                ok = zperrs<0.1

            # Ensure these are regular floats for JSON (thanks to PH)
            zps = [Measurement(float(zp), float(zperr), 1) for zp, zperr
                   in zip(zps[ok], zperrs[ok])]

            if len(zps) == 0:
                print("No good photometric sources found")
                all_zp.append(Measurement(None, None, 0))
                continue

            # Collapse all the Measurements to a single value + error
            if len(zps) > 2:
                # TODO: need better than the 1-sigma clip!
                stats = _stats(zps)
                m, s = stats.value, stats.std
                zps = [z for z in zps if abs(z.value - m) < s]

            ext_zp = _stats(zps, weights='variance') if len(zps)>1 else zps[0]

            # Report average extinction measurement
            ext_cloud = _arith(_arith(ext_zp, 'sub', npz), 'mul', -1)

            # Individual extinction measurements for all sources
            all_cloud.extend([_arith(_arith(zp, 'sub', npz), 'mul', -1)
                               for zp in zps])
            all_zp.append(ext_zp)


            # Only if we've managed to measure at least one zeropoint
        if any(zp.value for zp in all_zp):
            avg_cloud = _stats(all_cloud, weights='variance') # Extinction
            qastatus = _get_qa_band('cc', ad, avg_cloud, qa.ccBands, simple=False)


        else:
            print("    Could not measure zeropoint - no catalog sources associated")

        return adinputs



    def measureIQ(self, img):
        
               
            iq_overlays = []
            measure_iq = True

            # Check that the data is not an image with non-square binning
            xbin = get_x_bin(img)
            ybin = get_y_bin(img)
            if xbin != ybin:
                print("No IQ measurement possible, image {} is {} x {} binned data".format(filename, xbin, ybin))
                measure_iq = False

            good_source = gt.clip_sources(adiq)


            # Check for no sources found: good_source is a list of Tables
            if len(good_source) <1:
                measure_iq=false

            if measure_iq:
                # Descriptors and other things will be the same for ad and adiq
                try:
                    zcorr = get_airmass()**(-0.6)
                except:
                    zcorr = None


                info_list = []

                for src in good_source:
                    ellip = Measurement(None, None, 0)
                        
                    # Weighted mean of clipped FWHM and ellipticity
                    if "weight" in src.columns:
                        mean_fwhm = np.average(src["fwhm_arcsec"],
                                                  weights=src["weight"])
                        std_fwhm = np.sqrt(np.average((src["fwhm_arcsec"] -
                                   mean_fwhm)**2, weights=src["weight"]))
                    else:
                        mean_fwhm = np.mean(src["fwhm_arcsec"])
                        std_fwhm = np.std(src["fwhm_arcsec"])
                    fwhm = Measurement(float(mean_fwhm), float(std_fwhm),
                                           len(src))
                    ellip = Measurement(float(np.mean(src['ellipticity'])),
                    float(np.std(src['ellipticity'])), len(src))

                # Find the corrected FWHM. For AO observations, the IQ
                # constraint band is taken from the AO-estimated seeing
                # except for GSAOI, which has some magic formula that kind of works
                iq = fwhm

                if zcorr:
                    zfwhm = _arith(iq, 'mul', zcorr)
                else:
                    Print('Airmass not found, not correcting to zenith')
                    zfwhm = Measurement(None, None, 0)

        return zfwhm







##
#   Helper Functions
##

def _stats(stats_list, weights='sample'):
    """
    Estimates overall mean and standard deviation from measurements that have
    already been compressed, so the original data don't exist
    Parameters
    ----------
    stats_list: list of Measurements
        The input statistics
    weights: 'variance'/'sample'/None
        how to weight the measurements
    Returns
    -------
    Measurement: mean, standard deviation, total number of measurements
    """
    try:
        use_list = [m for m in stats_list if m.value is not None]
        if weights == 'variance':
            wt = [1.0 / (m.std * m.std) for m in use_list]
        elif weights == 'sample':
            wt = [m.samples for m in use_list]
        else:
            wt = [1.0] * len(use_list)
        total_samples = sum(m.samples for m in use_list)
        mean = np.average([m.value for m in use_list], weights=wt)
        var1 = np.average([(m.value - mean)**2 for m in use_list],
                                   weights = wt)
        var2 = sum(w*m.std*m.std for w, m in zip(wt, use_list))/sum(wt)
        sigma = np.sqrt(var1 + var2)
    except:
        return Measurement(None, None, 0)
    return Measurement(mean, sigma, total_samples)


def _arith(m, op, operand):
    """Performs an arithmetic operation on a value and its uncertainty"""
    if op in ['mul', 'div', 'truediv']:
        return Measurement(getattr(operator, op)(m.value, operand),
                getattr(operator, op)(m.std, abs(operand)) if m.std else m.std,
                m.samples)
    else:
        return Measurement(getattr(operator, op)(m.value, operand),
                           m.std, m.samples)
