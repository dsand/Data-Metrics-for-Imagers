# Real Time Data Metrics for Steward's Small Telescope Imagers

The goal of this code is to support queue and classical users and support staff in operating and maintaining Steward's small telescope imagers, specifically the Kuiper+Mont4k and Bok+90Prime, with possible expansion to the VATT at a later time.

Real time image quality metrics will be provided for each science and calibration image taken.  Necessary statistics include:
1. Sky brightness level.
2. Transparency, potentiall in the form of a delta-zeropoint from nominal 'clear skies' value to roughly dtermine image depth.
3. Image quality, in the form of image FWHM and ellipticity, to aid in observing decision-making.  This should be measured for stars across the field of view of the imager.
4. Pointing accuracy, as determined through a WCS solution of the image (or perhaps something simpler)

These metrics should be saved in a databased and presented in a clear visual format to the classical/queue observer.  This database should also have access to TCS and weather values for later analysis.  Data may also be placed in image headers, but perhpas that is for version 2 of the code.

Over time, tracking of these values will help to assess the overall health of the imager.  

Real time metrics of calibration data, specifically bias and flat field frames, will also be important.  Bias frames should be monitored for stability and noise properties.  Flat fields for count values, noise and broad changes over time.
