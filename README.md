# Real Time Data Metrics for Steward's Small Telescope Imagers

The goal of this code is to support queue and classical users and support staff in operating and maintaining Steward's small telescope imagers, specifically the Kuiper+Mont4k and Bok+90Prime, with possible expansion to the VATT at a later time.

Real time image quality metrics will be provided for each science and calibration image taken.  Necessary statistics include:
1. Sky brightness level.
2. Transparency, potentiall in the form of a delta-zeropoint from nominal 'clear skies' value to roughly dtermine image depth.
3. Image quality, in the form of image FWHM and ellipticity, to aid in observing decision-making.  This should be measured for stars across the field of view of the imager.
4. Pointing accuracy, as determined through a WCS solution of the image (or perhaps something simpler)

Will update more soon.
