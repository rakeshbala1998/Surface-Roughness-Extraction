This MATLAB code extracts various surface roughness features from the data file having x and y points representing the surface profile.
In this code, the sliding window feature extracts the properties in the window. And this script can find the valley accurately using adaptive index searching and forming equal size negative and positive slopes to get a more accurate valley root radius. This code can help us to measure the root radius from CT- images of the surface instead of conventional NDT such as profilometer measurement.

The following features can be calculated using this MATLAB script:
-Average Surface Roughness
-Difference between the highest peak and lowest valley
-Ten-point Roughness
-Root radius
