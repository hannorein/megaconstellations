# Visibility Predictions for Near-Future Satellite Megaconstellations: Latitudes near 50 Degrees will Experience the Worst Light Pollution

This repository contains the code to reproduce the figures of [S. M. Lawler, A. C. Boley, H. Rein (2021)](https://arxiv.org/abs/2109.04328) and the code for the webapp currently running at http://megaconstellations.hanno-rein.de.


## Directories

### `/data`
The `/data` directory contains two data files. One lists 23 observations of starlink satellites. The other contains a list of all bright stars in the Hipparcos catalogue.

### `/model`
The `mega.py` file in this directory includes the code for out model. It contains data for all satellite constallations, load the data into [REBOUND](https://github.com/hannorein/rebound), and then calculates the location of the satellites in the night sky given a location, date, and time. 

### `/paper`
This directory includes Jupyter notebooks to recreate the figures of [S. M. Lawler, A. C. Boley, H. Rein (2021)](https://arxiv.org/abs/2109.04328).
The Jupyter notebook make use of the data files in `/data` and the model in `/model/mega.py`.

### `/webapp`
This directory contains the code required to run the interactive webapp running at http://megaconstellations.hanno-rein.de.
It requires flask, rebound, numpy, and matplotlib as well as the data files in `/data` and the model in `model/mega.py` on the server side. It uses jquery on the client side.
Note that when run for the first time, the app creates several binary REBOUND files to later reload the satellite data more quickly. These binary files are only a cache and can be deleted at any time.

## iOS App
You might also be interested in the iOS app [Mega Constellations](https://apps.apple.com/us/app/mega-constellations/id1598820453).

## Authors
- [Hanno Rein](mailto:hanno.rein@utoronto.ca)
- Samantha M. Lawler
- Aaron C. Boley

## License
This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 2.0 Generic License](http://creativecommons.org/licenses/by-nc-sa/2.0/).

