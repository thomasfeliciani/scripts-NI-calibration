_Latest edits: 15 November 2021._

# Calibrated repulsive influence model
This repository concerns an agent-based model (ABM) of opinion dynamics. It runs an empirically calibrated, spatially explicit version of the repulsive influence model (aka “negative influence” model- or simply _RI-_ or _NI-_ model). The ABM and this repository are used in a research paper currently in preparation.

This code runs in [R 4.2.0](https://www.r-project.org/) (R Core Team, 2022).

Follows a description of the contents of this repository.


## ODD protocol.docx
This document provides a description of the ABM and its implementation. It follows the standardized “ODD”(Overview, Design concepts, and Details) protocol by Grimm et al. (2010).


## util.r
This script contains utility functions that are ancillary to the simulation model.
_Requires libraries:_ “ggplot2”_,_ “ggpubr”_,_ “ggmap”_, and_ “ggsn”_._


## simulation.r
This script is the core of the ABM, as it defines the function that executes a simulation run. An usage example is provided at the bottom of the script. 

_Requires library:_ “ggplot2”_._

_Sources_ “util.r”_._

_The script also attempts to load the file containing the calibration data necessary to run the ABM:_ “./cityData/geodata_Rotterdam.RData”_._ Note: due to weight limitations, the file is not included in the repositoy, so the script will offer to download via download link.


## plots.r
This script reproduces all figures used in the research paper, and so it also plots the simulation results.

_Requires libraries:_ “ggplot2”_,_ “gridExtra”_,_ “reshape2” _and_ “scales”_._

_The script also attempts to load the file containing the simulation results produced for the research paper,_ “./simOutput/completeDataset.RDATA”_._ Again due to weight limitations, this file is also not included in the repositoy but is available to download separately via public link. This will be provided by the script.


## ./battery/
This folder contains all R scripts and shell scripts used to run simulations in parallel on a compute cluster and to process its results. In essence, we used these scripts to execute parallel simulations of the ABM using various parameter configurations, and then to summarize the simulation results in one results data file (“./simOutput/completeDataset.RDATA”).


## ./cityData/
This folder is to store three kinds of contents, all pertaining to the empirical calibration of the ABM.

**First**, two publicly available dataset published by Statistics Netherlands (aka “CBS”):
* [Statistische gegevens per vierkant 2000-2014](https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data/kaart-van-100-meter-bij-100-meter-met-statistieken) (CBS, 2021).
* [Wijk- en Buurtkaart 2014 versie 3](https://www.nationaalgeoregister.nl/geonetwork/srv/api/records/8eee71f9-9365-490d-9484-00263c6fa35c/formatters/xsl-view?view=advanced&portalLink=) (NGR, 2017).

These dataset are not included in the repository.

**Second**, a script “calibration.r”, that we used to extract the necessary (spatial) data from the two dataset and to save them into the **third** kind of file in this folder: “./geodata_Rotterdam.RData”.
The file “geodata_Rotterdam.RData”, is the only content of this folder strictly needed for running the simulation model (see script “simulation.r” for download directions).


## References
* CBS (2021). _Kaart van 100 meter bij 100 meter met statistieken_. https://www.cbs.nl/nl-nl/dossier/nederland-regionaal/geografische-data/kaart-van-100-meter-bij-100-meter-met-statistieken
* Grimm, V., Berger, U., DeAngelis, D. L., Polhill, J. G., Giske, J., & Railsback, S. F. (2010). The ODD protocol: a review and first update. _Ecological modelling_, 221(23), 2760-2768.
* NGR (2017). _Wijk- en Buurtkaart 2014 versie 3_. https://www.nationaalgeoregister.nl/geonetwork/srv/api/records/8eee71f9-9365-490d-9484-00263c6fa35c/formatters/xsl-view?view=advanced&portalLink=
* R Core Team (2022). _R: A language and environment for statistical computing_. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.

