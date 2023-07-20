# Apollo

**APOLLO** is a gravitational wave (GW) alert broker system. Apollo consists of several Python functions that compose both a listener and scheduler. The listener employs the GCN Kafka Client for Python to receive JSON notices from the LIGO-Virgo-KAGRA (LVK) Collaboration. Once received, an alert is parsed for the GW event parameters and the attached probability skymap. The scheduler generates observable targets from this skymap for following up a GW alert using the latest Galaxy List for the Advanced Detector Era (GLADE) catalog.

---

## Usage

From the terminal, run the main script `apollo.py`:

`$ python apollo.py`

All parameters are read from the `config.ini` file. One can engage/disengage the various functions within the configuration file.

## Structure

The listener uses the wrapper function `retrieve_notice` to receive incoming LVK alerts. The JSON notice is read by the function `parse_notice` which displays the event and skymap parameters to the terminal. For testing purposes, a collection of mock JSON notices and skymaps can be read using the `test_notice`, `read_flatres_skymap`, and `read_multires_skymap` functions.

The scheduler is primarily composed of the functions `generate_targets` and `graph_targets`. Galaxy targets are generated using the GLADE catalog and a series of cuts to the catalog. The top 50 likely follow up candidates are saved to a CSV file. The skymap with these targets superimposed on the image are also saved.

## Installation

Apollo's scheduler uses the latest GLADE catalog (https://glade.elte.hu/) for selecting galaxy targets to follow up. Download the catalog as an ASCII text file (i.e. GLADE+.txt), and save it to the same working directory as the main Apollo module.

One requires an account on the GCN platform (https://gcn.nasa.gov/). The parameters `client_id` and `client_secret` are uniquely generated per user, and must be plugged into the configuration file in order to use the listener.

I have developed Apollo on the Ubuntu 22.04 LTS operating system using Python 3.10.6.

## References

- G. Dálya et al. (2022), GLADE+: An Extended Galaxy Catalogue for Multimessenger Searches with Advanced Gravitational-Wave Detectors, MNRAS, Volume 514, Issue 1, https://doi.org/10.1093/mnras/stac1443

- LVK Collaboration (2022), LIGO/Virgo/KAGRA Public Alerts User Guide, https://emfollow.docs.ligo.org/userguide/

---

19 Jul 2023<br>
Last update: 19 Jul 2023

Richard Camuccio<br>
rcamuccio@gmail.com