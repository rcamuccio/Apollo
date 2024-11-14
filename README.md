# Apollo

**APOLLO** is an astronomical alert and broker system. APOLLO consists of several Python functions that compose both a listener and scheduler. The listener employs the GCN Kafka Client for Python to receive classic and JSON notices from the NASA General Coordinates Network (GCN).

---

## Overview

### Code Structure

The code is run from the main APOLLO directory through a set of scripts that pull functions from module libraries. The scripts currently available include:

- GCN listener (`gcn_listener.py`)
- Survey simulator (`survey_simulator.py`)

### Running APOLLO

To run the GCN listener, navigate to the main APOLLO directory, and run the following command:

```
$ python -m scripts.gcn_listener
```

## Installation

### Setup APOLLO

First, download a copy of APOLLO and navigate to the main directory. Open the configuration file (`config.py`) and change the following parameters to suit your specific environment:

- `APOLLO_DIRECTORY`
- `DELPHI_DIRECTORY`

Then, in a terminal and from the main APOLLO directory, run the following command:

```
$ python -m scripts.setup_apollo
```

This will create the `Delphi` directory at the specified location in the configuration file. The `Delphi` directory is the place for saving alert contents, analysis files, catalogs, and logs.

### Select Observatory

There are currently a few observatories supported in the code. One may freely add their own observatories to the code (in the `observatory.py` library). To select an observatory, choose one that is successfully encoded by setting the `OBSERVATORY` parameter in the configuration file.

### Download Catalogs

Apollo's scheduler uses the latest GLADE catalog (https://glade.elte.hu/) for selecting galaxy targets to follow up. Download the catalog as an ASCII text file (i.e. GLADE+.txt), and save it to the directory `~/Delphi/analysis/survey`.

### Create GCN Account

One requires an account on the GCN platform (https://gcn.nasa.gov/). The parameters `CLIENT_ID` and `CLIENT_SECRET` are uniquely generated per user, and must be plugged into the configuration file in order to use the listener.

### Create Email Account

The GCN listener sends alerts to a mailing list from a specified email client. Set the email that will send the alerts via the parameter `EMAIL` in the configuration file. You can add email addresses to the `MAILING_LIST` parameter in order to send alerts to willing participants.

## Development

### Platform

I have developed APOLLO on the Ubuntu 24.04.1 LTS operating system using a Conda environment and Python 3.12.2.

### Testing

I have not written any tests of APOLLO yet. Suggestions are always welcome!

## References

- G. Dálya et al. (2022), GLADE+: An Extended Galaxy Catalogue for Multimessenger Searches with Advanced Gravitational-Wave Detectors, MNRAS, Volume 514, Issue 1, https://doi.org/10.1093/mnras/stac1443

- LVK Collaboration (2022), LIGO/Virgo/KAGRA Public Alerts User Guide, https://emfollow.docs.ligo.org/userguide/

- NASA, General Coordinates Network (GCN): NASA's Time-Domain and Multimessenger Alert System, https://gcn.nasa.gov/

---

19 Jul 2023<br>
Last update: 14 Nov 2024

Richard Camuccio<br>
rcamuccio@gmail.com