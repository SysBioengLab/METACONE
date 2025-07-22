# METACONE
Identification of interaction networks based on metabolic conversions.
All the functions and scripts need to be run with the [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html).
All the functons and scripts were created and tested with Matlab, vresion R2020a, Update 8 (9.8.0.1873465).

## Contents
* `example1_metaCone.m`: Simple script to demonstrate the use of the `metaCone()` algorithm.
* `example2_QModel.m`:  Example script to use the usage of the `QModelCreation()` algorithm.
* `models` folder: contains the models used to obtain the results.
* `functions` folder: that hold the main algorithm developed.

### Models

1. `toymodel1.bio`: Toy model with 15 internal reactions (including transport reactions) and 8 exchange reaction, with a mock-biomass reaction.
2. `bmodel`       : E. coli "core" model with the `EX_biomass(e)` reaction added and its respective metabolite.
3. `LS`           : Slightly modified `Clostridium_clostridioforme_CM201` model found in VMH.
4. `PD`           : Slightly modified `Bacteroides_dorei_DSM_17855` model found in VMH.

### Functions

1. `metaCone`: Runs the METACONE algorithm.
2. `QMatrixCreation`: Builds the stoichiometric matrix describing a consortium and its community exchanges.
3. `QModelCreation`:  Builds a model structure for the entire consortium.
