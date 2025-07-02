# METACONE
Identification of interaction networks based on metabolic conversions.
All the functions and scripts need to be run with the [COBRA Tool Box](https://opencobra.github.io/cobratoolbox/stable/index.html).
All the functons and scripts were created and tested with Matlab, vresion R2020a, Update 8 (9.8.0.1873465).

## Contents
* `example1_metaCone.m`: Simple script to demonstrate the use of the `metaCone()` algorithm.
* `example2_QModel.m`:  Example script to use the usage of the `QModelCreation()` algorithm.
* `models` folder: contains the models used to obtain the results.
* `functions` folder: that hold the main algorithm developed.

### Models

1. `toymodel1.bio`                  : Toy model with 15 internal reactions (including transport reactions) and 8 exchange reaction, with a mock-biomass reaction.
2. `bmodel`                         : E. coli "core" model with the `EX_biomass(e)` reaction added and its respective metabolite.
3. `Clostridium_symbiosum_WAL_14673`: Slightly modified `Clostridium_clostridioforme_CM201` model found in VMH.
4. `Bacteroides_dorei_DSM_17855`    : Slightly modified `Bacteroides_dorei_DSM_17855` model found in VMH.
5. `iML1515`                        : E. coli model, downloaded from BiGG ('http://bigg.ucsd.edu/models/iML1515')

### Functions

1. `metaCone`       : greedy, LP-based algorithm for computing a linearly independent basis of the Conversion Cone of a COBRA metabolic model, capturing diverse substrate-product conversions under growth constraints.
2. `QMatrixCreation`: function that generates a community-level stoichiometric matrix by aligning and merging multiple conversion bases (C_ext) from `metaCone()` across species, normalizing by biomass flux and ordering exchange metabolites consistently.
3. `QModelCreation` : function that builds a community-level COBRA model using METACONE-generated conversion bases, producing a stoichiometric matrix internally that represents metabolite exchanges and species interactions.

### Usage

All presented functions can be used independently, but a few remarks are necessary. 
The `metaCone()` function will take one model as input, with certain parameters, and the main output will be a matrix (C_ext).
The usage of `metaCone()` is explained in the `example1_metaCone.m`.
The `QMatrixCreation()` is used internally and automatically by `QModelCreation()`, as well as `metaCone()`, but previously
calculated matrices, one for each model to be included, can be provided as argument.
The usage of `QModelCreation()` is explained in the `example2_QModel.m`.