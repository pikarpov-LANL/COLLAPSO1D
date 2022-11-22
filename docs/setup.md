# Setup Parameters

## CCSN Code

Parameters in the `project/1dmlmix/setup`

???+ Quote "setup"
    ```html
    --8<-- "project/1dmlmix/setup"
    ```

### Extra options
#### [Input PyTorch Model](#__codelineno-0-6)
|Name|Description|
|:---|:---|
|`None` | sets turbulent pressure to 0 and doesn't use any ML subroutines
| `ModelName.pt` | uses ML subroutines (model name can be anything)

#### [EOS options](#__codelineno-0-22)

|EOS|Subroutine|Description|
|:---|:---|:---|
| `1` | eospg | perfect gas equation of state
| `2` | eospgr | equation of state that includes gas and radiation pressure
| `3` | eos3 | Ocean eos assuming NSE for low density and Swesty-Lattimer eos for high density above `rhoswe`. Uses *energy* as the primary variable.
| `4` | eos3 | Ocean eos assuming NSE for low density and Swesty-Lattimer eos for high density above `rhoswe`. Uses *entropy* as the primary variable.
| `5` | eos5 | SFHo EOS Tables with *entropy* is primary variable.
| `6` | eos5 | SFHo EOS Tables with *energy* is primary variable.

## Data Preparation

Edit `prep_data/setup_prep` for initial parameters, then run `make data` from root directory.

???+ Quote "setup_prep"
    ```html
    --8<-- "prep_data/setup_prep"
    ```

After setting the target grid size and the convection region resolution, the last thing we need to actually vary is the enclosed mass of the convection region, since this is a Lagrangian code. We wanted target a convective region resolution of ~1000 to cover the region between 20 and 200 km, so here are the ecnlosed masses:

??? Abstract "conv_region=1000"
    | Progenitor Mass | Enclosed Mass |
    | --------------- | ------------- |
    | 9               | 1.            |
    | 10              | 1.            |
    | 11              | 1.31          |
    | 12              | 1.31          |
    | 13              | 1.33          |
    | 14              | 1.35          |
    | 15              | 1.35          |
    | 16              | 1.37          |
    | 17              | 1.37          |
    | 18              | 1.35          |
    | 19              | 1.37          |
    | 20              | 1.39          |
    | 21              | 1.            |
    | 22              | 1.            |
    | 23              | 1.            |
    | 24              | 1.            |
    | 25              | 1.            |

??? Abstract "conv_region=2000, pns=800, grid=4000"
    | Progenitor Mass | Enclosed Mass |
    | --------------- | ------------- |
    | 11              | 1.4           |
    | 12              | 1.42          |
    | 13              | 1.46          |
    | 14              | 1.46          |
    | 15              | 1.46          |
    | 16              | 1.46          |
    | 17              | 1.47          |
    | 18              | 1.48          |
    | 19              | 1.49          |
    | 20              | 1.5           |

pns_cutoff = Enclosed Mass-0.16 for the 4k grid above. Maxrad was 5e9 cm for those runs.

Note: all of the above checkout ~200km radius ~10ms after bounce

### Outdated
Parameters for near ~4k resolution based on Sukhbold 2016 data. Maximum radius is 1e10 cm, i.e. 1e10 in code units.

??? Abstract "OUTDATED: maxrad=1e1"
    | Progenitor Mass | Initial Cell Mass | Result Grid Size |
    | --------------- | ----------------- | ---------------- |
    | 9               | 1.9e-5            | 3985             |
    | 10              | 2.5e-5            | 3960             |
    | 11              | 3.1e-5            | 3955             |
    | 12              | 4e-5              | 3975             |
    | 13              | 4.75e-5           | 3958             |
    | 14              | 5.7e-5            | 3982             |
    | 15              | 6.9e-5            | 3980             |
    | 16              | 7.5e-5            | 3962             |
    | 17              | 8e-5              | 3982             |
    | 18              | 9.3e-5            | 3958             |
    | 19              | 9.2e-5            | 3966             |
    | 20              | 1.05e-4           | 3955             |
    | 21              | 1.15e-4           | 4000             |
    | 22              | 1.2e-4            | 3955             |
    | 23              | 1.31e-4           | 3984             |
    | 24              | 1.43e-4           | 3991             |
    | 25              | 1.6e-4            | 3943             |

??? Abstract "OUTDATED: variable maxrad"
    | Progenitor Mass | Initial Cell Mass | maxrad | Result Grid Size |
    | --------------- | ----------------- | ------ | ---------------- |
    | 9               | 1.9e-5            | 1e1    | 3985             |
    | 10              | 0.000063          | 3e4    | 3989             |
    | 12              | 0.000085          | 6.5e4  | 3991             |
    | 13              | 0.00009           | 6.5e4  | 3965             |
    | 15              | 0.00012           | 6.5e4  | 3976             |
    | 17              | 0.00014           | 6.5e4  | 3921             |
    | 19              | 0.00015           | 6.5e4  | 3956             |


## Binary to Readable

Parameters in the `project/1dmlmix/setup_read`

???+ Quote "setup_readout"
    ```html
    --8<-- "project/1dmlmix/setup_readout"
    ```