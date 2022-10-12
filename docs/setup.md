# Setup Parameters

## CCSN Code

Parameters in the `project/1dmlmix/setup`

??? Quote "setup"
    ```html
    --8<-- "project/1dmlmix/setup"
    ```

### Extra options
#### [Input PyTorch Model](#__codelineno-0-6)

* `None` - sets turbulent pressure to 0 and doesn't use any ML subroutines
* `ModelName.pt` - uses ML subroutines (model name can be anything)

#### [EOS options](#__codelineno-0-22)

* `1` - uses eospg - perfect gas equation of state
* `2` - uses eospgr - equation of state that includes gas and radiation pressure
* `3` - uses eos3 - Ocean eos assuming NSE and Swesty-Lattimer eos below `rhoswe`. Switches back to the internal energy variable of state.
* `4` - uses eos3 - Ocean eos assuming NSE and Swesty-Lattimer eos below `rhoswe`

## Data Preparation

Edit `prep_data/setup_prep` for initial parameters, then run `make data` from root directory.

??? Quote "setup_prep"
    ```html
    --8<-- "prep_data/setup_prep"
    ```

Parameters for near ~4k resolution based on Sukhbold 2016 data. Maximum radius is 1e10 cm, i.e. 1e10 in code units.

??? Abstract "maxrad=1e1"
    | Progenitor Mass | Initial Cell Mass | Result Grid Size |
    | --------------- | ----------------- | ---------------- |
    | 9               | 1.9e-5            | 3985             |
    | 10              | 2.5e-5            | 3960             |
    | 11              | 3.1e-5            | 3955             |
    | 12              | 4e-5              | 3975             |
    | 13              | 4.75e-5           | 3958             |
    | 14              |
    | 15              |                   |                  |
    | 16              |
    | 17              |                   |                  |
    | 18              |
    | 19              | 9.2e-5            | 3966             |
    | 20              | 1.05e-4           | 3955             |
    | 21              |
    | 22              |
    | 23              |
    | 24              |
    | 25              |

??? Abstract "variable maxrad"
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

??? Quote "setup_readout"
    ```html
    --8<-- "project/1dmlmix/setup_readout"
    ```