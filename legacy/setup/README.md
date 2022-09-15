## Compile

```bash
gfortran read_woosley.f -o a.out
```

## Run

```bash
./a.out
```

It will prompt you for the *innermost zone mass* and then evolves the masses based on prescriptions that worked well a long time ago. Zone mass should be <1. Lastly, it will ask for the final *binary file name* that will be used as an input for the simulation. The name of the input and output files should not contain <ins>underscores</ins>.

## Defaults

Input file, i.e. progenitor mass, should be edited directly in the `read_woosley.f`. There are 3 available from Heger & Woosley, 2000:

1. s15presn 
2. s20presn (Default)
3. s25presn

For the prompts:

| Prompt              | Value  |
| ------------------- | ------ |
| Innermost zone mass | 0.0003 |
| Binary file name    | Data   |

This will result in `grid_szie=1703`

