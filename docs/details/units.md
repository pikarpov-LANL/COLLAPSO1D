# Code Units

Unit conversion from numerical to CGS can be found in the `#!fortran subroutine unit` in [project/1dmlmix/1dmlmix.f90](https://github.com/pikarpov-LANL/1dccsn/blob/master/project/1dmlmix/1dmlmix.f90). Here is a summary table of basic units:

| Quantity    |  Units   |
| ----------- | :------: |
| Mass        | 2e33 $g$ |
| Distance    | 1e9 $cm$ |
| Temperature | 1e9 $K$  |
| Time        | 1e1 $s$  |

## Conversion Constants

| Name     | Value     | Comment       |
| -------- | --------- | ------------- |
| ggcgs    | 6.67e-8   |
| avo      | 6.02e23   |
| aradcgs  | 7.565e-15 |
| boltzk   | 1.381e-16 |
| hbar     | 1.055e-27 |
| ccgs     | 3e10      |
| emssmev  | 0.511e0   |
| boltzmev | 8.617e-11 |
| ergmev   | 6.2422e5  |
| sigma1   | 9d-44     |
| sigma2   | 5.6d-45   |
| c2cgs    | 6.15e-4   |
| c3cgs    | 5.04e-10  |
| fermi    | 1d-13     |
| uergg    | 1d16      | ergs per gram |