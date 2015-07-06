### Author
[Sean Holden](http://www.cl.cam.ac.uk/~sbh11/)  
Computer Laboratory  
University of Cambridge

### Quick guide

- You'll need to go to the `c` folder and do an `R CMD SHLIB
outer.c`. The C code is only needed to speed up computing of the
kernel matrices. At present this only implements RBFs. `kernels-c.R`
encapsulates this.

- If you do `source("demo.R")` then `demo4P()` you should get a simple
example of a problem with 4 classes learned using probabilistic
outputs.  `easyData.R` is a bunch of stuff for making trivial data
sets for testing.

- `cv.R` is for stratified cross-validation and is needed because it
forms part of the method for turning outputs into probabilities.

- `multiShared.R` is mostly stuff for transforming multi-class
problems into sets of binary problems.

- You then have a whole bunch of stuff named `lpSVM*-c.R` If the name
includes `Multi` it's for multiple classes, otherwise it's binary. `P`
denotes probabilistic outputs and `2K` indicates it's the 2-kernel
version I developed to combine LOPIT and GO.

- In the non-`2K` code data is supplied as a `list(X, y)` with one
example per row in `X`. `y` is +/-1 for binary problems and positive
integers (preferably `1:n`) for `n` classes.

- In the `2K` versions there is in addition `X2`, which is the
features for the auxilliary data-only the primary data should be in
`X`. Also in the routines for classifying new stuff you have `x` for
the primary data and `x2` for the auxilliary data. Thus `c(x, x2)` is
the actual new vector you're classifying.

### LICENSE

This code is licensed under GNU General Public License, version 2
(GPL-2.0). A copy of the [license](./LICENSE) if distributed with the
code.

```
    Linear Programming SVM Transfer Learning
    https://github.com/ComputationalProteomicsUnit/lpsvm-tl-code
    Copyright (C) 2015  Sean Holden

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
```
