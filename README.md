VNLB | Non-local Bayesian Video Denoising
=========================================

* Author    : Pablo Arias <pariasm@gmail.com>, see `AUTHORS`
* Copyright : (C) 2019, Pablo Arias <pariasm@gmail.com>

OVERVIEW
--------

This code provides an implementation of the video denoising method VNLB-H
described in:

[P. Arias, J.-M. Morel. "Video denoising via empirical Bayesian estimation of
space-time patches", Journal of Mathematical Imaging and Vision, 60(1),
January 2018.](https://link.springer.com/article/10.1007%2Fs10851-017-0742-4)

Please cite the publication if you use results obtained with this code in your
research.

The following libraries are also included as part of the code:
* For computing the optical flow, it includes [the IPOL
implementation](http://www.ipol.im/pub/art/2013/26/) of
the [TV-L1 optical flow method of Zack and Pock and
Bischof](https://link.springer.com/chapter/10.1007/978-3-540-74936-3_22).
* For image I/O, we use [Enric Meinhardt's iio](https://github.com/mnhrdt/iio).

**Dependencies:** The code depends on the following packages:
* [CBLAS](http://www.netlib.org/blas/#_cblas),
[LAPACKE](https://www.netlib.org/lapack/lapacke.html): operations with matrices
* OpenMP: parallelization [optional, but recommended]
* libpng, libtiff and libjpeg: image i/o


**Compilation:** 
Compilation was tested on Ubuntu Linux 16.04 and 18.04.
Configure and compile the source code using cmake and make.
It is recommended that you create a folder for building:
```
$ mkdir build; cd build
$ cmake ..
$ make
```

Binaries will be created in `build/bin folder`.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). Use the
`OMP_NUM_THREADS` enviroment variable to control the number of threads
used.

USAGE
-----

The compilation creates the following binaries:
* `awgn` add additive white Gaussian noise to an image sequence
* `psnr` compute PSNR (and RMSE) between two image sequences
* `tvl1flow` compute TVL1 optical flow between two frames
* `vnlbayes` run VNLB-H on an image sequence

In addition, the following helper scripts will be installed in `bin/`
* `tvl1flow-seq.sh` compute forward/backward optical flow for an image sequence
* `vnlb.sh` computes optical flow and run `vnlbayes`
* `vnlb-gt.sh` given a clean sequence: adds noise, computes optical flow, runs
denoising and computes PSNR

-----
**Usage of `vnlbayes`**

This is the main program. The simplest use is shown via the following example:

```
vnlbayes -i /my/video/frame-%03d.png -f 1 -l 10 -sigma 10
```

The method reads the video as a sequence of images. The sequence of images is passed
as a pattern in printf format, thus `frame-%03d.png` means that frames have the following
filenames: `frame-001.png`, `frame-002.png`, etc. The first and last frame are
given with '-f' and '-l' (1 and 10 in the example). The denoising method performs
two steps or iterations. The outputs will be stored by default
in the folder from where the program is invoked:
* `bsic_%03d.png`: first iteration (or basic estimate)
* `deno_%03d.png`: second (and final) iteration

If an optical flow (forward and backward) has been computed, it can be given to the
method as
```
vnlbayes -i /my/video/frame-%03d.png -f 1 -l 10 -sigma 10\
         -fof /fwd/flow/frame-%03d.flo \
         -bof /bwd/flow/frame-%03d.flo
```

Several options can be given. For example, to use `5x5x2` patches
in the first step, run:
```
vnlbayes -i /my/video/frame-%03d.png -f 1 -l 10 -sigma 10\
         -fof /fwd/flow/frame-%03d.flo \
         -bof /bwd/flow/frame-%03d.flo \
         -px1 5 -pt1 2
```
For a list of all parameters run `vnlbayes --help`, and for more information about
the default parameters, see [this section.](#default-parameters)


*NOTE: The optical flow files given to the method should be numbered according to the
following convention:*
* *Forward flow at frame *i*: optical flow from frame *i* to *i+1
* *Backward flow at frame *i*: optical flow from frame *i* to *i-1

*For the first frame, there is no backward flow, and for the last frame there is
no forward flow. **However, the code of vnlbayes expects forward and backward
flow files for each frame, even for the first and the last frames** (this is due
to the author being lazy). The backward flow for the first frame and the forward
flow for the last frame will not be used, so they can be anything (an image full
of zeros, for example). The script `tvl1flow-seq.sh` computes the optical flows for
a sequence following the naming convention used by VNLB.*

-----
**Usage of `vnlb.sh`**

This is just a helper scripts that computes the optical flow and then runs the
`vnlbayes` binary. Its usage is the following:

```
vnlb.sh /my/video/frame-%03d.png first-frame last-frame sigma out-folder ["vnlb-params"]
```

 `sigma` refers to the noise
standard deviation. Outputs are going to be left in `out-folder/`. The `"vnlb-params"` is an 
optional string to override any of the default parameters of the denoising. Note that the
parameters given in this string *need to be between quotes*. For example, to use `5x5x2` patches
in the first step, run:
```
vnlb.sh /my/video/frame-%03d.png first-frame last-frame sigma out-folder "-px1 5 -pt1 2"
```

For a list of all parameters and options run `vnlbayes --help`.

Outputs include:
* `out-folder/tvl1_%03d_f.flo` forward optical flows
* `out-folder/tvl1_%03d_b.flo` backward optical flows
* `out-folder/bsic_%03d.tif`: first iteration (or basic estimate)
* `out-folder/deno_%03d.tif`: second (and final) iteration

-----
**Usage of `vnlb-gt.sh`**

This helper script is useful for evaluating the performance of the denoising method by
comparing its output to the ground truth clean video. It assumes that the input sequence
has no noise, then adds noise of the given sigma, denoises it by calling the `vnlb.sh` script,
and then computes the RMSR/PSNR with respect to the clean video. Its usage is exactly the
same as `vnlb.sh`:

```
vnlb-gt.sh /my/video/frame-%03d.png first-frame last-frame sigma out-folder ["vnlb-params"]
```

`sigma` is the standard deviation of the noise that *will be added to the
sequence*. Outputs are going to be left in `out-folder/`. The `"vnlb-params"`
is an optional string to override any of the default parameters of the
denoising. Note that the parameters given in this string need to be between
quotes. For example, to use `5x5x2` patches in the first step, run:
```
vnlb-gt.sh /my/video/frame-%03d.png first-frame last-frame sigma out-folder "-px1 5 -pt1 2"
```

For a list of all parameters and options run `vnlbayes --help`.

The outputs are those of the `vnlb.sh` script (optical flows, denoised videos),
plus the following ones:
* `out-folder/%03d.tif`: frames with noise added (as tif floating point images)
* `out-folder/measures-bsic`: global and per-frame RMSE/PSNR of first iteration
* `out-folder/measures-deno`: same for final iteration

-----
**Usage of `awgn`**

Adds additive white Gaussian noise to an image sequence
```
awgn sigma input-path output-path [first-frame last-frame]
```

Paths for input/output frames have to be given in printf format: e.g.
```
awgn 10 /input/frame-%02d.png /output/frame-%02.tiff 1 20
```
It can also be applied to a single image:
```
awgn 10 /input/image.png /output/image.tiff
```

-----
**Usage of `psnr`**

Computes PSNR and RMSE between two image sequences

```
psnr input1-path input2-path [first-frame last-frame [output-file]]
```

Paths for input/output frames have to be given in printf format.
Prints to standard output the global PSNR and RMSE and, if an output filename
is given, it appends to the file the global and per-frame PSNR and RMSE.

-----
**Usage of `tvl1flow-seq.sh`**

```tvl1flow-seq.sh /my/video/frame-%03d.png first-frame last-frame out-folder [dir]```

Computes the optical flow in direction `dir` (either `fwd`, `bwd` or `both`),
and stores the output files in `out-folder/tvl1_%03d_f.flo` and
`out-folder/tvl1_%03d_b.flo` (forward and backward flows respectively).


DEFAULT PARAMETERS
------------------

The method has several paramaters (patch size, search window size, number of
similar patches, variance threshold, etc.).
The default parameters are set as a function of the patch size and the noise level sigma.
If no patch size is given, one is assigned by default (see below). This is
the one that maximized the PSNR on a small training set. Note that larger patches
take longer to run.

We have tuned the parameters for the following patch sizes:
* Grayscale:
     - 7x7x4
     - 8x8x3
     - **10x10x2** (default)
     - 14x14x1
     - 8x8x2
     - 5x5x4
     - 6x6x3
     - 7x7x2
     - 10x10x1
* RGB:
     - 5x5x4
     - **7x7x2** (default)
     - 10x10x1

If the patch size given by the user is not one of these, the program will choose
the default parameters from another patch size and a warning message will be printed.


FILES
-----

This project contains the following source files:
```
root/
├── cmake/                      used by cmake to find dependencies
│   └── modules/
│       ├── FindCBLAS.cmake
│       ├── FindLAPACKE.cmake
│       └── FindMKL.cmake
├── lib/
│   ├── iio/                    image i/o 3rd party lib
│   └── tvl1flow/               optical flow 3rd party lib
├── scripts/                    helper scripts
│   ├── CMakeLists.txt
│   ├── tvl1flow-seq.sh         compute optical flow on sequence
│   ├── vnlb-gt.sh              run VNLB with ground truth and compute PSNR
│   └── vnlb.sh                 run VNLB
└── src/
    ├── awgn.cpp                main function for awgn
    ├── CMakeLists.txt
    ├── cmd_option.h            library command line options
    ├── main_vnlb.cpp           main function for VNLB
    ├── psnr.cpp                main function for psnr
    ├── VidUtils/               small library with video tools
    │   ├── CMakeLists.txt
    │   ├── LibVideoT.hpp/cpp   template video class with some tools
    │   └── mt19937ar.h/c       random number generator
    └── VNLBayes/               VNLB code
        ├── CMakeLists.txt
        ├── LibMatrix.h/cpp     basic matrix computations
        └── VideoNLBayes.h/cpp  VNLB code
```


LICENSE
-------

Apart from a few exceptional files (see below), the the code of VNLB is
licensed under the GNU Affero General Public License v3.0, see `LICENSE`.

The following files are copied verbatim (or with trivial changes) from their
original sources. They are included for convenience and they are not essential
for VNLB. They are distributed under their own licences specified inside each file.

```
lib/iio/*                : written by Enric Meinhardt-Llopis
lib/tvl1flow/*           : written by Javier Sánchez
src/cmd_option.h         : written by David Tshumperlé
src/VidUtils/mt19937ar.* : written by Makoto Matsumoto and Takuji Nishimura
```

