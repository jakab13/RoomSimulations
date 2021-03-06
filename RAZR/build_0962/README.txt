===========================================
RAZR SOFTWARE PACKAGE FOR MATHWORK'S MATLAB
===========================================

RAZR - Room acoustics simulator

Version 0.96

Copyright (c) 2014-2020, Torben Wendt, Steven van de Par, Stephan Ewert,
University of Oldenburg, Germany.
All rights reserved. 

Authors: Torben Wendt, Stephan Ewert, Christoph Kirsch, Henning Steffens, 
Josef Poppitz, Oliver Buttler

Contact: Stephan Ewert
e-mail: stephan.ewert(at)uni-oldenburg.de

www.razrengine.com


===============================================================================
Algorithm reference
===============================================================================

Wendt, T., van de Par, S., Ewert, S.D. (2014): A Computationally-Efficient and
Perceptually-Plausible Algorithm for Binaural Room Impulse Response Simulation.
J. Audio Eng. Soc., Vol. 62, No. 11, pp. 748-766.

Please note: The default options contain some optimizations compared to the 
setup described in the publication above. You can however load the original
setup by using GET_DEFAULT_OPTIONS_2014.M. To see how to use options, see also
the subsection "Getting started/Data structures", below in this file.


===============================================================================
What's in the distribution?
===============================================================================

This distribution adds a toolbox for synthesizing and analyzing binaural
room impulse responses (BRIRs) for headphone auralization and multi-channel
room impulse responses (MRIRs) for loudspeaker-array rendering to your 
Matlab installation. Example simulations are included.


Included Files
--------------

./*                     - RAZR top-level functions, README, license, disclaimer
ANALYSIS_TOOLS/
    *.M                 - Collection of tools to analyze BRIRs
    HEADPHONE_EQ/*.mat  - Equalization files for headphones
EXAMPLES/*.M            - Some examples to create BRIRs
BASE/**                 - RAZR core files
    ABSCOEFF/*          - Collection of absorption coefficients and access tools
    COUPLED/*           - Functions for coupled-rooms simulation
    EXTERNAL/**         - Files created by other Authors, see Credits below
    HRTF/**             - Supplementary files for HRTF application
    MEXFDN/*            - Compiled FDN mainloop
    TOOLS/*.M           - Supplementary functions for simulation and analysis

All these files must not be changed unless not stated otherwise (see also 
RAZR's license). Changes might moreover harm proper function of this 
distribution. For exceptions see the following section "License and permissions".


===============================================================================
License and permissions
===============================================================================

Unless otherwise stated, the RAZR distribution, including all files is licensed
under Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International
(CC BY-NC-ND 4.0).
In short, this means that you are free to use and share (copy, distribute and
transmit) the RAZR distribution under the following conditions:

Attribution - You must attribute the RAZR distribution by acknowledgement of
              the author if it appears or if was used in any form in your work.
              The attribution must not in any way that suggests that the author
              endorse you or your use of the work.

Noncommercial - You may not use RAZR for commercial purposes.
 
No Derivative Works - You may not alter, transform, or build upon RAZR.


Exceptions are:
- The files ./BASE/EXTERNAL/**, see their respective license.

Users are explicitly granted to modify, derive their own modifications and to
distribute their modifications for the configurations, options and examples:
- The files RAZR_CFG_USER.M and SELECT_RAZR_CFG.M in the main folder
- All files in the folder ./EXAMPLES

For further exceptions, see the license headers in the respecive files and the
section "Credits" in this file.


===============================================================================
Installation
===============================================================================

Requirements
------------
- Windows 95/98, Windows NT/2000/XP or higher, or MAC OSX, or Linux
- Tested under Mathwork's MATLAB 2019b.


Recommended
-----------
- Sound device
- Headphones


Installation
------------
- Unzip all included files to your preferred directory.


===============================================================================
Getting Started
===============================================================================

Demo and Example scripts
------------------------
To start the demo, type "razr demo".

Besides, in the directory "EXAMPLES", you can find some example scripts, that 
demonstrate the usage of this program package. See the respective files for
details. After having read and tried out all examples, you are familiar with the 
most important features of RAZR.
To start any of the examples, type "razr addpath" to temporally add all required
folders to MATLAB's path.
Then type, e.g., "example_default" or "example_options". 

A formal usage description is provided in the documentation of the top-level
function razr(). Type "razr help" to display it.

Please see also the changelog if you have already used earlier versions of RAZR.


Data structures
---------------
RAZR uses three essential data structures, that occur in the probably most
common usage of RAZR:

ir = razr(room, op);

(1) A room is defined as a Matlab structure, called "room" in most scripts and
functions. It contains information about the room size, wall absorption,
positions of source and receiver, and other specifications. For the required
fields of a room structure see the example rooms defined in
EXAMPLES/GET_ROOM_*.M. Both, required and optional fields are explained in the
help entry of RAZR.M. Depending on your taste and needs, you can
generate and store own rooms either in functions like the examples, as *.mat
files, directly within scripts, or by any other means. For example, for
automated RIR creation for many rooms in a loop, .mat files might be the most
convenient way to store rooms as they can easily be loaded by filename.

(2) Options for the BRIR synthesis are defined in a structure, called "op" in
most scripts and functions. All possible fields and their default values are
stored in BASE/GET_DEFAULT_OPTIONS.M. To display options, check out also the
function "sop".
Some of them are experimental, refer to features in development (and thus not
contained in the current version), or might read unclear. In doubt, stick to
the default or contact the author. See EXAMPLES/EXAMPLE_OPTIONS.M to see how
to use options.

(3) Room impulse responses are stored as a structure, called "ir" in most
scripts and functions. It contains the time signal itself and some metadata,
depending on your chosen options. For instance, ir can contain also the signals
for direct sound, early and late reflections in separate fields, if the option
"return_rir_parts" is set to "true".
Most of the functions in ANALYSIS_TOOLS expect a BRIR in the format of such a
structure as input parameter.

For further information, type "razr help".


Configuration
-------------
In contrast to the options structure, which specifies how a BRIR is created, in
the configuration it is mainly defined, where RAZR finds certain files on your
computer. Hence, on different computers, the same options can be used, whereas
the configuration has to be adapted.
Configuration structures are defined in RAZR_CFG_DEFAULT.M and RAZR_CFG_USER.M.
(The latter one is automatically created when you run razr for the first time.)
Internally, RAZR_CFG_USER overloads RAZR_CFG_DEFAULT. You can also create and
store several configuration files in parallel. In the file SELECT_RAZR_CFG.M,
you can specify what configuration is actually used to override
RAZR_CFG_DEFAULT. RAZR_CFG_USER is the predefined one.


===============================================================================
Using HRTFs/HRIRs
===============================================================================

This program package does not contain any HRTF data. In the default settings, a
spherical head model (SHM) is used, implemented after [2]. However, although the
SHM yields reasonably good results (and is computationally very efficient),
using HRTFs is an essential part of BRIR simulation. If you would like to use
HRTFs you need to obtain a database from external sources. 
To link an HRTF database to RAZR, download the database to the preferred 
directory on your harddrive. Then, this database must be specified in 
./RAZR_CFG.M as described in the following:

SOFA format
-----------
RAZR supports HRIRs that are stored in the SOFA format [3], using the SOFA con-
vention "SimpleFreeFieldHRIR" [3]. If you would like to use SOFA-formatted HRIRs
you need to download the SOFA API for Matlab/Octave from

    https://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs

Store and unzip the API in your preferred directory. In RAZR, open the file
RAZR_CFG_USER and set the field "cfg.sofa_api_path" to the path where you stored
the API.
Then, there are two ways to use a SOFA database in RAZR:
(1) Specify the name of the .sofa file directly as option, e.g.,
    op.hrtf_database = path/to/database.sofa;
(2) Use a shortcut. In ./RAZR_CFG.M, a shortcut to a .sofa file can be linked
    as, e.g,.:
    cfg.sofa_file__<shortcut> = path/to/database.sofa;
    (Note the double underscore before <shortcut>.) Then, the option to use this
    Database can be set as:
    op.hrtf_database = 'shortcut.sofa';
    (The suffix ".sofa" is required although it is not part of the shortcut.)

Other formats
-------------
For an HRTF database that is stored in any other format, its path on your hard-
drive must be specified in ./RAZR_CFG.M:

cfg.hrtf_path__<shortcut> = path/to/database;

where <shortcut> is the unique key for the database. For already supported data-
bases, configuration fields have already been prepared there. (Note: Some of
them are not open access databases or not published, yet.)

If you would like to use an unsupported (non-SOFA) HRTF database, you need
create two Matlab functions, which will be executed internally in RAZR:
  - hrtf_params_<shortcut>  to set database-specific parameters.
  - pick_hrir_<shortcut>    to pick a set of HRIRs from the database according
                            to a set of desired azimuth and elevation angles.
For the syntax of these functions and required parameters please see the
examples provided in ./BASE/HRTF.

If you create these functions for an open-access database, be cordially invited
to share them with the author, such that they can be included in a future
release.

All formats (SOFA and non-SOFA)
-------------------------------
To actually use an HRTF database (SOFA or non-SOFA), the following option must
be set and passed to razr:
op.spat_mode = 'hrtf';

Note that the sampling rate of the HRIRs from a database must match the sampling
rate of the BRIR generated by RAZR (44100 Hz). If necessary, it can be adjusted
using the option op.fs, e.g.,
op.fs = 48000;

HRTF usage is also demonstrated in EXAMPLES/EXAMPLE_HRTF.M which is precon-
figured for the FABIAN database [4]. To run this example with the FABIAN data-
base, download the database from [3] and store it in a preferred directory on
your harddrive. Then replace 'path_to_database' in ./RAZR_CFG.M in the SOFA file
definition for the FABIAN database with the path to the directory where you
stored the database.


===============================================================================
References
===============================================================================

[1] Wendt, T., van de Par, S., Ewert, S.D. (2014): A Computationally-Efficient
    and Perceptually-Plausible Algorithm for Binaural Room Impulse Response
    Simulation. J. Audio Eng. Soc., Vol. 62, No. 11
[2] C. P. Brown and R.O. Duda (1998): "A Structural Model for Binaural  Sound
    Synthesis," IEEE Trans. Speech Audio Processing, vol. 6, no. 5
[3] https://www.sofaconventions.org
[4] Brinkmann, F., Lindau, A., Weinzierl, S., Geissler, G., & van de Par, S.
    (2013). A high resolution head-related transfer function database including
    different orientations of head above the torso. AIA-DAGA 2013.

===============================================================================
Credits
===============================================================================

This work was supported by
- DFG-FOR 1732 "Individualisierte Hoerakustik",
- DFG SFB 1330 352015383 - C5
- Cluster of Excellence EXC 1077/1 Hearing4all.

We'd like to thank Torben Wendt for the implementation of the original version of
this software and Stephan D. Ewert and Steven van de Par for continuous support 
on the algorithm development and evaluation.

We'd like to thank all the people who tested earlier versions of the Software,
contributed to it, reported on bugs, or made suggestions for improvements:
Thomas Biberger, Oliver Buttler, Stephan D. Ewert, Jan-Hendrik Flessner, Nico
Goessling, Julian Grosse, Andreas Haeussler, Stefan Klockgether, Steven van de
Par, Josef Poppitz, Joachim Thiemann, Henning Steffens, Christoph Kirsch 
(Sorry if we forgot someone here.)

All files in the directory ./BASE/EXTERNAL and all subdirectories are written
by other authors. Please see the file headers or license.txt files further
informations on licenses. In some cases I did code modifications. For details,
see the comments in the code. The authors are:

Daniyal Amir        OVERWRITE_MERGE.M
Robert Baumgartner  HOR2GEO.M
Piotr Majdak        GEO2HORPOLAR.M
Ofek Shilon         RANDORTHMAT.M (modified Version for RAZR: RANDORTHMAT3.M)
Todd Welti          SMOOTHNEW3.M
Archontis Politis   Vector Base Amplitude Panning Library

Furthermore, some supplementary functions for HRTF databases, located in
respective subfolders of ./BASE/HRTF, are written by the following authors:

Hendrik Kayser (kayser-database), Gunnar Geissler (mk2-database).

The unreverberated sound samples in ./BASE/EXTERNAL/SAMPLES are
part of the Oldenburg Sentence Test for AFC Software Package by Stephan Ewert,
Daniel Berg, Hoertech gGmbH. See also
BASE/EXTERNAL/SAMPLES/OLSA_SENTENCES_README.TXT.
