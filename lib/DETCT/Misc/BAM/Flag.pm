## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::BAM::Flag;
## VERSION
## use critic

# ABSTRACT: BAM flag field bits

## Author         : is1
## Maintainer     : is1
## Created        : 2014-04-03
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Readonly;

# Constants
Readonly our $READ_PAIRED         => 1;
Readonly our $PROPER_PAIR         => 2;
Readonly our $READ_UNMAPPED       => 4;
Readonly our $MATE_UNMAPPED       => 8;
Readonly our $READ_REVERSE_STRAND => 16;
Readonly our $MATE_REVERSE_STRAND => 32;
Readonly our $FIRST_IN_PAIR       => 64;
Readonly our $SECOND_IN_PAIR      => 128;
Readonly our $SECONDARY           => 256;
Readonly our $QC_FAILED           => 512;
Readonly our $DUPLICATE           => 1024;
Readonly our $SUPPLEMENTARY       => 2048;

1;
