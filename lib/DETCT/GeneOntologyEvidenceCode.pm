## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::GeneOntologyEvidenceCode;
## VERSION
## use critic

# ABSTRACT: GO evidence codes

## Author         : is1
## Maintainer     : is1
## Created        : 2014-03-21
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
Readonly our %DESCRIPTION_FOR_EVIDENCE_CODE => (
    EXP => 'Inferred from Experiment',
    IDA => 'Inferred from Direct Assay',
    IPI => 'Inferred from Physical Interaction',
    IMP => 'Inferred from Mutant Phenotype',
    IGI => 'Inferred from Genetic Interaction',
    IEP => 'Inferred from Expression Pattern',
    ISS => 'Inferred from Sequence or Structural Similarity',
    ISO => 'Inferred from Sequence Orthology',
    ISA => 'Inferred from Sequence Alignment',
    ISM => 'Inferred from Sequence Model',
    IGC => 'Inferred from Genomic Context',
    IBA => 'Inferred from Biological aspect of Ancestor',
    IBD => 'Inferred from Biological aspect of Descendant',
    IKR => 'Inferred from Key Residues',
    IRD => 'Inferred from Rapid Divergence',
    RCA => 'inferred from Reviewed Computational Analysis',
    TAS => 'Traceable Author Statement',
    NAS => 'Non-traceable Author Statement',
    IC  => 'Inferred by Curator',
    ND  => 'No biological Data available',
    IEA => 'Inferred from Electronic Annotation',
    NR  => 'Not Recorded',
);

1;
