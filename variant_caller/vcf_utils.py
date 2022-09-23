
from typing import List

import pysam
import functools
import operator
import numpy as np


from config import minTotalDepth, minCandidatesDepth, minMappingQuality, minBaseQuality
from structs import Position
from utils import genotype_likelihood, error_probability, to_error_probability, to_phred_scale

def create_vcf_header() -> pysam.VariantHeader:
    vcfHeader = pysam.VariantHeader()

    vcfHeader.add_meta('INFO', items=[
        ('ID', 'DP'),
        ('Number', 1),
        ('Type', 'Integer'),
        ('Description', 'Raw Depth')
    ])

    vcfHeader.add_meta('INFO', items=[
        ('ID', 'AD'),
        ('Number', 1),
        ('Type', 'Integer'),
        ('Description', 'Allele Depth')
    ])

    vcfHeader.add_meta('INFO', items=[
        ('ID', 'MBQ'),
        ('Number', 4),
        ('Type', 'Float'),
        ('Description', 'Mean Allele Base Qualities (A, T, C, G)')
    ])

    vcfHeader.add_meta('INFO', items=[
        ('ID', 'MMQ'),
        ('Number', 4),
        ('Type', 'Float'),
        ('Description', 'Mean Allele Mapping Qualities (A, T, C, G)')
    ])

    vcfHeader.add_meta('INFO', items=[
        ('ID', 'AF'),
        ('Number', 4),
        ('Type', 'Integer'),
        ('Description', 'Allele Frequency (A, T, C, G)')
    ])

    return vcfHeader
