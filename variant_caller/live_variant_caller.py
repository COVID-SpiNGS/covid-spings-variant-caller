import functools
import operator
import pickle
import math
from typing import List

import pysam
import numpy as np

from tqdm import tqdm
from pysam import AlignedSegment
from time import strftime, localtime

import logging
import config_util.logging as log

from .structs import Site, Variant
from .utils import genotype_likelihood, from_phred_scale, to_phred_scale

import constants as c

class LiveVariantCaller:
    def __init__(self, referenceFasta: str, minBaseQuality: int, minMappingQuality: int, minTotalDepth: int,
                 minAlleleDepth: int, minEvidenceRatio: float, maxVariants: int):
        self.minBaseQuality = minBaseQuality
        self.minMappingQuality = minMappingQuality
        self.minTotalDepth = minTotalDepth
        self.minAlleleDepth = minAlleleDepth
        self.minEvidenceRatio = minEvidenceRatio
        self.maxVariants = maxVariants
        self.fastaFile = pysam.FastaFile(referenceFasta)
        self.memory = {}
        self.reset_memory()

    def __del__(self):
        self.fastaFile.close()

    def reset_memory(self):
        self.memory: dict[int, Site] = {}


    def process_bam(self, inputBam: str, referenceIndex=0):
        bamFile = pysam.AlignmentFile(inputBam, 'rb')
        pileupColumns = bamFile.pileup(
            min_mapping_quality=self.minMappingQuality,
            min_base_quality=self.minBaseQuality,
            reference=self.fastaFile.references[referenceIndex]
        )

        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        progressBar = tqdm(
            pileupColumns,
            desc=f'{timestamp} Processing {inputBam}',
            total=bamFile.get_reference_length(self.fastaFile.references[referenceIndex])
        )

        for pileupColumn in progressBar:
            self.process_pileup_column(pileupColumn)

        bamFile.close()

    def process_pileup_column(self, pileupColumn: AlignedSegment):
        totalDepth = len(pileupColumn.pileups)

        if pileupColumn.reference_pos not in self.memory:
            reference = self.fastaFile.fetch(reference=pileupColumn.reference_name)

            self.memory[pileupColumn.reference_pos] = {
                'reference': reference[pileupColumn.reference_pos],
                'totalDepth': totalDepth,
                'snvs': {},
                'indels': {}
            }
        else:
            self.memory[pileupColumn.reference_pos]['totalDepth'] += totalDepth

        for pileup in pileupColumn.pileups:
            self.process_pileup_at_position(pileupColumn.reference_pos, pileup)

    def process_pileup_at_position(self, position: int, pileup):
        self.process_svn(position, pileup)
        # self.process_indel(position, pileup)

    def process_svn(self, position, pileup):
        if not pileup.is_del and not pileup.is_refskip:
            svn = pileup.alignment.query_sequence[pileup.query_position]

            if svn not in self.memory[position]['snvs'].keys():
                self.memory[position]['snvs'][svn] = []

            self.memory[position]['snvs'][svn].append(pileup.alignment.query_qualities[pileup.query_position])

    def process_indel(self, position, pileup):
        if pileup.is_del or pileup.is_refskip:
            indel = '-' if pileup.is_del else f'+{pileup.alignment.query_sequence[pileup.query_position]}'

            if indel not in self.memory[position]['indels'].keys():
                # We could store more information here in the memory. But as the Base Qualty is the the only information that matters for further calculation we 
                # save memory and only store them

                self.memory[position]['indels'][indel] = []

            if pileup.is_refskip:
                self.memory[position]['indels'][indel].append(pileup.alignment.query_qualities[pileup.query_position])
            else:
                self.memory[position]['indels'][inde/].append(None)

    def prepare_variants(self):
        timestamp = strftime('[%Y-%m-%d %H:%M:%S]', localtime())
        progressBar = tqdm(
            self.memory,
            desc=f'{timestamp} Calculating statistics',
            total=len(self.memory.keys())
        )

        variants: List[Variant] = []

        for position in progressBar:
            if self.memory[position]['totalDepth'] >= self.minTotalDepth:
                snvs = {
                    allele: [
                        from_phred_scale(quality)
                        for quality in self.memory[position]['snvs'][allele]
                    ]
                    for allele in self.memory[position]['snvs'].keys()
                }

                genotypeLikelihoods = {
                    allele: genotype_likelihood(allele, snvs)
                    for allele in snvs.keys()
                }

                sumGenotypeLikelihoods = functools.reduce(operator.add, genotypeLikelihoods.values(), 0.0)
                sumGenotypeLikelihoods = sumGenotypeLikelihoods if sumGenotypeLikelihoods != 0 else 1.0

                for allele in snvs.keys():
                    alleleDepth = len(snvs[allele])

                    filterConstrains = [
                        self.memory[position]['reference'] != allele,
                        alleleDepth >= self.minAlleleDepth,
                        alleleDepth / self.memory[position]['totalDepth'] >= self.minEvidenceRatio
                    ]

                    if all(filterConstrains):
                        genotypeLikelihood = genotypeLikelihoods[allele]

                        if genotypeLikelihood != 0:
                            gl = math.log10(genotypeLikelihood)
                            pl = round(-10.0 * gl)
                        else:
                            gl = 0
                            pl = 0

                        score = to_phred_scale(1.0 - (genotypeLikelihoods[allele] / sumGenotypeLikelihoods))
                        qual = np.mean(snvs[allele])

                        variants.append({
                            c.VCF_START: position,
                            c.VCF_STOP: position + 1,
                            c.VCF_ALLELES: (
                                self.memory[position]['reference'],
                                allele
                            ),
                            c.VCF_QUAL: qual,
                            c.VCF_INFO: {
                                c.VCF_DP: self.memory[position]['totalDepth'],
                                c.VCF_AD: alleleDepth,
                                c.VCF_GL: gl,
                                c.VCF_PL: pl,
                                c.VCF_SCORE: score
                            }
                        })

                for indel in self.memory[position]['indels'].keys():
                    alleleDepth = len(self.memory[position]['indels'][indel])

                    filterConstrains = [
                        alleleDepth >= self.minAlleleDepth,
                        alleleDepth / self.memory[position]['totalDepth'] >= self.minEvidenceRatio
                    ]

                    if all(filterConstrains):
                        if indel == '-':
                            variants.append({
                                c.VCF_START: position,
                                c.VCF_STOP: position + 1,
                                c.VCF_ALLELES: (
                                    self.memory[position]['reference'],
                                    '*'
                                ),
                                c.VCF_QUAL: 0,
                                c.VCF_INFO: {
                                    c.VCF_DP: self.memory[position]['totalDepth'],
                                    c.VCF_AD: alleleDepth,
                                    c.VCF_GL: 0,
                                    c.VCF_PL: 0,
                                    c.VCF_SCORE: 0
                                }
                            })
                        else:
                            variants.append({
                                c.VCF_START: position,
                                c.VCF_STOP: position + 1,
                                c.VCF_ALLELES: (
                                    '*',
                                    indel[1:]
                                ),
                                c.VCF_QUAL: 0,
                                c.VCF_INFO: {
                                    c.VCF_DP: self.memory[position]['totalDepth'],
                                    c.VCF_ED: alleleDepth,
                                    c.VCF_GL: 0,
                                    c.VCF_PL: 0,
                                    c.VCF_SCORE: 0
                                }
                            })

        return variants


    def prev_variant(self, variants: List[Variant], variant: Variant):
        return next(
            (
                v for v in variants
                if v[c.VCF_START] == variant[c.VCF_START] - 1
            ),
            None
        )

    def next_variant(self, variants: List[Variant], variant: Variant):
        return next(
            (
                v for v in variants
                if v[c.VCF_START] == variant[c.VCF_START] + 1
            ),
            None
        )

    def concat_deletions(self, variants: List[Variant]):
        concatinatedVariants: List[Variant] = []
        currentVariant: Variant = None

        for variant in variants:
            if variant[c.VCF_ALLELES][1] == '*':
                nextVariant = self.next_variant(variants, variant)

                if nextVariant:
                    if not currentVariant:
                        currentVariant = variant
                    else:
                        currentVariant = {
                            c.VCF_START: currentVariant[c.VCF_START],
                            c.VCF_STOP: variant[c.VCF_STOP],
                            c.VCF_ALLELES: (
                                f'{currentVariant["alleles"][0]}{variant["alleles"][0]}',
                                '*'
                            ),
                            c.VCF_QUAL: variant[c.VCF_QUAL],  # must be combined
                            c.VCF_INFO: variant[c.VCF_INFO]  # must be combined

                        }
                else:
                    if currentVariant:
                        concatinatedVariants.append(currentVariant)
                        currentVariant = None


            else:
                concatinatedVariants.append(variant)

        return concatinatedVariants

    def concat_insertions(self, variants: List[Variant]):
        return variants
