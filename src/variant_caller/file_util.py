import pysam
from . import vcf_file_constants as c
import pickle
import config_util.logging as log

def write_vcf(self, output_vcf: str):
        print("VCF output", output_vcf)
        vcf_header = pysam.VariantHeader()

        vcf_header.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, c.VCF_DP),
            (c.VCF_NUMBER, 1),
            (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
            (c.VCF_DESCRIPTION, c.VCF_TOTAL_DEPTH_STR)
        ])

        vcf_header.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, c.VCF_AD),
            (c.VCF_NUMBER, 1),
            (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
            (c.VCF_DESCRIPTION, c.VCF_ALLELE_DEPTH)
        ])

        vcf_header.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, c.VCF_GL),
            (c.VCF_NUMBER, 1),
            (c.VCF_TYPE, c.VCF_TYPE_FLOAT),
            (c.VCF_DESCRIPTION,
             'Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields')
        ])

        vcf_header.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, c.VCF_PL),
            (c.VCF_NUMBER, 1),
            (c.VCF_TYPE, c.VCF_TYPE_INTEGER),
            (c.VCF_DESCRIPTION,
             'The phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)')
        ])

        vcf_header.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, c.VCF_SCORE),
            (c.VCF_NUMBER, 1),
            (c.VCF_TYPE, c.VCF_TYPE_FLOAT),
            (c.VCF_DESCRIPTION, 'Custom scoring function')
        ])

        for reference in self.fasta_file.references:
            vcf_header.contigs.add(
                reference,
                self.fasta_file.get_reference_length(reference)
            )

        vcf_file = pysam.VariantFile(output_vcf, mode='w', header=vcf_header)

        variants = self.prepare_variants()
        # gvariants = self.concat_deletions(variants)

        for index, variant in enumerate(
                sorted(variants, key=lambda variant: (variant[c.VCF_START], variant[c.VCF_INFO][c.VCF_SCORE]))):
            vcf_file.write(
                vcf_file.new_record(
                    start=variant[c.VCF_START],
                    stop=variant[c.VCF_STOP],
                    alleles=variant[c.VCF_ALLELES],
                    qual=variant[c.VCF_QUAL],
                    info=variant[c.VCF_INFO],
                )
            )

        vcf_file.close()


def create_checkpoint(self, filename):
    log.print_and_log(f'Creating checkpoint {filename}', log.INFO)
    print('SELF.MEMORY', type(self.memory))
    file = open(filename, 'wb')
    pickle.dump(self.memory, file)
    file.close()

def load_checkpoint(self, filename):
    log.print_and_log(f'Loading checkpoint {filename}', log.INFO)

    file = open(filename, 'rb')
    self.memory = pickle.load(file)
    file.close()