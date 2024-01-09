def write_vcf(self, outputVfc: str):
        print("VFC output", outputVfc)
        vcfHeader = pysam.VariantHeader()

        vcfHeader.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, 'DP'),
            ('Number', 1),
            ('Type', 'Integer'),
            (c.VCF_DESCRIPTION, 'Total Depth')
        ])

        vcfHeader.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, 'AD'),
            ('Number', 1),
            ('Type', 'Integer'),
            (c.VCF_DESCRIPTION, 'Allele Depth')
        ])

        vcfHeader.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, 'GL'),
            ('Number', 1),
            ('Type', 'Float'),
            (c.VCF_DESCRIPTION,
             'Genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields')
        ])

        vcfHeader.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, 'PL'),
            ('Number', 1),
            ('Type', 'Integer'),
            (c.VCF_DESCRIPTION,
             'The phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)')
        ])

        vcfHeader.add_meta(c.VCF_INFO, items=[
            (c.VCF_ID, 'SCORE'),
            ('Number', 1),
            ('Type', 'Float'),
            (c.VCF_DESCRIPTION, 'Custom scoring function')
        ])

        for reference in self.fastaFile.references:
            vcfHeader.contigs.add(
                reference,
                self.fastaFile.get_reference_length(reference)
            )

        vcfFile = pysam.VariantFile(outputVfc, mode='w', header=vcfHeader)

        variants = self.prepare_variants()
        # gvariants = self.concat_deletions(variants)

        for index, variant in enumerate(
                sorted(variants, key=lambda variant: (variant[c.VCF_START], variant[c.VCF_INFO]['SCORE']))):
            vcfFile.write(
                vcfFile.new_record(
                    start=variant[c.VCF_START],
                    stop=variant[c.VCF_STOP],
                    alleles=variant[c.VCF_ALLELES],
                    qual=variant[c.VCF_QUAL],
                    info=variant[c.VCF_INFO],
                )
            )

        vcfFile.close()


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