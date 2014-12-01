mafs
----
Contains MAF files named [sample].snvs.maf and [sample].indels.maf.


bams
----
Contains folders of BAM files for each sample. Folder structure is
    bams/
        [sample]/
            tumour/
                a.bam
                a.bam.bai
                b.bam
                b.bam.bai
            normal/
                a.bam
                a.bam.bai
                b.bam
                b.bam.bai

by-sample
---------
After making, will contain one MAF file for each sample, with reference and
normal counts summed over each of the input bams.

by-patient
----------
Also contains one MAF file for each sample after making, but is augmented with
additional variants not found in the sample's original MAF file, but which were
found in another MAF file for the same patient.
