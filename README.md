# Analysis of SARS-CoV-2 reads in sequencing of 2018-2019 Antarctica samples in PRJNA692319

The samples analyzed here are described in [this preprint](https://assets.researchsquare.com/files/rs-1177047/v1_covered.pdf), which is a pre-print by Istvan Csabai and co-workers that describes SARS-CoV-2 reads in samples from Antarctica sequencing in China.
I was originally alerted to the pre-print by Carl Zimmer on Dec-23-2021.
Istvan Csabai and coworkers subsequently posted a [second pre-print](https://www.researchsquare.com/article/rs-1330800/v1) that also analyzes the host reads. 

## Repeating key parts of the analysis
The code in this repo independently repeats some of the analyses.

To run the analysis, build the `conda` environment in [environment.yml](environment.yml) and then run the analysis using [Snakefile](Snakefile).
To do this on the Hutch cluster, using [run.bash](run.bash):

    sbatch -c 16 run.bash

The results are placed in the [./results/](results) subdirectory.
Most of the results files are not tracked due to file-size limitations, but the following key files are tracked:

 - [results/alignment_counts.csv](results/alignment_counts.csv) gives the number of reads aligning to SARS-CoV-2 for each sample. This confirms that three accessions (`SRR13441704`, `SRR13441705`, and `SRR13441708`) have most of the SARS-CoV-2 reads, although a few other samples also have some.
 - [results/variant_analysis.csv](results/variant_analysis.csv) reports all variants found in the samples relative to Wuhan-Hu-1.
 - [results/variant_analysis_to_outgroup.csv](results/variant_analysis_to_outgroup.csv) reports the variants found in the samples that represent mutations from Wuhan-Hu-1 **towards** the two closest bat coronavirus relatives, RaTG13 and BANAL-20-52. Note that some of the reads contain three key mutations relative to Wuhan-Hu-1 (C8782T, C18060T, and T28144C) that move the sequence closer to the bat coronavirus relatives. These mutations define one of the two plausible progenitors for all currently known human SARS-CoV-2 sequences (see [Kumar et al (2021)](https://academic.oup.com/mbe/article/38/8/3046/6257226) and [Bloom (2021)](https://academic.oup.com/mbe/article/38/12/5211/6353034)).

## Archived links after initially hearing about pre-print
I archived the following links on Dec-23-2021 after hearing about the pre-print from Carl Zimmer:

  - [https://assets.researchsquare.com/files/rs-1177047/v1_covered.pdf](https://web.archive.org/web/20211223025110/https://assets.researchsquare.com/files/rs-1177047/v1_covered.pdf)
  - [https://www.frontiersin.org/articles/10.3389/fmicb.2020.573302/full](https://web.archive.org/web/20210512222144/https://www.frontiersin.org/articles/10.3389/fmicb.2020.573302/full)
  - [https://bg.copernicus.org/articles/16/4113/2019/](https://web.archive.org/web/20211223060912/https://bg.copernicus.org/articles/16/4113/2019/)
  - [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA692319/](https://web.archive.org/web/20211223055435/https://www.ncbi.nlm.nih.gov/bioproject/PRJNA692319/)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441700](https://archive.md/971aq)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441701](https://archive.md/flkdT)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441702](https://archive.md/EYEDf)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441703](https://archive.md/5ED3i)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441704](https://archive.md/vhYsE)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441705](https://archive.md/OORe2)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441706](https://archive.md/X3hT5)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441707](https://archive.md/1bgWa)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441708](https://archive.md/yGU5Z)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441709](https://archive.md/3gAXy)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441710](https://archive.md/g7VaV)
  - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441711](https://archive.md/TzUMX)

## Deletion of some samples from SRA
On Jan-3-2022, I received an e-mail one of the pre-print authors, Istvan Csabai, saying that three of the samples (appearing to be the ones with the most SARS-CoV-2 reads) had been removed from the SRA.
He also noted that _bioRxiv_ had refused to publish their pre-print without explanation; the file he attached indicates the submission ID was `BIORXIV-2021-472446v1`.
I confirmed that three of the accessions had indeed been removed from the SRA as shown in the following archived links:
 - [https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=692319](https://web.archive.org/web/20220104013243/https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=692319)
 - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441704](https://web.archive.org/web/20220104013347/https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441704)
 - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441705](https://web.archive.org/web/20220104013624/https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441705)
 - [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441708](https://web.archive.org/web/20220104013634/https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441708)

I also e-mailed Richard Sever at _bioRxiv_ to ask why the pre-print was rejected, and explained I had repeated and validated the key findings.
Richard Sever said he could not give details about the pre-print review process, but that in the future the authors could appeal if they thought the rejection was unfounded.

## Details from Istvan Csabai
On Jan-4-2022, I chatted with Istvan Csabai.
He had contacted the authors of the pre-print, and shared their reply to him.
The authors had prepped the samples in early 2019, and submitted to Sangon BioTech for sequencing in December, getting the results back in early January.

## Second pre-print from Csabai and restoration of deleted files
Istvan Csabai then worked on a second pre-print that analyzed host reads and made various findings, including co-contamination with African green monkey (Vero?) and human DNA.
He sent me pre-print drafts on Jan-16-2022 and on Jan-24-2022, and I provided comments on both drafts and agreed to be listed in the Acknowledgments.

On Feb-3-2022, Istvan Csabai told me that the second pre-print had also been rejected from _bioRxiv_.
Because I had previously contacted Richard Sever when I heard the first pre-print was rejected, I suggested Istvan could CC me on an e-mail to Richard Sever appealing the rejection, which he did.
Unfortunately, Richard Sever declined the appeal, so instead Istvan posted the pre-print on [Resarch Square](https://www.researchsquare.com/article/rs-1330800/v1).

At that point on Feb-3-2022, I also re-checked the three deletion accessions (SRR13441704, SRR13441705, and SRR13441708).
To my surprise, all three were now again available by public access.
Here are archived links demonstrating that they were again available:
 - still missing from this overview page: [https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP301869&o=acc_s%3Aa](https://archive.is/AuOj2)
 - again active: [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441704](https://archive.is/2o7WK)
 - again active [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441705](https://archive.is/OORe2)
 - again active [https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13441708](https://archive.is/yGU5Z)

I confirmed that the replaced accessions were identical to the deleted ones.

## Inquiry to authors of PRJNA692319
On Feb-8-2022, I e-mailed the Chinese authors of the paper to ask about the sample deletion and restoration.
They e-mailed back almost immediately.
They confirmed what they had told Istvan: they had sequenced the samples with Sangon Biotech (Shanghai) after extracting the DNA in December 2019 from their samples.
The suspect that contamination of the samples happened at Sangon Biotech.
They deleted the three most contaminated samples from the Sequence Read Archive.
They do not know why the samples were then "un-deleted."

