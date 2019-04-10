# Chem-seq Read Enrichment Discovery: CRED

**C**hem-seq **R**ead **E**nrichment **D**iscovery (__CRED__) is a simple peak caller written in C for identifying non-canonical feature enrichments in paired Chem-seq data.

Committed April 10, 2019

### Installation

CRED adheres requires a compatible compiler (e.g. GCC), and utilizes [HTSlib](http://www.htslib.org/) to process BAM files. Please clone the repo inside CRED before compiling. Additionally, `make` on certain macOS systems may not automatically compile HTSLib with CRED, an additional call of `make` within `htslib/` may be required in some cases.

Use `PREFIX` to specify a directory to install CRED with `make install`:

Sample installation guide:
```
export INSTALL_PREFIX=/usr/local/bin # or a different location
git clone https://github.com/jlincbio/cred.git
cd cred
git clone https://github.com/samtools/htslib.git
cd htslib && make # required on macOS 10.14 (GNU Make 3.81)
cd ..
make
make PREFIX=$INSTALL_PREFIX install
```


### Inputs and Parameters
```
CRED: Chem-seq Read Enrichment Discovery (Version 0.1, Apr 2019 Initial Release)
Command  :  cred [options] -t TREATMENT.BAM -c CONTROL.BAM > OUTPUT.BED
Required :
  -t  TREATMENT.BAM  Path to the treatment ("pulldown") track [BAM]
  -c  CONTROL.BAM    Path to the control ("input") Chem-seq track [BAM]
Optional :
  -p  [P-VALUE]      Significance level [default 0.0001]
  -q  [SCORE]        Minimum MAPQ quality for reads to count [default 30]
  -w  [INTEGER]      Size of differential windows [default 1200 bp]
  -k                 Evaluate site significance with Kolmogorov-Smirnov
Reminders:
   1. BAM files must be sorted and indexed.
   2. Use a pipe (">") to capture CRED output.
```

* BAM files for the treatment and control BAM's, respectively: use `-t` and `-c` to specify the pair. BAM's should be aligned, coordinate-sorted AND indexed (both .bam and .bam.bai should be present).
* Quality score cutoff (MAPQ, option `-q`): this is the minimum required mapping quality score as defined in the SAM format specification. CRED defaults to 30. 
* Significance (alpha) level (option `-p`): please specify this either as a decimal (e.g. `0.0001`) or a fraction (e.g. `1/10000`). CRED will check all regions against this cutoff and output only features more significant than this predefined alpha level.
* Size of sliding windows may also be adjusted with option `-w`; defaults to 1200 bp (maximum size used in Lin et al. PLoS ONE 2016).
* Method of evaluating the significance of enrichment may also be modified from Welch's t-test (default) to Kolmogorov-Smirnov by the `-k` toggle.

A helper program, "Batch CRED" or `BCRED`, is also supplied here for running CRED on more hardware-limited systems (e.g., those without access to a lot of RAM). `BCRED` is written in Perl, and requires `samtools` as well as `Parallel::ForkManager` (a Perl module available on CPAN) for operation; this will split the input BAM pairs per chromosome, dispatch `CRED` calls, and merge the result. Multithreading support is also made possible via `Parallel::ForkManager`.

To launch `BCRED`, look for `bcred` in the same folder as `cred`. The inputs and parameters are largely the same as `cred`:

```
bcred: a batch assistant for CRED Ver. 0.1 (Apr 2019 Initial Release)
Usage: bcred [options] -t TREATMENT.BAM -c CONTROL.BAM > OUTPUT.BED
Required :
  -t  TREATMENT.BAM  Path to the treatment ("pulldown") track [BAM]
  -c  CONTROL.BAM    Path to the control ("input") Chem-seq track [BAM]
Optional :
  -p  [P-VALUE]      Significance level [default 0.0001]
  -q  [SCORE]        Minimum MAPQ quality for reads to count [default 30]
  -w  [INTEGER]      Size of differential windows [default 1200 bp]
  -n  [INTEGER]      Number of threads to utilize [default 1]
  -k                 Evaluate site significance with Kolmogorov-Smirnov
Reminders:
   1. BAM files must be sorted.
   2. Use a pipe (">") to capture CRED output.
```


### Output
The current version of CRED writes to STDOUT so the results can be streamed in-line for subsequent tasks, e.g. checking for motif intersects with BEDTools and immediate compressing the results with GZip. To store the output to a file, use a pipe (">"). The output is presented in a BED-like format directly interpretable in genome browsers such as IGV. The columns are as follows:

1. Chromosome ID
2. Start of the site
3. End position of the site
4. Peak ID (numerically ordered)
5. log ratio of relative enrichment in the pulldown vs. control track
6. (Unused)
7. Significance level of the relative enrichment


### Example
Under `samples/` there is a set of simulated treatment and control BAM's one can use to test CRED:
* `simReads_hg19_treatment-chr20.bam`
* `simReads_hg19_control-chr20.bam`

Those files were created by using [DWGSIM](https://github.com/nh13/DWGSIM) with a BED file containing a list of reference features ("sites") as true positives (N = 3M), followed by random regions as control (N = 50M). The treatment track was compiled by merging the two to ensure genomic enrichment (see `simReads-generate.pl` for an example script). Following processing and alignment by [LAST](http://last.cbrc.jp/), reads located in chr20 were extracted so that the file sizes will be under 100M per GitHub rules.

To try out CRED with these two tracks this way with default settings:
```
cred -t samples/simReads_hg19_treatment-chr20.bam -c samples/simReads_hg19_control-chr20.bam > simReads_chr20-cred.bed
```

If BEDTools is installed, intersecting the resultant BEDs with the true positives would reveal that the CRED results includes more sites containing the positive spikes (700+) compared to MACS (~500). These regions can also be confirmed by IGV (see `results_sample_igv_snapshot.png` for an example). On a 3.5GHz 6-core Mac Pro with 64GB of RAM running MacOS 10.14.4, the CRED run completed ~20 seconds for CRED and about a minute for MACS.
