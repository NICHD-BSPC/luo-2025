import pybedtools
from Bio import Seq
from Bio import SeqIO
import pybedtools
import pyfaidx
import primer3
import itertools
import tempfile
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_dna_complement
import numpy as np
import seaborn as sns

from . import constants
from . import utils


class CandidateFinder(object):

    def __init__(
        self,
        fasta,
        annotations,
        pam,
        guide_length,
        min_gc=0,
        max_gc=100,
        het_tm=40,
        end_tm=30,
        bowtie_index_prefix="tmp/index",
        outfn='candidates',
        plots=True,
    ):
        """
        Represents a collection of sgRNA candidates.

        Parameters
        ----------

        fasta : str
            Path to genome fasta file

        annotations : str or pybedtools.BedTool
            If str, should be path to BED or GTF file that will be converted to
            pybedtools.BedTool file. This conversion will automatically handle
            1-based GTF and 0-based BED annotations to avoid off-by-one errors.

        pam : str
            PAM sequence, e.g. "NGG" for wildtype Cas9

        guide_length : int
            How long of a guide to select

        min_gc, max_gc : int
            Range of acceptable GC percentage (inclusive)

        het_tm : int
            Degrees C threshold for Primer3's heterodimer routine

        end_tm : int
            Degrees C threshold for Primer3's end stability routine

        bowtie_index_prefix : str
            Path to directory containing existing bowtie2 index. If missing,
            index will be created.

        outfn : str
            Output files path and basename. Default is 'candidates'
        plots : bool
            Generate coverage plots
        """
        fasta = str(fasta)
        self.fasta_file = fasta

        self.fasta = pyfaidx.Fasta(fasta)

        self.genome = fasta + ".genome"
        with open(fasta + ".fai") as fin:
            with open(self.genome, "w") as fout:
                for line in fin:
                    toks = line.strip().split("\t")
                    fout.write(f"{toks[0]}\t{toks[1]}\n")

        self.annotations = pybedtools.BedTool(annotations)
        self.pam = pam
        self.pam_length = len(pam)
        self.pam_seqs = list(utils.expand_ambiguous(pam))
        self.rev_pam_seqs = list(utils.expand_ambiguous(Seq.Seq(pam).reverse_complement()))
        self.guide_length = guide_length
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.het_tm = het_tm
        self.end_tm = end_tm
        self.bowtie_index_prefix = bowtie_index_prefix
        self.outputtmp = outfn
        self.outputfiles = self.filtered_candidates_to_files(outfn)
        self.plotfile = self.plot_grnas_coverage(outfn, plots)

    def populate_candidates(self, force=True):
        """
        Stores a (cached) dict on the object. Keys are integers.
        """
        if hasattr(self, "_candidates") and not force:
            return self._candidates
        self._candidates = {}
        i = 0
        for annotation in self.annotations:
            pyfaidx_seq_obj = self.fasta.get_seq(
                annotation.chrom, annotation.start, annotation.stop
            )
            seq = Seq.Seq(pyfaidx_seq_obj.seq)
            for guide, start, stop, strand in find_guides(
                seq, self.pam_seqs, self.guide_length
            ):
                i += 1
                cg = CandidateGuide(
                    guide,
                    annotation.chrom,
                    annotation.start,
                    annotation.stop,
                    start,
                    stop,
                    id=i,
                    strand=strand,
                )
                self._candidates[i] = cg
        return self._candidates

    @property
    def candidates(self):
        if not hasattr(self, "_candidates"):
            self.populate_candidates()
        return self._candidates

    def candidates_to_file(self, output):
        with open(output, "w") as fout:
            for candidate in self.candidates.values():
                fout.write(f">{candidate.id}\n{candidate.seq}\n")
        return output



    def off_target_filter(self):
        """
        Returns a list of keys corresponding to off-target guides.
        """
        # Build an index if it doesn't already exist
        utils.bowtie_index(self.fasta_file, self.bowtie_index_prefix)

        # Save current candidates to file and align
        candidates = self.candidates_to_file(self.outputtmp + "_intermediate.fa")
        utils.bowtie(candidates, self.bowtie_index_prefix, self.outputtmp + "_output.tsv")

        offtarget = []
        hits = (
            pybedtools.BedTool(self.outputtmp + "_output.tsv")
            .bamtobed()

            # only return those that do NOT map inside provided annotations
            .intersect(self.annotations, v=True)

            # of those, get the flanking PAM (in a strand-specific manner)
            .flank(l=0, r=self.pam_length, g=self.genome, s=True)
            .saveas(self.outputtmp + "_flank.bed")
        )

        pam_seqs_dict = {'+': self.pam_seqs, '-': self.rev_pam_seqs}
        for hit in hits:
            flanking = self.fasta.get_seq(hit.chrom, hit.start + 1, hit.stop)
            pam_seqs = pam_seqs_dict[hit.strand]
            if flanking.seq in pam_seqs:
                offtarget.append(int(hit.name))
                self.candidates[int(hit.name)].filters["offtarget"] = flanking.seq
        return offtarget

    def gc_filter(self):
        """
        Returns a list of keys corresponding to out-of-range GC content"
        """
        gc_fail = []
        for i, candidate in self.candidates.items():
            if not candidate.validate_gc(self.min_gc, self.max_gc):
                gc_fail.append(i)
        return gc_fail

    def heterodimer_filter(self):
        """
        Returns a list of keys corresponding to possible heterodimers
        """
        heterodimer_fail = []
        for i, candidate in self.candidates.items():
            if not candidate.validate_heterodimer(self.het_tm):
                heterodimer_fail.append(i)
        return heterodimer_fail

    def end_stability_filter(self):
        """
        Returns a list of keys corresponding to possible end_stabilitys
        """
        end_stability_fail = []
        for i, candidate in self.candidates.items():
            if not candidate.validate_end_stability(self.end_tm):
                end_stability_fail.append(i)
        return end_stability_fail

    def filtered_candidates(self):
        """
        Returns a dictionary of keys = id, values = candidate object
        """
        gc = self.gc_filter()
        hd = self.heterodimer_filter()
        es = self.end_stability_filter()
        ot = self.off_target_filter()
        remove = set(gc + hd + es + ot)
        return {k:v for k, v in self.candidates.items() if k not in remove}

    def filtered_candidates_to_files(self, outfn):
        """
        Side effect of creating FASTA and BED of filtered guides
        and CSV of oligos
        Return dict of filenames
        """
        filtered = self.filtered_candidates()
        seqs = []
        with open(outfn + '_grnas.fa', "w") as fout, \
                open(outfn + '_oligos.csv', "w") as csvout, \
                open(outfn + '_grnas.bed', "w") as bedout:
            for filtcandidate in filtered.values():
                bedstart = filtcandidate.start -1
                bedstop = filtcandidate.stop -1
                bedout.write(
                    f"{filtcandidate.contig}\t" +
                    f"{bedstart}\t" +
                    f"{bedstop}\t" +
                    f"{filtcandidate.id}_{filtcandidate.seq}\t.\t" +
                    f"{filtcandidate.strand}\n")
                if filtcandidate.seq not in seqs:
                    seqs.append(filtcandidate.seq)
                    fout.write(f">{filtcandidate.id}\n{filtcandidate.seq}\n")
                    csvout.write(f"{filtcandidate.id},{filtcandidate.oligo}{constants.PRIMER_REVCOMP}\n")
        return {'fa':outfn + '_grnas.fa',
                'csv':outfn + '_oligos.csv',
                'bed':outfn + '_grnas.bed'}

    def plot_grnas_coverage(self, outfn, plots):
        """
        if plots=True, generate coverage plots
        for each interval in self.annotation
        of output bed stored in self.outputfiles
        """
        if plots == False:
            return(False)
        ann = self.annotations
        grnas = pybedtools.BedTool(self.outputfiles['bed'])
        utils.plot_coverage(ann, grnas, outfn)
        return(True)


class CandidateGuide(object):
    """
    Represents a single candidate guide.

    Has validation routines.
    """

    def __init__(
        self,
        seq,
        contig,
        annotation_start,
        annotation_stop,
        guide_start,
        guide_stop,
        id,
        strand,
    ):
        """
        Parameters
        ----------

        seq : str
            Represents the guide RNA

        contig, start, stop : str, int, int
            Coordinates of the parent transcript
        """
        self.seq = seq
        self.id = id
        self.strand = strand
        self.contig = contig
        self.annotation_start = annotation_start
        self.annotation_stop = annotation_stop
        self.guide_start = guide_start
        self.guide_stop = guide_stop
        self.start = annotation_start + guide_start
        self.stop = annotation_start + guide_stop
        self.oligo = self._compose_oligo()
        self.filters = {}

    def _gc(self):
        return (self.seq.count("C") + self.seq.count("G")) / len(self.seq) * 100

    def _compose_oligo(self):
        if self.seq.startswith("GG"):
            spacer = ""
        elif self.seq.startswith("G"):
            spacer = "G"
        else:
            spacer = "GG"
        return constants.T7 + spacer + self.seq + constants.SCAFFOLD #+ constants.PRIMER_REVCOMP

    def validate_gc(self, min_gc, max_gc):
        """
        Returns False if GC content is outside requested range.

        If validation fails, update self.filters.
        """
        gc = self._gc()
        gc_pass = min_gc <= gc <= max_gc
        if not gc_pass:
            self.filters["GC"] = gc
        return gc_pass

    def calc_heterodimer(self):
        return primer3.bindings.calcHeterodimer(
            constants.PRIMER_REVCOMP, str(self.oligo)
        ).tm

    def calc_end_stability(self):
        return primer3.bindings.calcEndStability(
            constants.PRIMER_REVCOMP, str(self.oligo)
        ).tm

    def validate_heterodimer(self, het_tm):
        """
        Returns False if the melting temperature of a heterodimer between
        constants.PRIMER_REVCOMP and the guide sequence is over `het_tm`
        (indicating high predicted internal binding).

        If validation fails, update self.filters.
        """
        tm = self.calc_heterodimer()
        valid = tm <= het_tm
        if not valid:
            self.filters["het"] = tm
        return valid

    def validate_end_stability(self, end_tm):
        """
        Returns False if the melting temperature of the end stability between
        constants.PRIMER_REVCOMP and this guide is over `end_tm`.

        If validation fails, update self.filters.
        """
        tm = self.calc_end_stability()
        valid = tm <= end_tm
        if not valid:
            self.filters["end"] = tm
        return valid

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.contig}:{self.start}-{self.stop} [{self.strand}]>"


def find_guides(seq, pams, guide_length):
    """
    Given a sequence (e.g., of an rRNA gene), a list of all possible expansions
    of a PAM sequence, and a length, yields tuples of (guide, start, stop).
    Start and stop are always in the coordinate space of the provided seq, even
    if found on the reverse complement.
    """

    def _candidates_in_range(seq, thispam, strand):
        seqlen = len(seq)

        # Original implementation seems to have had an off-by-one error (or
        # maybe only truly should be given 1-based GTF rather than 0-based
        # BED files). BED files are 0-based start, 1-based end, (-) and (+)
        # strands have different offset.
        # This can be controlled with the constants.LEGACY flag.
        if constants.LEGACY:
            FIRST_BASE = 1
            if strand == '+':
                OFFSET = 0
            else:
                OFFSET = 1
        else:
            FIRST_BASE = 0
            OFFSET = 1

        for i in utils.find_all_substrings(seq, thispam):
            start = i - guide_length

            if start <= (FIRST_BASE - OFFSET):
                continue

            candidate = seq[start:i]
            if strand == "-":
                start = seqlen - i
                end = start + guide_length
            else:
                end = i
            yield candidate, start, end, strand

    for pam in pams:
        for result in _candidates_in_range(seq, pam, strand="+"):
            yield result
        for result in _candidates_in_range(seq.reverse_complement(), pam, strand="-"):
            yield result
