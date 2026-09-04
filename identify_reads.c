#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "htslib/kseq.h"

KSEQ_INIT(gzFile, gzread)

#define BARCODE_LEN   20
#define MAX_SEQ_LEN   65536
#define MAX_ID_LEN    1024
#define MAX_REF_SEQS  16

typedef struct {
    char id[MAX_ID_LEN];
    char seq[MAX_SEQ_LEN];
    int  len;
} RefSeq;

static RefSeq refs[MAX_REF_SEQS];
static int    n_refs = 0;

static const char *get_ref(const char *id) {
    for (int i = 0; i < n_refs; i++)
        if (strcmp(refs[i].id, id) == 0)
            return refs[i].seq;
    return NULL;
}

static int get_ref_len(const char *id) {
    for (int i = 0; i < n_refs; i++)
        if (strcmp(refs[i].id, id) == 0)
            return refs[i].len;
    return 0;
}

static void load_references(const char *fasta_path) {
    gzFile fp = gzopen(fasta_path, "r");
    if (!fp) { fprintf(stderr, "Cannot open reference: %s\n", fasta_path); exit(1); }
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0 && n_refs < MAX_REF_SEQS) {
        strncpy(refs[n_refs].id,  seq->name.s, MAX_ID_LEN  - 1);
        strncpy(refs[n_refs].seq, seq->seq.s,  MAX_SEQ_LEN - 1);
        refs[n_refs].len = (int)seq->seq.l;
        n_refs++;
    }
    kseq_destroy(seq);
    gzclose(fp);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <lib> <readfile> <suffix>\n", argv[0]);
        return 1;
    }
    const char *lib      = argv[1];
    const char *readfile = argv[2];
    const char *suff     = argv[3];

    load_references("./alignments/reference.fasta");

    const char *left_flank = get_ref("pKMW7_U:L");
    const char *tn5_seq    = get_ref("Tn5_recognition");

    if (!left_flank || !tn5_seq) {
        fprintf(stderr, "Missing required reference sequences in FASTA.\n");
        return 1;
    }
    int left_len = get_ref_len("pKMW7_U:L");
    int tn5_len  = get_ref_len("Tn5_recognition");

    /* Output CSV */
    char outpath[2048];
    snprintf(outpath, sizeof(outpath),
             "processed/%s%s%s.csv", lib, readfile, suff);
    FILE *outfile = fopen(outpath, "w");
    if (!outfile) { fprintf(stderr, "Cannot open output: %s\n", outpath); return 1; }
    fprintf(outfile, "read_id,umi,barcode,genomic_region,genomic_phredscore,prunedGs\n");

    /* Input FASTQ (gzipped) */
    char inpath[2048];
    snprintf(inpath, sizeof(inpath),
             "/n/groups/springer/amogh/data/2025-12-05-xa65-transposon-screen/Reads/%s%s_001.fastq.gz",
             lib, readfile);
    gzFile fp = gzopen(inpath, "r");
    if (!fp) { fprintf(stderr, "Cannot open input: %s\n", inpath); return 1; }
    kseq_t *seq = kseq_init(fp);

    long count = 0, count_pass = 0;

    while (kseq_read(seq) >= 0) {
        count++;
        const char *read = seq->seq.s;
        const char *qual = seq->qual.s;
        int read_len     = (int)seq->seq.l;

        /* 1. Find left flank; UMI is everything before it */
        const char *left_pos = strstr(read, left_flank);
        if (!left_pos) continue;

        int left_start = (int)(left_pos - read);
        int left_end   = left_start + left_len;
        int umi_len    = left_start;
        int bc_start   = left_end;
        int bc_end     = bc_start + BARCODE_LEN;
        if (bc_end > read_len) continue;

        /* 2. Find Tn5 recognition sequence */
        const char *tn5_pos = strstr(read, tn5_seq);
        if (!tn5_pos) continue;

        int tn5_end = (int)(tn5_pos - read) + tn5_len;
        if (tn5_end >= read_len) continue;

        /* 3. Genomic region runs to end of read */
        int genome_start = tn5_end;
        int genome_end   = read_len;
        int glen         = genome_end - genome_start;

        /* 4. Mean Phred quality score over genomic region */
        double phred_sum = 0.0;
        for (int i = genome_start; i < genome_end; i++)
            phred_sum += (unsigned char)qual[i] - 33;
        double phredscore = phred_sum / glen;

        /* 5. Write record: read_id, umi, barcode, genomic_region, phredscore, prunedGs */
        fprintf(outfile, "%s,%.*s,%.*s,%.*s,%.4f,0\n",
                seq->name.s,
                umi_len,    read,
                BARCODE_LEN, read + bc_start,
                glen,        read + genome_start,
                phredscore);
        count_pass++;
    }

    kseq_destroy(seq);
    gzclose(fp);
    fclose(outfile);

    /* Append to log */
    FILE *logfile = fopen("log.txt", "a");
    if (logfile) {
        fprintf(logfile, "Processed %s%s%s: %ld reads, %ld passed filters\n",
                lib, readfile, suff, count, count_pass);
        fclose(logfile);
    }
    printf("Done: %ld reads processed, %ld passed filters\n", count, count_pass);
    return 0;
}
