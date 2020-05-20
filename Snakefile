

from glob import glob 
import pandas as pd 
import os 

       




def find_pair(folder):
    files = [i for i in glob(f"{folder}/*.fastq.gz") if "Undetermined" not in i ]
    R1 = next(filter(lambda x : "R1" in x, files))
    R2 = next(filter(lambda x : "R2" in x, files ))

    return R1,R2


rule all:
    input:
        "191/stat.demux.txt",
        "192/stat.demux.txt",
        "190-701/stat.demux.txt",
        "190-702/stat.demux.txt",
        "190-703/stat.demux.txt",
        "190-705/stat.demux.txt",

        "191/stat.txt",
        "192/stat.txt",
        "190-701/stat.txt",
        "190-702/stat.txt",
        "190-703/stat.txt",
        "190-705/stat.txt"


rule merge:
    input: lambda wildcards: find_pair(wildcards.folder)

    output:
        "{folder}/merged.fastq"

    log:
        "{folder}/merged.log"

    shell:
        "flash --max-overlap 200 {input} -d {wildcards.folder}  > {wildcards.folder}/out.flash.log ;"
        "mv {wildcards.folder}/out.extendedFrags.fastq {output}"


rule create_sabre_file:
    input:
        "{folder}/barcode.txt"
    output:
        "{folder}/sabre.barcode.txt"
    run:
        dirname = os.path.dirname(input[0])
        df = pd.read_csv(input[0], sep="\t", header=None)
        #df.columns = ["name", "barcode",""]
        df = df.rename({0:"name", 1:"barcode"}, axis=1)
        df["file"] = dirname +f"/demux/" + df["name"] + ".fastq"
        df[["barcode", "file"]].to_csv(output[0], index=False, sep="\t", header=False)


rule reverse_and_trim:
    input:
        "{folder}/merged.fastq"
    output:
        "{folder}/merged.trim.fastq"
    shell:
        "seqkit seq -p -r -t dna {input} | seqkit subseq -r 23:-1> {output}"

rule demux:
    input:
        fastq = "{folder}/merged.trim.fastq",
        barcode = "{folder}/sabre.barcode.txt"
    params:
        mismatch = 2
    output:
        directory("{folder}/demux"),
        unknown = "{folder}/demux.unknown.fastq"

    log:
        "{folder}/demux.log"

    shell:
        "mkdir -p {wildcards.folder}/demux;"
        "sabre se -m {params.mismatch} -f {input.fastq} -b {input.barcode} -u {output.unknown}  > {log}"


rule stats:
    input:
        unmerged = lambda wildcards: find_pair(wildcards.folder)[0],
        merged = "{folder}/merged.fastq",
        trim = "{folder}/merged.trim.fastq",
        undemux = "{folder}/demux.unknown.fastq"
    output:
        "{folder}/stat.txt"
    run:
        shell(f"touch {output}")
        data = []
        for i in input:
            generator = shell(f"seqkit stat -b {i}|tail -n+2|sed -E 's/[[:space:]]+/;/g'", iterable=True)
            data.append(next(generator).strip().split(";"))
        df = pd.DataFrame(data)
        df.columns = ["name", "type","dna","count","total count","min","mean","max"]
        df["count"] = df["count"].str.replace(",","").astype(int)

        df["percent"] = df["count"]  /df["count"][0] * 100.0
        df["percent"] = df["percent"].round(2)
        df[["name","count","min","mean","max", "percent"]].to_csv(output[0], index=False, sep="\t")


rule demux_stats:
    input:
        "{folder}/demux"
    output:
        "{folder}/stat.demux.txt"

    run:
        seqs = {
            "FORWARD_PRIMER_N1": "CATTTCGCTGATTTTGGGGTC",
            "FORWARD_PRIMER_N2": "GCCAATGTTTGTAATCAGTTCCTT",
            "FORWARD_PRIMER_HUMAN": "CGCTCGCAGGTCCAAATC",

            "N1": "GCCAGTTGAATCTGAGGGTCCACCAAACGTAATGCGGGGTG",
            "N2": "CGACATTCCGAAGAACGCTGAAGCGCTGGGGGCAAATTGTGCAATTTGCG",
            "HUMAIN": "AAGTCCGCGCAGAGCCTTCAGGTCAGAACCCGCTCGCAGGTCCAAATC",

            "SPIKE_N1": "CCTAGCCCTATCTCCGTGGCCGTTCAGCGAAAGTCTGGGACCGACAAATATCTCGCCGACAACATTTATTCTGATGCCTGGTGTCATCACCAATGGCGCC",
            "SPIKE_N2": "AACGCCAAGGGGCGTAATGTGCAACTACAAGCGTTGTATACAGGAACGTAGGGCCGACTGCGGACTTTTAGGCAATGTAGTTGGAGTGAACTTTTTGCAA",
            "SPIKE_HUMAN":"GACCTTAGGAATGCGTGAACATATCCCCATTGATGATTCCAGATGTCCAACAAGCGAGACCCAACATTGGAACCGTCCAATGCTTAGCCATGGAGTCTTG"
        }

        all_df = []
        for i in glob(input[0]+"/*.fastq"):
            count = {}
            basename = os.path.basename(i)

            count["TOTAL"] = next(shell(f"seqkit seq -s {i}|wc -l", iterable=True))

            for key in seqs.keys():
                subseq = seqs[key]
                result = shell(f"seqkit grep -s -p {subseq} {i}|seqkit seq -s |wc -l", iterable=True)
                count[key] = next(result)
            
            all_df.append(pd.DataFrame(count, index = [basename]))

        all_df = pd.concat(all_df)
        all_df = all_df.rename({0 : "file"}, axis=1)
        all_df = all_df.astype(int)

        all_df["pN1"] = all_df["N1"] / all_df["TOTAL"]
        all_df["pN2"] = all_df["N2"] / all_df["TOTAL"]
        all_df["pH"] = all_df["HUMAIN"] / all_df["TOTAL"]

        all_df["pSN1"] = all_df["SPIKE_N1"] / all_df["TOTAL"]
        all_df["pSN2"] = all_df["SPIKE_N2"] / all_df["TOTAL"]
        all_df["pSH"] = all_df["SPIKE_HUMAN"] / all_df["TOTAL"]

        all_df.to_csv(output[0])
