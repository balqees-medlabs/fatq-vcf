/*
 * hap.py: benchmark a query VCF against a GIAB truth set inside the high-confidence BED.
 * Emits the summary.csv (SNV/indel precision/recall/F1), parsed into MultiQC custom content.
 */
process HAPPY {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.summary.csv"),      emit: summary
    tuple val(meta), path("*.extended.csv"),     emit: extended, optional: true
    path "*_happy_mqc.tsv",                       emit: mqc
    path "versions.yml",                          emit: versions

    script:
    def prefix = "${meta.id}.happy"
    def args   = task.ext.args ?: ''
    """
    export HGREF=${fasta}
    hap.py ${truth_vcf} ${query_vcf} \\
        -f ${truth_bed} -r ${fasta} -o ${prefix} \\
        --threads ${task.cpus} ${args}

    # Summary → MultiQC custom table (SNP/INDEL PASS rows)
    python3 - <<'PY'
import csv
rows = list(csv.DictReader(open("${prefix}.summary.csv")))
keep = [r for r in rows if r.get("Filter") == "PASS"]
with open("${meta.id}_happy_mqc.tsv", "w") as fh:
    fh.write("# id: happy_${meta.id}\\n# section_name: 'hap.py ${meta.id}'\\n# plot_type: 'table'\\n")
    fh.write("Type\\tRecall\\tPrecision\\tF1\\tTRUTH.TOTAL\\tFP\\tFN\\n")
    for r in keep:
        fh.write("\\t".join([r["Type"], r.get("METRIC.Recall",""), r.get("METRIC.Precision",""),
                             r.get("METRIC.F1_Score",""), r.get("TRUTH.TOTAL",""),
                             r.get("QUERY.FP",""), r.get("TRUTH.FN","")]) + "\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: \$(hap.py --version 2>&1 | grep -oP '[0-9.]+' | head -1 || echo 0.3.15)
    END_VERSIONS
    """
}
