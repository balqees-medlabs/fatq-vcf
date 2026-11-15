#!/usr/bin/env nextflow
/*
 * fatq-vcf — germline WGS/WES pipeline on GRCh38.
 *
 *   nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
 *   nextflow run main.nf -profile docker --prepare_references          # one-off reference bundle
 *
 * Design document: docs/PIPELINE_PLAN.md
 */

nextflow.enable.dsl = 2

include { FATQVCF            } from './workflows/fatqvcf'
include { PREPARE_REFERENCES } from './workflows/prepare_references'
include { VALIDATE           } from './workflows/validate'
include { PON_BUILD          } from './workflows/pon_build'

workflow {
    if (params.help) {
        log.info helpMessage()
        exit 0
    }
    if (params.prepare_references) {
        PREPARE_REFERENCES()
    }
    else if (params.validate) {
        VALIDATE()
    }
    else if (params.pon_build) {
        PON_BUILD()
    }
    else {
        if (!params.input) {
            error "ERROR: --input samplesheet.csv is required. See assets/samplesheet.example.csv"
        }
        FATQVCF( file(params.input, checkIfExists: true) )
    }
}

def helpMessage() {
    """
    fatq-vcf ${workflow.manifest.version}

    Usage:
      nextflow run main.nf -profile docker --input samplesheet.csv [--outdir results]

    Samplesheet columns: sample,lane,fastq_1,fastq_2,seq_type(wgs|wes),capture_bed,sex(male|female|unknown)

    Key options:
      --aligner        bwamem3 (default) | bwamem2
      --skip_bqsr      skip base quality score recalibration
      --skip_mito      skip mitochondrial calling
      --skip_sv_cnv    skip Manta / CNVpytor / gCNV
      --scatter_count  interval shards for BQSR + HaplotypeCaller (default 20)
      --gcnv_model_dir WES panel-of-normals model (enables gCNV for wes samples)

    One-off:
      nextflow run main.nf -profile docker --prepare_references [--test_data false] [--validation_data]
    """.stripIndent()
}
