#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format

params.seqs = "$baseDir/data/*.fa"

params.refs = "$baseDir/data/*.ref"

// input guide trees in Newick format. Or `false` to generate trees
//params.trees = "/users/cn/egarriga/datasets/homfam/trees/*.{FAMSA,CLUSTALO,CLUSTALO-RANDOM,MAFFT_PARTTREE}.dnd"
params.trees = false

// which alignment methods to run
params.align_method = "FAMSA"

// which tree methods to run if `trees` == `false`
params.tree_method = "FAMSA"

// generate regressive alignments ?
params.regressive_align = true

// create progressive alignments ?
params.progressive_align = true

// create default alignments ? 
params.default_align = false

// evaluate alignments ?
params.evaluate = true

// bucket sizes for regressive algorithm
params.buckets= '1000'

// output directory
params.output = "$baseDir/resultsPROGRESSIVE"


log.info """\
         R E G R E S S I V E   M S A   A n a l y s i s  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA)               : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Output directory (DIRECTORY)                   : ${params.output}
         Alignment methods                              : ${params.align_method}
         Tree methods                                   : ${params.tree_method}
         Generate default alignments                    : ${params.default_align}
         Generate progressive alignments                : ${params.progressive_align}
         Generate regressive alignments (DPA)           : ${params.regressive_align}
         Bucket Sizes for regressive alignments         : ${params.buckets}
         Perform evaluation? Requires reference         : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.output}
         """
         .stripIndent()


// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  //.view()
  .into { seqsCh; seqs2 }
}

if ( params.refs ) {
  Channel
  .fromPath(params.refs)
  .map { item -> [ item.baseName, item] }
  .set { refs }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
    .set { trees }
}
else { 
  Channel
    .empty()
    .set { trees }
}

tree_methods = params.tree_method
align_methods = params.align_method


/*
 * GENERATE GUIDE TREES USING MEHTODS DEFINED WITH "--tree_method"
 *
 * NOTE: THIS IS ONLY IF GUIDE TREES ARE NOT PROVIDED BY THE USER
 * BY USING THE `--trees` PARAMETER
 */

process generate_trees {

    tag "${id}.${tree_method}"
    publishDir "${params.output}/guide_trees", mode: 'copy', overwrite: true
   
    input:
    set val(id), \
         file(seqs) \
         from seqsCh
     each tree_method from tree_methods.tokenize(',') 

   output:
     set val(id), \
       val(tree_method), \
       file("${id}.${tree_method}.dnd") \
       into treesGenerated

   when:
     !params.trees

   script:
   """
   famsa -gt_export ${id}.${tree_method}.dnd ${seqs} ${id}.aln
   """
}


treesGenerated
  .mix ( trees )
  .combine ( seqs2, by:0 )
  .into { seqsAndTreesForStandardAlignment; seqsAndTreesForRegressiveAlignment; seqsAndTreesForProgressiveAlignment }


process regressive_alignment {

    tag "${id}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input:
        set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForRegressiveAlignment

      each bucket_size from params.buckets.tokenize(',')

      each align_method from align_methods.tokenize(',')

    when:
      params.regressive_align

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("reg_align"), \
        val(bucket_size), \
        file("${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.aln") \
        into regressiveOut

    script:
    """
        t_coffee -reg -reg_method famsa_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.aln
    """
}

process progressive_alignment {

    tag "${id}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input:
        set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForProgressiveAlignment

    each align_method from align_methods.tokenize(',')

    when:
      params.progressive_align

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("prog_align"), \
        val("NA"), \
        file("${id}.prog_align.NA.${align_method}.with.${tree_method}.tree.aln") \
        into progressiveOut

    script:
    """
      famsa -gt_import ${guide_tree} ${seqs} ${id}.prog_align.NA.${align_method}.with.${tree_method}.tree.aln
    """
}

progressiveOut
  .mix ( regressiveOut )
  .set { all_alignments }

refs
  .cross (all_alignments )
  .map { it -> [it[0][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5], it[0][1]] }
  .set { toEvaluate }

process evaluate {
    
    tag "${id}.${align_method}.${tree_method}.${align_type}.${bucket_size}"
    publishDir "${params.output}/individual_scores", mode: 'copy', overwrite: true

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(align_type), \
          val(bucket_size), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set val(id), \
          val(tree_method), \
          val(align_method), \
          val(align_type), \
          val(bucket_size), \
          file("*.sp"), \
          file("*.tc"), \
          file("*.col") \
          into scores

    when:
      params.evaluate

     script:
     """
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"

    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
