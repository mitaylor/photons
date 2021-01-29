#!/usr/bin/env bash

histogram() {
    tag=$1

    ./bin/populate configs/populate_${tag}_raw.conf \
        data/populate_${tag}_raw.root &
    ./bin/populate configs/populate_${tag}_bkg.conf \
        data/populate_${tag}_bkg.root &
}

arithmetic() {
    tag=$1

    ./bin/manipulate configs/manipulate_${tag}.conf \
        data/manipulate_${tag}.root
    ./bin/accumulate configs/accumulate_${tag}.conf \
        data/accumulate_${tag}.root
}

nominal() {
    sample=$1

    histogram ${sample}
    wait

    arithmetic ${sample}
    arithmetic ${sample}_loose_purity
    arithmetic ${sample}_tight_purity
}

systematic() {
    sample=$1
    syst=$2

    histogram ${sample}_${syst}
    wait

    arithmetic ${sample}_${syst}
}

samples=(pp aa pp_smear_0_10 pp_smear_10_30 pp_smear_30_50 pp_smear_50_90)

for sample in ${samples[@]}; do
    nominal $sample

    for syst in jeu_down jeu_up wo_ele_rej qcd qcd_gen_iso; do
        systematic $sample $syst
    done

    # manipulation of indirect systematics
    ./bin/granulate configs/granulate_${sample}.conf \
        data/granulate_${sample}.root

    # final systematic uncertainties
    ./bin/obnubilate configs/obnubilate_${sample}.conf \
        data/obnubilate_${sample}.root
done
