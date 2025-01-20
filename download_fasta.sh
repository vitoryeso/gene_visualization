#!/bin/bash

# Diretório para salvar os arquivos FASTA
output_dir="bio_data"
mkdir -p "$output_dir"

# Lista de IDs de genomas ou proteínas para baixar
ids=(
    # Bactérias
    "NC_000913.3"  # Escherichia coli K-12
    "NC_002516.2"  # Pseudomonas aeruginosa PAO1
    "NC_000962.3"  # Mycobacterium tuberculosis

    # Vírus
    "NC_001416.1"  # Vírus lambda
    "NC_001802.1"  # HIV-1
    "NC_010210.1"  # Influenza A

    # Eucariotos
    "NC_001301.1"  # Homo sapiens (humano)
    "NC_000913.3"  # Saccharomyces cerevisiae (levedura)

    # Arqueias
    "NC_002693.2"  # Methanobacterium thermoautotrophicum
    "NC_003879.4"  # Haloferax volcanii
)

# Base URL do NCBI para baixar genomas
base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Loop para baixar cada genoma/proteína
for id in "${ids[@]}"; do
    echo "Baixando $id..."
    curl -o "$output_dir/${id}.fasta" \
        "$base_url?db=nuccore&id=$id&rettype=fasta&retmode=text"
done

echo "Downloads concluídos! Arquivos salvos em $output_dir/"
