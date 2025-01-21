#!/bin/bash

# Diretório para salvar os arquivos FASTA
output_dir="bio_data"
mkdir -p "$output_dir"

# Arrays para diferentes grupos de organismos
declare -A genome_groups=(
    ["bacteria"]="
        NC_000913.3   # E. coli K-12
        NC_002516.2   # P. aeruginosa PAO1
        NC_000962.3   # M. tuberculosis H37Rv
        NC_003198.1   # S. enterica Typhi
        NC_002695.2   # E. coli O157:H7
        NC_003112.2   # Neisseria meningitidis
        NC_002929.2   # Bordetella pertussis
        NC_003197.2   # Salmonella typhimurium
        NC_002505.1   # Vibrio cholerae
        NC_002737.2   # Streptococcus pyogenes"
        
    ["archaea"]="
        NC_000917.1   # Archaeoglobus fulgidus
        NC_000916.1   # Methanothermobacter thermautotrophicus
        NC_000961.1   # Thermoplasma acidophilum
        NC_009634.1   # Methanococcus maripaludis
        NC_015474.1   # Pyrococcus sp."
        
    ["viruses"]="
        NC_001802.1   # HIV-1
        NC_001416.1   # Bacteriophage lambda
        NC_001348.1   # HSV-1
        NC_001806.2   # HHV-6
        NC_001498.1   # HTLV-1
        NC_001722.1   # DENV-2
        NC_045512.2   # SARS-CoV-2
        NC_001477.1   # Dengue virus 1
        NC_001563.2   # West Nile virus
        NC_002549.1   # Zaire ebolavirus"
    
    ["eukaryotes"]="
        NC_001133.9   # S. cerevisiae chr 1
        NC_001134.8   # S. cerevisiae chr 2
        NC_001135.5   # S. cerevisiae chr 3
        NC_001224.1   # S. cerevisiae mitochondrion
        NC_037301.1   # D. melanogaster chr X"
)

# Base URL do NCBI
base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Função para baixar um genoma
download_genome() {
    local id=$1
    local group=$2
    echo "Baixando $id (Grupo: $group)..."
    
    # Criar subdiretório para o grupo
    mkdir -p "$output_dir/$group"
    
    # Baixar o arquivo
    curl -s -o "$output_dir/$group/${id}.fasta" \
        "$base_url?db=nuccore&id=$id&rettype=fasta&retmode=text"
    
    # Verificar o tamanho do arquivo
    local size=$(stat -f%z "$output_dir/$group/${id}.fasta" 2>/dev/null || stat -c%s "$output_dir/$group/${id}.fasta")
    echo "Tamanho do arquivo $id: $size bytes"
}

# Loop através dos grupos e seus genomas
for group in "${!genome_groups[@]}"; do
    echo "Processando grupo: $group"
    
    # Loop através dos IDs no grupo
    for id in ${genome_groups[$group]}; do
        # Ignorar linhas de comentário
        [[ $id =~ ^#.* ]] && continue
        
        # Remover espaços em branco
        id=$(echo $id | tr -d '[:space:]')
        
        # Se o ID não estiver vazio, fazer o download
        if [ ! -z "$id" ]; then
            download_genome "$id" "$group"
            # Adicionar um pequeno delay para não sobrecarregar o servidor
            sleep 2
        fi
    done
done

echo "Downloads concluídos! Arquivos salvos em $output_dir/"