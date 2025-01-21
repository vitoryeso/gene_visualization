# gene_visualization
Visualização de dados em fasta, para disciplina de bioinfo.

# Visualização de Genomas

Este projeto tem como objetivo visualizar dados de genomas em formato FASTA. Ele inclui scripts para baixar os dados, processá-los e gerar uma visualização interativa em 3D.

## Requisitos

- Python 3.6+
- Bibliotecas Python: numpy, biopython, transformers, torch, scikit-learn, plotly, tqdm

## Passos para Executar o Projeto

### 1. Clonar o Repositório

Clone o repositório para sua máquina local:

```sh
git clone https://github.com/vitoryeso/gene_visualization.git
cd gene_visualization
```

### 2. Instalar Dependências
```sh
pip install numpy biopython transformers torch scikit-learn plotly tqdm
```

### 3. Baixar dados
```sh
chmod +x download_fasta.sh
./download_fasta.sh
```

### 4. Gerar embeddings e criar visualização
```sh
python plot_fasta_data.py
```

### 5. Visualizar os resultados
Abra o arquivo sequences_3d.html gerado em seu navegador.