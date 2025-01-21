import numpy as np
from Bio import SeqIO
from transformers import T5EncoderModel, T5Tokenizer
import torch
from sklearn.manifold import TSNE
import plotly.graph_objects as go
import os
from tqdm import tqdm
import colorsys

def print_sequence_statistics(sequences, lengths, organism_types):
    print("\n=== Estatísticas das Sequências ===")
    
    # Estatísticas gerais
    total_seqs = len(sequences)
    total_bases = sum(lengths)
    avg_length = np.mean(lengths)
    median_length = np.median(lengths)
    std_length = np.std(lengths)
    
    print(f"\nEstatísticas Gerais:")
    print(f"Total de bases: {total_bases:,}")
    print(f"Comprimento médio: {avg_length:,.2f} bases")
    print(f"Comprimento mediano: {median_length:,.2f} bases")
    print(f"Desvio padrão: {std_length:,.2f} bases")
    print(f"Menor sequência: {min(lengths):,} bases")
    print(f"Maior sequência: {max(lengths):,} bases")
    
    # Estatísticas por tipo de organismo
    print("\nEstatísticas por Tipo de Organismo:")
    for org_type in sorted(set(organism_types)):
        org_lengths = [l for l, t in zip(lengths, organism_types) if t == org_type]
        
        print(f"\n{org_type.upper()}:")
        print(f"  Número de sequências: {len(org_lengths)}")
        print(f"  Total de bases: {sum(org_lengths):,}")
        print(f"  Comprimento médio: {np.mean(org_lengths):,.2f} bases")
        print(f"  Comprimento mediano: {np.median(org_lengths):,.2f} bases")
        print(f"  Desvio padrão: {np.std(org_lengths):,.2f} bases")
        print(f"  Menor sequência: {min(org_lengths):,} bases")
        print(f"  Maior sequência: {max(org_lengths):,} bases")

def load_fasta_sequences(fasta_dir):
    sequences = []
    ids = []
    lengths = []
    organism_types = []  # Para guardar o tipo de organismo (bacteria, virus, etc)
    
    print("\nCarregando sequências FASTA...")
    
    # Percorre cada subdiretório
    for organism_type in os.listdir(fasta_dir):
        organism_path = os.path.join(fasta_dir, organism_type)
        
        # Verifica se é um diretório
        if os.path.isdir(organism_path):
            print(f"\nProcessando {organism_type}...")
            
            # Lista todos os arquivos FASTA no subdiretório
            fasta_files = [f for f in os.listdir(organism_path) 
                          if f.endswith(('.fasta', '.fa'))]
            
            for file in tqdm(fasta_files, desc=f"Arquivos de {organism_type}"):
                file_path = os.path.join(organism_path, file)
                try:
                    for record in SeqIO.parse(file_path, 'fasta'):
                        sequences.append(str(record.seq))
                        # Adiciona o tipo de organismo ao ID para melhor visualização
                        modified_id = f"{organism_type}_{record.id}"
                        ids.append(modified_id)
                        lengths.append(len(record.seq))
                        organism_types.append(organism_type)
                        #print(f"Sequência {modified_id}: {len(record.seq)} bases")
                except Exception as e:
                    print(f"Erro ao processar {file_path}: {str(e)}")

    print(f"\nTotal de sequências carregadas: {len(sequences)}")
    for org_type in set(organism_types):
        count = organism_types.count(org_type)
        print(f"{org_type}: {count} sequências")
    
    return sequences, ids, lengths, organism_types

def generate_professional_colors(organism_types):
    # Definir cores fixas para cada tipo de organismo
    color_map = {
        'bacteria': 'rgb(66, 133, 244)',   # Google Blue
        'viruses': 'rgb(219, 68, 55)',     # Google Red
        'archaea': 'rgb(244, 180, 0)',     # Google Yellow
        'eukaryotes': 'rgb(15, 157, 88)'   # Google Green
    }
    
    return [color_map[org_type] for org_type in organism_types]

def plot_3d(reduced_embeddings, ids, lengths, organism_types):
    # Normalize lengths for marker scaling
    min_size = 15
    max_size = 40
    normalized_sizes = np.interp(lengths, (min(lengths), max(lengths)), (min_size, max_size))
    
    # Generate professional color palette based on organism types
    colors = generate_professional_colors(organism_types)
    
    # Criar uma figura para cada tipo de organismo
    traces = []
    
    # Pegar tipos únicos de organismos
    unique_organisms = list(set(organism_types))
    
    for org_type in unique_organisms:
        # Criar máscara para este tipo de organismo
        mask = [t == org_type for t in organism_types]
        
        # Criar scatter plot para este tipo
        traces.append(go.Scatter3d(
            x=reduced_embeddings[mask, 0],
            y=reduced_embeddings[mask, 1],
            z=reduced_embeddings[mask, 2],
            mode='markers+text',
            name=org_type,
            text=[id for id, m in zip(ids, mask) if m],
            hovertext=[f'ID: {id}<br>Type: {org_type}<br>Length: {length} bases' 
                      for id, length, m in zip(ids, lengths, mask) if m],
            hoverinfo='text',
            marker=dict(
                size=[s for s, m in zip(normalized_sizes, mask) if m],
                color=generate_professional_colors([org_type] * sum(mask))[0],
                opacity=0.8,
                symbol='circle',
                line=dict(color='white', width=1)
            )
        ))
    
    fig = go.Figure(data=traces)
    
    fig.update_layout(
        title=dict(
            text='Visualização 3D de Genomas por Tipo de Organismo',
            y=0.95,
            x=0.5,
            xanchor='center',
            yanchor='top',
            font=dict(size=24, color='rgb(50, 50, 50)')
        ),
        scene=dict(
            xaxis_title=dict(text='PC1', font=dict(size=14)),
            yaxis_title=dict(text='PC2', font=dict(size=14)),
            zaxis_title=dict(text='PC3', font=dict(size=14)),
            xaxis=dict(backgroundcolor="rgb(240, 240, 240)", 
                      gridcolor="white", 
                      showbackground=True,
                      gridwidth=2),
            yaxis=dict(backgroundcolor="rgb(240, 240, 240)", 
                      gridcolor="white", 
                      showbackground=True,
                      gridwidth=2),
            zaxis=dict(backgroundcolor="rgb(240, 240, 240)", 
                      gridcolor="white", 
                      showbackground=True,
                      gridwidth=2),
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        paper_bgcolor="rgb(255, 255, 255)",
        plot_bgcolor="rgb(255, 255, 255)",
        font=dict(family="Arial", size=12),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        )
    )
    
    fig.write_html('sequences_3d.html')

def get_embeddings(sequences, model_name='Rostlab/prot_t5_xl_half_uniref50-enc'):
    print("\nGerando embeddings...")
    print(f"Utilizando modelo: {model_name}")
    
    # Inicializar tokenizer com legacy=False para usar o novo comportamento
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False, legacy=False)
    
    print("Carregando modelo...")
    model = T5EncoderModel.from_pretrained(model_name)
    model = model.eval()
    
    embeddings = []
    batch_size = 4
    
    print(f"\nProcessando {len(sequences)} sequências em lotes de {batch_size}...")
    
    for i in tqdm(range(0, len(sequences), batch_size), desc="Processando lotes"):
        batch = sequences[i:i + batch_size]
        
        # Tokenização com padding e truncation explícitos
        ids = tokenizer.batch_encode_plus(
            batch,
            add_special_tokens=True,
            padding=True,
            truncation=True,
            max_length=512,  # Definindo um tamanho máximo para evitar problemas de memória
            return_tensors="pt"  # Retorna tensores PyTorch diretamente
        )
        
        input_ids = ids['input_ids']
        attention_mask = ids['attention_mask']
        
        try:
            with torch.no_grad():
                embedding = model(input_ids=input_ids, attention_mask=attention_mask)
                embedding = embedding.last_hidden_state.mean(dim=1)
                embeddings.append(embedding.numpy())
        except Exception as e:
            print(f"\nErro ao processar lote {i//batch_size + 1}: {str(e)}")
            print(f"Tamanhos das sequências no lote: {[len(seq) for seq in batch]}")
            continue
            
    if not embeddings:
        raise ValueError("Nenhum embedding foi gerado com sucesso!")
        
    final_embeddings = np.vstack(embeddings)
    print(f"\nGerados embeddings com formato: {final_embeddings.shape}")
    
    return final_embeddings

def reduce_dimensions(embeddings, n_components=3, perplexity=30, n_iter=300):
    tsne = TSNE(n_components=n_components, perplexity=perplexity, n_iter=n_iter, random_state=42)
    reduced = tsne.fit_transform(embeddings)
    return reduced

def main(fasta_dir, output_file='embeddings.npy'):
    # Load sequences
    sequences, ids, lengths, organism_types = load_fasta_sequences(fasta_dir)

    print_sequence_statistics(sequences, lengths, organism_types)
    
    # Generate embeddings
    embeddings = get_embeddings(sequences)
    np.save(output_file, embeddings)
    
    # Reduce dimensions and plot
    print("\nReduzindo dimensionalidade...")
    reduced = reduce_dimensions(embeddings, perplexity=30, n_iter=300)
    plot_3d(reduced, ids, lengths, organism_types)
    print("\nVisualização gerada com sucesso!")

if __name__ == '__main__':
    main('bio_data/')