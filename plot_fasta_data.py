import numpy as np
from Bio import SeqIO
from transformers import T5EncoderModel, T5Tokenizer
import torch
from sklearn.manifold import TSNE
import plotly.graph_objects as go
import os
import random

def load_fasta_sequences(fasta_dir):
    sequences = []
    ids = []
    for file in os.listdir(fasta_dir):
        if file.endswith('.fasta') or file.endswith('.fa'):
            for record in SeqIO.parse(os.path.join(fasta_dir, file), 'fasta'):
                sequences.append(str(record.seq))
                ids.append(record.id)
    return sequences, ids

def get_embeddings(sequences, model_name='Rostlab/prot_t5_xl_half_uniref50-enc'):
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_name)
    model = model.eval()
    
    embeddings = []
    batch_size = 4
    
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        ids = tokenizer.batch_encode_plus(batch, add_special_tokens=True, padding=True)
        input_ids = torch.tensor(ids['input_ids'])
        attention_mask = torch.tensor(ids['attention_mask'])
        
        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask)
            embedding = embedding.last_hidden_state.mean(dim=1)
            embeddings.append(embedding.numpy())
            
    return np.vstack(embeddings)

def reduce_dimensions(embeddings, n_components=3, perplexity=30, n_iter=300):
    tsne = TSNE(n_components=n_components, perplexity=perplexity, n_iter=n_iter, random_state=42)
    reduced = tsne.fit_transform(embeddings)
    return reduced

def plot_3d(reduced_embeddings, ids):
    # Gerar uma cor única para cada ponto
    colors = [f'rgb({random.randint(0, 255)}, {random.randint(0, 255)}, {random.randint(0, 255)})' for _ in ids]
    
    fig = go.Figure(data=[go.Scatter3d(
        x=reduced_embeddings[:, 0],
        y=reduced_embeddings[:, 1],
        z=reduced_embeddings[:, 2],
        mode='markers+text',
        text=ids,
        hoverinfo='text',
        marker=dict(
            size=5,  # Ajustar o tamanho dos pontos
            color=colors,  # Atribuir cores únicas
            opacity=0.8
        )
    )])
    
    fig.update_layout(
        scene=dict(
            xaxis_title='PC1',
            yaxis_title='PC2',
            zaxis_title='PC3',
            xaxis=dict(backgroundcolor="rgb(200, 200, 200)", gridcolor="white", showbackground=True),
            yaxis=dict(backgroundcolor="rgb(200, 200, 200)", gridcolor="white", showbackground=True),
            zaxis=dict(backgroundcolor="rgb(200, 200, 200)", gridcolor="white", showbackground=True)
        ),
        margin=dict(l=0, r=0, b=0, t=0),  # Ajustar margens
        paper_bgcolor="rgb(240, 240, 240)",  # Cor de fundo do gráfico
        font=dict(size=12)  # Tamanho da fonte
    )
    
    fig.write_html('sequences_3d.html')

def main(fasta_dir, output_file='embeddings.npy'):
    # Load sequences
    sequences, ids = load_fasta_sequences(fasta_dir)
    
    # Generate embeddings
    embeddings = get_embeddings(sequences)
    np.save(output_file, embeddings)
    
    # Reduce dimensions and plot
    reduced = reduce_dimensions(embeddings, perplexity=30, n_iter=300)
    plot_3d(reduced, ids)

if __name__ == '__main__':
    main('bio_data/')