import tkinter as tk
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.MeltingTemp import Tm_Wallace

#FONTE
import pyglet,tkinter
pyglet.font.add_file(r"C:\Users\joaov\Documents\Biopython Resources\MYRIADPRO-REGULAR.ttf")

#PROGRAMA
Entrez.email = "your_email@example.com"

def fetch_sequence_info(accession, email):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        sequence = record.seq

        seq_length = len(sequence)
        mol_weight = ProteinAnalysis(str(sequence)).molecular_weight()

        gene_name = None
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                break

        nucleotides_counts = {aa: sequence.count(aa) for aa in set(sequence)}
        total_aa = sum(nucleotides_counts.values())
        nucleotides_percent = {aa: count / total_aa * 100 for aa, count in nucleotides_counts.items()}

        result_text = f"Identificação do gene: {gene_name}\nTamanho da sequência: {seq_length}\nPeso molecular (g/mol): {mol_weight}\nPorcentagem de bases: {nucleotides_percent}"
        result_label.config(text=result_text)

    except Exception as e:
        result_label.config(text=f"Error: {str(e)}")


def design_primers(sequence, primer_length):
    primer_forward = sequence[:primer_length]
    primer_reverse = sequence[-primer_length:].reverse_complement()

    tm_forward = Tm_Wallace(primer_forward)
    tm_reverse = Tm_Wallace(primer_reverse)

    return primer_forward, primer_reverse, tm_forward, tm_reverse

def fetch_and_design_primers(accession, email, primer_length):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        sequence = record.seq

        primer_forward, primer_reverse, tm_forward, tm_reverse = design_primers(sequence, primer_length)

        result_text = f"Forward Primer: {primer_forward}\nReverse Primer: {primer_reverse}\nPonto de fusão do forward: {tm_forward}°C\nPonto de fusão do reverse: {tm_reverse}°C"
        result_label.config(text=result_text)

    except Exception as e:
        result_label.config(text=f"Error: {str(e)}")

def on_primer_length_change(entry):
    try:
        primer_length = int(entry.get())
        if primer_length <= 0:
            raise ValueError("O tamanho do primer deve ser maior que zero!")
        primer_length_var.set(primer_length)
    except ValueError:
        entry.delete(0, tk.END)
        entry.insert(0, primer_length_var.get())

#GUI VISUAL
app = tk.Tk()
app_icon = tk.PhotoImage(file=r"C:\Users\joaov\Documents\Biopython Resources\dna_icon.ico")
app.iconphoto(False, app_icon)
app.title("Sequenciamento e montagem de primers por Accession do GenBank")
app.configure(bg='lightgreen')
app.option_add("*Font", "MyriadPro-Regular 14 bold")

#GUI INTERFACE
accession_label = tk.Label(app, text="Accession:", bg='#6DB192')
accession_entry = tk.Entry(app)
email_label = tk.Label(app, text="E-mail", bg='#6DB192')
email_entry = tk.Entry(app)
primer_length_label = tk.Label(app, text="Tamanho do primer em nucleotídeos", bg='#6DB192')
primer_length_var = tk.IntVar(value=20)
primer_length_entry = tk.Entry(app, textvariable=primer_length_var)
primer_length_entry.bind("<FocusOut>", lambda event: on_primer_length_change(primer_length_entry))
fetch_button = tk.Button(app, text="Analisar", command=lambda: fetch_sequence_info(accession_entry.get(), email_entry.get()), bg='#ADD8E6')
primer_button = tk.Button(app, text="Montar primers", command=lambda: fetch_and_design_primers(accession_entry.get(), email_entry.get(), primer_length_var.get()), bg='#ADD8E6')
result_label = tk.Label(app, text="Resultados aparecerão aqui! ^^", bg='#6DB192')

#GUI LAYOUT
accession_label.pack(pady=10)
accession_entry.pack(pady=5)
email_label.pack(pady=5)
email_entry.pack(pady=5)
primer_length_label.pack(pady=5)
primer_length_entry.pack(pady=5)
fetch_button.pack(pady=5)
primer_button.pack(pady=5)
result_label.pack(pady=10)

#RODAR
app.mainloop()
