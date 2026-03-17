import streamlit as st
import requests
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import random
import collections

st.title("🧬 DNA CSI LAB")

st.write("Laboratorio genético educativo")

sequence = st.text_area("Introduce una secuencia ADN", "ATGGATTTATCTGCTCTTCGATCGATCGATCG")


# =============================
# DETECCIÓN SNP
# =============================

snp_markers = {
    "rs12913832": ("eyes","AA","AG","GG"),
    "rs16891982": ("hair","GG","CG","CC"),
    "rs1426654": ("skin","AA","AG","GG"),
    "rs6152": ("bald","AA","AG","GG")
}

def detect_snps(seq):

    detected = {}

    for snp in snp_markers:
        detected[snp] = random.choice(["AA","AG","GG"])

    return detected


# =============================
# PREDICCIÓN RASGOS
# =============================

def predict_traits(snps):

    traits={}

    # ojos
    if snps["rs12913832"]=="AA":
        traits["eyes"]="blue"
    elif snps["rs12913832"]=="AG":
        traits["eyes"]="green"
    else:
        traits["eyes"]="brown"

    # piel
    if snps["rs1426654"]=="AA":
        traits["skin"]="light"
    elif snps["rs1426654"]=="AG":
        traits["skin"]="medium"
    else:
        traits["skin"]="dark"

    # pelo
    if snps["rs16891982"]=="GG":
        traits["hair"]="blonde"
    elif snps["rs16891982"]=="CG":
        traits["hair"]="brown"
    else:
        traits["hair"]="black"

    # calvicie
    traits["bald"]= snps["rs6152"]=="GG"

    return traits


# =============================
# AVATAR
# =============================

def draw_avatar(traits):

    img=Image.new("RGB",(400,400),(255,255,255))
    d=ImageDraw.Draw(img)

    skin_colors={
        "light":(255,224,189),
        "medium":(205,133,63),
        "dark":(92,64,51)
    }

    d.ellipse((100,50,300,300),fill=skin_colors[traits["skin"]])

    eye_colors={
        "blue":"blue",
        "green":"green",
        "brown":"brown"
    }

    d.ellipse((150,150,170,170),fill=eye_colors[traits["eyes"]])
    d.ellipse((230,150,250,170),fill=eye_colors[traits["eyes"]])

    d.line((200,170,200,210),fill="black",width=3)
    d.arc((170,220,230,260),0,180,fill="red",width=3)

    if not traits["bald"]:

        hair_colors={
            "blonde":"yellow",
            "brown":"brown",
            "black":"black"
        }

        d.rectangle((120,40,280,90),fill=hair_colors[traits["hair"]])

    return img


# =============================
# ANCESTRÍA SIMULADA
# =============================

def ancestry_simulation():

    europe=random.randint(40,80)
    asia=random.randint(0,30)
    africa=random.randint(0,30)

    total=europe+asia+africa

    return {
        "Europe":round(europe/total*100),
        "Asia":round(asia/total*100),
        "Africa":round(africa/total*100)
    }


def draw_map(data):

    labels=list(data.keys())
    values=list(data.values())

    fig,ax=plt.subplots()
    ax.pie(values,labels=labels,autopct='%1.0f%%')
    ax.set_title("Ancestría estimada")

    return fig


# =============================
# CONSULTA API ENSEMBL
# =============================

def query_snp(rsid):

    url=f"https://rest.ensembl.org/variation/human/{rsid}"
    headers={"Content-Type":"application/json"}

    r=requests.get(url,headers=headers)

    if r.status_code==200:
        return r.json()

    return None


# =============================
# ESTADÍSTICAS ADN
# =============================

def dna_statistics(seq):

    seq=seq.upper()
    length=len(seq)

    counts=collections.Counter(seq)

    gc=(counts["G"]+counts["C"])/length*100

    return {
        "length":length,
        "A":counts["A"],
        "T":counts["T"],
        "G":counts["G"],
        "C":counts["C"],
        "GC_percent":round(gc,2)
    }


# =============================
# TRADUCCIÓN ADN → PROTEÍNA
# =============================

genetic_code={
'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
'TGC':'C','TGT':'C','TGA':'*','TGG':'W'
}

def translate_dna(seq):

    seq=seq.upper()
    protein=""

    for i in range(0,len(seq)-2,3):

        codon=seq[i:i+3]
        aa=genetic_code.get(codon,"?")
        protein+=aa

    return protein


# =============================
# MOTIVOS GENÉTICOS
# =============================

def detect_motifs(seq):

    motifs={
        "TATA_box":"TATAAA",
        "start_codon":"ATG",
        "stop_codon":["TAA","TAG","TGA"]
    }

    found={}

    for m in motifs:

        if isinstance(motifs[m],list):

            for codon in motifs[m]:
                if codon in seq:
                    found[m]=codon

        else:
            if motifs[m] in seq:
                found[m]=motifs[m]

    return found


# =============================
# BOTÓN ANALIZAR
# =============================

if st.button("Analizar ADN"):

    snps=detect_snps(sequence)

    st.subheader("SNP detectados")
    st.write(snps)

    traits=predict_traits(snps)

    st.subheader("Rasgos estimados")
    st.write(traits)

    avatar=draw_avatar(traits)

    st.subheader("Avatar generado")
    st.image(avatar)

    ancestry=ancestry_simulation()

    st.subheader("Origen poblacional estimado")
    fig=draw_map(ancestry)
    st.pyplot(fig)

    st.subheader("Consulta SNP en Ensembl")

    for snp in snps:

        data=query_snp(snp)

        if data:
            st.write(snp,"→",data["mappings"][0]["location"])

    stats=dna_statistics(sequence)

    st.subheader("Estadísticas ADN")
    st.write(stats)

    protein=translate_dna(sequence)

    st.subheader("Traducción proteína")
    st.write(protein)

    motifs=detect_motifs(sequence)

    st.subheader("Motivos genéticos detectados")
    st.write(motifs)
