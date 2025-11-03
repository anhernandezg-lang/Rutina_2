import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import statistics

##Luego procedo a definir dos variables para mis archivos
fasta = "bartonella_dataset/ncbi_dataset/data/GCF_019930925.1/GCF_019930925.1_ASM1993092v1_genomic.fna"
gff ="bartonella_dataset/ncbi_dataset/data/GCF_019930925.1/genomic.gff"

# Luego se procede a leer el genoma completo ponemos el condicional if not para eliminar ">" del encabezado
genome = ""
with open(fasta) as file:
    for line in file:
        if not line.startswith(">"):
            genome += line.strip().upper()
genes = []
with open(gff) as file:
    for line in file:
        if not line.startswith("#"):
            parts = line.strip().split("\t")
            if len(parts) > 8 and parts[2] == "gene":
                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                info = parts[8]
                genes.append([seqid, start, end, strand, info])
# Aquí transformo la lista de genes en un DataFrame y calculo la longitud de cada uno
df = pd.DataFrame(genes, columns=["seqid", "start", "end", "strand", "info"])
df["length"] = df["end"] - df["start"]
# Aqui utilizó el codigo de GC content realizado en clase
def GC_content(sequence):
    count_A = sequence.count('A')
    count_C = sequence.count('C')
    count_G = sequence.count('G')
    count_T = sequence.count('T')
    if len(sequence) == 0:
        return 0
    GC_content = (count_G + count_C) / len(sequence)
    GC_content = round(100 * GC_content, 2)
    return GC_content

gc_porcentajes = []

for i, fila in df.iterrows():
    secuencia = genome[fila["start"]:fila["end"]]
    gc_value = GC_content(secuencia)
    gc_porcentajes.append(gc_value)

df["GC_percent"] = gc_porcentajes

## Funcion de significancia con p-value tarea puesta en clase para incluir.
def analisis_significancia(df):
    media = statistics.mean(df["length"]) ##se calcula el promedio de las longitudes
    desvio = statistics.stdev(df["length"]) ## calcula de desviacion estandar
    valor_p = 0.05 ## valor de referencia equivalente al 95% confianza
    resultados = [] # aqui voy  a guardar los resultados
    for valor in df["length"]:
        if valor > media + 2 * desvio or valor < media - 2 * desvio: ## Si el valor está a más de 2 desviaciones de la media → significativo

            resultados.append("Significativo")
        else:
            resultados.append("No significativo") #si no se cumple la condicion anterior entonces --> No significativo

    df["significativo"] = resultados
    df["media_longitud"] = media
    df["desvio"] = desvio
    df["valor_p"] = valor_p
    return df
df = analisis_significancia(df)
### GRAFICAS ###

sns.set_theme(style="whitegrid")
# Histograma del %GC
plt.figure()
sns.histplot(df["GC_percent"], kde=True, color="skyblue")
plt.title("Distribución del contenido GC por gen")
plt.xlabel("%GC")
plt.ylabel("Frecuencia")
plt.show()

# 2 Histograma de longitudes
plt.figure()
sns.histplot(df["length"], kde=True, color="lightgreen")
plt.title("Distribución de longitudes génicas")
plt.xlabel("Longitud (bp)")
plt.ylabel("Frecuencia")
plt.show()

# 3 Boxplot de longitud según el rango de GC
df["rango_GC"] = pd.cut(df["GC_percent"], bins=10)
plt.figure(figsize=(9,4))
sns.boxplot(data=df, x="rango_GC", y="length", color="salmon")
plt.xticks(rotation=45)
plt.title("Longitud por rango de %GC")
plt.xlabel("Rango de GC (%)")
plt.ylabel("Longitud (bp)")
plt.show()

# 4 Violinplot del %GC
plt.figure()
sns.violinplot(y="GC_percent", data=df, color="orchid")
plt.title("Distribución del %GC por gen")
plt.ylabel("%GC")
plt.show()

# 5️ Scatterplot longitud vs GC con significancia
plt.figure()
sns.scatterplot(data=df, x="GC_percent", y="length",
                hue="significativo", palette="Set1", alpha=0.7)
plt.title("Longitud vs GC% (genes significativos)")
plt.xlabel("%GC")
plt.ylabel("Longitud (bp)")
plt.show()

# 6 Gráfico resumen: conteo de genes significativos
plt.figure(figsize=(8,5))
sns.countplot(data=df, x="rango_GC", hue="significativo", palette="coolwarm")
plt.title("Cantidad de genes significativos por rango de GC")
plt.xlabel("Rango de %GC")
plt.ylabel("Cantidad de genes")
plt.xticks(rotation=45)
plt.legend(title="Significancia")
plt.show()

media = round(statistics.mean(df["length"]), 2)
desvio = round(statistics.stdev(df["length"]), 2)
num_sig = df["significativo"].value_counts().get("Significativo", 0)
total = len(df)
porcentaje = round(num_sig * 100 / total, 2)

print("\n=== RESUMEN DEL ANÁLISIS ===")
print(f"Media de longitud: {media} bp")
print(f"Desviación estándar: {desvio} bp")
print(f"Genes totales: {total}")
print(f"Genes significativos: {num_sig} ({porcentaje}%)")

