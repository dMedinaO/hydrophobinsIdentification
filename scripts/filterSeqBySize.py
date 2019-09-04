import sys
from Bio import SeqIO

#funcion que permite exportar el contenido a un archivo fasta
def createExportFile(arraySequence, fileExport):

    fileOpen = open(fileExport, 'w')
    index = 0
    for sequence in arraySequence:

        sequenceData = ""
        for i in range(len(sequence)):
            sequenceData+=sequence[i]

        fileOpen.write(">"+str(index)+"\n")
        fileOpen.write(sequenceData)
        fileOpen.write("\n")
        index+=1
    fileOpen.close()

#filtro que permite evaluar la cantidad de cys en la secuencia
def prefilterNumberCys(sequence):

    cont=0
    for i in range(len(sequence)):
        if sequence[i] == 'C':
            cont+=1
    if cont>=8:
        return 0#cumple con la cantidad minima
    else:
        return 1#no cumple con la cantidad minima

#hacemos el filtro con respecto al tamano de secuencia normalmente son 100 pero se consideraran hasta 200 secuencias
seqFilterBySize = []

for record in SeqIO.parse(sys.argv[1], "fasta"):

    if len(record.seq)<=250:
        if prefilterNumberCys(record.seq) == 0:
            seqFilterBySize.append(record.seq)

createExportFile(seqFilterBySize, 'sequenceBySize.fasta')
