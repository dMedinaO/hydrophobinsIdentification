from Bio import SeqIO
import sys

#funcion que permite generar un archivo con las secuencias existentes en un arreglo
def createExportFile(arraySequence, idSequences, descriptionSequence, fileExport):

    fileOpen = open(fileExport, 'w')
    for i in range(len(arraySequence)):

        sequenceData = ""
        for j in range(len(arraySequence[i])):
            sequenceData+=arraySequence[i][j]

        header = idSequences[i]+" "+ descriptionSequence[i]
        fileOpen.write(">"+header+"\n")
        fileOpen.write(sequenceData)
        fileOpen.write("\n")
    fileOpen.close()

#recibimos las secuencias
db1 = '/home/dmedina/Escritorio/projects/hydrophobinsIdentification/dataBase/dataBase.fasta'

sequenceDB1 = []
idSequence1 = []
description1 = []

#hacemos la lectura del archivo de secuencias db1
for record in SeqIO.parse(db1, "fasta"):

    if record.seq not in sequenceDB1:
        sequenceDB1.append(record.seq)
        idSequence1.append(record.id)
        description1.append(record.description)

createExportFile(sequenceDB1, idSequence1, description1, 'sequenceDump.fasta')
