import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#funcion que permite exportar el archivo al directorio correspondiente
def exportDataFile(sequencesExport, id_seqExport, outputName):
    listRecord = []

    for i in range (len(id_seqExport)):

        record = SeqRecord(sequencesExport[i], id=id_seqExport[i])

        listRecord.append(record)

    #export file
    with open(outputName, "w") as output_handle:
        SeqIO.write(listRecord, output_handle, "fasta")


#recibimos los argumentos desde terminal
inputNameFile = sys.argv[1]
outputPaht = sys.argv[2]

id_seq = []
seqData = []

for record in SeqIO.parse(inputNameFile, "fasta"):

    id_seq.append(record.id)
    seqData.append(record.seq)

numberFiles = len(id_seq)/1000
restFiles = len(id_seq)%1000

index=0


for i in range(numberFiles):#creamos los primeros elementos

    sequencesExport = []
    id_seqExport = []

    for j in range(index, index+1000):
        sequencesExport.append(seqData[j])
        id_seqExport.append(id_seq[j])

    nameFile = "%sexport%d.fasta" % (outputPaht, i+1)
    exportDataFile(sequencesExport, id_seqExport, nameFile)

    index=index+1000

#los pendientes
sequencesExport = []
id_seqExport = []
for i in range(index, len(id_seq)):#creamos los primeros elementos

    sequencesExport.append(seqData[i])
    id_seqExport.append(id_seq[i])

nameFile = "%sexport_pending.fasta" % (outputPaht)
exportDataFile(sequencesExport, id_seqExport, nameFile)
