import sys

fileOpen = open(sys.argv[1], 'r')

lineData = []

line = fileOpen.readline()

while line:
    lineData.append(line.replace("\n", ""))
    line = fileOpen.readline()
fileOpen.close()

i=0
proteinsData = []

while i <len(lineData):

    if 'protein sequence' in lineData[i]:#encuentro donde inicia una proteina
        protein = []
        protein.append(lineData[i].replace("#", "").replace("protein sequence = [", "").replace(" ", ""))
        #desde esta posicion busco el final de la data de la protein
        i+=1
        for j in range(i, len(lineData)):
            i+=1
            protein.append(lineData[j].replace("#", "").replace(" ", "").replace("]", ""))
            if ']' in lineData[j]:
                break

        proteinValue = ""

        for element in protein:
            proteinValue+=element
        proteinsData.append(proteinValue)
    else:
        i+=1

#exportamos el doc a fasta
fileExport = open(sys.argv[2], 'w')

index=1
for element in proteinsData:
    line = ">Sequence: "+str(index)
    index+=1

    fileExport.write(line+"\n"+element+"\n")
fileExport.close()
