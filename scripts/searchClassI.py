import sys
from Bio import SeqIO
import json

#hacemos la busqueda de los dobles enlaces
def searchDobleEnlaces(sequence):

    cont=0
    for i in range(len(sequence)-1):
        if sequence[i] == 'C' and sequence[i+1] == 'C':
            cont+=1
    if cont>=2:
        return 0#cumple con la cantidad minima
    else:
        return 1#no cumple con la cantidad minima

#obtener posiciones de cys que no sean dobles enlaces
def getPosicionesCys(sequence):

    posiciones = []

    for i in range(len(sequence)):

        try:
            if sequence[i] == 'C' and sequence[i+1] != 'C':
                posiciones.append(i)
        except:
            pass
    return posiciones

#obtener posiciones de los dobles enlaces de cys
def getPosicionesDobleCys(sequence):

    posiciones = []

    for i in range(len(sequence)-1):

        if sequence[i] == 'C' and sequence[i+1] == 'C':
            posiciones.append(i)
    return posiciones


#obtener numero de cys
def getNumberCys(sequence):
    cont=0
    for i in range(len(sequence)):
        if sequence[i] == 'C':
            cont+=1
    if cont>=8:
        return 0#cumple con la cantidad minima
    else:
        return 1#no cumple con la cantidad minima

#funcion que permite obtener las distancias asociadas
def getDistanceResidues(sequence, idSequence):
    posicionesCys = getPosicionesCys(sequence)
    posicionesDobleCys = getPosicionesDobleCys(sequence)

    #formamos un arreglo con la distancia entre los elementos de cysteina
    distanceCys = {}
    distanceDobleCys = {}

    #obtenemos la diferencia entre las cisteinas y los dobles enlaces de cys
    for element in posicionesDobleCys:
        arrayDiff = []
        for values in posicionesCys:
            dif = element-values
            arrayDiff.append(dif)
        distanceDobleCys.update({str(element): arrayDiff})

    #obtenemos la diferencia entre los residuos de cisteinas
    for element in posicionesCys:
        arrayDiff = []
        for values in posicionesCys:
            if element != values:
                dif = element-values
                arrayDiff.append(dif)
        distanceCys.update({str(element):arrayDiff})

    #formamos la estructura de datos para generar la salida
    arrayResponse = {}
    arrayResponse.update({'IDSeq': idSequence})
    arrayResponse.update({'diffDobleCys':distanceDobleCys})
    arrayResponse.update({'diffCys':distanceCys})

    nameFile = "responseSeq"+idSequence+".json"
    with open(nameFile, 'w') as fp:
        json.dump(arrayResponse, fp)

#hacemos la lectura de datos
for record in SeqIO.parse(sys.argv[1], "fasta"):
    if searchDobleEnlaces(record.seq) == 0:#secuencia presenta minimo dos dobles enlaces
        if getNumberCys(record.seq) == 0:#secuencia presenta la minima cantidad de residuos
            getDistanceResidues(record.seq, record.id)
            
