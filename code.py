from mp_api.client import MPRester
import math
import csv
import networkx as nx
import random

def descargar_datos():
    
    with MPRester("BYA3Yb9XkXhGb6uOHdEGzKlsTi5Jb4zX") as mpr:
        diccionario, atomos = {},[]
        idt=str(random.randint(1, 20000))
        red_crist = mpr.materials.search(material_ids="mp-"+idt,fields=["structure"])[0]
        nom_atom = red_crist.structure.species 
        pos_atom = red_crist.structure.cart_coords
        print(red_crist.structure)
        
        
        for i in range(len(pos_atom)):
            atomos.append(str(nom_atom[i]))
            c = 0
            while True:
                indice = str(nom_atom[i])+str(c)
                if indice in diccionario:
                    c += 1
                    pass
                else:
                    diccionario[indice] = list(pos_atom[i])
                    break
                
        cel_coords = red_crist.structure.lattice.abc
        
        print(diccionario)
        print(atomos)
        return diccionario, atomos, cel_coords, idt

def celdas_circundantes(molecula, identificadores, cel_coords ):
    circun_cel = []
    for x in range(-1,2):
        circun_cel.append([])
        for y in range(-1,2):
            circun_cel[x+1].append([])
            for z in range(-1,2):
                circun_cel[x+1][y+1].append({})
                for identificador in identificadores:
                    print(identificador)
                    circun_cel[x+1][y+1][z+1][identificador] = (molecula[identificador][0]+(cel_coords[0]*x),
                                                                molecula[identificador][1]+(cel_coords[1]*y),
                                                                molecula[identificador][2]+(cel_coords[2]*z))
     
    return circun_cel

def generador_matrices_adyacencia(atomos,arreglo):
    # Abre una base de datos y compara ambos átomos junto con su longitud para comprobar si hay un enlace
    with open(r"C:\Users\Krozz\Code\Quimica-Matriz\Longitudes de enlace.csv","r") as archivo:
        Tabla_lista=[]
        Matriz=[]
        Tabla = csv.DictReader(archivo)
        for datos_enlace in Tabla:
            Tabla_lista.append(datos_enlace)
        for atomo_x in range(len(atomos)):
            Matriz.append([])
            for atomo_y in range(len(atomos)):
                for enlace_actual in Tabla_lista:
                    if (atomos[atomo_x]+"-"+atomos[atomo_y])==enlace_actual["Enlace"] or (atomos[atomo_y]+"-"+atomos[atomo_x])==enlace_actual["Enlace"]:
                        if math.isclose(arreglo[atomo_x][atomo_y],float(enlace_actual["Media"]), rel_tol = float(enlace_actual["Desviacion"])) or arreglo[atomo_x][atomo_y] == 0:
                            Matriz[atomo_x].append(1)
                try:
                    print(Matriz[atomo_x][atomo_y])
                except:
                    Matriz[atomo_x].append(0)
    
    
    with open(r"C:/Users/Krozz/Code/Quimica-Matriz/Radios covalentes.csv") as archivo:
        Matriz_covalente = []
        Datos_dicc = {}
        for datos_covalente in csv.DictReader(archivo):
            Datos_dicc[datos_covalente["Simbolo"]] = datos_covalente["Radio covalente"]
        for atomo_x in range(len(atomos)):
            Matriz_covalente.append([])
            for atomo_y in range(len(atomos)):
                try:
                    if (float(Datos_dicc[atomos[atomo_x]]) + float(Datos_dicc[atomos[atomo_y]])) > arreglo[atomo_x][atomo_y]:
                        Matriz_covalente[atomo_x].append(1)
                    else:
                        Matriz_covalente[atomo_x].append(0)
                except:
                    pass
                
    for x in range(len(atomos)):
        for y in range(len(atomos)):
            if Matriz[x][y] == 0 and Matriz_covalente[x][y] == 1:
                Matriz[x][y] = 1
    return Matriz

def enlaces_circundantes(molecula,circun_cel,atomos):
    Base = []
    for a in range(len(molecula)):
        Base.append([])
        for b in range(len(molecula)):
            Base[a].append(0)

    for x in range(3):
        for y in range(3):
            for z in range(3):
                arreglo = []
                contador = -1
                # Calcula la distancia por medio de un vector tridimensional y genera un arreglo.
                for atomo in molecula:
                    arreglo.append([])
                    contador +=1
                    for interaccion in molecula:
                        valor = 0
                        for eje in range(3):
                            valor += pow(circun_cel[x][y][z][interaccion][eje]-molecula[atomo][eje],2)
                        arreglo[contador].append(round(math.sqrt(valor),4))
                if x==y==z:
                    Matriz = generador_matrices_adyacencia(atomos, arreglo)
                    for a in range(len(atomos)):
                        for b in range(len(atomos)):
                            if Base[a][b] == 0 and Matriz[a][b] == 1:
                                Base[a][b] = 1
    return Base

def generador_matrices(molecula):
    
    arreglo = []
    contador = -1        
    for atomo in molecula:
        arreglo.append([])
        contador +=1
        for interaccion in molecula:
            valor = 0
            for eje in range(3):
                valor += pow(molecula[interaccion][eje]-molecula[atomo][eje],2)
            arreglo[contador].append(round(math.sqrt(valor),4))
    return arreglo



       
   
    

molecula, atomos, cel_coords, idt = descargar_datos()
identificadores = []
for i in molecula:
    identificadores.append(i)
circun_cel=celdas_circundantes(molecula,identificadores,cel_coords)
Matriz = enlaces_circundantes(molecula,circun_cel,atomos)


print("--------------------------------------------------------")
contador = 0
for y in Matriz:
    print("%s:%s" % (identificadores[contador],y))
    contador += 1
print(identificadores)

                    
                    
   
G = nx.Graph()
for identificador in identificadores:
    G.add_node(identificador)
for x in range(len(Matriz)):
    for y in range(len(Matriz[x])):
        if Matriz[x][y] == 1 and (not x == y):
            G.add_edge(identificadores[x],identificadores[y])
nx.draw(G,with_labels=(True))
print(idt)
