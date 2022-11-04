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
