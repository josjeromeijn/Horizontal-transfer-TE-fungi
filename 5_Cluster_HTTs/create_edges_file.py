#written by: Josje Romeijn Nov '24

#script to make a list of unique nodes from multiple edges files
import argparse, os, csv, re

#1------------------------GET INPUT FROM COMMAND--------------------------------
parser = argparse.ArgumentParser(description="Creates list of unique nodes from multiple files containing edges")
parser.add_argument('-i','--input',help='Directory of files candidate HTTs',required=True)
parser.add_argument('-o','--output',help='Output file (node file name)', required=True)
parser.add_argument('-c','--classifications',help='Directory of files with classifications')
parser.add_argument('-cc','--classifications_curlib',help='File with classifications of curated library Repeatmasker')

args=parser.parse_args()

#2----------------------------DEFINE FUNCTIONS----------------------------------
def load_classifications(dir):
    class_data = {}
    for file_name in os.listdir(dir):
        genome = file_name.replace('.txt','')
        genome = genome.replace('pastec_earlgray_annotation_','')
        with open(os.path.join(dir, file_name), 'r') as file:
            next(file)
            for line in file:
                id, _, _, _, _, _, order, superfam, fam, _ = line.split("\t")
                id = re.match(r'(rnd-\d+_family-\d+)', id).group(1)
                class_data[f"{genome}__{id.upper()}"] = f"{order}__{superfam}__{fam}"
    return class_data

def load_class_curlib(file, class_data):
    with open(file, 'r') as f:
        for line in f:
            id, _, _, _, _, _, order, superfam, fam, _ = line.split("\t")
            class_data[id.replace('>', '').upper()] = f"{order}__{superfam}__{fam}"
    return class_data

def load_can_HTT_data(file_name):
    te_pairs = {}
    with open(file_name, 'r') as file:
        next(file)
        for line in file:
            q, s, pid = line.strip().split("\t")
            te_pairs[f"{q}___{s}"] = pid
    return te_pairs

#3--------------------LOAD CLASSIFICATIONS AND EDGES----------------------------
#laod classifications
class_data = load_classifications(args.classifications)
print("Done loading classifications")

#load classifications of curated library of repeatmasker
class_data = load_class_curlib(args.classifications_curlib, class_data)

#load edges and connect with classification
te_pairs = load_can_HTT_data(args.input)
print("Done loading edges data")

#4---------------------WRITE DATA TO OUTPUT FILE--------------------------------

with open(args.output, 'w') as out:
    out.write("Node\tpid\tclass_o_s_f\tagreement_o_s_f\tclass_o_s\tagreement_o_s\tclass_o\tagreement_o\n")
    for key in te_pairs:
        #get classification of both TEs
        te1, te2 = key.split("___")
        pid = te_pairs.get(key)
        #classification on order, superfamily and family level
        classif1 = class_data.get(te1.split("::")[0].rsplit("_",1)[0]) or \
            class_data.get(te1.split("::")[0].rsplit("_", 1)[0].split("__", 1)[1]) #TE is from curated library
        classif2 = class_data.get(te2.split("::")[0].rsplit("_",1)[0]) or \
            class_data.get(te2.split("::")[0].rsplit("_", 1)[0].split("__", 1)[1]) # TE is from curated library
        code1 = 'same' if classif1 == classif2 else 'diff'

        if not classif1:
            classif1 = "None"

        if not classif2:
            classif2 = "None"

        #classification on order, superfamily level
        classif11 = classif1.rsplit("__",1)[0]
        classif22 = classif2.rsplit("__",1)[0]
        code11 = 'same' if classif11 == classif22 else 'diff'

        #classification on order level
        classif111 = classif11.rsplit("__",1)[0]
        classif222 = classif22.rsplit("__",1)[0]
        code111 = 'same' if classif111 == classif222 else 'diff'
        out.write(f"{key}\t{pid}\t{classif1}___{classif2}\t{code1}\t{classif11}___{classif22}\t{code11}\t{classif111}___{classif222}\t{code111}\n")

print("Done writing nodes file")
