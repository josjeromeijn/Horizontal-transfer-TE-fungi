import re, copy, subprocess, os, argparse, shutil
from itertools import combinations

parser = argparse.ArgumentParser(description="Calculate Ks for candidate TE-TE blastn hits")
parser.add_argument('-b','--blast',help='Blast output file',required=True)
parser.add_argument('-hm','--hmm_loc', help='file location of hmm output files', required=True)
parser.add_argument('-r','--rps_loc', help='file location of rps output files', required=True)
parser.add_argument('-t','--translated_te', help='file location of fasta files containing translated TEs', required=True)
parser.add_argument('-f','--fasta', help='file location of fasta files containing TEs', required=True)
parser.add_argument('-o','--output', help='output directory (script creates subdirectories with intermediate files in here)', required=True)
parser.add_argument('-k','--output_file_ks', help='output file with KaKs for TE-TE hits', required=True)


args = parser.parse_args()

#file locations
prot_domain_loc_hmm = args.hmm_loc
prot_domain_loc_rps = args.rps_loc
te_translated_orf_loc = args.translated_te
fasta_loc = args.fasta
output_dir = args.output
output_file_ks = args.output_file_ks

#functions
def translate_protein_to_dna(data_hmm):
    for item in set(data_hmm["query_name"]):
        # get indices of this query_name
        index = [i for i, n in enumerate(data_hmm["query_name"]) if n == item]
        # get transposon ID and ORF (number 1 through 6)
        tpid, orf = item.rsplit("_", 1)
        orf = int(orf)  # convert string to integer

        # for each index, convert coordinates
        for idx in index:
            #get length of TE
            end_coord = item.split(":")[-1].split("_")[0].split("-")
            end_coord = int(end_coord[1]) - int(end_coord[0])
            if orf in [1, 2, 3]:  # all the forward ORFs
                #print(f"{data_hmm["accession"][idx]}, start = {data_hmm["from_acc"][idx]}, end = {data_hmm["to_acc"][idx]}")

                data_hmm["from_acc"][idx] = 3 * data_hmm["from_acc"][idx] - 3 + orf
                data_hmm["to_acc"][idx] = 3 * data_hmm["to_acc"][idx] - 1 + orf
                #print(f"{data_hmm["accession"][idx]}, start = {data_hmm["from_acc"][idx]}, end = {data_hmm["to_acc"][idx]}")

                #if end coordinates exceed the length of the TE, remove last AA
                if data_hmm["to_acc"][idx] > end_coord:
                    data_hmm["to_acc"][idx] = data_hmm["to_acc"][idx] - 3
                    data_hmm["to_hmm"][idx] = int(data_hmm["to_hmm"][idx]) - 1

            elif orf in [4, 5, 6]:  # all the reverse ORFs

                l = end_coord % 3
                DIF = [0,1,2]
                dif = DIF[l - orf + 4]

                #print(f"{data_hmm["accession"][idx]}, orf = {orf}, l = {l}, DIF = {dif}, start = {data_hmm["from_acc"][idx]}, end = {data_hmm["to_acc"][idx]}")
                #print(f"{end_coord} - {3} * {data_hmm["from_acc"][idx]} - {dif} + 3 = {end_coord - 3 * data_hmm["from_acc"][idx]} - dif + 3 = {end_coord - 3 * data_hmm["from_acc"][idx] - dif + 3}")
                #calculate ORF based on DIF and l (ORF determination depends on total length of the TE, not where the protein starts compared to reverse compliment unfortunately...
                data_hmm["from_acc"][idx] = end_coord - 3 * data_hmm["from_acc"][idx] - dif + 3
                data_hmm["to_acc"][idx] = end_coord - 3 * data_hmm["to_acc"][idx] - dif + 1

                #if the coordinates start before the TE starts (e.g. start coordinate -1 or 0), start from 2nd AA of hit
                if data_hmm["to_acc"][idx] < 1:
                    data_hmm["to_acc"][idx] = data_hmm["to_acc"][idx] + 3
                    data_hmm["to_hmm"][idx] = int(data_hmm["to_hmm"][idx]) + 1

def get_overlapping_prot_domains(qry, sbj):
    #create dictionaries for query and subject hmm data
    data_hmm = {"target_name": [], "accession": [], "tlen": [], "query_name": [], "qlen": [], "e_value": [], "score": [], "bias": [],
                    "from_hmm": [], "to_hmm": [], "from_acc": [], "to_acc": [], "from_env": [], "to_env": [], "acc": []}
    column_indices = [0, 1, 2, 3, 5, 6, 7, 8, 15, 16, 17, 18, 19, 20, 21] #columns in hmm data output
    data_hmm_qry = {key: [] for key in data_hmm}
    data_hmm_sbj = {key: [] for key in data_hmm}

    with open(prot_domain_loc_hmm + qry.split("__")[0] + "_hmm.txt", "r") as hmm:
        for line in hmm:
            if not line.startswith("#"): #skip lines starting with #
                columns = re.split(r'\s+', line.strip())
                #only store lines that match query
                if columns[3].rsplit("_", 1)[0] == qry.split("__")[1]:
                    for key, idx in zip(data_hmm_qry.keys(), column_indices):
                        data_hmm_qry[key].append(columns[idx])

    with open(prot_domain_loc_hmm + sbj.split("__")[0] + "_hmm.txt", "r") as hmm:
        for line in hmm:
            if not line.startswith("#"): #skip lines starting with #
                columns = re.split(r'\s+', line.strip())
                #only store lines that match subject
                if columns[3].rsplit("_", 1)[0] == sbj.split("__")[1]:
                    for key, idx in zip(data_hmm_sbj.keys(), column_indices):
                        data_hmm_sbj[key].append(columns[idx])

    #translate protein coordinates of TE (here: query) to DNA coordinates
    for key in ["qlen", "from_acc", "to_acc"]:
        data_hmm_qry[key] = [int(value) for value in data_hmm_qry[key]]
        data_hmm_sbj[key] = [int(value) for value in data_hmm_sbj[key]]
    translate_protein_to_dna(data_hmm_qry)
    translate_protein_to_dna(data_hmm_sbj)


    #get overlap of hmm's between query & subject
    overlap_hmm = list(set(data_hmm_qry["accession"]).intersection(data_hmm_sbj["accession"]))

    #only keep indeces of query and subject dataframe if it's overlapping
    qry_idx = [idx for idx, i in enumerate(data_hmm_qry["accession"]) if i in overlap_hmm]
    sbj_idx = [idx for idx, i in enumerate(data_hmm_sbj["accession"]) if i in overlap_hmm]


    for key in data_hmm_qry.keys():
        data_hmm_qry[key] = [data_hmm_qry[key][idx] for idx in qry_idx]

    for key in data_hmm_sbj.keys():
        data_hmm_sbj[key] = [data_hmm_sbj[key][idx] for idx in sbj_idx]

    #create dictionaries for query and subject rps data
    data_rps = {"query_name": [], "subject_name": [], "pident": [], "length": [], "mismatch": [], "gapopen": [],
                "qstart": [], "qend": [], "sstart": [], "send": [], "e_value": [], "bitscore": []}
    data_rps_qry = {key: [] for key in data_rps}
    data_rps_sbj = {key: [] for key in data_rps}

    #do the same for rpstblastn
    with open(prot_domain_loc_rps + qry.split("__")[0], "r") as rps:
        for line in rps:
            if not line.startswith("#"): #skip lines starting with #
                columns = re.split(r'\s+', line.strip())
                if columns[0] == qry.split("__")[1]:
                    for key, value in zip(data_rps_qry.keys(), columns):
                        data_rps_qry[key].append(value)


    with open(prot_domain_loc_rps + sbj.split("__")[0], "r") as rps:
        for line in rps:
            if not line.startswith("#"): #skip lines starting with #
                columns = re.split(r'\s+', line.strip())
                if columns[0] == sbj.split("__")[1]:
                    for key, value in zip(data_rps_sbj.keys(), columns):
                        data_rps_sbj[key].append(value)

    #get overlap of cdd profiles (rpstblastn output)
    overlap_rps = list(set(data_rps_qry["subject_name"]).intersection(data_rps_sbj["subject_name"]))

    #only keep indeces of query and subject dataframe if it's overlapping
    qry_idx = [idx for idx, i in enumerate(data_rps_qry["subject_name"]) if i in overlap_rps]
    sbj_idx = [idx for idx, i in enumerate(data_rps_sbj["subject_name"]) if i in overlap_rps]

    for key in data_rps_qry.keys():
        data_rps_qry[key] = [data_rps_qry[key][idx] for idx in qry_idx]

    for key in data_rps_sbj.keys():
        data_rps_sbj[key] = [data_rps_sbj[key][idx] for idx in sbj_idx]

    return data_hmm_qry, data_hmm_sbj, data_rps_qry, data_rps_sbj

def within(number_to_check,range1, range2):
    #if range2 is larger than range1, swap the numbers
    if range2 < range1:
        tmp = range1
        range1 = range2
        range2 = tmp
    #if number is outside range1 and range2:
    if (number_to_check < range1 and number_to_check < range2) or (number_to_check > range1 and number_to_check > range2):
        return False
    else: #if number is within range1 and range2
        return True

#function checking if there is overlap between to ranges (for example: range 2:10 and range 8:12 have 10, 11 and 12 in common in one-based system)
#returns TRUE/FALSE
def overlap_ranges(rangeA_1, rangeA_2, rangeB_1, rangeB_2):
    if rangeA_1 > rangeA_2:
        tmp = rangeA_1
        rangeA_1 = rangeA_2
        rangeA_2 = tmp
    if rangeB_1 > rangeB_2:
        tmp = rangeB_1
        rangeB_1 = rangeB_2
        rangeB_2 = tmp
    #python is zerobased but we need onebased, so add +1 to range*_2 numbers
    #if there is no overlap:
    if len(list(set(range(rangeA_1, rangeA_2 + 1)).intersection(range(rangeB_1, rangeB_2 + 1)))) == 0:
        return False
    else: #if there is overlap:
        return True

#function to determine overlap between to ranges (for example: range 2:10 and range 8:12 have 10, 11 and 12 in common in one-based system)
#returns overlap coordinates
def find_overlap_ranges(rangeA_1, rangeA_2, rangeB_1, rangeB_2):
    if rangeA_1 > rangeA_2:
        tmp = rangeA_1
        rangeA_1 = rangeA_2
        rangeA_2 = tmp
    if rangeB_1 > rangeB_2:
        tmp = rangeB_1
        rangeB_1 = rangeB_2
        rangeB_2 = tmp
    #python is zerobased but we need onebased, so add +1 to range*_2 numbers
    #if there is no overlap:
    return list(sorted(set(range(rangeA_1, rangeA_2 + 1)).intersection(range(rangeB_1, rangeB_2 + 1))))

def prot_dom_coordinates(dict_prot_domain, query, subject, qstart, qend, sstart, send, SUBJECT_REVERSE):
    #if dict_prot_domain is empty
    if not dict_prot_domain["query_name"]:
        return dict_prot_domain

    #if df contains hmm output
    if "from_hmm" in dict_prot_domain.keys():
        #find out if this dict is about query or subject TE (remember: HMM query name has "_ORF" at end, remove this)
        if query.endswith(dict_prot_domain["query_name"][0].rsplit("_", 1)[0]):
            #dict is hmm and about query TE, so always forward.
            idx_to_keep = [i for i in range(0,len(dict_prot_domain["e_value"]))
                           if within(dict_prot_domain["from_acc"][i], qstart, qend) or within(dict_prot_domain["to_acc"][i], qstart, qend)]

            #remove indices that are outside the alignment
            for key in dict_prot_domain.keys():
                dict_prot_domain[key] = [dict_prot_domain[key][idx] for idx in idx_to_keep]

            #for remaining protein domains, convert coordinates to alignment coordinates
            for key in ["from_acc", "to_acc"]:
                dict_prot_domain[key] = [i - qstart for i in dict_prot_domain[key]]


        elif subject.endswith(dict_prot_domain["query_name"][0].rsplit("_",1)[0]):
            #dict is hmm and about subject TE.
            idx_to_keep = [i for i in range(0,len(dict_prot_domain["e_value"]))
                          if within(dict_prot_domain["from_acc"][i], send, sstart) or within(dict_prot_domain["to_acc"][i], send, sstart)]

            #remove indices that are outside the alignment
            for key in dict_prot_domain.keys():
                dict_prot_domain[key] = [dict_prot_domain[key][idx] for idx in idx_to_keep]

            #for remaining protein domains, convert coordinates to alignment coordinates
            if not SUBJECT_REVERSE:
                for key in ["from_acc", "to_acc"]:
                    dict_prot_domain[key] = [i - sstart for i in dict_prot_domain[key]]
            elif SUBJECT_REVERSE:
                for key in ["from_acc", "to_acc"]:
                    dict_prot_domain[key] = [sstart - i for i in dict_prot_domain[key]]

        else:
            print(f"something went wrong with recognizing whether hmm dict is query or subject TE inside the prot_dom_coordinates function")
            return {}



    #if dict contains rps output
    elif "qstart" in dict_prot_domain.keys():
        #find out if this dict is about query or subject TE
        if query.endswith(dict_prot_domain["query_name"][0]):
            #dict is hmm and about query TE, so always forward. Keep in mind: qstart is alignment coord on query TE, qstart in dict is protein domain coords on query TE
            idx_to_keep = [i for i in range(0,len(dict_prot_domain["e_value"]))
                           if within(int(dict_prot_domain["qstart"][i]), qstart, qend) or within(int(dict_prot_domain["qend"][i]), qstart, qend)]


            #remove indices that are outside the alignment
            for key in dict_prot_domain.keys():
                dict_prot_domain[key] = [dict_prot_domain[key][idx] for idx in idx_to_keep]

            #for remaining protein domains, convert coordinates to alignment coordinates
            for key in ["qstart", "qend"]:
                dict_prot_domain[key] = [int(i) - qstart for i in dict_prot_domain[key]]


        elif subject.endswith(dict_prot_domain["query_name"][0]):
            #dict is rps and about subject TE.
            idx_to_keep = [i for i in range(0,len(dict_prot_domain["e_value"]))
                               if within(int(dict_prot_domain["qstart"][i]), sstart, send) or within(int(dict_prot_domain["qend"][i]), sstart, send)]


            #remove indices that are outside the alignment
            for key in dict_prot_domain.keys():
                dict_prot_domain[key] = [dict_prot_domain[key][idx] for idx in idx_to_keep]

            #for remaining protein domains, convert coordinates to alignment coordinates
            if not SUBJECT_REVERSE:
                for key in ["qstart", "qend"]:
                    dict_prot_domain[key] = [int(i) - sstart for i in dict_prot_domain[key]]
            elif SUBJECT_REVERSE:
                for key in ["qstart", "qend"]:
                    dict_prot_domain[key] = [sstart - int(i) for i in dict_prot_domain[key]]


    return dict_prot_domain

#function to merge hmm and rps dictionaries with protein domain info
def merge_hmm_rps(dict_hmm_prot, dict_rps_prot):
    #first create empty dict to store merged results
    merged = {"TE_name": [], "Prot_name": [], "TE_start": [], "TE_end": [], "prot_start": [], "prot_end": []}

    #merge TE_name column
    merged["TE_name"] = copy.deepcopy(dict_hmm_prot["query_name"])
    merged["TE_name"] = [i.rsplit("_", 1)[0] for i in merged["TE_name"]]
    merged["TE_name"].extend(dict_rps_prot["query_name"])

    #merge TE and protein start & end coordinates
    hmm_names = ["accession", "from_acc", "to_acc", "from_hmm", "to_hmm"]
    rps_names = ["subject_name", "qstart", "qend", "sstart", "send"]
    idx = 0
    for key in list(merged.keys())[1:]:
        merged[key] = copy.deepcopy(dict_hmm_prot[hmm_names[idx]])
        merged[key].extend(dict_rps_prot[rps_names[idx]])
        idx += 1

    #transform coordinate keys from string to integers
    for key in list(merged.keys())[2:]:
        merged[key] = [int(value) for value in merged[key]]

    return merged

def filter_out_non_overlapping_prot_dom_coords(qry_prot, sbj_prot, overlap):
        #create empty dict to store indices to remove of query protein domain dict and subject protein domain dict
    idx_to_remove = {"qry": [], "sbj": []}

    for prot_dom in overlap:
        idx_q = [i for i, x in enumerate(qry_prot["Prot_name"]) if x == prot_dom]
        idx_s = [i for i, x in enumerate(sbj_prot["Prot_name"]) if x == prot_dom]

        #if both protein domain occur once per TE
        if len(idx_q) == 1 and len(idx_s) == 1:
            #if protein domains do NOT overlap (on the protein domain coordinates)
            if not overlap_ranges(qry_prot["prot_start"][idx_q[0]], qry_prot["prot_end"][idx_q[0]], sbj_prot["prot_start"][idx_s[0]], sbj_prot["prot_end"][idx_s[0]]):
                #print(prot_dom)
                idx_to_remove["qry"].append(idx_q)
                idx_to_remove["sbj"].append(idx_s)

        #if protein domain occurs once on one TE and multiple times on other
        elif len(idx_q) == 1 or len(idx_s) == 1:
            if len(idx_q) == 1:
                for i in idx_s:
                    if not overlap_ranges(qry_prot["prot_start"][idx_q[0]], qry_prot["prot_end"][idx_q[0]], sbj_prot["prot_start"][i], sbj_prot["prot_end"][i]):
                        idx_to_remove["sbj"].append(str(i))
            elif len(idx_s) == 1:
                for i in idx_q:
                    if not overlap_ranges(qry_prot["prot_start"][i], qry_prot["prot_end"][i], sbj_prot["prot_start"][idx_s[0]], sbj_prot["prot_end"][idx_s[0]]):
                        idx_to_remove["qry"].append(str(i))
            #if all elements of the TE with multiple copies hits of this protein domain are getting removed, remove also
            #the protein domain on the TE with 1 copy. Make sure idx_to_remove is not empty (otherwise all() will give True).
            if idx_to_remove["qry"] and all(elem in idx_q for elem in idx_to_remove["qry"]):
                idx_to_remove["sbj"].append(idx_s)
            if idx_to_remove["sbj"] and all(elem in idx_s for elem in idx_to_remove["sbj"]):
                idx_to_remove["qry"].append(idx_q)

        #if protein domain occurs more than once on both TEs
        else:
            #1 check for all query TE protein domains if they have overlap with subject TE domains
            for i in idx_q:
                if not any(overlap_ranges(qry_prot["prot_start"][i], qry_prot["prot_end"][i], sbj_prot["prot_start"][s], sbj_prot["prot_end"][s]) for s in idx_s):
                    #so there is no overlap of this query TE prot domain with subject TE domains
                    idx_to_remove["qry"].append(i)

            #2 check the other way around: check subject TE protein domains if they have overlap with query TE domains
            for i in idx_s:
                if not any(overlap_ranges(qry_prot["prot_start"][q], qry_prot["prot_end"][q], sbj_prot["prot_start"][i], sbj_prot["prot_end"][i]) for q in idx_q):
                    #so there is no overlap of this query TE prot domain with subject TE domains
                    idx_to_remove["sbj"].append(i)

    #remove indices stored in idx_to_remove
    for key in qry_prot.keys():
        qry_prot[key] = [item for i, item in enumerate(qry_prot[key]) if i not in idx_to_remove["qry"]]
    for key in sbj_prot.keys():
        sbj_prot[key] = [item for i, item in enumerate(sbj_prot[key]) if i not in idx_to_remove["sbj"]]

#function to eliminate protein domains that overlap with SAME domain on SAME TE
#example: PF00075 domain on query TE with alignment coordinates 1-50 and another one on 30-330
#         check with subject TE which of the two domains (^) has biggest overlap with PF00075 domain on subject TE
def filter_out_overlapping_alignment_coords_per_dom(qry_prot, sbj_prot):
    #find overlapping protein domains
    overlap = set(qry_prot["Prot_name"]).intersection(sbj_prot["Prot_name"])

    for prot_dom in overlap:
        idx_q = [i for i, x in enumerate(qry_prot["Prot_name"]) if x == prot_dom]
        idx_s = [i for i, x in enumerate(sbj_prot["Prot_name"]) if x == prot_dom]

        #if protein domain only occurs once on both TEs, continue
        if len(idx_q) == 1 and len(idx_s) == 1:
            continue
        #protein domain occurs more than once on one TE, or both TEs.
        else:
            #1 start with query TE
            to_remove = []
            if len(idx_q) != 1:
                #go over all possible pairs within these protein domains
                for pair in list(combinations(idx_q, 2)):
                    #if one of the two indexes in this pair is already listed to be removed, continue
                    if pair[0] in to_remove or pair[1] in to_remove:
                        continue

                    #if they overlap, only keep domain which has the largest overlap in the protein domain coords with the other TE (subject in this case)
                    if overlap_ranges(qry_prot["TE_start"][pair[0]], qry_prot["TE_end"][pair[0]], qry_prot["TE_start"][pair[1]], qry_prot["TE_end"][pair[1]]):
                        overlap_coords_0 = [len(find_overlap_ranges(qry_prot["TE_start"][pair[0]], qry_prot["TE_end"][pair[0]], sbj_prot["TE_start"][i], sbj_prot["TE_end"][i]))
                                          for i in idx_s]
                        overlap_coords_1 = [len(find_overlap_ranges(qry_prot["TE_start"][pair[1]], qry_prot["TE_end"][pair[1]], sbj_prot["TE_start"][i], sbj_prot["TE_end"][i]))
                                          for i in idx_s]

                        #add index of query protein domain that has shortest overlap with subject protein domain
                        to_remove.append(pair[1] if max(overlap_coords_0) >= max(overlap_coords_1) else pair[0])

                #remove query protein domains
                for key in qry_prot.keys():
                    qry_prot[key] = [item for i, item in enumerate(qry_prot[key]) if i not in to_remove]

                #after removing query prot domains, update idx_q
                idx_q = [i for i, x in enumerate(qry_prot["Prot_name"]) if x == prot_dom]


            #2 do the same for subject TE
            to_remove = []
            if len(idx_s) != 1:
                #go over all possible pairs within these protein domains
                for pair in list(combinations(idx_s, 2)):
                    #if they overlap, only keep domain which has the largest overlap in the protein domain coords with the other TE (subject in this case)
                    #if one of the two indexes in this pair is already listed to be removed, continue
                    if pair[0] in to_remove or pair[1] in to_remove:
                        continue

                    #if they overlap, only keep domain which has the largest overlap in the protein domain coords with the other TE (subject in this case)
                    if overlap_ranges(sbj_prot["TE_start"][pair[0]], sbj_prot["TE_end"][pair[0]], sbj_prot["TE_start"][pair[1]], sbj_prot["TE_end"][pair[1]]):
                        overlap_coords_0 = [len(find_overlap_ranges(sbj_prot["TE_start"][pair[0]], sbj_prot["TE_end"][pair[0]], qry_prot["TE_start"][i], qry_prot["TE_end"][i]))
                                          for i in idx_q]
                        overlap_coords_1 = [len(find_overlap_ranges(sbj_prot["TE_start"][pair[1]], sbj_prot["TE_end"][pair[1]], qry_prot["TE_start"][i], qry_prot["TE_end"][i]))
                                          for i in idx_q]

                        #add index of sbj protein domain that has shortest overlap with query protein domain
                        to_remove.append(pair[1] if max(overlap_coords_0) >= max(overlap_coords_1) else pair[0])

                #remove query protein domains
                for key in sbj_prot.keys():
                    sbj_prot[key] = [item for i, item in enumerate(sbj_prot[key]) if i not in to_remove]

# Calculate differences and translate coordinates
def calculate_and_append_coordinates(prot, idx, coord_key, prot_coord_overlap, is_forward, alncoord_to_extract):
    diff_start = prot_coord_overlap[0] - prot["prot_start"][idx]
    diff_end = prot["prot_end"][idx] - prot_coord_overlap[-1]
    if is_forward:
        alncoord_to_extract[f"{coord_key}_start"].append(prot["TE_start"][idx] + (diff_start * 3))
        alncoord_to_extract[f"{coord_key}_end"].append(prot["TE_end"][idx] - (diff_end * 3))
    else:
        alncoord_to_extract[f"{coord_key}_start"].append(prot["TE_start"][idx] - (diff_start * 3))
        alncoord_to_extract[f"{coord_key}_end"].append(prot["TE_end"][idx] + (diff_end * 3))

#filter out protein domains that are shared between TEs, but are on different alignment coords.
#(for example: PF00075 is on sbj on 15-200 and 400-500, and on qry on 15-200, max difference is 15. Sbj 400-500 must be filtered out,
#it is not on query TE.)
def get_align_coord_overlap_prot_domains(qry_prot, sbj_prot, MAX_DIFFERENCE):
    #find overlapping protein domains
    overlap = set(qry_prot["Prot_name"]).intersection(sbj_prot["Prot_name"])

    #create dictionary with align_coords to extract
    alncoord_to_extract = {"qry_start": [], "qry_end": [], "sbj_start": [], "sbj_end": [], "domain": []}

    for prot_dom in overlap:
        idx_q = [i for i, x in enumerate(qry_prot["Prot_name"]) if x == prot_dom]
        idx_s = [i for i, x in enumerate(sbj_prot["Prot_name"]) if x == prot_dom]

        #1 go over query indices of this prot domain
        to_remove = []
        for i in idx_q:
            #if it has overlap with subject protein domain, keep it
            if any(overlap_ranges(qry_prot["TE_start"][i], qry_prot["TE_end"][i], sbj_prot["TE_start"][j] - MAX_DIFFERENCE, sbj_prot["TE_end"][j] - MAX_DIFFERENCE) for j in idx_s) or \
            any(overlap_ranges(qry_prot["TE_start"][i], qry_prot["TE_end"][i], sbj_prot["TE_start"][j] + MAX_DIFFERENCE, sbj_prot["TE_end"][j] + MAX_DIFFERENCE) for j in idx_s):
                for j in idx_s: #for this specific index, double check it has overlapping align coordinates with sbj TE
                    if overlap_ranges(qry_prot["TE_start"][i], qry_prot["TE_end"][i], sbj_prot["TE_start"][j] - MAX_DIFFERENCE, sbj_prot["TE_end"][j] - MAX_DIFFERENCE) or \
                    overlap_ranges(qry_prot["TE_start"][i], qry_prot["TE_end"][i], sbj_prot["TE_start"][j] + MAX_DIFFERENCE, sbj_prot["TE_end"][j] + MAX_DIFFERENCE):
                        #now find the overlap between PROTEIN DOMAIN coords, and translate these back to alignment coords to extract DNA from TEs
                        prot_coord_overlap = find_overlap_ranges(qry_prot["prot_start"][i], qry_prot["prot_end"][i], sbj_prot["prot_start"][j], sbj_prot["prot_end"][j])

                        #double check it did not return an empty list
                        if len(prot_coord_overlap) > 4:


                            # Determine direction and calculate differences
                            is_forward_qry = qry_prot["TE_start"][i] < qry_prot["TE_end"][i]
                            is_forward_sbj = sbj_prot["TE_start"][j] < sbj_prot["TE_end"][j]

                            #both domains must be in same direction (this is ALIGNMENT coords, not TE coords). In case they are not, skip this domain.
                            if is_forward_qry is not is_forward_sbj:
                                continue

                            #add this domain already to final dataframe
                            alncoord_to_extract["domain"].append(prot_dom)

                            # Process query and subject protein domains
                            calculate_and_append_coordinates(qry_prot, i, "qry", prot_coord_overlap, is_forward_qry, alncoord_to_extract)
                            calculate_and_append_coordinates(sbj_prot, j, "sbj", prot_coord_overlap, is_forward_qry, alncoord_to_extract)

                            #in rare cases, some coordinates are <0, meaning that there likely are many gaps in the alignment of the protein domain profile and the TE.
                            #there is no easy solution, therefore remove these protein domains. Additionally, sometimes because of the gaps the directionality
                            #changes. Check this, and if directionality has changed -> remove domain.
                            is_forward_qry_te = alncoord_to_extract["qry_start"][-1] < alncoord_to_extract["qry_end"][-1]
                            is_forward_sbj_te = alncoord_to_extract["sbj_start"][-1] < alncoord_to_extract["sbj_end"][-1]

                            if any(alncoord_to_extract[k][-1] < 0 for k in ["qry_start", "qry_end", "sbj_start", "sbj_end"]) or \
                            is_forward_qry_te is not is_forward_sbj_te or is_forward_qry is not is_forward_qry_te:
                                for key in alncoord_to_extract:
                                    alncoord_to_extract[key] = alncoord_to_extract[key][:-1]

            #this protein domain on query does not overlap with protein domains on subject TE
            else:
                to_remove.append(i)

        #remove query protein domains
        for key in qry_prot.keys():
            qry_prot[key] = [item for i, item in enumerate(qry_prot[key]) if i not in to_remove]

        #after removing query prot domains, update idx_q
        idx_q = [i for i, x in enumerate(qry_prot["Prot_name"]) if x == prot_dom]

        #2 go over subject indices of this prot domain (we already have all the align coords we need, but for completeness we also remove
        #protein domains on the subject that do not share alignment coordinates with the query protein domain.
        to_remove = []
        for i in idx_s:
            if not any(overlap_ranges(qry_prot["TE_start"][j], qry_prot["TE_end"][j], sbj_prot["TE_start"][i] - MAX_DIFFERENCE, sbj_prot["TE_end"][i] - MAX_DIFFERENCE) for j in idx_q) or not \
            any(overlap_ranges(qry_prot["TE_start"][j], qry_prot["TE_end"][j], sbj_prot["TE_start"][i] + MAX_DIFFERENCE, sbj_prot["TE_end"][i] + MAX_DIFFERENCE) for j in idx_q):
                to_remove.append(i)

        for key in sbj_prot.keys():
            sbj_prot[key] = [item for i, item in enumerate(sbj_prot[key]) if i not in to_remove]

    return alncoord_to_extract

#function to translate alignment coordinates back to TE coordinates
def aligncoord_to_tecoord(aln_coord_to_extract, SUBJECT_REVERSE, qstart, sstart, send):
    #for query TE, add start coordinates of TE alignment to get TE coordinates of protein domains
    for key in ["qry_start", "qry_end"]:
        aln_coord_to_extract[key] = [num + qstart for num in aln_coord_to_extract[key]]

    #do the same for subject TE, but use send for if subject TE is in reverse direction, and sstart if it is forward.
    #if subject is reverse:
    if SUBJECT_REVERSE:
        for key in ["sbj_start", "sbj_end"]:
            aln_coord_to_extract[key] = [abs(num - sstart) for num in aln_coord_to_extract[key]]
    #if subject is forward:
    else:
        for key in ["sbj_start", "sbj_end"]:
            aln_coord_to_extract[key] = [num + sstart for num in aln_coord_to_extract[key]]

    to_remove = set()
    #remove overlapping protein domains (so different protein domains, but are located on same location on TE)
    for i in range(0,len(aln_coord_to_extract["qry_start"]) - 1):
        for j in range(i + 1, len(aln_coord_to_extract["qry_start"])):
            if overlap_ranges(aln_coord_to_extract["qry_start"][i],
                               aln_coord_to_extract["qry_end"][i],
                               aln_coord_to_extract["qry_start"][j],
                               aln_coord_to_extract["qry_end"][j]):
                if abs(aln_coord_to_extract["qry_start"][i] - aln_coord_to_extract["qry_end"][i]) > abs(aln_coord_to_extract["qry_start"][j] - aln_coord_to_extract["qry_end"][j]):
                    to_remove.add(j)
                else:
                    to_remove.add(i)

    for key in aln_coord_to_extract.keys():
        aln_coord_to_extract[key] = [item for k, item in enumerate(aln_coord_to_extract[key]) if k not in to_remove]

def suff_length_aa(aln_coord_to_extract):
    qry_sum = sum([abs(aln_coord_to_extract["qry_start"][i] - aln_coord_to_extract["qry_end"][i]) + 1 for i in range(0,len(aln_coord_to_extract["qry_start"]))])
    sbj_sum = sum([abs(aln_coord_to_extract["sbj_start"][i] - aln_coord_to_extract["sbj_end"][i]) + 1 for i in range(0,len(aln_coord_to_extract["sbj_start"]))])
    return qry_sum, sbj_sum

def protein_alignment(aln_coord_to_extract, query, subject, qry_end, sbj_end):
    #store ORF and AA coordinates
    for key in ["ORF_qry","ORF_sbj","AA_start_qry","AA_end_qry","AA_start_sbj","AA_end_sbj"]:
        aln_coord_to_extract[key] = []


    #loop over length of lists in dictionary
    for i in range(0, len(aln_coord_to_extract["domain"])):
        #determine ORF for query and subject
        for key in ["qry","sbj"]:
            te = query if key == "qry" else subject
            te = te.split("__")[1]

            te1, te2 = (query, subject) if key == "qry" else (subject, query)


            if aln_coord_to_extract[f"{key}_start"][i] < aln_coord_to_extract[f"{key}_end"][i]:
                #it's forward, so ORF1, 2, or 3
                ORF = ((aln_coord_to_extract[f"{key}_start"][i] - 1) % 3) + 1
                #print(f"{aln_coord_to_extract[f"{key}_start"][i]} % 3 + 1 = {ORF}")
                #calculate nucleotide coordinates of AA start and end
                aln_coord_to_extract[f"AA_start_{key}"].append((aln_coord_to_extract[f"{key}_start"][i] - ORF + 3) // 3)
                aln_coord_to_extract[f"AA_end_{key}"].append((aln_coord_to_extract[f"{key}_end"][i] - ORF + 1) // 3)


                #create bed file with DNA coords
                with open(f"{output_dir}bed/{te1}_{te2}_{aln_coord_to_extract["domain"][i]}_{i}_DNA.bed", 'w') as bed_dna:
                    bed_dna.write(f"{te}\t{aln_coord_to_extract[f"{key}_start"][i] - 1}\t{aln_coord_to_extract[f"{key}_end"][i]}\t{aln_coord_to_extract["domain"][i]}\t0\t+")

            else:
                #start = qry_start if key == "qry" else sbj_start, otherwise it is
                end = qry_end if key == "qry" else sbj_end
                #print(end)
                #it's reverse, so ORF4, 5, 6
                #calculate at which nucleotide the ORF starts from reverse compliment (DIF = 0, 1, or 2) and calculate remainder of te length (l = 0, 1 or 2)
                DIF = (end - aln_coord_to_extract[f"{key}_start"][i]) % 3
                l = end % 3

                #calculate ORF based on DIF and l (ORF determination depends on total length of the TE, not where the protein starts compared to reverse compliment unfortunately...
                orf = [4,6,5]
                ORF = orf[DIF - l]

                #calculate nucleotide coordinates of AA start and end
                aln_coord_to_extract[f"AA_start_{key}"].append((end - aln_coord_to_extract[f"{key}_start"][i] - DIF + 3) // 3)
                aln_coord_to_extract[f"AA_end_{key}"].append((end - aln_coord_to_extract[f"{key}_end"][i] - DIF + 1) // 3)

                #create bed file with DNA coords (these are swapped)
                with open(f"{output_dir}bed/{te1}_{te2}_{aln_coord_to_extract["domain"][i]}_{i}_DNA.bed", 'w') as bed_dna:
                    bed_dna.write(f"{te}\t{aln_coord_to_extract[f"{key}_end"][i] - 1}\t{aln_coord_to_extract[f"{key}_start"][i]}\t{aln_coord_to_extract["domain"][i]}\t0\t-")

            #create bed file with AA coords (these do not need to be swapped for reverse ORFs)
            with open(f"{output_dir}bed/{te1}_{te2}_{aln_coord_to_extract["domain"][i]}_{i}_AA.bed", 'w') as bed_aa:
                bed_aa.write(f"{te}_{ORF}\t{aln_coord_to_extract[f"AA_start_{key}"][i] - 1}\t{aln_coord_to_extract[f"AA_end_{key}"][i]}\t{aln_coord_to_extract["domain"][i]}\t0\t+")
            aln_coord_to_extract[f"ORF_{key}"].append(ORF)
            #print(aln_coord_to_extract)






    #create bed file with coordinates, then extract sequences with bedtools.
    #put qry and sbj sequences together
    #get dna sequences of qry and sbj

def run_bedtools_getfasta(input_fasta, bed_file, output_file):
    command = ["bedtools", "getfasta", "-s", "-fi", input_fasta, "-bed", bed_file, "-fo", output_file]
    c = subprocess.Popen(command)
    c.wait()
    if not os.path.exists(output_file):
        raise RuntimeError(f"Output file {output_file} was not created correctly.")


def merge_fasta(input_files, merged_fasta):
    sequences = {}

    for file in input_files:
        with open(file) as f:
            header = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:]
                    if len(header.split(":")) > 4:
                        header = line[1:].rsplit(":",1)[0].rsplit("_",1)[0]
                    #print(header)
                else:
                    if header in sequences:
                        sequences[header] += line
                    else:
                        sequences[header] = line

    with open(merged_fasta, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n{seq}\n")

# 1 get list of overlapping prot domains

not_suff = 0
suff = 0
no_prot = 0

#with open("blastn_Mycafr1_chained_suff_align.txt", "r") as blast_out:
Ks_files = []
with open(args.blast, "r") as blast_out:
#with open("blastn_Blugra3_16_test.txt", "r") as blast_out:
    lines = blast_out.readlines()
    #loop over hits
    for hit in lines:
        #print(hit)
        #get list of overlapping prot domains, first extract query/subject/coordinates of alignment
        query, subject, _, length, _, _, qstart, qend, sstart, send, *_ = hit.split("\t")
        length, qstart, qend, sstart, send = map(int, (length, qstart, qend, sstart, send)) #convert to integers
        hmm_query, hmm_subject, rps_query, rps_subject = get_overlapping_prot_domains(query, subject)

        #2 find out if subject is in reverse direction (query is always forward direction)
        SUBJECT_REVERSE = False
        if sstart > send:
            SUBJECT_REVERSE = True

        #3 make alignment coordinates (translate TE coordinates to "shared" alignment coordinates
        #the difference between coordinates comes from potential gaps in the alignment, introduced by either blastn or the chaining.
        #the difference between the align_qend and align_send coordinates is the exact max difference that is allowed between protein
        #domains on the query and subject TE.
        align_qstart = 1
        align_qend = qend - qstart
        align_sstart = 1
        align_send = abs(send - sstart)

        MAX_DIFFERENCE = abs(align_qend - align_send)

        #if MAX_DIFFERENCE > 100:
            #print(f"max difference = {MAX_DIFFERENCE}")
            #print(f"align_qstart: {align_qstart}, align_qend: {align_qend}, align_sstart: {align_sstart}, align_send: {align_send}")



        #3 filter out protein domains outside of TE alignment and convert protein domain coordinates to alignment coordinates
        hmm_query, hmm_subject, rps_query, rps_subject =  [prot_dom_coordinates(dct, query, subject, qstart, qend, sstart, send, SUBJECT_REVERSE)
                                                           for dct in (hmm_query, hmm_subject, rps_query, rps_subject)]

        #4 merge protein domain dictionaries (for example: merge query_rps_hmm and subject_rps_hmm)
        query_prot = merge_hmm_rps(hmm_query, rps_query)
        subject_prot = merge_hmm_rps(hmm_subject, rps_subject)


        #5 check again that protein domains are shared between query_prot and subject_prot
        overlap = set(query_prot["Prot_name"]).intersection(subject_prot["Prot_name"])
        for prot_dict in [query_prot, subject_prot]:
            idx = [i for i in range(0,len(prot_dict["TE_name"])) if prot_dict["Prot_name"][i] in overlap]
            for key in prot_dict.keys():
                prot_dict[key] = [prot_dict[key][i] for i in idx]

        #6 filter out protein domains not referring to same part of protein domain
        filter_out_non_overlapping_prot_dom_coords(query_prot, subject_prot, overlap)

        #7 filter out protein domains that are overlapping with each other on the SAME TE. (so 2 SAME domains overlap for example on query TE
        #with alignment coords -> this makes biologically no sense. Two different domains sharing alignment coordinates will be filtered out later).
        filter_out_overlapping_alignment_coords_per_dom(query_prot, subject_prot)

        #8 get TE alignment coordinates of overlapping parts of protein domains. They must refer to same TE alignment coordinates, with
        #an error-marge of MAX_difference (the difference in length between the subject and query hit).
        aln_coord_to_extract = get_align_coord_overlap_prot_domains(query_prot, subject_prot, MAX_DIFFERENCE)

        #9 turn alignment coords back into TE coords, and remove (different) protein domains that overlap on the same TE. Only keep the longest protein domain.
        aligncoord_to_tecoord(aln_coord_to_extract, SUBJECT_REVERSE, qstart, sstart, send)

        #10 check that protein alignments together are longer than 100 aa, and extract sequences and translate to protein
        query_prot_bp, subject_prot_bp = suff_length_aa(aln_coord_to_extract)

        #create directories for temp files
        for dir in ["bed", "AA", "DNA", "AA_sep", "DNA_sep", "DNA_glued", "codon", "Ks_tmp", "insuf"]:
            os.makedirs(output_dir + dir + "/", exist_ok=True)

        if query_prot_bp > 300 and subject_prot_bp > 300: #if prot alignment is longer than 100 aa
            suff += 1
        #if prot alignment is shorter than 300 bp, write hit to insufficient aa lenght file
        elif query_prot_bp == 0 and subject_prot_bp == 0:
            no_prot += 1
            with open(f"{output_dir}insuf/{query.split("__")[0]}_no_AA_length.txt", "a") as f:
                f.write(hit)
            continue

        else:
            not_suff += 1
            with open(f"{output_dir}insuf/{query.split("__")[0]}_AA_insuf_length.txt", "a") as f:
                f.write(hit)
            continue

        #print(f"number of TE hits with sufficient protein length in alignment vs insufficient vs no protein: {suff}, {not_suff}, {no_prot}")




        #11 make protein alignment
        #11.1 retrieve ORF and calculate AA coordinates. Create BED files with AA coordinates, and BED files with DNA coordinates.
        #get total length of query and subject (coordinates are 0-based) (needed for calculation in case of reverse ORF)
        length_qry = int(query.split(":")[3].split("-")[1]) - int(query.split(":")[3].split("-")[0])
        length_sbj = int(subject.split(":")[3].split("-")[1]) - int(subject.split(":")[3].split("-")[0])

        protein_alignment(aln_coord_to_extract, query, subject, length_qry, length_sbj)


        qry_sp = query.split("__")[0]
        sbj_sp = subject.split("__")[0]

        #11.2 per domain, extract AA and DNA seqs from bed file
        DNA_files = []
        for i in range(0, len(aln_coord_to_extract["domain"])):
            bed_qry_AA = f"{output_dir}bed/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}_AA.bed"
            bed_qry_DNA = f"{output_dir}bed/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}_DNA.bed"

            bed_sbj_AA = f"{output_dir}bed/{subject}_{query}_{aln_coord_to_extract["domain"][i]}_{i}_AA.bed"
            bed_sbj_DNA = f"{output_dir}bed/{subject}_{query}_{aln_coord_to_extract["domain"][i]}_{i}_DNA.bed"

            output_qry_AA = f"{output_dir}AA/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}.faa"
            output_qry_DNA = f"{output_dir}DNA/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}.fna"

            output_sbj_AA = f"{output_dir}AA/{subject}_{query}_{aln_coord_to_extract["domain"][i]}_{i}.faa"
            output_sbj_DNA = f"{output_dir}DNA/{subject}_{query}_{aln_coord_to_extract["domain"][i]}_{i}.fna"

            #final AA and DNA fasta files:
            output_AA = f"{output_dir}AA_sep/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}.faa.aln"
            output_DNA = f"{output_dir}DNA_sep/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}.fna"

            #codon alignment file:
            codon = f"{output_dir}codon/{query}_{subject}_{aln_coord_to_extract["domain"][i]}_{i}.codon.aln"

            DNA_files.append(codon)

            #11.3 use bedtools to extract DNA and AA from both query and subject
            run_bedtools_getfasta(f"{te_translated_orf_loc}{qry_sp}_translated.fasta", bed_qry_AA, output_qry_AA)
            run_bedtools_getfasta(f"{te_translated_orf_loc}{sbj_sp}_translated.fasta", bed_sbj_AA, output_sbj_AA)
            run_bedtools_getfasta(f"{fasta_loc}{qry_sp}_extracted_TE.fasta", bed_qry_DNA, output_qry_DNA)
            run_bedtools_getfasta(f"{fasta_loc}{sbj_sp}_extracted_TE.fasta", bed_sbj_DNA, output_sbj_DNA)

            #11.4 combine DNA and AA sequences of query and subject, and make AA alignment.
            subprocess.run("cat {} {} | sed 's/*/X/g' | mafft --thread 1 --auto --quiet - | awk '/^>/ {{if (seq) print seq; print; seq=\"\"}} !/^>/ {{seq=seq$0}} END {{if (seq) print seq}}' > {}".format(output_qry_AA, output_sbj_AA, output_AA), shell=True, check=True)
            subprocess.run(f"cat {output_qry_DNA} {output_sbj_DNA} > {output_DNA}", shell=True, check=True)

            #11.5 make codon alignment
            subprocess.run(f"pal2nal.pl {output_AA} {output_DNA} -output fasta -nostderr -nogap > {codon}", shell=True, check=True)


            #11.6 delete intermediate files
            for file in [bed_qry_AA, bed_qry_DNA, bed_sbj_AA, bed_sbj_DNA, output_qry_AA, output_qry_DNA, output_sbj_AA, output_sbj_DNA]:
                os.remove(file)


        #12 glue together DNA sequences and AA alignment
        output_DNA_glued = f"{output_dir}DNA_glued/{query}_{subject}_DNA.fna"

        #13 glue together
        merge_fasta(DNA_files, output_DNA_glued)

        #remove temporary files
        for file in DNA_files:
            os.remove(file)

        #14 calculate Ks.
        output_Ks = f"{output_dir}Ks_tmp/{query}_{subject}_Ks.txt"
        subprocess.run(f"Rscript KaKs.R {output_DNA_glued} {output_Ks} {query} {subject}", shell = True, check = True)

        with open(output_Ks, 'r') as Ks_tmp, open(output_file_ks + "_tmp", 'a') as Ks:
                Ks.write(Ks_tmp.read())

        os.remove(output_Ks)

#after completion, copy temporary Ks file into final Ks file
#if there were no hits with enough overlapping protein domains, create empty Ks file
if os.path.isfile(output_file_ks + "_tmp"):
    os.rename(output_file_ks + "_tmp", output_file_ks)
else:
    with open(output_file_ks, 'a') as file:
        file.write('')
