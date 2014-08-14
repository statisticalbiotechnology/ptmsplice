
import os
import lxml.etree as etree


def parse_percolator_xml(filename, q_value_cutoff=0.01):
    """
    parse_percolator_xml - takes a file name and q value cutoff and returns 
        a dictionary keyed on identified peptides (give the specified cutoff) 
        with a list of potential protein matches as values.
    """
    NS = '{http://per-colator.com/percolator_out/15}'
    with open(filename, 'r') as f:
        perc_output = etree.parse(f).getroot()
    peptides = {}
    for pep in perc_output.iter('{}peptide'.format(NS)):
        pep_seq = pep.attrib.get('{}peptide_id'.format(NS))
        
        if pep_seq.find("UNIMOD:1") == -1:
            continue

        q_value_str = pep.find('{}q_value'.format(NS)).text
        try:
            q_value = float(q_value_str)
        except ValueError:
            raise ValueError("Invalid q value {}".format(q_value_str))
        if pep_seq in peptides:
            raise RuntimeError("Duplicate peptide in percolator xml output")
        if q_value > q_value_cutoff:
            continue
        protein_ids = [p.text for p in pep.iter('{}protein_id'.format(NS))]
        peptides[pep_seq] = protein_ids
    return peptides


def main():
    inFolder = "../input/Acetylation_data/XML_FILES/"
    outFolder = "../output/acetylation/"

    for inFile in os.listdir(inFolder):
        dict_peptides = parse_percolator_xml(inFolder + inFile)

        with open(outFolder + inFile.replace(".xml", ".txt"), "w") as outFile:
            for key in dict_peptides:
                outFile.write(key.replace("[UNIMOD:1]", "") + ":" + ("".join([";" + entry for entry in dict_peptides[key]]))[1:] + "\n")




if __name__ == '__main__':
    main()