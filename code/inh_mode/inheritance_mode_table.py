"""
The goal of inh_mode is to provide a list of inheritance modes and genotype names based on:
- mother genotype, father genotype, child genotype
- chromosome
- sex of child
- novoPP (novocaller posterior probability)

This module provides some examples for incfheritance_mode.
- based on combinations of ["0/0", "0/1", "1/1"]
- or the real variants on cgpaptest for NA12879, provided in "<DATA_DIR>/variants_NA12879.json"
"""

from os import path
import json

import pandas as pd

import inheritance_mode

DATA_DIR = path.join(base.ROOT_DIR, "data")

def all_scenarios_table():
    """
    loop through scenarios for
    - mother genotype, father genotype, child genotype [0/0, 0/1, 1/1]
    - chromosome
    - sex of child
    - novoPP (novocaller posterior probability)

    and create a table for output of inheritance_mode, output as tsv:
    "<DATA_DIR>all_scenarios_table.tsv"
    """

    variant = {
        "samplegeno": [
            {
                "samplegeno_role": "self",
                "samplegeno_numgt": "0/1",
                "samplegeno_sex": "male"
            },
            {
                "samplegeno_role": "mother",
                "samplegeno_numgt": "0/1",
                "samplegeno_sex": "female"
            },
            {
                "samplegeno_role": "father",
                "samplegeno_numgt": "0/1",
                "samplegeno_sex": "male"
            }],
        "variant": {"CHROM": "1"},
        "novoPP": 0.5
    }

    genotypes = ["0/0", "0/1", "1/1"]
    rows = []
    for genotype_mother in genotypes:
        for genotype_father in genotypes:
            for genotype_self in genotypes:
                if genotype_mother == "0/0" and genotype_father == "0/0" and genotype_self == "0/0":
                    continue
                conditions = ["autosome", "chrX - male child", "chrY - male child",
                              "chrX - female child", "chrY - female child"]
                if genotype_mother == "0/0" and genotype_father == "0/0" and genotype_self == "0/1":
                    conditions = ["autosome - novocaller high", "autosome - novocaller low",
                                  "chrX - male child", "chrY - male child",
                                  "chrX - female child - novoPP=0",
                                  "chrX - female child - novoPP=None",
                                  "chrY - female child"]
                if genotype_mother == "0/0" and genotype_father == "0/0" and genotype_self == "1/1":
                    conditions = ["autosome",
                                  "chrX - male child - novoPP=0", "chrX - male child - novoPP=None",
                                  "chrY - male child - novoPP=0", "chrY - male child - novoPP=None",
                                  "chrX - female child", "chrY - female child"]
                for condition in conditions:
                    chrom = condition.split(" ")[0]
                    variant["variant"]["CHROM"] = "1" if chrom == "autosome" else chrom[3] #chrXY
                    variant["samplegeno"][0]["samplegeno_sex"] = "female" if "female" in condition else "male"
                    variant["novoPP"] = 1 if "novocaller high" in condition else 0
                    if "novoPP=None" in condition:
                        variant.pop("novoPP")
                    variant["samplegeno"][0]["samplegeno_numgt"] = genotype_self
                    variant["samplegeno"][1]["samplegeno_numgt"] = genotype_mother
                    variant["samplegeno"][2]["samplegeno_numgt"] = genotype_father

                    result = inheritance_mode.inheritance_mode(variant)
                    result["genotype_mother"] = genotype_mother
                    result["genotype_father"] = genotype_father
                    result["genotype_self"] = genotype_self
                    result["condition"] = condition
                    for role in ["mother", "father", "self"]:
                        result["genotype_label_" + role] = result["genotype_label"][role]
                    result.pop("genotype_label")
                    rows.append(result)

    dataframe = pd.DataFrame(rows, columns=["genotype_mother", "genotype_father",
                                            "genotype_self", "condition",
                                            "genotype_label_mother", "genotype_label_father",
                                            "genotype_label_self", "inheritance_modes"])
    dataframe.to_csv(path.join(DATA_DIR, "all_scenarios_table.tsv"), sep="\t", index=False)

def NA12879_table(for_excel=False):
    """
    call inheritance_mode for each variant in "<DATA_DIR>/variants_NA12879.json"
    and create a table for output of inheritance_mode, output as tsv:
    "<DATA_DIR>NA12879_genotype_table.tsv"

    1/1 becomes Jan 1 in excel.
    Let's just add a "'" in front of gt columns if for_excel is true.
    """
    sample = "NA12879"
    filename_variant = path.join(DATA_DIR, "variants_%s.json" % sample)
    file_variant = open(filename_variant) ## let it raise an error if file is missing.
    variants = json.load(file_variant)

    role_of_id = {
        'NA12878_sample': "mother",
        'NA12877_sample': "father",
        'NA12879_sample': "self"
    }

    column_order = ["GT_mother", "GT_father", "GT_self", "chrom", "novoPP",
                    "AD_mother", "AD_father", "AD_self", "title",
                    "GT_label_mother", "GT_label_father",
                    "GT_label_self", "inheritance_modes"]

    rows = []
    for variant in variants:
        sample_geno = variant.get("samplegeno")

        for i in range(len(sample_geno)):
            role = role_of_id[sample_geno[i]['samplegeno_sampleid']]
            sample_geno[i]["samplegeno_role"] = role
            sample_geno[i]["samplegeno_sex"] = "male" if role == "father" else "female"
        inh_mod_result = inheritance_mode.inheritance_mode(variant)

        output = {}
        for isample_geno in sample_geno:
            role = role_of_id[isample_geno['samplegeno_sampleid']]
            output["GT_" + role] = isample_geno["samplegeno_numgt"]
            output["AD_" + role] = isample_geno["samplegeno_ad"]
            output["GT_label_" + role] = inh_mod_result["genotype_label"][role]
        chrom = "chr" + variant.get("variant", {}).get("CHROM")
        if not chrom in ["chrX", "chrY"]:
            chrom = "autosome"
        output["chrom"] = chrom
        output["novoPP"] = variant.get("novoPP")
        output["title"] = variant["variant"]["display_title"]
        output["inheritance_modes"] = inh_mod_result["inheritance_modes"]

        rows.append(output)

    dataframe = pd.DataFrame(rows, columns=column_order)
    sort_by = ['chrom', 'GT_mother', 'GT_father', 'GT_self', 'novoPP']
    sort_ascending = [field != 'novoPP' for field in sort_by]
    dataframe.sort_values(by=sort_by, ascending=sort_ascending, inplace=True)

    if for_excel:
        for field in ['GT_mother', 'GT_father', 'GT_self', 'AD_mother', 'AD_father', 'AD_self']:
            dataframe[field] = "'" + dataframe[field]

    dataframe.to_csv(path.join(DATA_DIR, "NA12879_genotype_table.tsv"), sep="\t", index=False)

if __name__ == '__main__':
    all_scenarios_table()
    NA12879_table()
