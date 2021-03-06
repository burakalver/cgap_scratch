"""
Create a report for all variants found on cgap test.
"""

import json
from os import path, makedirs, remove
from urllib.parse import urlencode
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dcicutils import ff_utils, diff_utils 
from my_utils import nested_keys
from my_utils.Rutils import unique
from my_utils.pd_utils import print_full
import base


DATA_DIR = path.join(base.ROOT_DIR, "data")
FNAME_MAPPING_TABLE_VARIANT = "VCF Mapping Table - v0.4.8 variant table.tsv"
FNAME_MAPPING_TABLE_GENE = "VCF Mapping Table - GeneTable v0.4.6.tsv"

KEYNAME = "cgaptest"
VCF_FILE = "GAPFI2VBKGM7"


def search_result(params):
    '''
    Assuming the <KEYNAME> in the <keyfilename> is a valid admin key for cgapwolf.
    Perform a search based on params, e.g. {"type": "Gene"} and return result.
    '''
    keyfilename = path.expanduser("~") + '/keypairs.json'
    with open(keyfilename) as keyfile:
        keys = json.load(keyfile)
        key = keys[KEYNAME]
    base_url = "/search/"
    query = "%s?%s" % (base_url, urlencode(params))
    result = ff_utils.search_metadata(query, key=key)
    return result


def load_variants(sample="NA12879"):
    '''
    Search response for variants from VCF_FILE for <sample>
    should be in "DATA_DIR/variants_<sample>.json"

    If not, make it be.

    and return response
    '''
    filename_variant = path.join(DATA_DIR, "variants_%s.json" % sample) #variants file name
    if path.exists(filename_variant):
        with open(filename_variant) as file_variant:
            variants = json.load(file_variant)
    else:
        params = {"type": "VariantSample",
                  "file" : VCF_FILE}
        if sample != "all":
            params.update({"CALL_INFO" : '%s_sample' % sample})
        variants = search_result(params)
        with open(filename_variant, "w") as file_variant:
            json.dump(variants, file_variant)

    return variants


def load_genes():
    '''
    Search response for genes on cgapwolf should be in"DATA_DIR/genes.json"
    If not, make it be.
    and return response
    '''
    filename_gene = path.join(DATA_DIR, "genes.json") #genes file name
    if path.exists(filename_gene):
        with open(filename_gene) as file_gene:
            genes = json.load(file_gene)
    else:
        genes = search_result(params={"type": "Gene"})
        with open(filename_gene, "w") as file_gene:
            json.dump(genes, file_gene,indent=4)

    return genes


def get_values_pervar(search_results, field):
    '''
    if search_results=[{a1:1,a1:2},{a1:3}]
    and field = a1
    this function will return [[1,2],[3]]
    '''
    return [nested_keys.get_field_by_nested_key(search_result, field)
            for search_result in search_results]


def get_values_all(search_results, field):
    '''
    if search_results=[{a1:1,a1:2},{a1:3}]
    and field = a1
    this function will return [1,2,3]
    '''
    return np.array(nested_keys.get_field_by_nested_key(search_results, field))


def field_report(search_results, field, image_path=None):
    '''
    search_results is a search response for variantsamples or genes
    field is field name in snovault format, e.g. variant.ID
    if image_path is provided, a graph will be generated.
    returns a dict summary report.
    '''
    item_type = "variant" if "variant" in search_results[0].keys() else "gene"

    # How many search_results have how many entries for field.
    pervar_values = get_values_pervar(search_results, field)
    nvars = len(pervar_values)
    pervar_len = np.array([len(val) for val in pervar_values])
    pervar_lenmax = np.max(pervar_len)
    n_withvalue = np.count_nonzero(pervar_len)

    # If we collapse all search_results, how many entries are there? how many unique values.
    values = get_values_all(search_results, field)
    if values.dtype.name == "bool":
        values = values.astype("str")
    n_values = len(values)
    n_values_unique = len(set(values))
    subtext2 = "%d values for %d %ss" % (n_values, n_withvalue, item_type)

    stat_res_with_value_field = "%ss with Value" % item_type

    if n_values == 0:
        return {
            stat_res_with_value_field: 0,
            "Number of Values": 0,
            "Number of Unique Values": 0,
            "Example1": [],
            "Example2": []
        }

    #create summary statistics plot
    if image_path is not None:
        print("creating image for field: %s" % field)

        fig = plt.figure(figsize=(8, 10))
        fig.suptitle(field, fontsize=14, fontweight='bold')

        marleft = 0.1
        marright = 0.05
        marglobaltop = 0.05
        martop = 0.05
        marbottom = 0.05
        marleftstring = 0.3

        xsize = 1-marleft-marright
        ysize = (1 - 2*marbottom - 2*martop - marglobaltop)/2
        xsizestring = 1-marleftstring-marright

        #plot distribution of n_withvalue
        axes = plt.axes([marleft, 1-marglobaltop-martop-ysize, xsize, ysize])
        axes.hist(pervar_len, pervar_lenmax+1, (-0.5, pervar_lenmax+0.5), log=True)
        axes.set_title("Number of values per %s" % item_type)
        subtext1 = "%d / %d %ss have the field filled." % (n_withvalue, nvars, item_type)
        axes.text(0.01, 0.95, subtext1, transform=axes.transAxes)
        axes.set_ylim([1, axes.get_ylim()[1]*1.2])
        axes.set_xlabel('array length per %s' % item_type)
        axes.set_ylabel('number of %ss' % item_type)

        #plot distribution of value itself
        numeric = values.dtype.name == "float64" or values.dtype.name == "int64"

        if numeric:
            axes = plt.axes([marleft, marbottom, xsize, ysize])
            if values.dtype.name == "float64":
                axes.hist(values, 50, log=True)
            elif values.dtype.name == "int64":
                val_min = values.min()
                val_max = values.max()
                binsize = np.power(2, np.max([0, np.ceil(np.log2(val_max - val_min+1))-7]))
                bin_min = (round(np.round(val_min/binsize))-0.5)*binsize
                bin_max = (round(np.round(val_max/binsize))+0.5)*binsize
                nbins = round((bin_max-bin_min)/binsize)
                axes.hist(values, nbins, [bin_min, bin_max], log=True)

            axes.set_title("Distribution of values")
            axes.set_xlabel('field value')
            axes.set_ylabel('number of %ss' % item_type)
            axes.text(0.01, 0.95, subtext2, transform=axes.transAxes)
            axes.set_ylim([axes.get_ylim()[0], axes.get_ylim()[1]*1.2])

        else:
            marleft = 0.3
            axes = plt.axes([marleftstring, marbottom, xsizestring, ysize])

            uvalues, counts = np.unique(values, return_counts=True)
            counts_series = pd.Series(counts, uvalues)
            counts_series = counts_series.sort_values(ascending=False)

            axes.set_title('Distribution of values')
            axes.set_xlabel('counts_series')
            axes.text(0.01, 0.95, subtext2, transform=axes.transAxes)

            if len(counts_series) > 20:
                axes.set_title('Top 15 Value Examples')
                axes.text(0.01, 0.90,
                          "examples from %d unique values" % len(counts_series),
                          transform=axes.transAxes)
                counts_series = counts_series[0:14]

            y_pos = np.arange(len(counts_series), 0, -1)
            axes.barh(y_pos, counts_series, align='center', alpha=0.4)
            plt.yticks(y_pos, counts_series.index)
            axes.set_ylim(0, len(counts_series)+3)

        fig.savefig(image_path)
        plt.close()

    #gather some example values
    example1 = []
    example2 = []
    for val in pervar_values:
        if len(val) == 0:
            continue
        if len(example1) == 0:
            example1 = val
            continue
        example2 = val
        break

    return {
        stat_res_with_value_field: n_withvalue,
        "Number of Values" : n_values,
        "Number of Unique Values" : n_values_unique,
        "Example1" : example1,
        "Example2" : example2
    }


def load_clean_mapping_table(fname_mapping_table, do_import_Y_only = True):
    """
    Load mapping table, sort by no, filter down to do_import=Y,
    Add field_name, and set it as index.
    """
    map_table = pd.read_csv(path.join(DATA_DIR, fname_mapping_table), "\t", header=5)
    if do_import_Y_only:
        map_table = map_table[map_table["do_import"] == "Y"]
    map_table.sort_values(by="no", inplace=True)

    ## display_titles can appear twice, once because of link_to, once because they are already there.
    ## mark them for removal.
    map_table = map_table[~map_table["field_name"].str.endswith(".display_title")]
    
    addon3 = [".display_title"
              if pd.notna(x) else ""
              for x in map_table['links_to']]
    addon2 = [json.loads(x)["key"] + "."
              if pd.notna(x) else ""
              for x in map_table['sub_embedding_group']]
    addon1 = ["variant."
              if x == "variant" else ""
              for x in map_table["scope"]]
    map_table["field_name"] = ["".join(x) for x in zip(addon1, addon2,
                                                       list(map_table["field_name"]), addon3)]

    map_table.set_index("field_name", inplace=True)
    return map_table


def field_list_mapping_table(fname_mapping_table, fields=set()):
    """
    Get field list from mapping table.
    If fields is provided, append any of those that are missing in mapping table at the end.
    """
    map_table = load_clean_mapping_table(fname_mapping_table)
    fields_mapping_table = map_table.index.tolist()
    for field in fields:
        if field not in fields_mapping_table:
            fields_mapping_table.append(field)
    return fields_mapping_table


def create_all_reports(search_results, fields):
    """
    Gather stats about every field on variant or genes provided in search_result.
    Create an html.
    """
    report_dir = path.join(DATA_DIR, "report")
    image_dir = "images"

    image_dir_absolute = path.join(report_dir, image_dir)
    makedirs(image_dir_absolute, exist_ok=True)

    item_type = "variant" if "variant" in search_results[0].keys() else "gene"
    stat_res_with_value_field = "%ss with Value" % item_type

    stats = []
    for field in fields:
        print(field)
        image_path_relative = path.join(image_dir, "summary.%s.%s.png" % (item_type, field))
        image_path_absolute = path.join(report_dir, image_path_relative)
        if path.exists(image_path_absolute):
            image_path_absolute = None
        row = field_report(search_results, field, image_path=image_path_absolute)
        row["Stats"] = ""
        if row[stat_res_with_value_field] > 0:
            row["Stats"] = '<a href = "%s">link</a>' % image_path_relative
        row["field"] = field
        stats.append(row)

    column_order = ["field", stat_res_with_value_field,
                    "Number of Values",
                    "Number of Unique Values",
                    "Stats",
                    "Example1",
                    "Example2"]
    stats_table = pd.DataFrame(stats, columns=column_order)

    html = stats_table.to_html(render_links=True, escape=False, index=False)
    with open(path.join(report_dir, 'table_%s.html' % item_type), 'w') as outf:
        outf.write(html)


def delete_images_for_changed_fields(filename_old, filename_new, do_delete = False):
    differ = diff_utils.DiffManager()

    with open(filename_old) as file_variant:
        old = json.load(file_variant)

    with open(filename_new) as file_variant:
        new = json.load(file_variant)

    item_type = "variant" if "variant" in filename_new else "gene"
    print("comparing %s, %d from %s to %d from %s" % (item_type, len(old), filename_old, len(new), filename_new))
    comparison = differ.comparison(old, new)
    intermediate = unique([re.sub("\[(.*?)\]","",each_diff) for each_diff in comparison])
    intermediate = unique([re.sub(" :.*$","",each_diff) for each_diff in intermediate])
    fields = [re.sub("^.","",each_diff) for each_diff in intermediate]

    report_dir = path.join(DATA_DIR, "report")
    image_dir = "images"

    image_dir_absolute = path.join(report_dir, image_dir)
    makedirs(image_dir_absolute, exist_ok=True)


    for field in fields:
        image_path_relative = path.join(image_dir, "summary.%s.%s.png" % (item_type, field))
        image_path_absolute = path.join(report_dir, image_path_relative)
        if path.exists(image_path_absolute):
            if do_delete:
                print("deleting %s" % image_path_absolute)
                remove(image_path_absolute)
            else:
                print("will delete %s, if you rerun with do_delete=T" % image_path_absolute)
        else:
            print("the image file is already absent: %s" % image_path_absolute)


def create_links_report(search_results):
    """
    Gather links for fields with links.
    Create an html.
    """
    def html_link_for_value(link_base, value):
        link = link_base.replace("<ID>", value)
        return '<a href="%s">%s</a>' % (link, value)

    item_type = "variant" if "variant" in search_results[0].keys() else "gene"

    if item_type == "variant":
        fname_mapping_table = FNAME_MAPPING_TABLE_VARIANT
    else:
        fname_mapping_table = FNAME_MAPPING_TABLE_GENE

    map_table = load_clean_mapping_table(fname_mapping_table)
    map_table = map_table[pd.notnull(map_table["link"])]
    fields = map_table.index

    table_list = []

    n_examples = 5
    
    for field in fields:
        link_base = map_table.loc[field]["link"]
        pervar_values = get_values_pervar(search_results, field)
        examples = []
        for val in pervar_values:
            if len(val) == 0:
                continue
            examples.append(val)
            if len(examples) >= n_examples:
                break
        links = []
        for example in examples:
            links.append([html_link_for_value(link_base, value) for value in example])
        for i in range(len(links), n_examples):
            links.append([])
            
        row = {"field": field}
        for link_enum in enumerate(links):
            row["link%d" % link_enum[0]] = "<br>".join(link_enum[1])
        table_list.append(row)

    table = pd.DataFrame(table_list)
    html = table.to_html(render_links=True, escape=False, index=False)

    report_dir = path.join(DATA_DIR, "report")
    with open(path.join(report_dir, 'table_%s_link.html' % item_type), 'w') as outf:
        outf.write(html)


def identify_fields_with_no_values():
    """
    For every field in mapping table, check which ones have no values at the moment.
    output an updated do_import column with Not_yet values for missing fields.
    """

    for item_type in ["variant","gene"]:
        if item_type == "variant":
            sample = "NA12879"
            search_results = load_variants(sample)
            fname_mapping_table = FNAME_MAPPING_TABLE_VARIANT
        else:
            search_results = load_genes()
            fname_mapping_table = FNAME_MAPPING_TABLE_GENE

        map_table = load_clean_mapping_table(fname_mapping_table, do_import_Y_only=False)
        map_table = map_table["do_import"]
        for field in map_table.index:
            print(map_table.loc[field])
            if (map_table.loc[field]) != "Y":
                continue
            values = get_values_all(search_results, field)
            if(len(values)==0):
                map_table.loc[field] = "Not_yet"

        fname_csv = path.join(DATA_DIR, "do_import_%s.csv" % item_type)
        print("writing %s" % fname_csv)
        map_table.to_csv(fname_csv)


def sex_check():
    '''
    check that the father has lower coverage on chrX and too many hets on chrX
    '''
    variants = load_variants("NA12877")
    chrom = get_values_all(variants, 'variant.CHROM')
    allele_depth = get_values_all(variants, 'samplegeno.samplegeno_ad')
    sample = get_values_all(variants, 'samplegeno.samplegeno_sampleid')
    samples = ["NA12877_sample", "NA12878_sample", "NA12879_sample"]
    for i in range(len(chrom)):
        if samples[0] != sample[i*3] or samples[1] != sample[i*3+1] or samples[2] != sample[i*3+2]:
            raise Exception("sample order not good")
    coverage = [int(a.split("/")[0]) + int(a.split("/")[1]) for a in allele_depth]
    coverage_np = np.array(coverage).reshape(len(chrom), 3)
    coverage = pd.DataFrame(coverage_np)
    coverage.columns = samples
    coverage["chrom"] = chrom

    print(coverage.groupby('chrom').mean())

    fields_chosen = ['AD_ALT', 'AD_REF', 'variant.CHROM', 'GT']
    values = {field:get_values_all(variants, field) for field in fields_chosen}
    dataframe = pd.DataFrame(values)
    print_full(pd.crosstab(dataframe["GT"], dataframe["variant.CHROM"]))


def report_wrapper_variant():
    '''
    hardcoding variables for create_all_reports.
    '''
    sample = "NA12879"
    variants = load_variants(sample)
    fields = nested_keys.nested_keys(variants)

    fields = field_list_mapping_table(FNAME_MAPPING_TABLE_VARIANT, fields)
    create_all_reports(variants, fields)
    create_links_report(variants)


def report_wrapper_gene():
    '''
    hardcoding variables for create_all_reports.
    '''
    genes = load_genes()
    fields = nested_keys.nested_keys(genes)

    fields = field_list_mapping_table(FNAME_MAPPING_TABLE_GENE, fields)
    create_all_reports(genes, fields)
    create_links_report(genes)


def delete_images_wrapper(do_delete=False):
    sample="NA12879"
    filename_old = path.join(DATA_DIR, "variants_%s.json.bak" % sample)
    filename_new = path.join(DATA_DIR, "variants_%s.json" % sample)
    delete_images_for_changed_fields(filename_old, filename_new, do_delete)
        
    filename_old = path.join(DATA_DIR, "genes.json.bak")
    filename_new = path.join(DATA_DIR, "genes.json")
    delete_images_for_changed_fields(filename_old, filename_new, do_delete)


if __name__ == '__main__':
#    sex_check()
    VAL = 1
    delete_images_wrapper(False)
    report_wrapper_variant()
    report_wrapper_gene()
    identify_fields_with_no_values()
