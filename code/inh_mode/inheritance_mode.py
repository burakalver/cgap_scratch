"""
given family genotypes and other info from variant, calculate inheritance mode.
"""

# GENOTYPE_LABEL_DOT = "missing"
# GENOTYPE_LABEL_00 = "homozygus reference"
# GENOTYPE_LABEL_0M = "heterozygous"
# GENOTYPE_LABEL_MM = "homozygus alternate"
# GENOTYPE_LABEL_MN_KEYWORD = "multiallelic"
# GENOTYPE_LABEL_MN = "het alt/alt -  %s" % GENOTYPE_LABEL_MN_KEYWORD
# GENOTYPE_LABEL_MN_ADDON = " (%s in family)" % GENOTYPE_LABEL_MN_KEYWORD
# GENOTYPE_LABEL_0 = "hemizygous reference"
# GENOTYPE_LABEL_M = "hemizygous alternate"
# GENOTYPE_LABEL_FEMALE_CHRY = "-"
# GENOTYPE_LABEL_SEX_INCONSISTENT = "false"

GENOTYPE_LABEL_DOT = "missing"
GENOTYPE_LABEL_00 = "homo ref"
GENOTYPE_LABEL_0M = "het"
GENOTYPE_LABEL_MM = "homo alt"
GENOTYPE_LABEL_MN_KEYWORD = "multiallelic"
GENOTYPE_LABEL_MN = "het alt/alt -  %s" % GENOTYPE_LABEL_MN_KEYWORD
GENOTYPE_LABEL_MN_ADDON = " (%s in family)" % GENOTYPE_LABEL_MN_KEYWORD
GENOTYPE_LABEL_0 = "hemi ref"
GENOTYPE_LABEL_M = "hemi alt"
GENOTYPE_LABEL_FEMALE_CHRY = "-"
GENOTYPE_LABEL_SEX_INCONSISTENT = "false"

INHMODE_LABEL_DE_NOVO_STRONG = "de novo (strong)"
INHMODE_LABEL_DE_NOVO_MEDIUM = "de novo (medium)"
INHMODE_LABEL_DE_NOVO_WEAK = "de novo (weak)"
INHMODE_LABEL_DE_NOVO_CHRXY = "de novo (chrXY)"
INHMODE_DOMINANT_FATHER = "dominant (paternal)"
INHMODE_DOMINANT_MOTHER = "dominant (maternal)"
INHMODE_LABEL_RECESSIVE = "recessive"
INHMODE_LABEL_X_LINKED_RECESSIVE_MOTHER = "X-linked recessive (Maternal)"
INHMODE_LABEL_X_LINKED_DOMINANT_MOTHER = "X-linked dominant (Maternal)"
INHMODE_LABEL_X_LINKED_DOMINANT_FATHER = "X-linked dominant (Paternal)"
INHMODE_LABEL_Y_LINKED = "Y-linked dominant"
INHMODE_LABEL_LOH = "Loss of Heteozyogousity"

INHMODE_LABEL_NONE_DOT = "Low relevance, missing call(s) in family"
INHMODE_LABEL_NONE_MN = "Low relevance, multiallelic site family"
INHMODE_LABEL_NONE_SEX_INCONSISTENT = "Low relevance, mismatching chrXY genotype(s)"
INHMODE_LABEL_NONE_LOWDEPTH = "Low relevance, low depth"
INHMODE_LABEL_NONE_HOMOZYGOUS_PARENT = "Low relevance, homozygous in a parent"
INHMODE_LABEL_NONE_BOTH_PARENTS = "Low relevance, present in both parent(s)"
INHMODE_LABEL_NONE_OTHER = "Low relevance, other"

def multiallelic_site(genotype):
    '''
    boolean for "any genotype has allele number 2 or greater"
    '''
    for igenotype in genotype.values():
        genotype_allele1, genotype_allele2 = igenotype.split("/")
        if genotype_allele1 == ".":
            continue
        if int(genotype_allele1) > 1:
            return True
        if int(genotype_allele2) > 1:
            return True
    return False


def genotype_to_genotype_label_single(igenotype, isex, chrom):
    '''
    input is the genotype (e.g. 0/0) and sex (e.g. male) of one role
    output is the matching GENOTYPE_LABEL
    '''
    genotype_allele1, genotype_allele2 = igenotype.split("/")

    if genotype_allele1 == ".":
        if isex == "female" and chrom == "chrY":
            return GENOTYPE_LABEL_FEMALE_CHRY
        return GENOTYPE_LABEL_DOT

    genotype_allele1 = int(genotype_allele1)
    genotype_allele2 = int(genotype_allele2)

    if isex == "female" and chrom == "chrY":
        if genotype_allele1 == 0 and genotype_allele2 == 0:
            return GENOTYPE_LABEL_FEMALE_CHRY
        return GENOTYPE_LABEL_SEX_INCONSISTENT

    if isex == "male" and chrom in ["chrX", "chrY"]:
        if genotype_allele1 != genotype_allele2:
            return GENOTYPE_LABEL_SEX_INCONSISTENT
        if genotype_allele1 == 0:
            return GENOTYPE_LABEL_0
        return GENOTYPE_LABEL_M

    if genotype_allele1 == 0 and genotype_allele2 == 0:
        return GENOTYPE_LABEL_00
    if genotype_allele1 == 0:
        return GENOTYPE_LABEL_0M
    if genotype_allele1 == genotype_allele2:
        return GENOTYPE_LABEL_MM
    return GENOTYPE_LABEL_MN


def genotype_to_genotype_label_family(genotype, sex, chrom):
    '''
    input:
    genotype of family, e.g. genotype["self"]=0/1
    sex of family, e.g. sex["self"]="male"
    chrom, e.g. chrY

    output:
    genotype_label of family, e.g. genotype["self"]=GENOTYPE_LABEL_0M

    if the site is determined to be a multiallelic site,
    genotype_labels without GENOTYPE_LABEL_MN_KEYWORD label get GENOTYPE_LABEL_MN_ADDON
    '''

    roles = genotype.keys()
    genotype_label = {}
    for role in roles:
        genotype_label[role] = genotype_to_genotype_label_single(genotype[role], sex[role], chrom)

    if multiallelic_site(genotype):
        for role, igenotype_label in genotype_label.items():
            if not GENOTYPE_LABEL_MN_KEYWORD in igenotype_label:
                genotype_label[role] += GENOTYPE_LABEL_MN_ADDON
    return genotype_label


def inheritance_modes_trio(genotype, genotype_label, sex, chrom, novoPP):
    '''
    infer inheritance mode for trio based on genotypes, sexes, and chrom
    if novoPP has called a de novo, that takes precedence.
    '''
    roles = genotype.keys()

    if "father" not in roles or "mother" not in roles or "self" not in roles:
        return []

    if GENOTYPE_LABEL_DOT in genotype_label.values():
        return []
    if multiallelic_site(genotype):
        return []
    if GENOTYPE_LABEL_SEX_INCONSISTENT in genotype_label.values():
        return []

    if novoPP > 0.9:
        return [INHMODE_LABEL_DE_NOVO_STRONG]
    if novoPP > 0.1:
        return [INHMODE_LABEL_DE_NOVO_MEDIUM]

    if (genotype["mother"] == "0/0" and genotype["father"] == "0/0"
            and genotype["self"] == "0/1" and chrom == "autosome"):
        return [INHMODE_LABEL_DE_NOVO_WEAK]

    if (genotype["mother"] == "0/0" and genotype["father"] == "0/0"
            and ((genotype["self"] == "0/1" and sex["self"] == "female" and chrom == "chrX")
                 or (genotype["self"] == "1/1" and sex["self"] == "male" and chrom != "autosome"))):
        if novoPP == 0:
            return [INHMODE_LABEL_DE_NOVO_WEAK]
        if novoPP == -1:
            return [INHMODE_LABEL_DE_NOVO_CHRXY]
        raise ValueError("novoPP is different from 0 or -1 on sex chromosome: " + str(novoPP))

    if (genotype["mother"] == "0/0"
            and genotype_label["father"] == GENOTYPE_LABEL_0M
            and genotype["self"] == "0/1"):
        return [INHMODE_DOMINANT_FATHER]

    if (genotype["mother"] == "0/1" and genotype["father"] == "0/0"
            and genotype["self"] == "0/1"):
        return [INHMODE_DOMINANT_MOTHER]

    if (genotype["mother"] == "0/1" and genotype["father"] == "0/1"
            and genotype["self"] == "1/1"):
        return ["recessive"]

    if (genotype["mother"] == "0/1" and genotype["father"] == "0/0"
            and genotype["self"] == "1/1" and sex["self"] == "male" and chrom == "chrX"):
        return [INHMODE_LABEL_X_LINKED_RECESSIVE_MOTHER, INHMODE_LABEL_X_LINKED_DOMINANT_MOTHER]

    if (genotype["mother"] == "0/0" and genotype_label["father"] == GENOTYPE_LABEL_M and
            chrom == "chrX" and genotype_label["self"] in [GENOTYPE_LABEL_M, GENOTYPE_LABEL_0M]):
        return [INHMODE_LABEL_X_LINKED_DOMINANT_FATHER]

    if (genotype_label["father"] == GENOTYPE_LABEL_M and
            chrom == "chrY" and genotype_label["self"] == GENOTYPE_LABEL_M):
        return [INHMODE_LABEL_Y_LINKED]

    if (((genotype["mother"] == "0/1" and genotype["father"] == "0/0") or
         (genotype["mother"] == "0/0" and genotype["father"] == "0/1"))
            and genotype["self"] == "1/1"):
        return [INHMODE_LABEL_LOH]

    return []


def inheritance_modes_other_labels(genotype, genotype_label):
    '''
    we could leave inheritance mode field empty if none of the above fit.
    But maybe let's add an explainer value for why it is empty.
    '''
    roles = genotype.keys()

    #proband-only or duos can be added later
    if "father" not in roles or "mother" not in roles or "self" not in roles:
        return []

    if GENOTYPE_LABEL_DOT in genotype_label.values():
        return [INHMODE_LABEL_NONE_DOT]
    if multiallelic_site(genotype):
        return [INHMODE_LABEL_NONE_MN]
    if GENOTYPE_LABEL_SEX_INCONSISTENT in genotype_label.values():
        return [INHMODE_LABEL_NONE_SEX_INCONSISTENT]

    if genotype["mother"] == "1/1" or (
            genotype["father"] == "1/1" and genotype_label["father"] != GENOTYPE_LABEL_M):
        return [INHMODE_LABEL_NONE_HOMOZYGOUS_PARENT]
    if ((genotype["mother"] == "1/1" or genotype["mother"] == "0/1") and
            (genotype["father"] == "1/1" or genotype["father"] == "0/1")):
        return [INHMODE_LABEL_NONE_BOTH_PARENTS]

    return [INHMODE_LABEL_NONE_OTHER]


def inheritance_modes_cmphet(cmphet):
    '''
    summarize Comphound Het Caller results
    '''
    inheritance_modes_set = set()
    if cmphet is not None:
        for icmphet in cmphet:
            iphase = icmphet.get("comhet_phase")
            iimpact = icmphet.get("comhet_impact_gene")
            cmphet_str = "Compound Het (%s/%s)" % (iphase, iimpact.lower())
            inheritance_modes_set.add(cmphet_str)
    inheritance_modes = list(inheritance_modes_set)
    inheritance_modes.sort() # to make results reproducible
    return inheritance_modes


def inheritance_mode(variant):
    """
    variant is a sampleVariant item, including the samplegeno_role and samplegeno_sex fields.

    variantSample = {
        "samplegeno": [{
            "samplegeno_role": "self", # self must exist,
                                       # mother and father are also used if available
            "samplegeno_numgt": "0/1", # ./. or x/y where x and y are integers.
                                       # ref:0, alt1:1, alt2:2 etc.
            "samplegeno_sex": "male",  # male or female
            },{},{}],
        "variant": {
            "CHROM": 1
            },
        "novoPP": 0.5     # optional, must be float if exist.
        "cmphet": [{      # optional, but if values exist, they must include the keys:
            "comhet_phase": "Phased",           # "Phased" or Unphased.
            "comhet_impact_gene": "STRONG_PAIR" # STRONG, MEDIUM, OR WEAK
        }]

    returns genotypes for each role as a string and a list of possible inheritance_modes.
    {
     "genotype_label": {
     "self: "homo ref",
     "father: "homo ref",
     ... for each role
    },
    "inheritance_modes" : [de novo (weak)]}
        "genotype_label_father" : homo ref
        "genotype_label_proband" : het alt

    sex-mismatched genotypes get a "false" value, e.g. 0/1 on chrX for a male
    """

    sample_geno = variant.get("samplegeno")
    genotype = {s["samplegeno_role"]: s["samplegeno_numgt"] for s in sample_geno}
    sex = {s["samplegeno_role"]: s["samplegeno_sex"]  for s in sample_geno}
    roles = genotype.keys()

    chrom = "chr" + variant.get("variant", {}).get("CHROM")
    if not chrom in ["chrX", "chrY"]:
        chrom = "autosome"

    cmphet = variant.get("cmphet")
    novoPP = variant.get("novoPP", -1)

    if "self" not in roles:
        raise ValueError('variant["samplegeno"]["samplegeno_role"]="self" is missing')

    genotype_label = genotype_to_genotype_label_family(genotype, sex, chrom)
    inheritance_modes = inheritance_modes_trio(genotype, genotype_label, sex, chrom, novoPP)
    inheritance_modes += inheritance_modes_cmphet(cmphet)
    if len(inheritance_modes) == 0:
        inheritance_modes = inheritance_modes_other_labels(genotype, genotype_label)

    result = {
        "genotype_label": genotype_label,
        "inheritance_modes" : inheritance_modes
    }

    ## results are not valid for mitochondria, so wipe them.
    ## this could be moved up, but I wanted to partially get structure of results from above.
    if chrom == "chrM":
        for role in result["genotype_label"].keys():
            result["genotype_label"][role] = ""
        result["inheritance_modes"] = []

    return result
