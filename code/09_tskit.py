### ============================================================================
### phylo2
### Inferring genomic history based on ancestral recombination graph (tskit)
### ============================================================================

### Preparations
# python ts_samples.py


### ============================================================================
### 1. Create sample file (tsinfer)

#!/usr/bin/python

import tsinfer
import cyvcf2
import json


## Define functions
def add_diploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    # You may want to change the following line, e.g. here we allow * (a spanning
    # deletion) to be a valid allele state
    allele_chars = set("ATGCatgc*")
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        if pos == variant.POS:
            print(f"Duplicate entries at position {pos}, ignoring all but the first")
            continue
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF.upper()] + [v.upper() for v in variant.ALT]
        ancestral = variant.INFO.get("AA", ".")  # "." means unknown
        # some VCFs (e.g. from 1000G) have many values in the AA field: take the 1st
        ancestral = ancestral.split("|")[0].upper()
        if ancestral == "." or ancestral == "":
            # use the reference as ancestral, if unknown (NB: you may not want this)
            ancestral = variant.REF.upper()
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        # Check we have ATCG alleles
        for a in ordered_alleles:
            if len(set(a) - allele_chars) > 0:
                print(f"Ignoring site at pos {pos}: allele {a} not in {allele_chars}")
                continue
        allele_index = {
            old_index: ordered_alleles.index(a) for old_index, a in enumerate(alleles)
        }
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [
            allele_index[old_index]
            for row in variant.genotypes
            for old_index in row[0:2]
        ]
        samples.add_site(pos, genotypes=genotypes, alleles=ordered_alleles)

def chromosome_length(vcf):
    assert len(vcf.seqlens) == 1
    return vcf.seqlens[0]

def add_populations(vcf, samples):
    """
    Add tsinfer Population objects and returns a list of IDs corresponding to the VCF samples.
    """
    samples_pop = [sample_name[-6:] for sample_name in vcf.samples]
    pop_lookup = {}
    pop_lookup["abearc"] = samples.add_population(metadata={"population": "abearc"})
    pop_lookup["abeare"] = samples.add_population(metadata={"population": "abeare"})
    pop_lookup["abebel"] = samples.add_population(metadata={"population": "abebel"})
    pop_lookup["abeboc"] = samples.add_population(metadata={"population": "abeboc"})
    pop_lookup["abegun"] = samples.add_population(metadata={"population": "abegun"})
    pop_lookup["abehon"] = samples.add_population(metadata={"population": "abehon"})
    pop_lookup["abepri"] = samples.add_population(metadata={"population": "abepri"})
    pop_lookup["abequi"] = samples.add_population(metadata={"population": "abequi"})
    pop_lookup["abesan"] = samples.add_population(metadata={"population": "abesan"})
    pop_lookup["affboc"] = samples.add_population(metadata={"population": "affboc"})
    pop_lookup["affgun"] = samples.add_population(metadata={"population": "affgun"})
    pop_lookup["atlliz"] = samples.add_population(metadata={"population": "atlliz"})
    pop_lookup["atltam"] = samples.add_population(metadata={"population": "atltam"})
    pop_lookup["casliz"] = samples.add_population(metadata={"population": "casliz"})
    pop_lookup["chlbar"] = samples.add_population(metadata={"population": "chlbar"})
    pop_lookup["chlhon"] = samples.add_population(metadata={"population": "chlhon"})
    pop_lookup["chlpri"] = samples.add_population(metadata={"population": "chlpri"})
    pop_lookup["ecoarc"] = samples.add_population(metadata={"population": "ecoarc"})
    pop_lookup["esparc"] = samples.add_population(metadata={"population": "esparc"})
    pop_lookup["floarc"] = samples.add_population(metadata={"population": "floarc"})
    pop_lookup["floflk"] = samples.add_population(metadata={"population": "floflk"})
    pop_lookup["gemarc"] = samples.add_population(metadata={"population": "gemarc"})
    pop_lookup["gemare"] = samples.add_population(metadata={"population": "gemare"})
    pop_lookup["gemflk"] = samples.add_population(metadata={"population": "gemflk"})
    pop_lookup["gumboc"] = samples.add_population(metadata={"population": "gumboc"})
    pop_lookup["gumhon"] = samples.add_population(metadata={"population": "gumhon"})
    pop_lookup["gutbar"] = samples.add_population(metadata={"population": "gutbar"})
    pop_lookup["gutbel"] = samples.add_population(metadata={"population": "gutbel"})
    pop_lookup["gutboc"] = samples.add_population(metadata={"population": "gutboc"})
    pop_lookup["gutgun"] = samples.add_population(metadata={"population": "gutgun"})
    pop_lookup["guthon"] = samples.add_population(metadata={"population": "guthon"})
    pop_lookup["gutpri"] = samples.add_population(metadata={"population": "gutpri"})
    pop_lookup["gutsan"] = samples.add_population(metadata={"population": "gutsan"})
    pop_lookup["indbel"] = samples.add_population(metadata={"population": "indbel"})
    pop_lookup["indgun"] = samples.add_population(metadata={"population": "indgun"})
    pop_lookup["indhon"] = samples.add_population(metadata={"population": "indhon"})
    pop_lookup["indpri"] = samples.add_population(metadata={"population": "indpri"})
    pop_lookup["indqui"] = samples.add_population(metadata={"population": "indqui"})
    pop_lookup["libhai"] = samples.add_population(metadata={"population": "libhai"})
    pop_lookup["maybel"] = samples.add_population(metadata={"population": "maybel"})
    pop_lookup["nigbel"] = samples.add_population(metadata={"population": "nigbel"})
    pop_lookup["nigboc"] = samples.add_population(metadata={"population": "nigboc"})
    pop_lookup["nigflk"] = samples.add_population(metadata={"population": "nigflk"})
    pop_lookup["niggun"] = samples.add_population(metadata={"population": "niggun"})
    pop_lookup["nighon"] = samples.add_population(metadata={"population": "nighon"})
    pop_lookup["nigqui"] = samples.add_population(metadata={"population": "nigqui"})
    pop_lookup["prohon"] = samples.add_population(metadata={"population": "prohon"})
    pop_lookup["prosan"] = samples.add_population(metadata={"population": "prosan"})
    pop_lookup["pueala"] = samples.add_population(metadata={"population": "pueala"})
    pop_lookup["puearc"] = samples.add_population(metadata={"population": "puearc"})
    pop_lookup["pueare"] = samples.add_population(metadata={"population": "pueare"})
    pop_lookup["puebar"] = samples.add_population(metadata={"population": "puebar"})
    pop_lookup["puebel"] = samples.add_population(metadata={"population": "puebel"})
    pop_lookup["pueboc"] = samples.add_population(metadata={"population": "pueboc"})
    pop_lookup["pueflk"] = samples.add_population(metadata={"population": "pueflk"})
    pop_lookup["puegun"] = samples.add_population(metadata={"population": "puegun"})
    pop_lookup["puehon"] = samples.add_population(metadata={"population": "puehon"})
    pop_lookup["puepri"] = samples.add_population(metadata={"population": "puepri"})
    pop_lookup["puequi"] = samples.add_population(metadata={"population": "puequi"})
    pop_lookup["ranarc"] = samples.add_population(metadata={"population": "ranarc"})
    pop_lookup["ranbel"] = samples.add_population(metadata={"population": "ranbel"})
    pop_lookup["ranhon"] = samples.add_population(metadata={"population": "ranhon"})
    pop_lookup["ransan"] = samples.add_population(metadata={"population": "ransan"})
    pop_lookup["tabhon"] = samples.add_population(metadata={"population": "tabhon"})
    pop_lookup["tanbar"] = samples.add_population(metadata={"population": "tanbar"})
    pop_lookup["tanpri"] = samples.add_population(metadata={"population": "tanpri"})
    pop_lookup["tigboc"] = samples.add_population(metadata={"population": "tigboc"})
    pop_lookup["tighon"] = samples.add_population(metadata={"population": "tighon"})
    pop_lookup["torboc"] = samples.add_population(metadata={"population": "torboc"})
    pop_lookup["unibel"] = samples.add_population(metadata={"population": "unibel"})
    pop_lookup["uniboc"] = samples.add_population(metadata={"population": "uniboc"})
    pop_lookup["uniflk"] = samples.add_population(metadata={"population": "uniflk"})
    pop_lookup["unihon"] = samples.add_population(metadata={"population": "unihon"})
    pop_lookup["unipri"] = samples.add_population(metadata={"population": "unipri"})
    pop_lookup["uniqui"] = samples.add_population(metadata={"population": "uniqui"})
    return [pop_lookup[pop] for pop in samples_pop]

def add_species(vcf, samples):
    """
    Add tsinfer Population objects and returns a list of IDs corresponding to the VCF samples.
    """
    samples_pop = [sample_name[-6:-3] for sample_name in vcf.samples]
    pop_lookup = {}
    pop_lookup["abe"] = samples.add_population(metadata={"species": "aberrans"})
    pop_lookup["aff"] = samples.add_population(metadata={"species": "affinis"})
    pop_lookup["atl"] = samples.add_population(metadata={"species": "atlahua"})
    pop_lookup["cas"] = samples.add_population(metadata={"species": "castroaguirrei"})
    pop_lookup["chl"] = samples.add_population(metadata={"species": "chlorurus"})
    pop_lookup["eco"] = samples.add_population(metadata={"species": "ecosur"})
    pop_lookup["esp"] = samples.add_population(metadata={"species": "espinozaperezei"})
    pop_lookup["flo"] = samples.add_population(metadata={"species": "floridae"})
    pop_lookup["gem"] = samples.add_population(metadata={"species": "gemma"})
    pop_lookup["gum"] = samples.add_population(metadata={"species": "gummigutta"})
    pop_lookup["gut"] = samples.add_population(metadata={"species": "guttavarius"})
    pop_lookup["ind"] = samples.add_population(metadata={"species": "indigo"})
    pop_lookup["lib"] = samples.add_population(metadata={"species": "liberte"})
    pop_lookup["may"] = samples.add_population(metadata={"species": "maya"})
    pop_lookup["nig"] = samples.add_population(metadata={"species": "nigricans"})
    pop_lookup["pro"] = samples.add_population(metadata={"species": "providencianus"})
    pop_lookup["pue"] = samples.add_population(metadata={"species": "puella"})
    pop_lookup["ran"] = samples.add_population(metadata={"species": "randallorum"})
    pop_lookup["tan"] = samples.add_population(metadata={"species": "tan"})
    pop_lookup["uni"] = samples.add_population(metadata={"species": "unicolor"})
    pop_lookup["tab"] = samples.add_population(metadata={"species": "tabacarius"})
    pop_lookup["tig"] = samples.add_population(metadata={"species": "tigrinus"})
    pop_lookup["tor"] = samples.add_population(metadata={"species": "tortugarum"})
    return [pop_lookup[pop] for pop in samples_pop]

def add_clades(vcf, samples):
    """
    Add tsinfer Population objects and returns a list of IDs corresponding to the VCF samples.
    """
    populations = [sample_name[-6:] for sample_name in vcf.samples]
    for i in range(len(populations)):
        if populations[i] in ["atlliz", "atltam", "ecoarc", "floarc", "floflk"]:
            populations[i] = "small"
        elif populations[i] in ["abearc", "abeare", "abebel", "abeboc", "abegun", "abehon", "abepri", "abequi", \
            "abesan", "affboc", "affgun", "casliz", "chlbar", "chlhon", "chlpri", "esparc", "gemarc", "gemare", \
            "gemflk", "gumboc", "gumhon", "gutbar", "gutbel", "gutboc", "gutgun", "guthon", "gutpri", "gutsan", \
            "indbel", "indgun", "indhon", "indpri", "indqui", "libhai", "maybel", "nigbel", "nigboc", "nigflk", \
            "niggun", "nighon", "nigqui", "prohon", "prosan", "pueala", "puearc", "pueare", "puebar", "puebel", \
            "pueboc", "pueflk", "puegun", "puehon", "puepri", "puequi", "ranarc", "ranbel", "ranhon", "ransan", \
            "tanbar", "tanpri", "unibel", "uniboc", "uniflk", "unihon", "unipri", "uniqui"]:
            populations[i] = "large"
        elif populations[i] in ["tabhon", "tigboc", "tighon", "torboc"]:
            populations[i] = "out"
    pop_lookup = {}
    pop_lookup["small"] = samples.add_population(metadata = {"clade": "small"})
    pop_lookup["large"] = samples.add_population(metadata = {"clade": "large"})
    pop_lookup["out"] = samples.add_population(metadata = {"clade": "out"})
    return [pop_lookup[pop] for pop in populations]

def add_diploid_individuals(vcf, samples, populations):
    for name, population in zip(vcf.samples, populations):
        samples.add_individual(ploidy=2, metadata={"name": name}, population=population)


## ----------------------------------------------------------------------------
## Load data and add metadata (populations / species / clades)

## Specify dataset
dataset = "phylo2e-S"
level = "pop"  # pop | sp | cld

# l = ["LG03"]
l = list([f"LG{i:02}" for i in range(1, 4)])
for region in l:

    vcf = cyvcf2.VCF("../3_vcf/" + dataset + "_"+ region + "_n1.vcf")

    if level == "pop":
        with tsinfer.SampleData(
            path = dataset + "_" + region + "_" + level + "_n1.samples", sequence_length = chromosome_length(vcf)
        ) as samples:
            populations = add_populations(vcf, samples)
            add_diploid_individuals(vcf, samples, populations)
            add_diploid_sites(vcf, samples)
    elif level == "sp":
        with tsinfer.SampleData(
            path = dataset + "_" + region + "_" + level + "_n1.samples", sequence_length = chromosome_length(vcf)
        ) as samples:
            populations = add_species(vcf, samples)
            add_diploid_individuals(vcf, samples, populations)
            add_diploid_sites(vcf, samples)
    elif level == "cld":
        with tsinfer.SampleData(
            path = dataset + "_" + region + "_" + level + "_n1.samples", sequence_length = chromosome_length(vcf)
        ) as samples:
            populations = add_clades(vcf, samples)
            add_diploid_individuals(vcf, samples, populations)
            add_diploid_sites(vcf, samples)

    print(
        "Sample file created for {} samples ".format(samples.num_samples)
        + "({} individuals) ".format(samples.num_individuals)
        + "with {} variable sites.".format(samples.num_sites),
        flush = True,
    )


### ============================================================================
### 2. Infer and date tree sequence (tsinfer / tsdate)

#!/usr/bin/python

import tsinfer
import tsdate
import json

## Specify dataset
dataset = "phylo2e-S"
level = "pop"  # pop | sp | cld

# l = ["LG02"]
l = list([f"LG{i:02}" for i in range(1, 4)])
for region in l:

## Load ts samples file
    samples = tsinfer.load("../4_samples/" + dataset + "_" + region + "_" + level + "_n1.samples")


## Infer tree sequences
    ts = tsinfer.infer(samples)
    print(
        "Inferred tree sequence `ts`: {} trees over {} Mb".format(
            ts.num_trees, ts.sequence_length / 1e6
        )
    )
    ts.dump(dataset + "_" + region + "_" + level + "_n1.ts")

    # # Print metadata *** fix "clade"
    # for sample_node_id in ts.samples():
    #     individual_id = ts.node(sample_node_id).individual
    #     population_id = ts.node(sample_node_id).population
    #     print(
    #         "Node",
    #         sample_node_id,
    #         "is a",
    #         l,
    #         "segment sampled from individual",
    #         json.loads(ts.individual(individual_id).metadata),
    #         "of",
    #         json.loads(ts.population(population_id).metadata)["clade"],
    #         "clade"
    #     )


## ----------------------------------------------------------------------------
## Include dating

## Remove unary nodes
    tss = ts.simplify(keep_unary = False)
    print(
        "Removed unary nodes: `tss`, {} of {} nodes remaining".format(
            tss.num_nodes, ts.num_nodes
        )
    )
    #tss.dump(dataset + "_" + level + ".tss")


## Date trees
    tsd = tsdate.date(tss, Ne = 100000, mutation_rate = 1e-8)
    tsd.dump(dataset + "_" + region + "_" + level + "_n1.tsd")
    print(tsd)
    print(tsd.tables.sites)
    print("Tree sequence `tsd` dated and written to file\n")


## ----------------------------------------------------------------------------
## Tree stats

# LG02_pop: length x, mutations y, trees z
# LG03_pop: length x, mutations y, trees z
# LG12_pop: length x, mutations y, trees z


### ============================================================================
### 3. Calculate mean nucleotide diversity (tskit)

#!/usr/bin/python

import tskit


## Specify dataset
dataset = "phylo2e"
level = "pop"  # pop | sp | cld

# l = ["LG02"]
l = list([f"LG{i:02}" for i in range(1, 25)])
for region in l:


## Calculate nucleotide diversity
    tsd = tskit.load("../5_tsd/" + dataset + "_" + region + "_" + level + "_n1.tsd")
    n = len(tsd.tables.populations)


    pi = tsd.diversity(
        sample_sets = [
            tsd.samples(i) for i in range(n)
        ]
    )

    print(pi)


## Print array to file
    print(
        *pi, sep = '\n',
        file = open("pi_" + dataset + "_" + region + "_" + level + "_n1.txt", "w+")
    )
    print("Pi calculated for each pair of individuals and written to file")


### ============================================================================
### 4. Calculate TMRCA (tskit)

#!/usr/bin/python

import tskit
import numpy as np


## Specify dataset
dataset = "phylo2e"
level = "pop"  # pop | sp | cld

# l = ["LG02"]
l = list([f"LG{i:02}" for i in range(1, 7)])
for region in l:

    tsd = tskit.load("../5_tsd/" + dataset + "_" + region + "_" + level + "_n1.tsd")
    print(tsd)
    n = len(tsd.tables.populations)


## Calculate TMRCA (= average branch lengths between pairs of individuals / 2)
    np.set_printoptions(threshold = np.inf)

    tmrca = tsd.divergence(
        sample_sets = [tsd.samples(i) for i in range(n)],
        indexes = [(i, j) for i in range(n) for j in range(n)],
        mode = "branch",
        span_normalise = True
    ) / 2


## Print array to file
    print(
        *tmrca, sep = '\n',
        file = open("tmrca_" + dataset + "_" + region + "_" + level + "_n1.txt", "w+")
    )
    print("TMRCA calculated for each pair of individuals and written to file")


### ============================================================================
### 5. GNN proportions (tskit)

#!/usr/bin/python

import tskit
import json
import pandas as pd


## Load data
tsd = tskit.load("../5_tsd/phylo2e_LG02_cld_n1.tsd")

## Filter out samples with metadata clade = "out"
filtered_sample_nodes = [
    tsd.node(node_id)
    for node_id in tsd.samples()
    if json.loads(tsd.population(tsd.node(node_id).population).metadata)["clade"] != "out"
]


## Filter out samples with clade = "out" from samples_listed_by_population
samples_listed_by_population = []
for pop_id in range(tsd.num_populations):
    samples_for_population = tsd.samples(population=pop_id)
    filtered_samples_for_population = [
        node_id for node_id in samples_for_population
        if tsd.node(node_id) in filtered_sample_nodes
    ]
    samples_listed_by_population.append(filtered_samples_for_population)


## Calculate GNN proportions only for the remaining samples
gnn = tsd.genealogical_nearest_neighbours(
    [n.id for n in filtered_sample_nodes], samples_listed_by_population
)


## Tabulate GNN proportions using a Pandas dataframe
sample_ids = [n.id for n in filtered_sample_nodes]
sample_names = [json.loads(tsd.individual(n.individual).metadata)["name"] for n in filtered_sample_nodes]
sample_clades = [json.loads(tsd.population(n.population).metadata)["clade"] for n in filtered_sample_nodes]

gnn_table = pd.DataFrame(
    data=gnn,
    index=[
        pd.Index(sample_ids, name="Sample node"),
        pd.Index(sample_names, name="Individual"),
        pd.Index(sample_clades, name="clade"),
    ],
    columns=[json.loads(p.metadata)["clade"] for p in tsd.populations()],
)
print(gnn_table)


## Summarize GNN for all individuals of the same clade
print(gnn_table.groupby(level = "clade").mean())

## Print to file
tfile = open("../8_gnn/gnn_phylo2e_LG02_cld_n1.csv", "w")
tfile.write(gnn_table.to_csv())
tfile.close()


### ============================================================================
### x. Manipulate trees / further test code

# Return various data and metadata
print(ts.tables)

# Genetic diversity across genome
pi = tsd.diversity()
print(pi)
#> 0.0013

pi1 = tsd.diversity(sample_sets=tsd.samples(population=1))
print(pi1)
#> 0.0007

# Extract genetic data
for variant in tsd.variants():
    print(
        "Variable site", variant.site.id,
        "at genome position", variant.site.position,
        ":", [variant.alleles[g] for g in variant.genotypes],
    )

# Genealogical distance between tips
branch_diversity = tsd.diversity(mode="branch")
print("Av. genealogical dist. between pairs of tips is", branch_diversity,  tsd.time_units)
#> 124018 generations on average

# Returns no. of distinct trees in tree sequence
tsd.num_trees
#> 2330

# Returns no. of nodes (genomes)
tsd.num_nodes
#> 3065

tsd.samples(population=1)

# Iterate over trees
for tree in tsd.trees():
    print(f"Tree {tree.index} covers {tree.interval}")
    if tree.index >= 4:
        print("...")
        break
print(f"Tree {tsd.last().index} covers {tsd.last().interval}")

# Check coalescence
import time

elapsed = time.time()
for tree in tsd.trees():
    if tree.has_multiple_roots:
        print("Tree {tree.index} has not coalesced")
        break
else:
    elapsed = time.time() - elapsed
    print(f"All {tsd.num_trees} trees coalesced")
    print(f"Checked in {elapsed:.6g} secs")

# Metadata
print("Metadata for individual 0:", tsd.individual(0).metadata)

# Plot first tree
first_tree = tsd.first()
print("Total branch length in first tree is", first_tree.total_branch_length, tsd.time_units)
print("The first of", tsd.num_trees, "trees is plotted below")
first_tree.draw_svg(y_axis=True)  # plot the tree: only useful for small trees
