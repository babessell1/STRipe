import yaml
import itertools
import os
#from tabulate import tabulate


with open("config.yaml", 'r') as handle:
    try:
        config = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        print(exc)


def get_sample_dict(config, init=False):
    """
    Create sample dictionaries where samples are keys, and values are either haplotype, filetype, url, or number of files
    this should be done line by line with a generator assming the manifest is a csv at configfile["short_manifest"]
    """
    sample_dicts = {
        "short": {
            "haplotype": {}, # haplotype of sample
            "url": {}, # location to download
            "file_num": {}, # number of split bams
            "ext": {}, # bam or cram
            "iext": {} # index extension
        },
        "hifi": {
            "haplotype": {},
            "url": {},
            "file_num": {}, # bam or cram
            "ext": {}, # "bam" or "cram
            "iext": {} # index extension
        },
        "clr": {
            "haplotype": {},
            "url": {},
            "file_num": {}, # bam or cram
            "ext": {}, # "bam" or "cram
            "iext": {} # index extension
        },
        "assembly": {
            "haplotype": {},
            "url": {},
            "file_num": {},
            "ext": {}, # fasta
            "iext": {} # index extension
        },
    }
    with open(config["SHORT_MANIFEST"]) as handle:
        # line is: sample_name,haplotype,file_num,datatype,short_read_url
        for line in handle.readlines()[1:]:
            sample, haplotype, file_num, datatype, url = line.split(",")
            url = url.strip()
            if init:
                sample = sample + "." + file_num
                sample_dicts["short"]["file_num"][sample] = file_num
            else:
                sample_dicts["short"]["file_num"][sample] = 1
            sample_dicts["short"]["haplotype"][sample] = haplotype
            sample_dicts["short"]["url"][sample] = url
            sample_dicts["short"]["ext"][sample] = os.path.splitext(url)[1]
            iext = ".bai" if sample_dicts["short"]["ext"][sample] == ".bam" else "crai"
            sample_dicts["short"]["iext"][sample] = iext
            # touch file so it exist for snakemake
    with open(config["LONG_MANIFEST"]) as handle:
        # line is: sample_name,haplotype,file_num,long_read_url
        for line in handle.readlines()[1:]:
            sample, haplotype, file_num, datatype, url = line.split(",")
            url = url.strip()
            if datatype == "FASTA":
                dkey = "assembly"
                iext = "fai"
            elif datatype == "HIFI":
                dkey = "hifi"
                iext = None
            elif datatype == "CLR":
                dkey = "clr"
                iext = None
            else:
                raise ValueError(f"Datatype in the long manifest to either FASTA, HIFI, or CLR, not '{datatype}'")
            if init:
                sample = sample + "." + file_num
                sample_dicts[dkey]["file_num"][sample] = file_num
            else:
                sample_dicts[dkey]["file_num"][sample] = 1
            # check if value exists at sample key
            sample_dicts[dkey]["haplotype"][sample] = haplotype
            sample_dicts[dkey]["url"][sample] = url
            sample_dicts[dkey]["file_num"][sample] = file_num
            sample_dicts[dkey]["ext"][sample] = os.path.splitext(url)[1]
            if not iext:
                iext = "bai" if sample_dicts[dkey]["ext"][sample] == ".bam" else "crai"
            sample_dicts[dkey]["iext"][sample] = iext
    
    return sample_dicts

def get_samples(sample_dict, dtype):
    return list(sample_dict[dtype]["url"].keys())


def get_num(sample_dict, dtype):
    return [val for val in list(sample_dict[dtype]["file_num"].values())]


def get_ext(sample_dict, dtype):
    return [val for val in list(sample_dict[dtype]["ext"].values())]


def get_iext(sample_dict, dtype):
    return [val for val in list(sample_dict[dtype]["iext"].values())]

def string_to_list(stringed_list):  # Ex. "ABC, DEF, GHI"
    try:
        items = stringed_list.replace(" ", "").split(",")
    except:
        raise SyntaxError(
            'Please set multiple values for an item in config.yaml as a string in the format: '
            '"ABC, DEF, GHI, ..."')

    return items


def get_samp_id(manifest_path):
    with open(manifest_path) as handle:
        ids = [line.split(",")[0] for line in handle.readlines()]

    return ids


def get_url(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        for line in handle.readlines():
            print(line)
        seqtypes = [line.split(",")[-1] for line in handle.readlines()]

    return seqtypes


def get_haplo_type(long_manifest_path):
    with open(long_manifest_path) as handle:
        hp = [line.split(",")[1] for line in handle.readlines()]

    return hp


def get_chromosomes(stringed_list_chroms):
    chroms = string_to_list(stringed_list_chroms)

    return list(set(chroms))

