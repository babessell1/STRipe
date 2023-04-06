import yaml
import itertools
import os
#from tabulate import tabulate

with open("config/config.yaml", 'r') as handle:
    try:
        config = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        print(exc)


def string_to_list(stringed_list: str) -> list[str]:  # Ex. "ABC, DEF, GHI"
    try:
        items = stringed_list.replace(" ", "").split(",")
    except:
        raise SyntaxError(
            'Please set multiple values for an item in config.yaml as a string in the format: '
            '"ABC, DEF, GHI, ..."')

    return items


def get_samp_id(sample_info_filepath: str) -> list[str]:
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[0] for line in handle.readlines()]

    return ids


def get_seqtype(sample_info_filepath: str) -> list[str]:
    with open(sample_info_filepath) as handle:
        seqtypes = [line.split("\t")[1] for line in handle.readlines()]

    return seqtypes


def get_haplo_type(sample_info_filepath: str) -> list[str]:
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[2] for line in handle.readlines()]

    return ids


def get_chromosomes(stringed_list_chroms: str) -> list[str]:
    chroms = string_to_list(stringed_list_chroms)

    return list(set(chroms))


def exp_samp_id(samp_ids: list[str]) -> list[str]:  # expand sample ids to generate zipped wildcards with all permutes
    return samp_ids*\
        len(get_mei_type(config["MEI"]))*\
        len(get_chromosomes(config["CHROMOSOMES"]))


def exp_seqtype(barcodes: list[str]) -> list[str]:
    return barcodes*\
        len(get_mei_type(config["MEI"]))*\
        len(get_chromosomes(config["CHROMOSOMES"]))


def exp_haplo_type(haplo_types: list[str]) -> list[str]:
    return haplo_types*\
       len(get_chromosomes(config["CHROMOSOMES"]))*\
       len(get_mei_type(config["MEI"]))


def exp_chromosomes(chromosomes: list[str], germ=False, phased=False) -> list[str]:
    info_filepath_str = config["RAW_INFO_FILEPATH"] if not phased else config["PHASED_INFO_FILEPATH"]
    return list(
        itertools.chain.from_iterable(
            itertools.repeat(x, len(get_samp_id(info_filepath_str, germ=germ))
                             ) for x in chromosomes))*len(get_mei_type(config["MEI"]))


def exp_meis(meis: list[str], germ=False, phased=False) -> list[str]:
    info_filepath_str = config["RAW_INFO_FILEPATH"] if not phased else config["PHASED_INFO_FILEPATH"]
    return list(  # ['ALU', 'LINE'] -> ['ALU', 'ALU', 'ALU', 'LINE', 'LINE', 'LINE']
        itertools.chain.from_iterable(
            itertools.repeat(x, len(get_samp_id(info_filepath_str, germ=germ))
                             ) for x in meis))*len(get_chromosomes(config["CHROMOSOMES"]))


def get_bam_dir(haplo: str) -> str:
    return config["RAW_DIR"] if haplo == "raw" else config["PHASED_DIR"]



