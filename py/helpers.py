import yaml
import itertools
import os
#from tabulate import tabulate

with open("config.yaml", 'r') as handle:
    try:
        config = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        print(exc)


def string_to_list(stringed_list):  # Ex. "ABC, DEF, GHI"
    try:
        items = stringed_list.replace(" ", "").split(",")
    except:
        raise SyntaxError(
            'Please set multiple values for an item in config.yaml as a string in the format: '
            '"ABC, DEF, GHI, ..."')

    return items


def get_samp_id(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[0] for line in handle.readlines()]

    return ids


def get_seqtype(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        seqtypes = [line.split("\t")[1] for line in handle.readlines()]

    return seqtypes


def get_haplo_type(sample_info_filepath):
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[2] for line in handle.readlines()]

    return ids


def get_chromosomes(stringed_list_chroms):
    chroms = string_to_list(stringed_list_chroms)

    return list(set(chroms))

