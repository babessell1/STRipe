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

