#!/usr/bin/env python

import json

def parse_json(json_files):
    dict_parameters = {"CellMind":None,"CellBody":None,"World":None,"Lineage":None}
    if json_files:
        print('Imported override parameters from files:')
        for json_file_name in json_files:
            json_file = open(json_file_name)
            dict_parameters.update(json.load(json_file))
            json_file.close()
            print("\t" + json_file_name)
    return Parameters(dict_parameters)

class Parameters:
    def __init__(self, in_dict:dict):
        assert isinstance(in_dict, dict)
        for key, val in in_dict.items():
            if isinstance(val, (list, tuple)):
                setattr(self, key, [Parameters(x) if isinstance(x, dict) else x for x in val])
            else:
                setattr(self, key, Parameters(val) if isinstance(val, dict) else val)
