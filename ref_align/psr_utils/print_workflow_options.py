#!/usr/bin/env python

"""
Print all available workflow options
"""

import os
import sys
import glob
import yaml
import argparse


def read_yaml(filename):
    with open(filename, 'r') as stream:
        try:
            yamlD = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return yamlD


def pretty_name(targ):
    split = "  ####################  " 
    name = '\n\n' + split + targ + split + '\n'
    return name

def print_available_workflows_and_tools(paramsD,print_workflows=True,print_rules=False,only_rules=False):
    if only_rules:
        print_workflows = False
    if print_workflows:
        # print available workflows
        sys.stdout.write(pretty_name("Available Pseudoref Workflows") + '\n')
        psr_flows = paramsD.get('pseudoref_workflows')
        for workflow,tools in psr_flows.items():
            sys.stdout.write('\n  ' + workflow + ':\n\t')
            for t in tools.values():
                sys.stdout.write('\n\t'.join(t) + '\n\n')

    if print_rules:
        # print available rules
        sys.stdout.write(pretty_name("Advanced usage: all available rules") + '\n\t')
        rules_dir = paramsD['pseudoref_directories'].get('rules','rules')
        rules = glob.glob(os.path.join(rules_dir,'*','*.rule'))
        rule_names = [os.path.basename(r).split('.rule')[0] for r in rules]
        sys.stdout.write('\n\t'.join(rule_names) + '\n\n\n')


if __name__=="__main__":
    p = argparse.ArgumentParser()
    p.add_argument('paramsfile')
    p.add_argument('-w','--print_workflows',action='store_true',default=True)
    p.add_argement('-r','--print_rules',action='store_true',default=False)
    p.add_arguement('--only_rules',action='store_true',default=False)
    args=p.parse_args()
    # read params here, to enable passing in paramsD instead of paramsfile for run_pseudoref
    params=read_yaml(args.paramsfile)
    sys.exit(print_available_workflows_and_tools(params,args.print_workflows,args.print_rules,args.only_rules))
