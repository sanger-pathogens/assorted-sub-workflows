#!/usr/bin/env python3

import json
import argparse


def passorfail(fastqc_reports, pass_criteria, fail_criteria=[], passvals=['PASS'], failvals=['FAIL']):
    for nffqc in fastqc_reports:
        with open(nffqc) as ffqc:
            fqc = {}
            for line in ffqc:
                val, crit, filename = line.strip('\n').split('\t')
                if crit in pass_criteria:
                    if val not in passvals:
                        return 'fail'
                if crit in fail_criteria:
                    if val in failvals:
                        return 'fail'
    return 'pass'    

def main():
    parser = argparse.ArgumentParser(
                        prog='assorted-sub-workflows/qc/bin/pass_or_fail_fastqc.py',
                        description='parse JSON input to define which item of a FastQC report should have what value for the short reads set to be considered \"pass\" for the PaM informatics dataset generator pipeline')

    parser.add_argument('-p', '--pass_criteria', type=str, help="JSON file containing definition of an array specifying which item of a FastQC report (keys) are required to have the value PASS for the whole report to be considered a pass")
    parser.add_argument('-f', '--fail_criteria', type=str, help="JSON file containing definition of an array specifying which item of a FastQC report (keys) are required to NOT have the value FAIL for the whole report to be considered a pass")
    parser.add_argument('fastqc_reports', nargs='+', type=str)

    args = parser.parse_args()

    with open(args.pass_criteria) as fpcriteria:
        try:
            pass_criteria = tuple(json.load(fpcriteria))
        except TypeError:
            raise ValueError("bad formating of the JSON file given to '--pass_criteria' option; items should be listed in a iterable (tuple, list or set)")
    
    if args.fail_criteria:
        with open(args.fail_criteria) as ffcriteria:
            try:
                fail_criteria = tuple(json.load(ffcriteria))
            except TypeError:
                raise ValueError("bad formating of the JSON file given to '--fail_criteria' option; items should be listed in a iterable (tuple, list or set)")

    print( passorfail(args.fastqc_reports, pass_criteria=pass_criteria, fail_criteria=fail_criteria) )

