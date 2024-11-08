#!/usr/bin/env python3

import json
import argparse


def passorfail(fastqc_reports, pass_criteria, no_fail_criteria=[], passvals=['PASS'], failvals=['FAIL']):
    for fastqc_report in fastqc_reports:
        with open(fastqc_report) as report:
            for line in report:
                val, crit, filename = line.strip('\n').split('\t')
                if crit in pass_criteria:
                    if val not in passvals:
                        return 'fail'
                if crit in no_fail_criteria:
                    if val in failvals:
                        return 'fail'
    return 'pass'    

def main():
    parser = argparse.ArgumentParser(
                        prog='assorted-sub-workflows/qc/bin/pass_or_fail_fastqc.py',
                        description='parse JSON input to define which item of a FastQC report should have what value for the short reads set to be considered \"pass\"')

    parser.add_argument('-p', '--pass_criteria', type=str, help="JSON file containing definition of an array specifying which item of a FastQC report (keys) are required to have the value PASS for the whole report to be considered a pass")
    parser.add_argument('-f', '--no_fail_criteria', type=str, help="JSON file containing definition of an array specifying which item of a FastQC report (keys) are required to NOT have the value FAIL for the whole report to be considered a pass")
    parser.add_argument('fastqc_reports', nargs='+', type=str)

    args = parser.parse_args()

    with open(args.pass_criteria) as pass_criteria_file:
        try:
            pass_criteria = tuple(json.load(pass_criteria_file))
        except TypeError:
            raise ValueError("bad formating of the JSON file given to '--pass_criteria' option; items should be listed in a iterable (tuple, list or set)")
    
    if args.no_fail_criteria:
        with open(args.no_fail_criteria) as no_fail_criteria_file:
            try:
                no_fail_criteria = tuple(json.load(no_fail_criteria_file))
            except TypeError:
                raise ValueError("bad formating of the JSON file given to '--no_fail_criteria' option; items should be listed in a iterable (tuple, list or set)")

    print( passorfail(args.fastqc_reports, pass_criteria=pass_criteria, no_fail_criteria=no_fail_criteria) )

