#!/usr/bin/env python3

import json
import argparse


class PassOrFailFastqcError(Exception):
    pass


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


def get_criteria(criteria_filepath, option):
    with open(criteria_filepath) as criteria_file:
        try:
            return tuple(json.load(criteria_file))
        except json.JSONDecodeError:
            raise PassOrFailFastqcError(f"File given to '--{option}' option is not valid JSON format")
        except TypeError:
            raise PassOrFailFastqcError(f"Bad formatting of the JSON file given to '--{option}' option; items should be listed in an iterable (tuple, list or set)")



def main():
    parser = argparse.ArgumentParser(
                        prog='assorted-sub-workflows/qc/bin/pass_or_fail_fastqc.py',
                        description='parse JSON input to define which item of a FastQC report should have what value for the short reads set to be considered \"pass\"')

    parser.add_argument(
        '-p',
        '--pass_criteria',
        type=str,
        help="JSON file containing definition of an array specifying which item of a FastQC report (keys) are required to have the value PASS for the whole report to be considered a pass"
    )
    parser.add_argument(
        '-f',
        '--no_fail_criteria',
        type=str,
        help="JSON file containing definition of an array specifying which item of a FastQC report (keys) are required to NOT have the value FAIL for the whole report to be considered a pass"
    )
    parser.add_argument(
        'fastqc_reports',
        nargs='+',
        type=str,
        help="One or more summary.txt files output from FastQC"
    )

    args = parser.parse_args()

    pass_criteria = get_criteria(args.pass_criteria, 'pass_criteria')
    
    no_fail_criteria = []
    if args.no_fail_criteria:
        no_fail_criteria = get_criteria(args.no_fail_criteria, 'no_fail_criteria')

    print( passorfail(args.fastqc_reports, pass_criteria=pass_criteria, no_fail_criteria=no_fail_criteria) )

