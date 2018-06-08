#!/usr/bin/env python

import subprocess
import sys
import os

from pbcommand.engine import run_cmd


def run(args):
    RUN_SEPARATELY = set(["converters", "filters", "laa", "barcoding"])
    output_dir = os.getcwd()
    if len(args) == 1:
        output_dir = args[0]
        assert os.path.isdir(output_dir), "Not a directory: %s"%output_dir
    module_dir = os.path.join(os.path.dirname(__file__), "pbcoretools", "tasks")
    for file_name in os.listdir(module_dir):
        if file_name.endswith(".py") and not file_name.startswith("_"):
            if file_name[:-3] in RUN_SEPARATELY:
                continue
            module_name = "pbcoretools.tasks.{m}".format(m=file_name[:-3])
            json_file = os.path.join(output_dir,
                "{m}_tool_contract.json".format(m=module_name))
            cmd = "python -m {m} --emit-tool-contract > {j}".format(
                m=module_name, j=json_file)
            run_cmd(cmd, sys.stdout, sys.stderr)
    for module_name in list(RUN_SEPARATELY):
        cmd = "python -m pbcoretools.tasks.{m} emit-tool-contracts -o {d}".format(m=module_name, d=output_dir)
        run_cmd(cmd, sys.stdout, sys.stderr)


if __name__ == "__main__":
    run(sys.argv[1:])
