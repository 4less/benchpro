import subprocess as sb
from loguru import logger

def execute(params):
    cmd = "{} {} {} {} {}".format(
        params["program_name"],
        #params["params"],
        params["input"],
        params["input2"], #TODO: Change that, this is only a hackfix
        params["input_meta"], #TODO: Change that, this is only a hackfix
        params["output"]
    )
    print(cmd)
    logger.info(f"CMD: {cmd}")
    sb.run([cmd], shell=True)

def run_rbenchpro(input_file, input_detailed, input_meta, output_file):
    params = {
        "program_executer": "",
        "program_name": "benchpro.R",
        "input": "--input {}".format(input_file),
        "input2": "--input_detailed {}".format(input_detailed),
        "input_meta": "--meta {}".format(input_meta),
        "output": "--output {}".format(output_file)
    }
    return(params)
    