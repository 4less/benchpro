import subprocess as sb

def execute(params):
    cmd = "{} {} {}".format(
        params["program_name"],
        #params["params"],
        params["input"],
        params["output"]
    )
    sb.run([cmd], shell=True)

def run_rbenchpro(input_file, output_file):
    params = {
        "program_executer": "",
        "program_name": "benchpro.R",
        "input": "--input {}".format(input_file),
        "output": "--output {}".format(output_file)
    }
    return(params)
    