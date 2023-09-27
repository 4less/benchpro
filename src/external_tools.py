import subprocess as sb

def execute(params):
    cmd = "{} {} {} {}".format(
        params["program_name"],
        #params["params"],
        params["input"],
        params["input2"], #TODO: Change that, this is only a hackfix
        params["output"]
    )
    print(cmd)
    sb.run([cmd], shell=True)

def run_rbenchpro(input_file, input_detailed, output_file):
    params = {
        "program_executer": "",
        "program_name": "benchpro.R",
        "input": "--input {}".format(input_file),
        "input2": "--input_detailed {}".format(input_detailed),
        "output": "--output {}".format(output_file)
    }
    return(params)
    