import argparse
import os
import csv

def is_float(entry):
    try:
        float_entry = float(entry)
        return True
    except ValueError:
        return False

def main():

    parser = argparse.ArgumentParser(description = 'Fill blank spaces in times, value and iteration CSV Files')
    parser.add_argument(
        '-inputDirPath', type=str, default="./", help='Where are the csv-files? Default: ./'
    )
    parser.add_argument(
        '-referenceConf', type=str, default="T_LP_SI", help='If a conf violates by too much from the referenceConf, we consider it to be an Error. Default: D_BVI'
    )

    arguments = parser.parse_args()
    referenceConf = arguments.referenceConf
    input_dir = arguments.inputDirPath

    keys = []
    iters = []
    values = []
    times = []

    if(not os.path.exists(os.path.join(input_dir, "values.csv"))):
        raise Exception("PATH %s DOES NOT EXIST!" % (os.path.join(input_dir, "values.csv")))

    with open(os.path.join(input_dir, "times.csv"), 'r') as timesFile:
        i = 0
        times = list(csv.DictReader(timesFile))

        keys = times[0].keys()
        confs = list(filter(lambda x : (x != 'Model') and (x != '#States'), keys))
        if (referenceConf not in keys):
            raise Exception("Refrence Conf %s not in list of confs %s" % (referenceConf, str(confs)))
        for i in range(1, len(times)):
            model_name = times[i]['Model']

            if not is_float(times[i][referenceConf]):
                continue

            for conf in confs:
                if (conf != referenceConf):
                    continue
                conf_folder = os.path.join(input_dir, conf)
                model_log_path = os.path.join(conf_folder, model_name+".log")
                with open(model_log_path, 'r') as logFile:
                    time_to_solve_dtmc = 0
                    for log_line in logFile:
                        if ("Solve DTMC: " in log_line):
                            time_to_solve_dtmc = float(log_line.replace("Solve DTMC: ", ''))/1000.0
                    print(f'For Model {model_name} {referenceConf} needed {time_to_solve_dtmc} DTMC and {times[i][conf]} overall')
                    times[i][conf] = float(times[i][conf]) - time_to_solve_dtmc


    with open(os.path.join(input_dir, "times.csv"), 'w', encoding='utf8', newline='\n') as timesFile:
        timesFileContent = csv.DictWriter(timesFile, fieldnames=keys)
        timesFileContent.writeheader()
        timesFileContent.writerows(times)

main()