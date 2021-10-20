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
        '-referenceConf', type=str, default="D_BVI", help='If a conf violates by too much from the referenceConf, we consider it to be an Error. Default: D_BVI'
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

    with open(os.path.join(input_dir, "iters.csv"), 'r') as itersFile:
        with open(os.path.join(input_dir, "times.csv"), 'r') as timesFile:
            with open(os.path.join(input_dir, "values.csv"), 'r') as valuesFile:
                i = 0
                iters = list(csv.DictReader(itersFile))
                times = list(csv.DictReader(timesFile))
                values = list(csv.DictReader(valuesFile))

                keys = times[0].keys()
                confs = list(filter(lambda x : (x != 'Model') and (x != '#States'), keys))
                if (referenceConf not in keys):
                    raise Exception("Refrence Conf %s not in list of confs %s" % (referenceConf, str(confs)))
                for i in range(1, len(times)):
                    model_name = times[i]['Model']
                    for conf in confs:
                        if (conf == "WP"):
                            continue
                        conf_folder = os.path.join(input_dir, conf)
                        model_log_path = os.path.join(conf_folder, model_name+".log")
                        with open(model_log_path, 'r') as logFile:
                            logContent = logFile.read()
                            if ("OutOfMemory" in logContent):
                                times[i][conf] = 'M'
                                iters[i][conf] = 'M'
                                values[i][conf] = 'M'
                            elif ("Error" in logContent or "Exception" in logContent):
                                times[i][conf] = 'E'
                                iters[i][conf] = 'E'
                                values[i][conf] = 'E'
                            elif (not is_float(times[i][conf])): # If there is no time but also no error, it's a timeout
                                times[i][conf] = 'T'
                                iters[i][conf] = 'T'
                                values[i][conf] = 'T'
                            else:
                                if (not is_float(values[i][referenceConf])):
                                    print("Reference Conf %s cannot be used as reference for Model %s since it's result here is: %s" % (referenceConf, model_name, values[i][referenceConf]))
                                    continue
                                    
                                if (not is_float(values[i][conf])):
                                    if (values[i][conf] == 'V'):
                                        values[i][conf] = 'V'
                                        times[i][conf] = 'V'
                                        iters[i][conf] = 'V'
                                    print("Conf %s cannot be corrected for Model %s since it's result here is: %s" % (conf, model_name, values[i][conf]))
                                    continue

                                #If we got here, there should be no error and no timeout -> there must be result
                                refrenceValue = float(values[i][referenceConf])

                                confValue = float(values[i][conf])
                                #Theoretically we have precision 1E-6
                                if (abs(refrenceValue-confValue)>1E-6):
                                    values[i][conf] = 'V'
                                    times[i][conf] = 'V'
                                    iters[i][conf] = 'V'

    with open(os.path.join(input_dir, "values.csv"), 'w', encoding='utf8', newline='\n') as valuesFile:
        valuesFileContent = csv.DictWriter(valuesFile, fieldnames=keys)
        valuesFileContent.writeheader()
        valuesFileContent.writerows(values)

    with open(os.path.join(input_dir, "iters.csv"), 'w', encoding='utf8', newline='\n') as itersFile:
        itersFileContent = csv.DictWriter(itersFile, fieldnames=keys)
        itersFileContent.writeheader()
        itersFileContent.writerows(iters)

    with open(os.path.join(input_dir, "times.csv"), 'w', encoding='utf8', newline='\n') as timesFile:
        timesFileContent = csv.DictWriter(timesFile, fieldnames=keys)
        timesFileContent.writeheader()
        timesFileContent.writerows(times)

main()