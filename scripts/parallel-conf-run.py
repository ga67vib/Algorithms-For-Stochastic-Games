def parallelize_confs(conf):
    print(conf)
    os.system("mkdir -p " + output_dir + "/" + conf)
    for model_count, model in enumerate(sorted(models.keys())):
        counting_str = ("Conf: [%s], Model: [%d/%d] - " %(conf, model_count + 1, len(models)))
        print("\t"+counting_str + model)
        for i in range(1, reps+1):
            print("\t\t"+str(i))
            rep_string = "" if reps == 1 else "_rep" + str(i)
            if exists(output_dir + "/" + conf + "/" + model + rep_string + ".log"):
                print("\t\tAlready there, skipping")
                continue
            prismParams = "" # "-javamaxmem 32g -javastack 16g"  # Change this appropriately
            command = "timeout 1m " + configurations[conf][0] + " " + \
                models[model] + " " + configurations[conf][1] + " " + prismParams + \
                " > " + output_dir + "/" + conf + "/" + model + rep_string + ".log"
            try:
                os.system(command)
            except:
                e = sys.exc_info()[0]
                print(e)