import os
import json
import time

def main():
    out_type = 1  # 1: population final state; 2: time series
    habitat_avail = [0.7, 1.0]
    resource_avail = 0.7
    habitat_length = 20
    init_pop_dens = 0.6
    end_time = 2000
    init_replic = 0
    replic_num = 10
    span = 1000
    topt = 16.7
    ctmax = 23.0
    sociality = [0, 1]  # 0: non-social; 1: social
    coop_efficiency = 10
    half_const = 2
    repro_max = 2
    cost_rate = 0.1
    init_coop_prop = 1
    random_seed = 42  # note that negative random seed will not be used

    temperature = [17, 17.2, 17.4, 17.6, 17.8, 18, 18.2, 18.4, 18.6, 18.8, 19, 19.2, 19.4, 19.6,19.8, 20, 20.2, 20.4, 20.6, 20.8, 21, 22, 23]

    input_list  = [[], []]
    output_list = [[], []]
    
    out_type_abbr = {1: "pfs", 2: "ts"}
    sociality_abbr = {0: "n", 1: "s"}
    
    if not os.path.isdir(f"data/{out_type_abbr[out_type]}"):
        os.makedirs(f"data/{out_type_abbr[out_type]}")
        
    for s in sociality:
        input_file_name_format = "data/{out_type}/{out_type}_{s}_{:02}xx_i.txt"
        
        for a in habitat_avail:
            input_file_name = input_file_name_format.format(int(a*10), out_type=out_type_abbr[out_type], s=sociality_abbr[s])
            input_list[s].append(input_file_name)
            output_list[s].append(input_file_name[:-6]+"_o.csv")
            with open(input_file_name, "w") as f:
                f.write(f"{out_type}\t{a}\t{habitat_length}\t{resource_avail}\t{init_pop_dens}\t{end_time}\t{init_replic}\t{replic_num}\t{span}\t{topt}\t{ctmax}\t{s}\t{coop_efficiency}\t{half_const}\t{repro_max}\t{cost_rate}\t{init_coop_prop}\t{random_seed}\n")
                
            for t in temperature:
                with open(input_file_name, "a") as f:
                    f.write(f"{t}\n")
                    
        with open(f"{out_type_abbr[out_type]}_{sociality_abbr[s]}_input_file_list.txt", "w") as fin:
            json.dump(input_list[s], fin)
        with open(f"{out_type_abbr[out_type]}_{sociality_abbr[s]}_output_file_list.txt", "w") as fout:
            json.dump(output_list[s], fout)
        
    
    
        
if __name__ == "__main__":
    main()