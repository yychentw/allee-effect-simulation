import time
import json
import subprocess

def main():
    out_type = 1  # 1: population final state; 2: time series
    sociality = 1 # 0: non-social; 1: social
    
    out_type_abbr = {1: "pfs", 2: "ts"}
    sociality_abbr = {0: "n", 1: "s"}
    
    s_time = time.time()

    with open(f"{out_type_abbr[out_type]}_{sociality_abbr[sociality]}_input_file_list.txt", "r") as fin:
        input_list = json.load(fin)

    with open(f"{out_type_abbr[out_type]}_{sociality_abbr[sociality]}_output_file_list.txt", "r") as fout:
        output_list = json.load(fout)

    assert len(input_list) == len(output_list)

    for i in range(len(input_list)):
        ref_time = time.time()
        fi = open(input_list[i],  "r")
        fo = open(output_list[i], "w")
        subprocess.run(["./ca.out"], stdin=fi, stdout=fo)  # ca.out is the compiled program
        fi.close()
        fo.close()
        print(f"[{i+1}/{len(input_list)}] time: {time.strftime('%Mm%Ss', time.gmtime(time.time()-ref_time))}")

    e_time = time.time()
    print(f"Total consumed time: {time.strftime('%Hh%Mm%Ss', time.gmtime(e_time-s_time))}")


if __name__ == "__main__":
    main()