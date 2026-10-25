import csv
import sys


csv_list = list(sys.argv[1].split(","))
new_csv = sys.argv[2] 


def main():
    #print(csv_list)
    #reset content of new csv file
    with open(new_csv, "w") as csv_f:
        writer = csv.writer(csv_f)
        writer.writerow("")

    with open(new_csv, "a") as csv_f:
        writer = csv.writer(csv_f)
        for file in csv_list:
            #print(file)
            writer.writerow([file])
            with open(file, "r") as current_file:
                csvreader = csv.reader(current_file)
                for row in csvreader:
                    writer.writerow(row)

    
  


if __name__=="__main__":
    main()