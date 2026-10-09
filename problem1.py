def count_ATGC(dataset_adress):
    with open(dataset_adress, "r") as f:
        s = f.read()
        if len(s)>=1000:
            return
        print(f"{s.count("A")} {s.count("C")} {s.count("G")} {s.count("T")}")

count_ATGC("rosalind_dna.txt")
