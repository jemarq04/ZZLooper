def combine():
    skip = []
    with open("skip_vars.txt", "r") as infile:
        skip = [line.strip() for line in infile]

    for i,key in enumerate(["addr", "vars"]):
        eeee = []
        eemm = []
        mmmm = []
        eeee_eemm = []
        mmmm_eemm = []
        each = []

        f_eeee = []
        f_eemm = []
        f_mmmm = []

        with open("comp/eeee_%s.txt" % key, "r") as infile:
            f_eeee = infile.readlines()

        with open("comp/eemm_%s.txt" % key, "r") as infile:
            f_eemm = infile.readlines()

        with open("comp/mmmm_%s.txt" % key, "r") as infile:
            f_mmmm = infile.readlines()

        eeee = [x for x in f_eeee if x not in f_eemm and x not in f_mmmm]
        eemm = [x for x in f_eemm if x not in f_eeee and x not in f_mmmm]
        mmmm = [x for x in f_mmmm if x not in f_eeee and x not in f_eemm]
        eeee_eemm = [x for x in f_eeee if x in f_eemm and x not in f_mmmm]
        mmmm_eemm = [x for x in f_mmmm if x in f_eemm and x not in f_eeee]
        each = [x for x in f_eeee if x in f_eemm and x in f_mmmm]

        with open("comp/%s.txt" % key, "w") as outfile:
            outfile.write("//Shared branches\n")
            for line in each:
                if i==0 and any(var in line for var in skip):
                    outfile.write("if (_isT%iMC){\n%s}\n" % (1 if "_ntuple1" in line else 2, line))
                else:
                    outfile.write(line)
            if i==0:
                outfile.write("\nif (_channel == \"eeee\"){\n")
            else:
                outfile.write("\n//eeee\n")
            for line in eeee:
                if i==0 and any(var in line for var in skip):
                    outfile.write("if (_isT%iMC){\n%s}\n" % (1 if "_ntuple1" in line else 2, line))
                else:
                    outfile.write(line)
            if i==0:
                outfile.write('}\n\nelse if (_channel == "eemm"){\n')
            else:
                outfile.write("\n//eemm\n")
            for line in eemm:
                if i==0 and any(var in line for var in skip):
                    outfile.write("if (_isT%iMC){\n%s}\n" % (1 if "_ntuple1" in line else 2, line))
                else:
                    outfile.write(line)
            if i==0:
                outfile.write('}\n\nelse if (_channel == "mmmm"){\n')
            else:
                outfile.write("\n//mmmm\n")
            for line in mmmm:
                if i==0 and any(var in line for var in skip):
                    outfile.write("if (_isT%iMC){\n%s}\n" % (1 if "_ntuple1" in line else 2, line))
                else:
                    outfile.write(line)
            if i==0:
                outfile.write('}\n\nif (_channel == "eeee" || _channel == "eemm"){\n')
            else:
                outfile.write("\n//eeee or eemm\n")
            for line in eeee_eemm:
                if i==0 and any(var in line for var in skip):
                    outfile.write("if (_isT%iMC){\n%s}\n" % (1 if "_ntuple1" in line else 2, line))
                else:
                    outfile.write(line)
            if i==0:
                outfile.write('}\n\nif (_channel == "mmmm" || _channel == "eemm"){\n')
            else:
                outfile.write("\n//mmmm or eemm\n")
            for line in mmmm_eemm:
                if i==0 and any(var in line for var in skip):
                    outfile.write("if (_isT%iMC){\n%s}\n" % (1 if "_ntuple1" in line else 2, line))
                else:
                    outfile.write(line)
            if i==0:
                outfile.write("}\n")
        print("Combined channels into comp/%s.txt" % key)

if __name__ == "__main__":
    combine()
