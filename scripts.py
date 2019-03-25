# The code below was written by A. Oliva

#This function will create the Dataset (DS) needed for those experimentations.
def Create_DataSet_Cluster(Dir,Nb_Of_Nucl,Model_Of_Analize,nb_tree,IndelRate):

    #Creating files where trees will be generated
    for j in xrange(1,nb_tree+1):
        subprocess.call(['mkdir', 'Arbre'+str(j)])
        os.chdir('Arbre%d' %j)
        #Tree generated with a R package
        Create_Tree()
        #This function will delete all superscripts present in the Newick string
        Get_Rid_Of_Scientific_Number()
        subprocess.call('cd ..', shell=True)

    for a in xrange(1,nb_tree+1):
        os.chdir('%s/Arbre%d' % (Dir,a))
        #Will return Rooted & Unrooted trees
        Res=Get_Both_Tree()
        Unrooted_Tree=Res[0]
        Rooted_Tree=Res[1]

        #Indelrate is only needed if you are created a dataset with gaps
        Use_Indelible(Unrooted_Tree=Unrooted_Tree, nb_nucl=Nb_Of_Nucl, Model=Model_Of_Analize, Gap=False,Indelrate=IndelRate)

        #Recover two list: L1=[Tips name, Tips sequences] , L2[Nodes name, Nodes sequences]
        Res = Get_Tips_Nodes(File='GTRout_TRUE.phy')
        Tips_list = Res[0]
        Nodes_Indelible = Res[1]


        #Will run all 3 software with the generated tree and sequences
        Use_PhyML(Unrooted_Tree = Unrooted_Tree, Tips_list= Tips_list, Model=Model_Of_Analize,Path= Path_PhyML)
        Use_RAXML(Rooted_Tree= Rooted_Tree, Tips_list = Tips_list,Model=Model_Of_Analize,Path=Path_RAXML)
        Use_PAML(Unrooted_Tree= Unrooted_Tree, Tips_list = Tips_list,Model= Model_Of_Analize,Path=Path_PAML)


########################################################################################################################
#will call the R script Create_Tree.R to create a random rooted Tree and also generate an unrooted Tree, needed for RAxML
def Create_Tree(Nb_Tips=0):
    command = 'Rscript'
    path2script = 'Create_Tree.R'
    arg = [str(Nb_Tips)]
    cmd = [command, path2script] + arg
    subprocess.check_output(cmd, universal_newlines=True)
    # The R script create two file, one with  the unrooted tree and another with the rooted tree
    lines=[]
    with open("Unrooted_Tree.tree", "r") as f:
        lines = f.read().split('\n')
    Unrooted_Tree = lines[0]
    #Now the Rooted Tree
    lines=[]
    with open("Rooted_Tree.tree", "r") as f:
        lines = f.read().split('\n')
    Rooted_Tree = lines[0]
    return (Unrooted_Tree,Rooted_Tree)

########################################################################################################################
#Will delete the exponent in the branch length from the newick string and will rewrite it
def Get_Rid_Of_Scientific_Number():
    with open("Unrooted_Tree.tree", "r") as f:
        lines = f.read().split('\n')
    lines = re.sub(r'(\d+\.\d+e-\d+)', lambda m: '{:.10f}'.format(float(m.group(0))), lines[0])
    with open('Unrooted_Tree.tree', 'w') as fa:
        fa.write(lines)
        lines=""
    with open("Rooted_Tree.tree", "r") as f:
        lines = f.read().split('\n')
    lines = re.sub(r'(\d+\.\d+e-\d+)', lambda m: '{:.10f}'.format(float(m.group(0))), lines[0])
    with open('Rooted_Tree.tree', 'w') as fa:
        fa.write(lines)

########################################################################################################################
#Will recover and return the rooted and unrooted trees
def Get_Both_Tree():
    with open("Unrooted_Tree.tree", "r") as f:
        lines = f.read().split('\n')
    Unrooted_Tree = lines[0]
    with open("Rooted_Tree.tree", "r") as f:
        lines = f.read().split('\n')
    Rooted_Tree = lines[0]
    return(Unrooted_Tree,Rooted_Tree)

########################################################################################################################
#Model= JC / K80 / HKY / GTR
#Run Indelible and recover the tips and nodes sequences generated along with their corresponding name
def Use_Indelible(Unrooted_Tree,nb_nucl,Model,Gap,Indelrate=0.01):
    print("################################################################\n"
          "Start the Indelible Run\n"
          "################################################################\n")

    #Run Indelible
    # Gap=False if you dont want indel / Gap = True if you want to include indels simulation
    #Indelrate will determine the rate of indels in the generated sequence
    Just_Run_Indelible(Unrooted_Tree,nb_nucl,Model,Gap,Indelrate)
    ModelOut_True=Model+"out_TRUE.phy"

    #This will read and return all tips name & corresponding sequences AND nodes name & corresponding sequences
    Res=Get_Tips_Nodes(File=ModelOut_True)
    Tips_Seq=Res[0]
    Nodes_Seq=Res[1]

    lines=[]
    #Clean all site with only GAPS and then rewrite the file
    #This would only be usefull if Gaps=True.
    #Site where there are only gaps, cannot be estimated by RAxML so we are removing them
    Clean_Align(Tips_Seq,Nodes_Seq)

    #Write over the old data just if some sequences changed (because of the Clean_Align function)
    Write_In_File(file=ModelOut_True, Tips_Seq=Tips_Seq, Nodes_Seq=Nodes_Seq)

    #Recover the tree from indelible
    with open("trees.txt","r") as f:
        lines = f.read().split('\n')
    for i, line in enumerate(lines):
        if re.findall(r'.*out\s+\w+\s+\d+\s+\d+\s+\d+\s+\d+.\d+\s\d+.\d+\s+\d+.\d+\s+(\(.*;)',line):
            Tree_Indelible=(re.findall(r'.*out\s+\w+\s+\d+\s+\d+\s+\d+\s+\d+.\d+\s\d+.\d+\s+\d+.\d+\s+(\(.*;)',line)[0])
    return(Tips_Seq,Nodes_Seq,Tree_Indelible)


########################################################################################################################
#This function will write the control file to start an Indelible run
#Model = JC / GTR - Gap = True / False
def Just_Run_Indelible(Unrooted_Tree,nb_nucl,Model,Gap,Indelrate):
    ModelExample=Model+'example'
    pModel='p'+Model
    ModelOut=Model+'out'
    ModelOut_True=Model+"out_TRUE.phy"
    if Model!="JC" and Model!="K80" and Model != "GTR" and Model != "HKY":
        print("Model unknown")
        return(0);

    #start to write the control file to run the indelible with the parameter needed
    #this file name is states as "control.txt"
    f = open('control.txt', 'w')
    #If we want to work with Nucleotide uncomment the next next line and comment the following one
    #f.write('[TYPE] NUCLEOTIDE 1\n\n')
    f.write('[TYPE] Codon 1\n\n')

    f.write('[MODEL]    {}\n'.format(ModelExample))
    if Model=="JC":
        f.write('[submodel]  {}\n\n'.format(Model))
    elif Model=="K80" or Model=="HKY":
        Value= input('Value of your parameter: ')
        f.write('[submodel]  {}\t{}\n\n'.format(Model,Value))
    else:
        #If here Model=GTR

        #Here we want to write the probas of different mutations; Transversion & Transition are not treated the same way

        #Value1 & 3 & 4 are transversion
        #Uncomment the following lines when we work with Nucleotides
        #Value1 = random.uniform(low= 0.5, high= 1.5, size= None)
        #Value3 = random.uniform(low= 0.5, high= 1.5, size= None)
        #Value4 = random.uniform(low= 0.5, high= 1.5, size= None)

        #Value 2 & 5 are the transition
        #Value2=0
        #Value5=0
        #while Value2 <1 or Value2 >10:
        #    Value2 = random.normal(loc= 4, scale= 4, size= None)
        #while Value5 <1 or Value2 >10:
        #    Value5 = random.normal(loc= 4, scale= 4, size= None)

        #write those probas
        #f.write('[submodel]  {}\t{}\t{}\t{}\t{}\t{}\n\n'.format(Model,Value1,abs(Value2),Value3,Value4,abs(Value5)))
        #uncomment the last lines and comment the next one if you work with Nucleotides
        f.write('[submodel]     2.5  0.5')


        #Then for the prior frequency of each nucleotide, we are picking them randomly (but sum of them must be =1)
        #sum=0
        #stat1=random.uniform(low= 0, high= 1, size= None)
        #sum=sum+stat1
        #stat2=random.uniform(low= 0, high= 1, size= None)
        #sum=sum+stat2
        #stat3=random.uniform(low= 0, high= 1, size= None)
        #sum=sum+stat3
        #stat4=random.uniform(low= 0, high= 1, size= None)
        #sum=sum+stat4
        #stat1=stat1/sum
        #stat2=stat2/sum
        #stat3=stat3/sum
        #stat4=stat4/sum

        #However we want a realistic dataset and a prior frequency for a nucleotide < 0.1 does not seem realistic for us
        #so we keep choosing probas until they are all >0.1
        #while(stat1<0.1 or stat2<0.1 or stat3<0.1 or stat4<0.1):
        #    sum = 0
        #    stat1 = random.uniform(low=0, high=1, size=None)
        #    sum = sum + stat1
        #    stat2 = random.uniform(low=0, high=1, size=None)
        #    sum = sum + stat2
        #    stat3 = random.uniform(low=0, high=1, size=None)
        #    sum = sum + stat3
        #    stat4 = random.uniform(low=0, high=1, size=None)
        #    sum = sum + stat4
        #    stat1 = stat1 / sum
        #    stat2 = stat2 / sum
        #    stat3 = stat3 / sum
        #    stat4 = stat4 / sum

        #then we write those frequencies
        #f.write('[statefreq] {} {} {} {}\n'.format(stat1,stat2,stat3,stat4))

    #if gap (indels) are wanted, we are adding them in the analyze , indelrate is given as an input argument and see Indelibl documentation for more infos about indelmodel
    if Gap==True:
        f.write('[indelmodel] POW 1.5 5\n')
        f.write('[indelrate] {}\n'.format(Indelrate))

    #Starting to write the Settings segment
    f.write('[SETTINGS]\n')
    #prints ancestral sequences in the same file as the leaf sequences.
    f.write(' [ancestralprint] SAME\n\n')
    #will print out sequences to file in either fasta, nexus, or phylip format respectively.
    f.write(' [output] PHYLIP\n\n')
    #TRUE will output inserted bases/residues that have been subsequently been deleted as * rather than - for easy identification.
    #FALSE will output all deletions as "-"
    f.write(' [markdeletedinsertions] FALSE\n\n')
    #These blocks can be used to set the parameters that control random rooted and unrooted guide tree generation or to define user-specified guide trees
    f.write('[TREE] arbre1 {}\n'.format(Unrooted_Tree))
    #This block type is used to set root lengths (nb nucl) and to choose which models are used with which trees.
    #nb_nucl can be considerated as Codons or AA. In our case we used nb_nucl=300 to create a 900 long dataset 
    f.write('[PARTITIONS] {}   [arbre1 {} {}]\n\n'.format(pModel,ModelExample,nb_nucl))
    #f.write('[PARTITIONS] {}   [arbre1 {} 300]\n\n'.format(pModel,ModelExample))
    # This block tells INDELible how many replicates to generate from each [PARTITIONS] block and what the output files should be named.
    f.write('[EVOLVE]\n')
    f.write('  {} 1 {}\n'.format(pModel,ModelOut))
    f.close()
    os.system("indelible")

########################################################################################################################
#Will read the file given as an input (PHYLIP format ONLY) and will return two list
#Tips_Seq=[Tip_name1,Tip_name2,Tip_name3,...][Sequence1,Sequence2,Sequence3,...]
#Nodes_Seq=[Node_name1,Node_name2,Node_name3,...][Sequence1,Sequence2,Sequence3,...]
def Get_Tips_Nodes(File):
    # Now I will recover all tips sequences in Tips_Seq and all Nodes sequences in Nodes_Seq
    tips_list = []
    node_list = []
    seq_list1 = []
    seq_list2 = []
    # Recover all sequences and name of sequences from tips and nodes
    with open(File, "r") as f:
        lines = f.read().split('\n')
    for i, line in enumerate(lines):
        if re.findall(r't\d+\s+[ACGT-]*\s*', line):
            tips_new = re.findall(r'(t\d+)\s+[ACGT-]*\s*', line)
            tips_list.insert(len(tips_list), tips_new[0])
            seq_new = re.findall(r't\d+\s+([ACGT-]*)\s*', line)
            seq_list1.insert(len(seq_list1), seq_new[0])
        if re.findall(r'N\d+\s+[ACGT-]*\s+', line):
            node_new = re.findall(r'(N\d+)\s+[ACGT-]*\s*', line)
            node_list.insert(len(node_list), node_new[0])
            seq_new = re.findall(r'N\d+\s+([ACGT-]*)\s*', line)
            seq_list2.insert(len(seq_list2), seq_new[0])
        if re.findall(r'ROOT.*', line):
            node_new = re.findall(r'(ROOT)\s+[ACGT-]*\s*', line)
            node_list.insert(len(node_list), node_new[0])
            seq_new = re.findall(r'ROOT\s+([ACGT-]*)\s*', line)
            seq_list2.insert(len(seq_list2), seq_new[0])

    Tips_Seq = [tips_list, seq_list1]
    Nodes_Seq = [node_list, seq_list2]
    return(Tips_Seq,Nodes_Seq)


########################################################################################################################
#Here I clean all site with only gaps (bc RAXML does not work if there is only gaps at one site)
def Clean_Align(Tips_Seq,Nodes_Seq):
    j=0
    while j < len(Tips_Seq[1][0]):
        if Tips_Seq[1][0][j]=='-':
            i=1
            ToErase = True
            while i < len(Tips_Seq[0]):
                ToErase=True
                if Tips_Seq[1][i][j]=='-':
                    i=i+1
                else:
                    ToErase=False
                    break
            i2=0
            i3=0
            #If only gaps, rewrite the sequences without this site
            if ToErase==True:
                while i2 < len(Tips_Seq[0]):
                    NewSeq = str(Tips_Seq[1][i2][:j]) + str(Tips_Seq[1][i2][j+1:])
                    Tips_Seq[1][i2] = NewSeq
                    i2=i2+1
                while i3 < len(Nodes_Seq[0]):
                    NewSeq = str(Nodes_Seq[1][i3][:j]) + str(Nodes_Seq[1][i3][j+1:])
                    Nodes_Seq[1][i3] = NewSeq
                    i3=i3+1
                j=j-1
        j=j+1

#########################################################################################################################
#Write sequences from tips and nodes in a given file (PHYLIP FORMAT ONLY)
def Write_In_File(file,Tips_Seq,Nodes_Seq):
    # Rewrite the file with the new values
    f = open(file, 'w')
    f.write("{}    {}\n".format(len(Tips_Seq[0]), len(Tips_Seq[1][0])))
    i = 0
    while i < len(Tips_Seq[0]):
        f.write("{}     {}      \n".format(Tips_Seq[0][i], Tips_Seq[1][i]))
        i = i + 1
    i = 0
    while i < len(Nodes_Seq[0]):
        f.write("{}     {}      \n".format(Nodes_Seq[0][i], Nodes_Seq[1][i]))
        i = i + 1
    f.close

########################################################################################################################
#Start a PhyML run
#Model HKY85 / JC69 / K80 / GTR
def Use_PhyML(Unrooted_Tree,Tips_list,Model,Path):
    print("################################################################\nStart the PhyML analyze\n"
          "################################################################\n")
    if Model!="JC69" and Model!="K80" and Model != "GTR" and Model != "HKY85":
        print("Model unknown")
        return(0);

    #ReWrite the sequence file with only the tips
    Write_Seq_Input(Tips_list)

    #Create a file with the Unrooted tree corresponding
    f = open ('TreeInputUnrooted.trees','w')
    f.write("{}".format(Unrooted_Tree))
    f.close()

    #launch PhyML with all parameters
    os.system("{}./src/phyml -i SequencesInput.phy -m {} -c 4 -a e -v 0.0 -u TreeInputUnrooted.trees -o lr --ancestral -b 0 --no_memory_check".format(Path,Model))

########################################################################################################################
#This function will create the file with all Tips sequences in PHYLIP format, needed for all runs
def Write_Seq_Input(Tips_list):
    f = open('SequencesInput.phy', 'w')
    f.write("{}  {}\n".format(len(Tips_list[0]), len(Tips_list[1][1])))
    i=0
    while i<len(Tips_list[0]):
        f.write("{}     {}\n".format(str(Tips_list[0][i]), Tips_list[1][i]))
        i = i+1
    f.close()

########################################################################################################################
#Start a RAxML run
def Use_RAXML(Rooted_Tree,Tips_list,Model,Path):
    print("################################################################\nStart the RAxML analyze\n"
          "################################################################\n")
    if Model!="JC69" and Model!="K80" and Model != "GTR" and Model != "HKY85":
        print("Model unknown")
        return(0);

    #ReWrite the sequence file with only the tips
    Write_Seq_Input(Tips_list)

    #Create a file with the Rooted tree corresponding
    f = open('TreeInputRooted.txt', 'w')
    f.write("{}".format(Rooted_Tree))
    f.close()

    #Launch RAxML with
    if Model=="GTR":
        os.system("/usr/bin/time -o TimeRAxML.txt -p {}/raxmlHPC -s SequencesInput.phy -n outputfile.txt -m GTRGAMMA -p 5 -f Ae -t TreeInputRooted.txt".format(Path))
    else:
        os.system("/usr/bin/time -o TimeRAxML.txt -p {}/raxmlHPC -s SequencesInput.phy -n outputfile.txt -m GTRGAMMA -p 5 -f Ae -t TreeInputRooted.txt --{}".format(Path,Model))


########################################################################################################################
#Create the "baseml.ctl" file needed to start a PAML run
def Use_PAML(Unrooted_Tree,Tips_list,Model,Path):
    print("################################################################\nStart the PAML analyze\n"
          "################################################################\n")
    if Model!="JC69" and Model!="K80" and Model != "REV" and Model != "HKY85" and Model != "GTR":
        print("Model unknown")
        return(0);

    #Convert Model into the associated number for PAML analysis
    if Model=="JC69":
        Model_nb=0
    elif Model=="K80":
        Model_nb = 1
    elif Model=="HKY85":
        Model_nb = 4
    elif Model=="REV" or Model=="GTR":
        Model_nb = 7

    #Rewrite a file with all tips sequences
    Write_Seq_Input(Tips_list)

    #Create a file with Unrooted file
    f = open('TreeInputUnrooted.txt', 'w')
    f.write("{}".format(Unrooted_Tree))
    f.close()

    #Then create the baseml file, which contain all parameter for the PAML run.
    f = open('baseml.ctl','w')
    #Fill the file with all parameters
    f.write("\tseqfile = SequencesInput.phy\n\ttreefile = TreeInputUnrooted.txt\n\n\toutfile = OutputFile\n\tnoisy = 2\n"
            "\tverbose = 2\n\trunmode = 0\n\n\tmodel = {}\n\n\tMgene = 0\n\n\tclock = 0\n\tfix_kappa = 0\n\tkappa = 5\n"
            "\n\tfix_alpha = 0\n\talpha = 0.5\n\tMalpha = 0\n\tncatG = 4\n\tnparK = 0\n\n\tnhomo = 0\n\tgetSE = 0\n\t"
            "RateAncestor = 1\n\n\tSmall_Diff = 7e-6\n\tcleandata = 0\n\tfix_blength = 1\n\tmethod = 0".format(Model_nb))
    f.close()
    #Launch PAML run
    os.system("/usr/bin/time -o TimePAML.txt -p {}/baseml".format(Path))
########################################################################################################################
